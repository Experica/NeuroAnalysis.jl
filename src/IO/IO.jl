using FileIO,JLD2,MAT,YAML,Mmap,LoopVectorization
import FileIO: save

include("SpikeGLX.jl")
include("Imager.jl")

"""
Drop array dims according to `MATLAB` convention.

- scalar: scalar instead of 1x1 array
- rvector: vector instead of 1xN array
- cvector: vector instead of Nx1 array
"""
dropmatdim!(x;scalar = true, rvector = true, cvector = true) = x
function dropmatdim!(x::Dict;scalar = true, rvector = true, cvector = true)
    foreach(k->x[k]=dropmatdim!(x[k],scalar=scalar,rvector=rvector,cvector=cvector),keys(x))
    return x
 end
function dropmatdim!(x::Array;scalar = true, rvector = true, cvector = true)
    if ndims(x) == 2
        nr,nc = size(x)
        if nr==1 && nc==1 && scalar
            return dropmatdim!(x[1,1],scalar=scalar,rvector=rvector,cvector=cvector)
        elseif nr==1 && nc>1 && rvector
            x = dropdims(x,dims=1)
        elseif nr>1 && nc==1 && cvector
            x = dropdims(x,dims=2)
        end
    end
    if x isa Array{Any}
        foreach(i->x[i]=dropmatdim!(x[i],scalar=scalar,rvector=rvector,cvector=cvector),eachindex(x))
    end
    return x
end

"""
Read variables of a `MATLAB` MAT file into a Dictionary

1. f: MAT file path
2. vars: variable names in the varset to read

- varset: variable names to select in vars
- scalar: scalar instead of 1x1 matrix
- rvector: vector instead of 1xN matrix
- cvector: vector instead of Nx1 matrix
"""
function readmat(f::AbstractString,vars...;varset=["spike","lfp","digital","analog1k","image"],scalar = true, rvector = true, cvector = true)
    if isempty(vars)
        d=matread(f)
    else
        d=Dict()
        matopen(f) do file
            fvs = names(file)
            nonselectvars=setdiff(fvs,varset)
            selectvars = intersect(varset,vars)
            for i in setdiff(selectvars,fvs)
                @warn """Selected variable "$i" not found in "$f"."""
            end
            readvars=union(nonselectvars,intersect(fvs,selectvars))
            if !isempty(readvars)
                for v in readvars
                    d[v]=read(file,v)
                end
            end
        end
    end
    if scalar || rvector || cvector
        d = dropmatdim!(d,scalar=scalar,rvector=rvector,cvector=cvector)
    end
    return d
end

"Read `DataExport` Metadata MAT File"
readmeta(f::AbstractString) = DataFrame(readmat(f)["Tests"])

"Load images in directory"
function loadimageset(dir;name=[],n=typemax(Int),alpha=false)
    isdir(dir) || error("Invalid Directory")
    for (root, dirs, files) in walkdir(dir)
        if isempty(name)
            n = min(n,length(files))
            n<=0 && error("n <= 0")
            imgs=load.(joinpath.([root],string.(1:n).*splitext(files[1])[2]))
        else
            imgs=load.(joinpath.([root],string.(Int.(name)).*splitext(files[1])[2]))
        end
        if alpha
            imgs = map(i->coloralpha.(i),imgs)
        end
        return imgs
    end
end

"Default digital input channels"
const CondDCh=1
const MarkDCh=2
const StartDCh=3
const StopDCh=4
const Bits16DCh=5
"Get digital channel time and value."
function digitaldata(dataset::Dict,ch)
    chidx = findfirst(dataset["digital"]["channel"].==ch)
    if !isempty(chidx)
        return dataset["digital"]["time"][chidx],dataset["digital"]["data"][chidx]
    else
        return [],[]
    end
end
"Parse bit event time and value."
function digitalbit(dt,dv,bits...)
    n = length(dt)
    maxbit = maximum(bits)
    bn = length(bits)
    bt = [zeros(Float64,n+1) for _ in 1:bn]
    bv = [zeros(UInt8,n+1) for _ in 1:bn]
    j=ones(Int,bn)
    for i in 1:n
        ds = digits(UInt8,dv[i],2,maxbit)
        for b in 1:bn
            v=ds[bits[b]]
            if v != bv[b][j[b]]
                j[b]=j[b]+1
                bt[b][j[b]]=dt[i]
                bv[b][j[b]]=v
            end
        end
    end
    for b in 1:bn
        bt[b]=bt[b][2:j[b]]
        bv[b]=bv[b][2:j[b]]
    end
    if bn==1
        return bt[1],bv[1]
    end
    return bt,bv
end

"Read and Prepare Dataset in `MATLAB` MAT File"
prepare(f::AbstractString,vars...;spikesorter="kilosort")=prepare!(readmat(f,vars...);spikesorter)
function prepare!(d::Dict;spikesorter="kilosort")
    if haskey(d,"secondperunit")
        settimeunit(d["secondperunit"])
    end
    if haskey(d,"sourceformat")
        fmt = d["sourceformat"]
        if fmt=="Ripple"
            d=prepare_ripple!(d)
        elseif fmt=="OI"
            d=prepare_oi!(d)
        elseif fmt=="Experica"
            d=prepare_experica!(d)
        end
    end
    if haskey(d,"imecindex")
        imecspike=Dict()
        for i in d["imecindex"]
            s = "spike$(i)_$spikesorter"
            if haskey(d,s)
                fun = Symbol("unitspike_",spikesorter)
                imecspike[i] = @eval $fun($d[$s],syncindex=$i,sortspike=true)
            end
        end
        d = mergeimecspike!(imecspike,d)
    end
    return d
end
"Prepare `Ripple` Data"
function prepare_ripple!(d::Dict)
    if haskey(d,"spike")
        st=d["spike"]["time"]
        su=d["spike"]["unitid"]
        for i in 1:length(st)
            if ndims(st[i])==0
                st[i]=[st[i]]
                su[i]=[su[i]]
            end
        end
    end
    return d
end
"Prepare `Optical Imaging` Block Data"
function prepare_oi!(d::Dict)
    return d
end
"Prepare `Experica` Experiment Data"
function prepare_experica!(d::Dict)
    return d
end

"Extract state time of condtest"
function statetime(ct::Dict;statetype::AbstractString="CONDSTATE",state::AbstractString="COND")
    if haskey(ct,statetype)
        filter!(l->!isempty(l),map(i->begin
        t = filter!(k->!isempty(k),map(j->haskey(j,state) ? j[state] : [],i))
        length(t)==1 ? t[1] : t
    end,ct[statetype]))
else
    []
end
end

function getparam(d::Dict{String,Any},name,fromobject="")
    for k in keys(d)
        startswith(k,name) && endswith(k,fromobject) && return d[k]
    end
    nothing
end

"Get `CondTestCond` DataFrame"
function condtestcond(ctcd::Dict)
    ctc = DataFrame(ctcd)
    return ctc
end
"Get `CondTest` DataFrame"
function condtest(ctd::Dict)
    ct = DataFrame(ctd)
    return ct
end
"Get `CondTest` and `CondTestCond` DataFrame"
function ctctc(ctd::Dict,ctcd::Dict)
    ct = condtest(ctd)
    ctc = condtestcond(ctcd)
    return ct,ctc
end
function ctctc(ex::Dict)
    return ctctc(ex["CondTest"],ex["CondTestCond"])
end

"Map time to data timeline, optionally add latency."
function maptodatatime(x,ex::Dict;addlatency=true)
    t=x*(1+ex["TimerDriftSpeed"])+ex["t0"]
    if addlatency
        t+=ex["Latency"]
    end
    return t
end

"Map time between ref and target, according to the timing diff between target and ref of the same sync signal"
function ref2sync(t,dataset::Dict,si="0")
    haskey(dataset,"syncdt") || return t
    tt = "syncdiff$si"
    if haskey(dataset,tt)
        ref2sync(t,dataset[tt],dataset["syncdt"])
    else
        @warn "No syncdiff for Sync$si, skip timing correction.";return t
    end
end
function sync2ref(t,dataset::Dict,si="0")
    haskey(dataset,"syncdt") || return t
    tt = "syncdiff$si"
    if haskey(dataset,tt)
        sync2ref(t,dataset[tt],dataset["syncdt"])
    else
        @warn "No syncdiff for Sync$si, skip timing correction.";return t
    end
end
ref2sync(t,syncdiff,syncdt=500) = t .+ @view syncdiff[clamp!(round.(Int,t./syncdt),1,length(syncdiff))]
sync2ref(t,syncdiff,syncdt=500) = t .- @view syncdiff[clamp!(round.(Int,t./syncdt),1,length(syncdiff))]

"""
Epochs of `Neuropixels` data stream (Channels x Samples), 
optionally gain corrected(voltage), line noise(60,120,180Hz) removed, 
bandpass filtered and common average referenced
"""
function epochsamplenp(x,fs,epochs,chs;meta=nothing,bandpass=[1,100],whiten=nothing,car=nothing)
    if isnothing(meta)
        g = i -> copy!(similar(i,Float64),i)
    else
        g = i -> gaincorrectnp(i,meta)
    end
    if isnothing(bandpass)
        f = identity
    else
        if bandpass[1] <= 250
            f = i -> hlpass!(rmline!(i,fs),fs,high=bandpass[1],low=bandpass[2])
        else
            f = i -> hlpass!(i,fs,high=bandpass[1],low=bandpass[2])
        end
    end
    if isnothing(car)
        r = identity
    else
        r = i -> i .- car(i,dims=1)
    end
    fun = r ∘ f ∘ g
    epochsample(x,fs,epochs,chs;fun,whiten)
end

"""
Organize each spiking unit from `Kilosort` result.
"""
function unitspike_kilosort(data::Dict;syncindex="0",sortspike::Bool=true)
    # for each unit
    unitid = data["clusterid"]
    unitgood = data["clustergood"].==1
    unitindex = [data["cluster"].==i for i in unitid]
    unitspike = map(i->data["time"][i],unitindex) # spike train
    unittemplate = map(i->data["template"][i],unitindex) # template id
    unitamplitude = map(i->data["amplitude"][i],unitindex) # scaled template amplitude

    chmap = data["chanmap"]
    chposition = data["channelposition"]
    unittemplatesindex = unique.(unittemplate)
    # mean position of unit templates
    unitposition = vcat(map(i->mean(data["templatesposition"][i,:],dims=1),unittemplatesindex)...)
    # mean amplitude of unit templates
    unittemplateamplitude = map(i->mean(data["templatesamplitude"][i]),unittemplatesindex)
    # first found template waveform
    unittemplatewaveform = data["templateswaveform"][first.(unittemplatesindex),:]
    # first found template waveform feature
    unittemplatefeature = Dict(k=>data["templateswaveformfeature"][k][first.(unittemplatesindex)] for k in keys(data["templateswaveformfeature"]))
    unitwaveform = haskey(data,"clusterwaveform") ? data["clusterwaveform"] : unittemplatewaveform
    unitfeature = haskey(data,"clusterwaveformfeature") ? data["clusterwaveformfeature"] : unittemplatefeature


    sortspike && foreach(sort!,unitspike)
    return Dict("unitid"=>unitid,"unitgood"=>unitgood,"unitspike"=>unitspike,"chposition"=>chposition,"unitposition"=>unitposition,
                "unitwaveform"=>unitwaveform,"unitfeature"=>unitfeature,"unittemplatewaveform"=>unittemplatewaveform,"unittemplatefeature"=>unittemplatefeature,
                "unittemplateamplitude"=>unittemplateamplitude,"isspikesorted"=>sortspike,"unitsync"=>fill(syncindex,size(unitspike)))
end

"""
Organize each spiking unit from `Kilosort3` result.
"""
function unitspike_kilosort3(data::Dict;syncindex="0",sortspike::Bool=true)
    # for each unit
    unitid = data["clusterid"]
    unitgood = data["clustergood"].==1
    unitindex = [data["cluster"].==i for i in unitid]
    unitspike = map(i->data["time"][i],unitindex) # spike train
    unittemplate = map(i->data["template"][i],unitindex) # template id
    unitamplitude = map(i->data["amplitude"][i],unitindex) # scaled template amplitude

    datapath = data["datapath"] # concat whiten drift corrected binary file
    t0 = data["t0"] # start time in `datapath` file
    fs = data["fs"]
    nsample = Int(data["nsample"]) # number of samples in `datapath` file
    nch = Int(data["nch"]) # number of channels in `datapath` file
    chmapraw = data["chanmapraw"] # map each channel in `datapath` file to raw binary file
    chposition = data["channelposition"]
    unittemplatesindex = unique.(unittemplate)
    # mean position of templates with which unit spikes are extracted
    unittemplateposition = vcat(map(i->mean(data["templatesposition"][i,:],dims=1),unittemplatesindex)...)
    # mean amplitude of unit templates
    unittemplatemeanamplitude = map(i->mean(data["templatesmeanamplitude"][i]),unittemplatesindex)
    # first found template waveform
    unittemplatewaveform = data["templateswaveform"][first.(unittemplatesindex),:]
    # first found template waveform feature
    unittemplatefeature = Dict(k=>data["templateswaveformfeature"][k][first.(unittemplatesindex)] for k in keys(data["templateswaveformfeature"]))
    if haskey(data,"clusterwaveform")
        unitwaveforms = data["clusterwaveforms"]
        unitposition = data["clusterposition"]
        unitwaveform = data["clusterwaveform"]
        unitfeature = data["clusterwaveformfeature"]
    else
        unitwaveform=unitfeature=unitposition=unitwaveforms=nothing
    end
    w = data["whiteningmatrix"]
    winv = data["whiteningmatrixinv"]
    wraw = data["whiteningmatrixraw"]
    winvraw = data["whiteningmatrixinvraw"]
    qm = data["qm"]

    sortspike && foreach(sort!,unitspike)
    return Dict("unitid"=>unitid,"unitgood"=>unitgood,"unitspike"=>unitspike,"datapath"=>datapath,"t0"=>t0,"fs"=>fs,
        "nsample"=>nsample,"nch"=>nch,"chmapraw"=>chmapraw,"chposition"=>chposition,"unittemplateposition"=>unittemplateposition,
        "unittemplatemeanamplitude"=>unittemplatemeanamplitude,"unittemplatewaveform"=>unittemplatewaveform,"unittemplatefeature"=>unittemplatefeature,
        "unitwaveforms"=>unitwaveforms,"unitposition"=>unitposition,"unitwaveform"=>unitwaveform,"unitfeature"=>unitfeature,
        "w"=>w,"winv"=>winv,"wraw"=>wraw,"winvraw"=>winvraw,"qm"=>qm,"isspikesorted"=>sortspike,"unitsync"=>fill(syncindex,size(unitspike)))
end

function mergeimecspike!(imecspike,d)
    n = length(imecspike)
    if n == 0
        @warn "No IMEC spike data"
    elseif n==1
        d["spike"] = first(values(imecspike))
    else
        @info "Merging spike from multiple IMEC, not implemented"
    end
    return d
end

"convert `Matlab` struct of array to `DataFrame`"
matdictarray2df(d) = DataFrame(Any[squeeze(v,2) for v in values(d)],[Symbol(k) for k in keys(d)])

"Regular expression to match `VLab` data file names"
function vlabregex(testtype;subject="[0-9]",maxch=1,cell="[a-z]",maxrepeat=3,format="mat")
    mr = lpad(maxrepeat,2,0)
    Regex("^$subject+[A-Za-z]+[1-$maxch]$cell$testtype[0-$(mr[1])][0-$(mr[2])]\\.$format")
end

"Get matched file names in path"
function getdatafile(testtype;subject="[0-9]",path="./data")
    matchfile(vlabregex(testtype,subject=subject),path=path)
end

"Regular Expression to match Optical Imaging `VDAQ` block file name"
function oifileregex(;subject="[A-Za-z0-9]",session="[A-Za-z0-9]",experimentid="[0-9]",format="mat")
    Regex("^$subject+_$session*_?E$experimentid+B[0-9]+[.]$format")
end

"Get Optical Imaging `VDAQ` matched files in directory"
function getoifile(;subject="[A-Za-z0-9]",session="[A-Za-z0-9]",experimentid="[0-9]",format="mat",dir="",adddir=false)
    matchfile(oifileregex(subject=subject,session=session,experimentid=experimentid,format=format),dir=dir,adddir=adddir)
end

"Regular Expression to match `Experica` data file name"
function expericafileregex(;subject="[A-Za-z0-9]",session="[A-Za-z0-9]",site="[A-Za-z0-9]",test="[A-Za-z0-9]",maxrepeat=10,format="mat")
    mr = lpad(maxrepeat,2,"0");d2=mr[1];d1=parse(Int,d2)>0 ? 9 : mr[2]
    Regex("^$subject+_$session*_?$site*_?$test+_[0-$d2]?[0-$d1][.]$format")
end

"Get matched `Experica` files in directory"
function getexpericafile(;subject="[A-Za-z0-9]",session="[A-Za-z0-9]",site="[A-Za-z0-9]",test="[A-Za-z0-9]",maxrepeat=10,format="mat",dir="",adddir=false)
    matchfile(expericafileregex(subject=subject,session=session,site=site,test=test,maxrepeat=maxrepeat,format=format),dir=dir,adddir=adddir)
end

"Get matched file names in directory, optionally join directory to get file path"
function matchfile(pattern::Regex;dir="",join::Bool=false)
    fs = filter!(f->occursin(pattern,f),readdir(dir))
    if join
        fs=joinpath.(dir,fs)
    end
    return fs
end

function save(filepath::Union{AbstractString,IO},cm::Dict)
    ext = splitext(filepath)[2]
    if ext == ".yaml"
        data = deepcopy(cm)
        foreach(n->data[n]=Dict(k=>cm isa String ? cm : vec.(cm) for (k,cm) in pairs(data[n])),keys(data))
        YAML.write_file(filepath,data)
    elseif ext == ".mat"
        data = deepcopy(cm)
        foreach(n->data[n]=Dict(k=>cm isa String ? cm : mapreduce(c->vec(c)',vcat,cm) for (k,cm) in pairs(data[n])),keys(data))
        matwrite(filepath,data)
    end
end
