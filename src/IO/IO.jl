using MAT,Query,FileIO

include("SpikeGLX.jl")

export readmat,readmeta,mat2julia!,loadimageset,CondDCh,MarkDCh,StartDCh,StopDCh,Bits16DCh,digitaldata,digitalbit,
prepare,prepare!,prepare_ripple!,prepare_oi!,prepare_experica!,
statetime,getparam,condtestcond,condtest,ctctc,maptodatatime,
oifileregex,getoifile,expericafileregex,getexpericafile,matchfile,querymeta,
subrm,reshape2ref,unitfyspike

"Read variables in `Matlab` MAT format data"
function readmat(f::AbstractString,vars...;optvars=["spike","lfp","digital","analog1k","image"])
    if isempty(vars)
        d=matread(f)
    else
        d=Dict()
        matopen(f) do file
            fvs = names(file)
            reqvars=setdiff(fvs,optvars)
            ovs = intersect(optvars,vars)
            for ivov in setdiff(ovs,fvs)
                @warn """variable "$ivov" not found."""
            end
            vs=union(reqvars,intersect(fvs,ovs))
            if !isempty(vs)
                for v in vs
                    d[v]=read(file,v)
                end
            end
        end
    end
    return d
end

"Read Metadata MAT File"
function readmeta(f::AbstractString)
    d = readmat(f)["Tests"]
    mat2julia!(d)
    DataFrame(d)
end

"Convert `Matlab` type to `Julia` type with proper dimention"
function mat2julia!(x;isscaler = true,isvector = true)
    if x isa Dict
        for k in keys(x)
            x[k]=mat2julia!(x[k],isscaler=isscaler,isvector=isvector)
        end
    elseif x isa Array
        if ndims(x)==2
            s = size(x)
            if s[1]==1 && s[2]==1
                x = isscaler ? x[1,1] : isvector ? dropdims(x,dims=2) : x
            elseif s[1]==1 && s[2]>1
                x = isvector ? dropdims(x,dims=1) : x
            elseif s[1]>1 && s[2]==1
                x = isvector ? dropdims(x,dims=2) : x
            end
        end
        if x isa Array{Any} || x isa Array{Array} || x isa Array{Dict}
            for i in 1:length(x)
                x[i]=mat2julia!(x[i],isscaler=isscaler,isvector=isvector)
            end
        elseif x isa Dict
            for k in keys(x)
                x[k]=mat2julia!(x[k],isscaler=isscaler,isvector=isvector)
            end
        end
    end
    return x
end

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

"Prepare Dataset"
prepare(f::AbstractString,vars...)=prepare!(readmat(f,vars...),true)
function prepare!(d::Dict,ismat=true)
    ismat && mat2julia!(d)
    if haskey(d,"secondperunit")
        settimeunit(d["secondperunit"])
    end
    if haskey(d,"sourceformat")
        sf = d["sourceformat"]
        if sf=="Ripple"
            d=prepare_ripple!(d)
        elseif sf=="OI"
            d=prepare_oi!(d)
        elseif sf=="Experica"
            d=prepare_experica!(d)
        end
    end
    if !haskey(d,"spike")
        if haskey(d,"spike_kilosort")
            d["spike"] = unitfyspike(d["spike_kilosort"])
        end
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

function getparam(d::Dict,name,fromobject="")
    for k in keys(d)
        startswith(k,name) && endswith(k,fromobject) && return d[k]
    end
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

function subrm(rm,fs,epochs;meta=[],chs=1:size(rm,1),bandpass=[1,100])
    nepoch = size(epochs,1)
    epochis = floor.(Int,epochs.*fs)
    minepochlength = minimum(diff(epochis,dims=2))
    ys=Array{Float64}(undef,length(chs),minepochlength,nepoch)
    for i in 1:nepoch
        y = rm[chs,range(max(1,epochis[i,1]),length=minepochlength)]
        if !isempty(meta)
            y=gaincorrectim(y,meta)
            if !isempty(bandpass)
                rmline!(y,fs)
                y = hlpass(y,low=bandpass[2],high=bandpass[1],fs=fs)
            end
        end
        ys[:,:,i] = y
    end
    return nepoch==1 ? dropdims(ys,dims=3) : ys
end

function reshape2ref(ys,refmask;cleanref=true)
    nrow,ncol=size(refmask)
    yss=Array{Float64}(undef,nrow,ncol,size(ys)[2:end]...)
    for c in 1:ncol
        yss[:,c,:,:] = ys[c:ncol:end,:,:]
    end
    if cleanref
        for (r,c) in Tuple.(findall(refmask))
            yss[r,c,:,:] = (yss[r+1,c,:,:] .+ yss[r-1,c,:,:]) / 2 # Local Average
        end
    end
    return yss
end

"get spiking units info"
function unitfyspike(data::Dict)
    # kilosort results
    rawspiketime = data["time"]
    rawspiketemplate = data["template"]
    rawspikecluster = data["cluster"]
    rawamplitude = data["amplitude"]

    templates = data["templates"] # nTemplates x nTimePoints x nChannels
    chposition = data["channelposition"]
    whiteninv = data["whiteningmatrixinv"]

    templatesunwhiten = zeros(size(templates))
    for t in 1:size(templates,1)
        templatesunwhiten[t,:,:] = templates[t,:,:]*whiteninv
    end
    templatesunwhiten_height = dropdims(map(i->-(-(i...)),extrema(templatesunwhiten,dims=2)),dims=2) # unwhiten template height between trough to peak, nTemplates x nChannels

    # each unit
    unitid = data["clusterid"]
    unitgood = data["good"].==1
    unitindex = [rawspikecluster.==i for i in unitid]
    unitspike = map(i->rawspiketime[i],unitindex)
    unitamplitude = map(i->rawamplitude[i],unitindex)
    unittemplate = map(i->rawspiketemplate[i][1],unitindex)
    unittemplates = map(i->templates[i,:,:],unittemplate)
    unittemplatesunwhiten = map(i->templatesunwhiten[i,:,:],unittemplate)
    unittemplatesunwhiten_height = map(i->templatesunwhiten_height[i,:],unittemplate)

    unitposition = vcat(map(w->sum(w.*chposition,dims=1)/sum(w),unittemplatesunwhiten_height)...) # center of mass from all weighted unit template channel positions

    return Dict("unitid"=>unitid,"unitgood"=>unitgood,"unitspike"=>unitspike,"chposition"=>chposition,"unitposition"=>unitposition)
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

"Get matched files in directory, optionally add directory to get file path"
function matchfile(pattern::Regex;dir="",adddir::Bool=false)
    fs = filter!(f->occursin(pattern,f),readdir(dir))
    if adddir
        fs=joinpath.(dir,fs)
    end
    return fs
end

function querymeta(meta::DataFrame;test="",sourceformat="",subject="")
    @from i in meta begin
    @where startswith(get(i.Subject_ID),subject)
    @where i.ID==test
    @where i.sourceformat==sourceformat
    @select {i.ID,i.UUID,i.files}
    @collect DataFrame
    end
end
