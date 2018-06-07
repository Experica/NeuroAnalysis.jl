export readmat,readmeta,mat2julia!,loadimageset,CondDCh,MarkDCh,StartDCh,StopDCh,digitaldata,
prepare,prepare!,prepare_ripple!,prepare_oi!,prepare_vlab!,
statetime,getparam,condtestfactor,condtest,ctctc,maptodatatime,
oifileregex,getoifile,vlabfileregex,getvlabfile,matchfile,querymeta

using MAT,DataFrames,Query,FileIO,Colors

"Read `Matlab` MAT format data"
function readmat(f::AbstractString,vars...)
    datavars=["spike","lfp","digital","analog1k","image"]
    if isempty(vars)
        d=matread(f)
    else
        d=Dict()
        matopen(f) do file
            fvs = names(file)
            basevars=setdiff(fvs,datavars)
            vs=union(basevars,intersect(datavars,vars))
            if !isempty(vs)
                for v in vs
                    d[v]=read(file,v)
                end
            end
        end
    end
    return d
end

"Read Metadata MAT file"
function readmeta(f::AbstractString)
    d = readmat(f)["Tests"]
    mat2julia!(d)
    DataFrame(d)
end

"Convert `Matlab` variable to `Julia` type with proper dimention"
function mat2julia!(x;isscaler = true,isvector = true)
    if x isa Dict
        for k in keys(x)
            x[k]=mat2julia!(x[k],isscaler=isscaler,isvector=isvector)
        end
    elseif x isa Array
        if ndims(x)==2
            s = size(x)
            if s[1]==1 && s[2]==1
                x = isscaler?x[1,1]:isvector?squeeze(x,2):x
            elseif s[1]==1 && s[2]>1
                x = isvector?squeeze(x,1):x
            elseif s[1]>1 && s[2]==1
                x = isvector?squeeze(x,2):x
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
"Get digital channel time and value."
function digitaldata(dataset::Dict,ch)
    chidx = find(dataset["digital"]["channel"].==ch)
    if !isempty(chidx)
        return dataset["digital"]["time"][chidx[1]],dataset["digital"]["data"][chidx[1]]
    else
        return [],[]
    end
end
"Prepare exported `Matlab` dataset"
prepare(f::AbstractString,vars...)=prepare!(readmat(f,vars...))
function prepare!(d::Dict)
    if haskey(d,"sourceformat")
        sf = d["sourceformat"]
        if sf=="Ripple"
            d=prepare_ripple!(d)
        elseif sf=="OI"
            d=prepare_oi!(d)
        elseif sf=="VLab"
            d=prepare_vlab!(d)
        end
    end
    return d
end
"Prepare Ripple Data"
function prepare_ripple!(d::Dict)
    mat2julia!(d)
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
"Prepare Optical Imaging Block Data"
function prepare_oi!(d::Dict)
    mat2julia!(d)
    return d
end
"Prepare VLab Experiment Data"
function prepare_vlab!(d::Dict)
    mat2julia!(d)
    return d
end

"Extract state time of condtest"
function statetime(ct::Dict;statetype::AbstractString="CONDSTATE",state::AbstractString="COND")
    if haskey(ct,statetype)
        filter!(l->!isempty(l),map(i->begin
        t = filter!(k->!isempty(k),map(j->haskey(j,state)?j[state]:[],i))
        length(t)==1?t[1]:t
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

"Get condition test factor DataFrame"
function condtestfactor(ctcd::Dict)
    ctc = DataFrame(ctcd)
    return ctc
end
"Get `condtest` DataFrame"
function condtest(ctd::Dict)
    ct = DataFrame(ctd)
    return ct
end
"Get `condtest` and `condtestcond` DataFrame"
function ctctc(ctd::Dict,ctcd::Dict)
    return condtest(ctd),condtestfactor(ctcd)
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
function getoifile(;subject="[A-Za-z0-9]",session="[A-Za-z0-9]",experimentid="[0-9]",format="mat",dir="",adddir=true)
    matchfile(oifileregex(subject=subject,session=session,experimentid=experimentid,format=format),dir=dir,adddir=adddir)
end

"Regular Expression to match `VLab` data file name"
function vlabfileregex(;subject="[A-Za-z0-9]",session="[A-Za-z0-9]",site="[A-Za-z0-9]",test="[A-Za-z0-9]",maxrepeat=10,format="mat")
    mr = lpad(maxrepeat,2,0);d2=mr[1];d1=parse(Int,d2)>0?9:mr[2]
    Regex("^$subject+_$session*_?$site*_?$test+_[0-$d2]?[0-$d1][.]$format")
end

"Get matched `VLab` files in directory"
function getvlabfile(;subject="[A-Za-z0-9]",session="[A-Za-z0-9]",site="[A-Za-z0-9]",test="[A-Za-z0-9]",maxrepeat=10,format="mat",dir="",adddir=false)
    matchfile(vlabfileregex(subject=subject,session=session,site=site,test=test,maxrepeat=maxrepeat,format=format),dir=dir,adddir=adddir)
end

"Get matched files in directory, optionally add directory to get file path"
function matchfile(pattern::Regex;dir="",adddir::Bool=false)
    fs = filter!(f->ismatch(pattern,f),readdir(dir))
    if adddir
        fs=joinpath.(dir,fs)
    end
    return fs
end

function querymeta(meta::DataFrame;test="",sourceformat="Ripple",species="C")
    @from i in meta begin
    @where startswith(get(i.Subject_ID),species)
    @where i.ID==test
    @where i.sourceformat==sourceformat
    @select {i.UUID,i.files}
    @collect DataFrame
    end
end