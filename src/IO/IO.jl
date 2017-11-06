export readmat,mat2julia!,loadimageset,CondDCh,MarkDCh,StartDCh,StopDCh,digitaldata,
prepare,prepare!,prepare_ripple!,prepare_oi!,prepare_vlab!,
statetime,trim,condfactor,condtest,maptodatatime,
oifileregex,getoifile,vlabfileregex,getvlabfile,matchfile

using MAT,DataFrames,FileIO,Colors

"Read exported `Matlab` MAT format data"
function readmat(f::AbstractString,v="dataset")
  d = matread(f)[v]
end

mat2julia!(a;isscaler=false)=mat2julia!(nothing,a,isscaler=isscaler)
"Convert Matlab variable to Julia type with proper dimention"
function mat2julia!(t,a;isscaler=false)
  da = ndims(a)
  if da==0
  elseif da==1
    t!=nothing && (t=Array{t,1})
  else
    sa = size(a)
    if sa[1]==1 && sa[2]==1
      if isscaler
      a=a[1,1]
      else
      a=squeeze(a,1)
      t!=nothing &&  ( t=Array{t,1})
      end
    elseif sa[1]==1 && sa[2]>1
      a=squeeze(a,1)
      t!=nothing &&  ( t=Array{t,1})
    elseif sa[2]==1 && sa[1]>1
      a=squeeze(a,2)
      t!=nothing && ( t=Array{t,1})
    else
      t!=nothing &&( t = Array{t,da})
    end
  end
  t!=nothing &&( a=convert(t,a))
  return a
end

"Load images in path"
function loadimageset(path;name=[],n=Inf,alpha=false)
    isdir(path) || error("Invalid Directory Path")
    n<=0 && error("n <= 0")
    for (root, dirs, files) in walkdir(path)
        if isempty(name)
            n = Int(min(n,length(files)))
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
"Prepare exported Matlab dataset file"
prepare(f::AbstractString,v="dataset")=prepare!(readmat(f,v))
function prepare!(d::Dict)
  if haskey(d,"ex")
    d["ex"]=prepare_vlab!(d["ex"])
  end
  if haskey(d,"sourceformat")
    sf = d["sourceformat"]
    if(sf=="Ripple")
      d=prepare_ripple!(d)
    elseif (sf=="OI")
      d=prepare_oi!(d)
    end
  end
  return d
end
"Prepare Ripple data"
function prepare_ripple!(d::Dict)
  if haskey(d,"spike")
    d["spike"]["electrodeid"]=mat2julia!(Int,d["spike"]["electrodeid"])
    d["spike"]["time"]=map(i->mat2julia!(Float64,i),mat2julia!(d["spike"]["time"]))
    d["spike"]["unitid"]=map(i->mat2julia!(Int,i),mat2julia!(d["spike"]["unitid"]))
    if ndims(d["spike"]["electrodeid"])==0
      d["spike"]["electrodeid"]=[d["spike"]["electrodeid"]]
      d["spike"]["time"]=[d["spike"]["time"]]
      d["spike"]["unitid"]=[d["spike"]["unitid"]]
    end
    d["spike"]["uuid"] = map(x->sort(unique(x)),d["spike"]["unitid"])
  end
  if haskey(d,"digital")
    dc = i-> begin
    s = split(i)
    c = s[1]=="SMA"?parse(Int,s[2]):i
  end
  d["digital"]["channel"]=map(i->dc(mat2julia!(i,isscaler=true)),mat2julia!(d["digital"]["channel"]))
  d["digital"]["time"]=map(i->mat2julia!(Float64,i),mat2julia!(d["digital"]["time"]))
  d["digital"]["data"]=map(i->mat2julia!(Int,i),mat2julia!(d["digital"]["data"]))
  if haskey(d,"ex")
    startdchidx = find(d["digital"]["channel"].==StartDCh)
    if !isempty(startdchidx)
    d["ex"]["t0"]=d["digital"]["time"][startdchidx[1]]
    else
    d["ex"]["t0"]=0
    end
  end
end
if haskey(d,"analog1k")
  d["analog1k"]["electrodeid"]=mat2julia!(Int,d["analog1k"]["electrodeid"])
  d["analog1k"]["time"]=mat2julia!(Float64,d["analog1k"]["time"])
  d["analog1k"]["data"]=mat2julia!(Float64,d["analog1k"]["data"])
end
return d
end
"Prepare Optical Imaging Block Data"
function prepare_oi!(d::Dict)
  if haskey(d,"imagehead")
    d["imagehead"]["listofstimuli"]=mat2julia!(Int,d["imagehead"]["listofstimuli"])
    d["imagehead"]["nframesperstim"]=mat2julia!(Int,d["imagehead"]["nframesperstim"])
  end
  return d
end
"Prepare VLab Experiment Dict"
function prepare_vlab!(d::Dict)
  if haskey(d,"CondTest")
    if haskey(d["CondTest"],"CONDSTATE")
      d["CondTest"]["CONDSTATE"]= map(i->mat2julia!(i),mat2julia!(d["CondTest"]["CONDSTATE"]))
    end
    if haskey(d["CondTest"],"CondIndex")
      d["CondTest"]["CondIndex"] = mat2julia!(Int,d["CondTest"]["CondIndex"])+1
      d["CondTest"]["CondRepeat"] = mat2julia!(Int,d["CondTest"]["CondRepeat"])
    end
  end
  if (haskey(d,"Cond") && length(d["Cond"])>0)
    for f in keys(d["Cond"])
      d["Cond"][f] = mat2julia!(d["Cond"][f])
    end
  end
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

"Trim vector values of Dict"
function trim(d::Dict)
  l=map(length,values(d))
  lmin = minimum(l)
  for k in keys(d)
    d[k]=d[k][1:lmin]
  end
  return d
end

"Get condition factor/level DataFrame"
function condfactor(cond::Dict)
  df = DataFrame(cond)
  return df
end
"Get condtest DataFrame, extract state time"
function condtest(ctd::Dict,cond::DataFrame)
  ctd = trim(ctd)
  ctd["preiciontime"]= Array{Float64}(statetime(ctd,statetype="CONDSTATE",state="PREICI"))
  ctd["condontime"]= Array{Float64}(statetime(ctd,statetype="CONDSTATE",state="COND"))
  ctd["suficiontime"]= Array{Float64}(statetime(ctd,statetype="CONDSTATE",state="SUFICI"))
  ct = DataFrame(ctd)

  ctcond = cond[ct[:CondIndex],:]
  [ct ctcond]
end
condtest(ctd::Dict,cond::Dict) = condtest(ctd,condfactor(cond))
"Get condtest DataFrame, optionally map time to data timeline."
function condtest(ctd::Dict,cond::DataFrame,ex::Dict;maptime=true)
  ct = condtest(ctd,cond)
  if maptime
  ct[:preiciontime]=maptodatatime(ct[:preiciontime],ex,addlatency=true)
  ct[:condontime]=maptodatatime(ct[:condontime],ex,addlatency=true)
  ct[:suficiontime]=maptodatatime(ct[:suficiontime],ex,addlatency=true)
  end
  return ct
end
condtest(ex::Dict;maptime=true) = condtest(ex["CondTest"],condfactor(ex["Cond"]),ex,maptime=maptime)

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

"Get Optical Imaging `VDAQ` matched files path in path"
function getoifile(;subject="[A-Za-z0-9]",session="[A-Za-z0-9]",experimentid="[0-9]",format="mat",dir="",path=true)
  matchfile(oifileregex(subject=subject,session=session,experimentid=experimentid,format=format),dir=dir,path=path)
end

"Regular Expression to match `VLab` data file name"
function vlabfileregex(;subject="[A-Za-z0-9]",session="[A-Za-z0-9]",site="[A-Za-z0-9]",test="[A-Za-z0-9]",maxrepeat=5,format="mat")
  mr = lpad(maxrepeat,2,0);d2=mr[1];d1=parse(Int,d2)>0?9:mr[2]
  Regex("^$subject+_$session*_?$site*_?$test+_[0-$d2]?[0-$d1][.]$format")
end

"Get matched `VLab` files in directory"
function getvlabfile(;subject="[A-Za-z0-9]",session="[A-Za-z0-9]",site="[A-Za-z0-9]",test="[A-Za-z0-9]",maxrepeat=5,format="mat",dir="",path=false)
  matchfile(vlabfileregex(subject=subject,session=session,site=site,test=test,maxrepeat=maxrepeat,format=format),dir=dir,path=path)
end

"Get matched files in directory, optionally return file path"
function matchfile(pattern::Regex;dir="",path::Bool=false)
  fs = filter!(f->ismatch(pattern,f),readdir(dir))
  if path
      fs=joinpath.(dir,fs)
  end
  return fs
end
