export readmat,mat2julia!,prepare,prepare!,prepare_ripple!,prepare_vlab!,statetime,
condfactor,condtest,aligntime,matchfile

using MAT,DataFrames

"Read exported `Matlab` MAT format data"
function readmat(f::AbstractString,v="dataset")
  d = matread(f)[v]
end

mat2julia!(a;isscaler=false)=mat2julia!(nothing,a,isscaler=isscaler)
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

prepare(f::AbstractString,v="dataset")=prepare!(readmat(f,v))
function prepare!(d::Dict)
  if haskey(d,"sourceformat")
    sf = d["sourceformat"]
    if(sf=="Ripple")
      d=prepare_ripple!(d)
    end
  end
  if haskey(d,"ex")
    d["ex"]=prepare_vlab!(d["ex"])
  end
  return d
end
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
  end
  if haskey(d,"digital")
    dc = i-> begin
    s = split(i)
    c = s[1]=="SMA"?parse(Int,s[2]):i
  end
  d["digital"]["channel"]=map(i->dc(mat2julia!(i,isscaler=true)),mat2julia!(d["digital"]["channel"]))
  d["digital"]["time"]=map(i->mat2julia!(Float64,i),mat2julia!(d["digital"]["time"]))
  d["digital"]["data"]=map(i->mat2julia!(Int,i),mat2julia!(d["digital"]["data"]))
  d["ex"]["t0"]=d["digital"]["time"][2]
end
if haskey(d,"analog1k")
  d["analog1k"]["electrodeid"]=mat2julia!(Int,d["analog1k"]["electrodeid"])
  d["analog1k"]["time"]=mat2julia!(Float64,d["analog1k"]["time"])
  d["analog1k"]["data"]=mat2julia!(Float64,d["analog1k"]["data"])
end
return d
end
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

function condfactor(cond::Dict)
  df = DataFrame(cond)
  return df
end
function condtest(ctd::Dict,cond::DataFrame)
  ctd["preiciontime"]= Array{Float64}(statetime(ctd,statetype="CONDSTATE",state="PREICI"))
  ctd["condontime"]= Array{Float64}(statetime(ctd,statetype="CONDSTATE",state="COND"))
  ctd["suficiontime"]= Array{Float64}(statetime(ctd,statetype="CONDSTATE",state="SUFICI"))
  ct = DataFrame(ctd)

  ctcond = cond[ct[:CondIndex],:]
  [ct ctcond]
end
condtest(ctd::Dict,cond::Dict) = condtest(ctd,condfactor(cond))
function condtest(ctd::Dict,cond::DataFrame,ex::Dict;addlatency=true)
  ct = condtest(ctd,cond)
  ct[:preiciontime]=aligntime(ct[:preiciontime],ex,addlatency=addlatency)
  ct[:condontime]=aligntime(ct[:condontime],ex,addlatency=addlatency)
  ct[:suficiontime]=aligntime(ct[:suficiontime],ex,addlatency=addlatency)
  return ct
end
condtest(ex::Dict;addlatency=true) = condtest(ex["CondTest"],condfactor(ex["Cond"]),ex,addlatency=addlatency)

function aligntime(x,ex::Dict;addlatency=true)
  t=(x+ex["t0"])*(1+ex["TimerDriftSpeed"])
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

function matchfile(pattern::Regex;path="")
  fs = readdir(path)
  mi = map(f->ismatch(pattern,f),fs)
  mfs = fs[mi]
end
