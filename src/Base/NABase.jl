module NABase

using NeuroAnalysis.NACore

include("NeuroDataQuery.jl")
include("NeuroDataPrepare.jl")

import Base: convert
export drv,randvec3,anscombe,fconds,fcondi,fcondis,psthdf
using DataFrames, Distributions

function drv(p,n=1,isfreq=false)
  d=Categorical(p)
  if n==1
    return rand(d)
  end
  if isfreq
    pn = round(p[1:end-1]*n)
    pn = [pn,n-sum(pn)]

    pnc = zeros(Int,length(p))
    ps = zeros(Int,n)
    for i in 1:n
      while true
        v = rand(d)
        if pnc[v] < pn[v]
          pnc[v] += 1
          ps[i]=v
          break
        end
      end
    end
    return ps
  else
    rand(d,n)
  end
end

function randvec3(xr::Vector,yr::Vector,zr::Vector)
  dx = Uniform(xr[1],xr[2])
  dy = Uniform(yr[1],yr[2])
  dz = Uniform(zr[1],zr[2])
  Vec3(rand(dx),rand(dy),rand(dz))
end
randvec3(;xmin=0,xmax=1,ymin=0,ymax=1,zmin=0,zmax=1)=randvec3([xmin,xmax],[ymin,ymax],[zmin,zmax])
function randvec3(absmin::Vector)
  d = Uniform(-1,1)
  av3 = zeros(3)
  for i=1:3
    while true
      r= rand(d)
      if abs(r) >= absmin[i]
        av3[i]=r
        break;
      end
    end
  end
  return convert(Vec3,av3)
end
randvec3(xabsmin,yabsmin,zabsmin)=randvec3([xabsmin,yabsmin,zabsmin])
function randvec3(radmin;is2d=false)
  innerabsmin = radmin*sin(pi/4)
  absmin = fill(innerabsmin,3)
  while true
    cv3 = randvec3(absmin)
    if is2d;cv3.z=0.0;end
    if length(cv3) >= radmin
      return cv3
    end
  end
end

function convert(::Type{DataFrame},ct::CoefTable)
  df = convert(DataFrame,ct.mat)
  names!(df,map(symbol,ct.colnms))
  [DataFrame(coefname=ct.rownms) df]
end

function anscombe(x)
  2*sqrt(x+(3/8))
end

function fconds(fl,f)
  conds = [[(f,l)] for l in fl[f]]
end
function fconds(fl,f1,f2)
  conds = [[(f1,l1),(f2,l2)] for l1 in fl[f1],l2 in fl[f2]][:]
end
function fconds(fl,f1,f2,f3)
  conds = [[(f1,l1),(f2,l2),(f3,l3)] for l1 in fl[f1],l2 in fl[f2],l3 in fl[f3]][:]
end
function fcondi(df::DataFrame,fcond)
  fcidx = trues(size(df,1))
  fcstr = ""
  for fc in fcond
    f = fc[1]
    l = fc[2]
    fcidx &= df[symbol(f)].==l
    fcstr = "$fcstr, $f=$l"
  end
  return find(fcidx),fcstr[3:end]
end
function fcondis(df::DataFrame,fconds)
  nfcs = length(fconds)
  fcis = Array(Vector{Int},nfcs)
  fcss = Array(String,nfcs)
  for i in 1:nfcs
    fcis[i],fcss[i] = fcondi(df,fconds[i])
  end
  return fcis,fcss
end
function psthdf(tps::TPsVector,binedges::TimePoints,condstr="")
  m,sd,x,n = psth(tps,binedges)
  df = DataFrame(x=x,y=m,ysd=sd,n=n,Condition=condstr)
end
function psthdf(tpsv::TPVVector,binedges::TimePoints,condstrs::Vector{String})
  nconds = length(condstrs)
  dfs = [psthdf(tpsv[i],binedges,condstrs[i]) for i=1:nconds]
  df=DataFrame()
  for i=1:nconds
    df = [df,dfs[i]]
  end
  return df
end

end # module
