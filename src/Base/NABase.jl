module NABase

using NeuroAnalysis.NACore

include("NeuroDataQuery.jl")
include("NeuroDataPrepare.jl")

import Base: convert
export anscombe,fconds,fcondi,fcondis,psthdf
using DataFrames

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
