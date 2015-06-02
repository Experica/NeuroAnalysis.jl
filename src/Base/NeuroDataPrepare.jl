export isi,histrv,histmatrix,psth
using DataFrames


function isi(rv::RealVector)
  diff(sort(rv))
end

histrv(rv::RealVector,binedges::RealVector;isminzero::Bool=false,ismaxzero::Bool=false) = subrv(rv,binedges[1:end-1],binedges[2:end],isminzero=isminzero,ismaxzero=ismaxzero)
function histrv(rv::RealVector,min::Real,max::Real;nbins::Integer=10,binwidth::Real=0.0,isminzero::Bool=false,ismaxzero::Bool=false)
  if binwidth <= 0.0
    binwidth = (max-min)/nbins
  end
  histrv(rv,min:binwidth:max,isminzero=isminzero,ismaxzero=ismaxzero)
end
histrv(rv::RealVector;nbins::Integer=10,binwidth::Real=0.0,isminzero::Bool=false,ismaxzero::Bool=false) = histrv(rv,minimum(rv),maximum(rv),nbins=nbins,binwidth=binwidth,isminzero=isminzero,ismaxzero=ismaxzero)
function histrv(rvs::RVVector,binedges::RealVector;isminzero::Bool=false,ismaxzero::Bool=false)
  yn = length(rvs)
  yn == 0 && warn("Empty RVVector in histrv(rvs::RVVector, binedges::RealVector)")
  ys = Array(RVVector,yn)
  ns = Array(Vector{Int},yn)
  ws = []
  is = Array(Vector{Vector{Int}},yn)
  for i in 1:yn
    ys[i],ns[i],ws,is[i] = histrv(rvs[i],binedges,isminzero=isminzero,ismaxzero=ismaxzero)
  end
  return ys,ns,ws,is
end
function histrv(rvs::RVVector,min::Real,max::Real;nbins::Integer=10,binwidth::Real=0.0,isminzero::Bool=false,ismaxzero::Bool=false)
  if binwidth <= 0.0
    binwidth = (max-min)/nbins
  end
  histrv(rvs,min:binwidth:max,isminzero=isminzero,ismaxzero=ismaxzero)
end

function histmatrix(hv::Vector{Vector{Int}},ws::RVVector)
  hn = length(hv)
  nbins = length(ws)
  ((hn == 0) || (nbins == 0)) && error("Arguments Empty.")
  binwidth = ws[1][2]-ws[1][1]
  hm = Array(Int,hn,nbins)
  for i in 1:hn
    hm[i,:] = hv[i]
  end
  bincenters = [ws[i][1]+binwidth/2.0 for i=1:nbins]
  return hm,bincenters
end
function histmatrix(rvs::RVVector,binedges::RealVector)
  ys,ns,ws,is = histrv(rvs,binedges)
  hm,x = histmatrix(ns,ws)
end

function psth(hm::Matrix{Int},x)
  binwidth = x[2]-x[1]
  hmr = hm / (binwidth*SecondPerUnit)
  n = size(hmr,1)
  m = mean(hmr,1)[:]
  sd = std(hmr,1)[:]
  return m,sd,n,x
end
function psth(hv::Vector{Vector{Int}},ws::RVVector)
  hm,x = histmatrix(hv,ws)
  psth(hm,x)
end
function psth(rvs::RVVector,binedges::RealVector)
  hm,x = histmatrix(rvs,binedges)
  psth(hm,x)
end
function psth(rvs::RVVector,binedges::RealVector,condstr)
  m,sd,n,x = psth(rvs,binedges)
  df = DataFrame(x=x,y=m,ysd=sd,n=n,condition=condstr)
end
function psth(rvvs::RVVVector,binedges::RealVector,condstrs::Vector{String})
  n = length(condstrs)
  dfs = [psth(rvvs[i],binedges,condstrs[i]) for i=1:n]
  df=DataFrame()
  for i=1:n
    df = [df,dfs[i]]
  end
  return df
end
