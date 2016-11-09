export subrv,subrvr,isi,flatrvs,histrv,histmatrix,psth

function subrv(rv::RealVector,min::Real,max::Real;isminzero::Bool=false,ismaxzero::Bool=false)
  if ismaxzero && isminzero
    error("Zero Setting Conflicts.")
  end
  i = find(min .<= rv .< max)
  n = length(i)
  y = rv[i]
  if isminzero
    y -= min
  end
  if ismaxzero
    y -= max
  end
  w = [min, max]
  return y,n,w,i
end
function subrv(rv::RealVector,mins::RealVector,maxs::RealVector;isminzero::Bool=false,ismaxzero::Bool=false,israte::Bool=false)
  yn = length(mins)
  if yn != length(maxs)
    error("Length of mins and maxs do not match.")
  end
  ys = Array(RealVector,yn)
  ns = zeros(Int,yn)
  ws = Array(RealVector,yn)
  is = Array(Vector{Int},yn)
  for i in 1:yn
    ys[i],ns[i],ws[i],is[i] = subrv(rv,mins[i],maxs[i],isminzero=isminzero,ismaxzero=ismaxzero)
  end
  if israte
    ns = ns./((maxs-mins)*SecondPerUnit)
  end
  return ys,ns,ws,is
end
function subrvr(rv::RealVector,mins::RealVector,maxs::RealVector)
  ys,ns,ws,is = subrv(rv,mins,maxs,israte=true)
  return ns
end

function isi(rv::RealVector)
  diff(sort(rv))
end

function flatrvs(rvs::RVVector,sv=[])
  nrv = length(rvs)
  if isempty(sv)
    issort=false
  elseif length(sv)==nrv
    issort=true
  else
    warn("Length of rvs and sv do not match, sorting ignored.")
    issort=false
  end
  if issort
    srvs=rvs[sortperm(sv)]
    ssv=sort(sv)
  else
    srvs=rvs
    ssv=sv
  end
  x=[];y=[];s=[]
  for i in 1:nrv
    rv = srvs[i];n=length(rv)
    if n==0;continue;end
    x = [x;rv];y = [y;ones(n)*i]
    if issort;s=[s;fill(ssv[i],n)];end
  end
  return x,y,s
end

function histrv(rv::RealVector,binedges::RealVector;isminzero::Bool=false,ismaxzero::Bool=false)
  nbinedges = length(binedges)
  nbinedges<2 && error("Have $nbinedges binedges, Need at least two binedges.")
  subrv(rv,binedges[1:end-1],binedges[2:end],isminzero=isminzero,ismaxzero=ismaxzero)
end
function histrv(rv::RealVector,min::Real,max::Real;nbins::Integer=10,binwidth::Real=0.0,isminzero::Bool=false,ismaxzero::Bool=false)
  if binwidth <= 0.0
    binwidth = (max-min)/nbins
  end
  histrv(rv,min:binwidth:max,isminzero=isminzero,ismaxzero=ismaxzero)
end
histrv(rv::RealVector;nbins::Integer=10,binwidth::Real=0.0,isminzero::Bool=false,ismaxzero::Bool=false) = histrv(rv,minimum(rv),maximum(rv),nbins=nbins,binwidth=binwidth,isminzero=isminzero,ismaxzero=ismaxzero)
function histrv(rvs::RVVector,binedges::RealVector;isminzero::Bool=false,ismaxzero::Bool=false)
  yn = length(rvs)
  yn == 0 && error("Empty RVVector in histrv(rvs::RVVector, binedges::RealVector)")
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
  hn!=0 && nbins!=0 && nbins!=length(hv[1]) && error("nbins does not match.")
  binwidth = ws[1][2]-ws[1][1]
  hm = Array(Int,hn,nbins)
  for i in 1:hn
    hm[i,:] = hv[i]
  end
  bincenters = Float64[ws[i][1]+binwidth/2.0 for i=1:nbins]
  return hm,bincenters
end
function histmatrix(rvs::RVVector,binedges::RealVector)
  ys,ns,ws,is = histrv(rvs,binedges)
  hm,x = histmatrix(ns,ws)
end

function psth(hm::Matrix{Int},x::RealVector;normfun=nothing)
  binwidth = x[2]-x[1]
  hmr = hm / (binwidth*SecondPerUnit)
  n = size(hmr,1)
  if normfun==nothing
    m = mean(hmr,1)[:]
    sd = std(hmr,1)[:]
  else
    for i=1:n
      hmr[i,:]=normfun(hmr[i,:])
    end
  end
  return m,sd,n,x
end
function psth(hv::Vector{Vector{Int}},ws::RVVector;normfun=nothing)
  hm,x = histmatrix(hv,ws)
  psth(hm,x,normfun=normfun)
end
function psth(rvs::RVVector,binedges::RealVector;normfun=nothing)
  hm,x = histmatrix(rvs,binedges)
  psth(hm,x,normfun=normfun)
end
function psth(rvs::RVVector,binedges::RealVector,cond;normfun=nothing)
  m,sd,n,x = psth(rvs,binedges,normfun=normfun)
  df = DataFrame(x=x,y=m,ysd=sd,n=n,condition=cond)
end
function psth(rvs::RVVector,binedges::RealVector,rvsidx,condstr;normfun=nothing)
    nc = length(condstr)
    dfs=[psth(rvs[rvsidx[i]],binedges,condstr[i],normfun=normfun) for i=1:nc]
    return vcat(dfs)
end
function psth(rvvs::RVVVector,binedges::RealVector,conds;normfun=nothing)
  n = length(rvvs)
  n!=length(conds) && error("Length of rvvs and conds don't match.")
  dfs = [psth(rvvs[i],binedges,conds[i],normfun=normfun) for i=1:n]
  return cat(1,dfs)
end
