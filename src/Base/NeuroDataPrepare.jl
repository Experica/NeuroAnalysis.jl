export itpi,histtps,histmatrix,psth


function itpi(tps::TimePoints)
  diff(sort(tps))
end

function histtps(tps::TimePoints,binedges::TimePoints;isstarttimeorigin::Bool=false,isendtimeorigin::Bool=false)
  subs,tpns,sis,wins = subtps(tps,binedges[1:end-1],binedges[2:end],isstarttimeorigin=isstarttimeorigin,isendtimeorigin=isendtimeorigin)
  return tpns,wins,subs,sis
end
function histtps(tps::TimePoints,starttime::Real,endtime::Real;nbins::Integer=10,binwidth::Real=0.0,isstarttimeorigin::Bool=false,isendtimeorigin::Bool=false)
  if binwidth <= 0.0
    binwidth = (endtime-starttime)/nbins
  end
  histtps(tps,starttime:binwidth:endtime,isstarttimeorigin=isstarttimeorigin,isendtimeorigin=isendtimeorigin)
end
histtps(tps::TimePoints;nbins::Integer=10,binwidth::Real=0.0,isstarttimeorigin::Bool=false,isendtimeorigin::Bool=false) = histtps(tps,minimum(tps),maximum(tps),nbins=nbins,binwidth=binwidth,isstarttimeorigin=isstarttimeorigin,isendtimeorigin=isendtimeorigin)
function histtps(tps::TPsVector,binedges::TimePoints;isstarttimeorigin::Bool=false,isendtimeorigin::Bool=false)
  ntps = length(tps)
  vtpns = Array(Vector{Int},ntps)
  wins = []
  vsubs = Array(TPsVector,ntps)
  vsis = Array(Vector{Vector{Int}},ntps)
  for i in 1:ntps
    vtpns[i],wins,vsubs[i],vsis[i] = histtps(tps[i],binedges,isstarttimeorigin=isstarttimeorigin,isendtimeorigin=isendtimeorigin)
  end
  ntps == 0 && warn("Empty TPsVector in histtps(tps::TPsVector,binedges::TimePoints)")
  return vtpns,wins,vsubs,vsis
end
function histtps(tps::TPsVector,starttime::Real,endtime::Real;nbins::Integer=10,binwidth::Real=0.0,isstarttimeorigin::Bool=false,isendtimeorigin::Bool=false)
  if binwidth <= 0.0
    binwidth = (endtime-starttime)/nbins
  end
  histtps(tps,starttime:binwidth:endtime,isstarttimeorigin=isstarttimeorigin,isendtimeorigin=isendtimeorigin)
end

function histmatrix(vhist::Vector{Vector{Int}},wins::TPsVector)
  hn = length(vhist)
  nbins = length(wins)
  ((hn == 0) || (nbins == 0)) && error("Arguments Empty.")
  binwidth = wins[1][2]-wins[1][1]
  hm = Array(Int,hn,nbins)
  for i in 1:hn
    hm[i,:] = vhist[i]
  end
  bincenters = [wins[i][1]+binwidth/2.0 for i=1:nbins]
  return hm,bincenters
end
function histmatrix(tps::TPsVector,binedges::TimePoints)
  vhist,wins,vsubs,vsis = histtps(tps,binedges)
  hm,x = histmatrix(vhist,wins)
end

function psth(hm::Matrix{Int},x)
  binwidth = x[2]-x[1]
  hrm = hm / (binwidth*0.001)
  n = size(hrm,1)
  m = mean(hrm,1.0)[:]
  sd = std(hrm,1)[:]
  return m,sd,x,n
end
function psth(vhist::Vector{Vector{Int}},wins::TPsVector)
  hm,x = histmatrix(vhist,wins)
  psth(hm,x)
end
function psth(tps::TPsVector,binedges::TimePoints)
  vhist,wins,vsubs,vsis = histtps(tps,binedges)
  hm,x = histmatrix(vhist,wins)
  psth(hm,x)
end
