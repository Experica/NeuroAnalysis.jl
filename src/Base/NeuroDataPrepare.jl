export histtps,histarray

function prepare(block)

end

function histtps(tps::TimePoints,binedges::TimePoints;isstarttimeorigin::Bool=false)
  subs,subtpn,sis,wins = subtps(tps,binedges[1:end-1],binedges[2:end],isstarttimeorigin=isstarttimeorigin)
  return subtpn,wins,subs,sis
end
function histtps(tps::TimePoints,starttime::Real,endtime::Real;nbins::Integer=10,binwidth::Real=0.0,isstarttimeorigin::Bool=false)
  if binwidth <= 0.0
    binwidth = (endtime-starttime)/nbins
  end
  histtps(tps,starttime:binwidth:endtime,isstarttimeorigin=isstarttimeorigin)
end
histtps(tps::TimePoints;nbins::Integer=10,binwidth::Real=0.0,isstarttimeorigin::Bool=false) = histtps(tps,minimum(tps),maximum(tps),nbins=nbins,binwidth=binwidth,isstarttimeorigin=isstarttimeorigin)


function histtps(tps::TPsVector,binedges::TimePoints;isstarttimeorigin::Bool=false)
    starttimes = binedges[1:end-1]
    endtimes = binedges[2:end]
    ntps = length(tps)
    vcounts = Array(AbstractVector{Integer},ntps);
    vsubs = Array(TPsVector,ntps)
    vsis = Array(AbstractVector{AbstractVector{Integer}},ntps)
    wins = []
    for i=1:ntps
        vsubs[i],vcounts[i],vsis[i],wins = subtps(tps[i],starttimes,endtimes,isstarttimeorigin=isstarttimeorigin)
    end
    return vcounts,wins,vsubs,vsis
end
function histtps(tps::TPsVector,starttime::Real,endtime::Real;nbins::Integer=10,binwidth::Real=0.0,isstarttimeorigin::Bool=false)
  if binwidth <= 0.0
    binwidth = (endtime-starttime)/nbins
  end
  histtps(tps,starttime:binwidth:endtime,isstarttimeorigin=isstarttimeorigin)
end

function histarray(vhist::AbstractVector{AbstractVector{Integer}},bins::TPsVector)
  hn = length(vhist)
  nbins = length(bins)
  binwidth = bins[1][2]-bins[1][1]
  harray = Array(Integer,hn,nbins)
  for i=1:hn
    harray[i,:] = vhist[i]
  end
  bincenters = [bins[i][1]+binwidth/2.0 for i=1:nbins]
  return harray,bincenters
end


function binst(st::TimePoints,bin::TimePoints)
  bin,counts = hist(st,bin)
  return counts
end
