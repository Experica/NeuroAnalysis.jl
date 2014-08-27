export histtps

function prepare(block)

end

function histtps(tps::TimePoints,binedges::TimePoints)
  subs,subtpn,sis,wins = subtps(tps,binedges[1:end-1],binedges[2:end])
  return subtpn,wins,subs,sis
end
# binwidth = 0.001 second to get binary spike train
histtps(tps::TimePoints,starttime::Real,endtime::Real,binwidth::Real=0.001) = histtps(tps,starttime:binwidth:endtime)
histtps(tps::TimePoints,binwidth::Real=0.001) = histtps(tps,minimum(tps),maximum(tps),binwidth)


function binst(st::TimePoints,bin::TimePoints)
  bin,counts = hist(st,bin)
  return counts
end
