function prepare(block)

end

function binst(st::TimePoints,bin::TimePoints)
  bin,counts = hist(st,bin)
  return counts
end
# binwidth = 0.001 second to get binary spike train
binst(st::TimePoints,begintime::Real,endtime::Real,bw::Real=0.001) = binst(st,begintime:bw:endtime)
