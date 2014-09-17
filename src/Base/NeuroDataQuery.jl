export subtps

function blockinfo()
end

function unitinfo()
end

function eventdata()
end

function indexbytime()
end

function timebyindex()
end


function subtps(tps::TimePoints,starttime::Real,endtime::Real;isstarttimeorigin::Bool=false,isendtimeorigin::Bool=false)
  if isendtimeorigin && isstarttimeorigin
    error("Time Origin Setting Conflicts.")
  end
  si = find(starttime .<= tps .< endtime)
  tpn = length(si)
  sub = tps[si]
  if isstarttimeorigin
    sub -= starttime
  end
  if isendtimeorigin
    sub -= endtime
  end
  win = [starttime, endtime]
  return sub,tpn,si,win
end
function subtps(tps::TimePoints,starttimes::TimePoints,endtimes::TimePoints;isstarttimeorigin::Bool=false,isendtimeorigin::Bool=false)
  subn = length(starttimes)
  if subn != length(endtimes)
    error("Lengths of starttimes and endtimes do not match.")
  end
  subs = Array(TimePoints,subn)
  tpns = zeros(Int,subn)
  sis = Array(Vector{Int},subn)
  wins = Array(TimePoints,subn)
  for i in 1:subn
    subs[i],tpns[i],sis[i],wins[i] = subtps(tps,starttimes[i],endtimes[i],isstarttimeorigin=isstarttimeorigin,isendtimeorigin=isendtimeorigin)
  end
  return subs,tpns,sis,wins
end
