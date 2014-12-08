export subtps,subtpsr

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

function subtpsr(tps::TimePoints,starttimes::TimePoints,endtimes::TimePoints;winstart::Real=0.0,winend::Real=0.0
                ,isstarttimeorigin::Bool=false,isendtimeorigin::Bool=false)
  if isendtimeorigin && isstarttimeorigin
    error("Response Window Time Origin Conflicts.")
  end
  if isstarttimeorigin
    startt = starttimes + winstart
    endt = starttimes + winend
  end
  if isendtimeorigin
    startt = endtimes + winstart
    endt = endtimes + winend
  end
  if !isendtimeorigin && !isstarttimeorigin
    startt = starttimes
    endt = endtimes
  end
  sts,stn,sis,wins = subtps(tps,startt,endt,isstarttimeorigin=isstarttimeorigin,isendtimeorigin=isendtimeorigin)
  response = stn ./ ((endt-startt)*0.001)
end
