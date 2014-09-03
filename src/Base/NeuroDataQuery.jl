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


function subtps(tps::TimePoints,starttime::Real,endtime::Real;isstarttimeorigin::Bool=false)
  si = find((starttime .<= tps) & (tps .< endtime))
  tpn = length(si)
  sub = tps[si]
  if isstarttimeorigin
    sub = sub - starttime
  end
  win = [starttime, endtime]
  return sub,tpn,si,win
end
function subtps(tps::TimePoints,starttimes::TimePoints,endtimes::TimePoints;isstarttimeorigin::Bool=false)
  subn = length(starttimes)
  if subn != length(endtimes)
    error("Lengths of starttimes and endtimes do not match.")
  end
  subs = Array(TimePoints,subn)
  subtpn = zeros(Integer,subn)
  sis = Array(AbstractVector{Integer},subn)
  wins = Array(TimePoints,subn)
  for i in 1:subn
    subs[i],subtpn[i],sis[i],wins[i] = subtps(tps,starttimes[i],endtimes[i],isstarttimeorigin=isstarttimeorigin)
  end
  return subs,subtpn,sis,wins
end
