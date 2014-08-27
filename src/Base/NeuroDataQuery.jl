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


function subtps(tps::TimePoints,starttime::Real,endtime::Real)
  si = find((starttime .<= tps) & (tps .< endtime))
  tpn = length(si)
  sub = tps[si]
  win = [starttime, endtime]
  return sub,tpn,si,win
end
function subtps(tps::TimePoints,starttimes::TimePoints,endtimes::TimePoints)
  subn = length(starttimes)
  if subn != length(endtimes)
    error("Lengths of starttimes and endtimes do not match.")
  end
  subs = cell(subn)
  subtpn = zeros(subn)
  sis = cell(subn)
  wins = cell(subn)
  for i in 1:subn
    subs[i],subtpn[i],sis[i],wins[i] = subtps(tps,starttimes[i],endtimes[i])
  end
  return subs,subtpn,sis,wins
end
