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


function cutt(t::TimePoints,begintime::Real,endtime::Real)
  t[(begintime .<= t) & (t .< endtime)]
end
function cutt(t::TimePoints,begintime::TimePoints,endtime::TimePoints)
  btn = length(begintime)

end
