export subrv,subrvr


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
function subrv(rv::RealVector,mins::RealVector,maxs::RealVector;isminzero::Bool=false,ismaxzero::Bool=false)
  yn = length(mins)
  if yn != length(maxs)
    error("Lengths of mins and maxs do not match.")
  end
  ys = Array(RealVector,yn)
  ns = zeros(Int,yn)
  ws = Array(RealVector,yn)
  is = Array(Vector{Int},yn)
  for i in 1:yn
    ys[i],ns[i],ws[i],is[i] = subrv(rv,mins[i],maxs[i],isminzero=isminzero,ismaxzero=ismaxzero)
  end
  return ys,ns,ws,is
end

function subrvr(rv::RealVector,mins::RealVector,maxs::RealVector;winmin::Real=0.0,winmax::Real=0.0
                ,isminzero::Bool=false,ismaxzero::Bool=false,isfiringrate::Bool=false)
  if ismaxzero && isminzero
    error("Zero Setting For Response Window Conflicts.")
  end
  if isminzero
    winmins = mins + winmin
    winmaxs = mins + winmax
  end
  if ismaxzero
    winmins = maxs + winmin
    winmaxs = maxs + winmax
  end
  if !ismaxzero && !isminzero
    winmins = mins
    winmaxs = maxs
  end
  ys,ns,ws,is = subrv(rv,winmins,winmaxs)
  isfiringrate?ns ./ ((winmaxs-winmins)*SecondPerUnit):ns
end
