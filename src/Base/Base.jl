include("NeuroDataType.jl")
include("Spike.jl")

import Base: convert
export sem,anscombe,flcond,findcond,flfln,setfln,testfln,condmean
using Distributions,DataFrames

sem(v) = std(v)/sqrt(length(v))

# function drv(p,n=1,isfreq=false)
#   d=Categorical(p)
#   if n==1
#     return rand(d)
#   end
#   if isfreq
#     pn = round(p[1:end-1]*n)
#     pn = [pn,n-sum(pn)]
#
#     pnc = zeros(Int,length(p))
#     ps = zeros(Int,n)
#     for i in 1:n
#       while true
#         v = rand(d)
#         if pnc[v] < pn[v]
#           pnc[v] += 1
#           ps[i]=v
#           break
#         end
#       end
#     end
#     return ps
#   else
#     rand(d,n)
#   end
# end
#
# function randvec3(xr::Vector,yr::Vector,zr::Vector)
#   dx = Uniform(xr[1],xr[2])
#   dy = Uniform(yr[1],yr[2])
#   dz = Uniform(zr[1],zr[2])
#   Vec3(rand(dx),rand(dy),rand(dz))
# end
# randvec3(;xmin=0,xmax=1,ymin=0,ymax=1,zmin=0,zmax=1)=randvec3([xmin,xmax],[ymin,ymax],[zmin,zmax])
# function randvec3(absmin::Vector)
#   d = Uniform(-1,1)
#   av3 = zeros(3)
#   for i=1:3
#     while true
#       r= rand(d)
#       if abs(r) >= absmin[i]
#         av3[i]=r
#         break;
#       end
#     end
#   end
#   return convert(Vec3,av3)
# end
# randvec3(xabsmin,yabsmin,zabsmin)=randvec3([xabsmin,yabsmin,zabsmin])
# function randvec3(radmin;is2d=false)
#   innerabsmin = radmin*sin(pi/4)
#   absmin = fill(innerabsmin,3)
#   while true
#     cv3 = randvec3(absmin)
#     if is2d;cv3.z=0.0;end
#     if length(cv3) >= radmin
#       return cv3
#     end
#   end
# end

function convert(::Type{DataFrame},ct::CoefTable)
  df = convert(DataFrame,ct.mat)
  names!(df,map(symbol,ct.colnms))
  [DataFrame(coefname=ct.rownms) df]
end

function anscombe(x)
  2*sqrt(x+(3/8))
end

# function flcond(fl,fs...)
#   if length(fs)==0
#     fs=collect(keys(fl))
#   end
#   ex = :([{(fs[1],l1)} for l1=fl[fs[1]]][:])
#   for i=2:length(fs)
#     lx=symbol("l$i")
#     push!(ex.args[1].args,:($lx=fl[fs[$i]]))
#     push!(ex.args[1].args[1].args,:((fs[$i],$lx)))
#   end
#   eval(ex)
# end
function flcond(fl,f1)
    conds = [Any[(f1,l1)] for l1 in fl[f1]]
end
function flcond(fl,f1,f2)
    conds = [Any[(f1,l1),(f2,l2)] for l1 in fl[f1],l2 in fl[f2]][:]
end
function flcond(fl,f1,f2,f3)
    conds = [Any[(f1,l1),(f2,l2),(f3,l3)] for l1 in fl[f1],l2 in fl[f2],l3 in fl[f3]][:]
end

function findcond(df::DataFrame,cond::Vector{Any})
  i = trues(size(df,1))
  condstr = ""
  for fl in cond
    f = fl[1]
    l = fl[2]
    i &= df[Symbol(f)].==l
    condstr = "$condstr, $f=$l"
  end
  return find(i),condstr[3:end]
end
function findcond(df::DataFrame,conds::Vector{Vector{Any}})
  n = length(conds)
  is = Array(Vector{Int},n)
  ss = Array(String,n)
  for i in 1:n
    is[i],ss[i] = findcond(df,conds[i])
  end
  return is,ss
end

"""
find levels(except NA) for each factor and repetition for each level
"""
function flfln(df::DataFrame,factors)
    fl=Dict();fln=Dict()
  for f in factors
    vft = dropna(df[Symbol(f)])
        ls = sort(unique(vft))
        ln = Int[countnz(vft.==l) for l in ls]
        fl[f]=ls;fln[f]=ln
  end
    return fl,fln
end

function setfln(fl::Dict,n::Int)
  fln=Dict()
  for f in keys(fl)
    for l in fl[f]
      fln[(f,l)] = n
    end
  end
  return fln
end
function setfln(fl::Dict,fn::Dict)
  fln=Dict()
  for f in keys(fl)
    for i=1:length(fl[f])
      fln[(f,fl[f][i])] = fn[f][i]
    end
  end
  return fln
end

function testfln(fln::Dict,minfln::Dict;showmsg::Bool=true)
  r=true
  for mkey in keys(minfln)
    if haskey(fln,mkey)
      if fln[mkey] < minfln[mkey]
        if showmsg;print("Level $(mkey[2]) of factor \"$(mkey[1])\" repeats $(fln[mkey]) < $(minfln[mkey]).\n");end
        r=false
      end
    else
      if showmsg;print("Level $(mkey[2]) of factor \"$(mkey[1])\" is missing.\n");end
      r=false
    end
  end
  return r
end
function testfln(ds::Vector,minfln::Dict;showmsg::Bool=true)
  vi = find(map(d->begin
             if showmsg;print("Testing factor/level of \"$(d["datafile"])\" ...\n");end
             testfln(d["fln"],minfln,showmsg=showmsg)
           end,ds))
end

function condmean(rs,ridx,conds)
    nc = length(conds)
    m = Array(Float64,nc)
    sd = similar(m)
    n = Array(Int,nc)
    for i=1:nc
        r=rs[ridx[i]]
        m[i]=mean(r)
        sd[i]=std(r)
        n[i]=length(r)
    end
    return m,sd,n
end
function condmean(rs,us,ridx,conds)
    msdn = map(i->condmean(i,ridx,conds),rs)
    m= map(i->i[1],msdn)
    sd = map(i->i[2],msdn)
    n=map(i->i[3],msdn)
    return m,hcat(sd...),hcat(n...)
end

function psth(ds::DataFrame,binedges::RealVector,conds::Vector{Vector{Any}};normfun=nothing,spike=:spike,isse::Bool=true)
  is,ss = findcond(ds,conds)
  df = psth(map(x->ds[spike][x],is),binedges,ss,normfun=normfun)
  if isse
    df[:ymin] = df[:y]-df[:ysd]./sqrt(df[:n])
    df[:ymax] = df[:y]+df[:ysd]./sqrt(df[:n])
  end
  return df,ss
end
