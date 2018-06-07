import Base: convert
using Distributions,DataFrames,StatsBase,GLM,LsqFit,HypothesisTests

include("NeuroDataType.jl")
include("Spike.jl")
include("Image.jl")

export anscombe,isresponsive,vmf,gvmf,circvar,circr,statsori,flcond,subcond,findcond,flni,condni,condfactor,finalfactor,condstring,condresponse,
setfln,testfln,condmean

anscombe(x) = 2*sqrt(x+(3/8))

isresponsive(baseline,response;alpha=0.05) = pvalue(SignedRankTest(baseline,response)) < alpha
isresponsive(baseline,response,is;alpha=0.05) = any(map(i->isresponsive(baseline[i],response[i],alpha=alpha),is))

vmf(α,β=1,μ=0.0,κ=1.0,n=1) = β*exp(κ*(cos(n*(α-μ))-1))
gvmf(α,β=1,μ₁=0.0,κ₁=1.0,μ₂=0.0,κ₂=1.0) = β*exp(κ₁*(cos(α-μ₁)-1) + κ₂*(cos(2*(α-μ₂))-1))

circvar(α::AbstractVector,w=ones(length(α))) = 1-circr(α,w)
function circr(α::AbstractVector,w=ones(length(α)))
    abs(sum(w.*exp.(im*α)))/sum(w)
end

function statsori(m,ori)
    a = deg2rad.(ori)
    dcv = circvar(a,m)
    dtf = curve_fit((x,p)->gvmf.(x,p...),a,m,Float64[1,0,0,0,0])

    # 1deg = 0.017rad, 0.01rad = 0.57deg
    x = collect(0:0.01:2pi)
    pdir = rad2deg(x[indmax(gvmf.(x,dtf.param...))])
    
    o = a.%pi
    ocv = circvar(2o,m)
    otf = curve_fit((x,p)->vmf.(x,p...,2),o,m,Float64[1,0,0])
    x = collect(0:0.01:pi)
    pori = rad2deg(x[indmax(vmf.(x,otf.param...,2))])

    Dict(:dcv=>dcv,:pdir=>pdir,:ocv=>ocv,:pori=>pori)
end

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

# function convert(::Type{DataFrame},ct::CoefTable)
#    df = convert(DataFrame,ct.mat)
#    names!(df,map(symbol,ct.colnms))
#    [DataFrame(coefname=ct.rownms) df]
#end



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
function subcond(conds,sc...)
    sci = map(i->issubset(sc,i),conds)
    return conds[sci]
end

function findcond(df::DataFrame,cond::Vector{Any};roundingdigit=3)
    i = trues(size(df,1))
    condstr = ""
    for fl in cond
        f = fl[1]
        l = fl[2]
        i &= df[Symbol(f)].==l
        if typeof(l)<:Real
            lv=round(l,roundingdigit)
        else
            lv=l
        end
        condstr = "$condstr, $f=$lv"
    end
    return find(i),condstr[3:end]
end
function findcond(df::DataFrame,conds::Vector{Vector{Any}};roundingdigit=3)
    n = length(conds)
    is = Array(Vector{Int},n)
    ss = Array(String,n)
    for i in 1:n
        is[i],ss[i] = findcond(df,conds[i],roundingdigit=roundingdigit)
    end
    vi = map(i->!isempty(i),is)
    return is[vi],ss[vi],conds[vi]
end

flni(cond::Dict)=flni(DataFrame(cond))
"""
Find levels(except missing) for each factor and repetition, indices for each level
"""
function flni(ctc::DataFrame)
    fl=Dict();fln=Dict();fli=Dict()
    for f in names(ctc)
        fv = skipmissing(ctc[f])
        ul = sort(unique(fv))
        uls = [fv.==l for l in ul]
        fl[f]=ul;fln[f]=countnz.(uls);fli[f]=find.(uls)
    end
    return fl,fln,fli
end

condni(cond::Dict)=condni(DataFrame(cond))
"""
Find unique condition and repetition, indices for each
"""
function condni(ctc::DataFrame)
    t = [ctc DataFrame(i=1:size(ctc,1))]
    sort(by(t, names(ctc),g->DataFrame(n=size(g,1), i=[g[:i]])))
end

condfactor(cond::DataFrame)=setdiff(names(cond),[:n,:i])

function finalfactor(cond::DataFrame)
    fs = String.(condfactor(cond))
    for i in length(fs):-1:1
        if !endswith(fs[i],"_Final") && any(fs.=="$(fs[i])_Final")
            deleteat!(fs,i)
        end
    end
    Symbol.(fs)
end

function condstring(cond::DataFrameRow,fs)
    join(["$f=$(cond[f])" for f in fs],", ")
end

function condresponse(rs,is)
    grs = [rs[i] for i in is]
    DataFrame(m=mean.(grs),se=sem.(grs))
end
function condresponse(rs,cond::DataFrame,u=0)
    crs = [rs[r[:i]] for r in eachrow(cond)]
    df = [DataFrame(m=mean.(crs),se=sem.(crs),u=fill(u,length(crs))) cond[condfactor(cond)]]
end
function condresponse(urs::Dict,cond::DataFrame)
    vcat([condresponse(v,cond,k) for (k,v) in urs]...)
end

function psth(rvs::RVVector,binedges::RealVector,c;normfun=nothing)
    m,se,x = psth(rvs,binedges,normfun=normfun)
    df = DataFrame(x=x,m=m,se=se,c=[c for _ in 1:length(x)])
end
function psth(rvs::RVVector,binedges::RealVector,cond::DataFrame;normfun=nothing)
    fs = finalfactor(cond)
    vcat([psth(rvs[r[:i]],binedges,condstring(r,fs),normfun=normfun) for r in eachrow(cond)])
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

function psth(rvvs::RVVVector,binedges::RealVector,conds;normfun=nothing)
    n = length(rvvs)
    n!=length(conds) && error("Length of rvvs and conds don't match.")
    dfs = [psth(rvvs[i],binedges,conds[i],normfun=normfun) for i=1:n]
    return cat(1,dfs)
end
