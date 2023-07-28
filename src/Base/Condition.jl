
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

"Find condition index with kwargs: Factor = level"
function findcond(df::DataFrame;fl...)
    i = trues(size(df,1))
    for f in keys(fl)
        i .&= df[!,f].==fl[f]
    end
    return findall(i)
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

"Check if `response` is significently different from `baseline` by `Wilcoxon Signed Rank Test`"
isresponsive(baseline,response;alpha=0.05) = pvalue(SignedRankTest(baseline,response)) < alpha
"Check if any `sub group of response` is significently different from `baseline` by `Wilcoxon Signed Rank Test`"
isresponsive(baseline,response,gi;alpha=0.05) = any(map(i->isresponsive(baseline[i],response[i];alpha),gi))
isresponsive(baseline::Vector,response::Matrix;alpha=0.05) = any(isresponsive.(baseline,response;alpha))
"Check if `fun` of a spatial-temporal kernal within response time window significently higher than that of the baseline time window"
function isresponsive(st::AbstractMatrix;fun=rms,bi=[],ri=[],sdfactor=3)
    @views p = [fun(st[:,i]) for i in 1:size(st,2)]
    pmaxi = argmax(p); pmax = p[pmaxi]
    bpm = median(p[bi]);bpsd=mad(p[bi],center=bpm,normalize=true)

    (;r=!(pmaxi in bi) && (pmaxi in ri) && pmax > bpm + sdfactor*bpsd, p=pmax, zp=(pmax-bpm)/bpsd, d=pmaxi)
end

# "Check if any factors and their interactions significently modulate response by `ANOVA`"
# function ismodulative(df;alpha=0.05,interact=true)
#     xns = filter(i->i!=:Y,propertynames(df))
#     foreach(i->df[!,i]=categorical(df[!,i]),xns)
#     if interact
#         f = term(:Y) ~ reduce(+,map(i->reduce(&,term.(i)),combinations(xns)))
#     else
#         f = term(:Y) ~ reduce(+,term.(xns))
#     end
#     lmr = fit(LinearModel,f,df,contrasts = Dict(x=>EffectsCoding() for x in xns))
#     anovatype = length(xns) <= 1 ? 2 : 3
#     any(Anova(lmr,anovatype = anovatype).p[1:end-1] .< alpha)
# end

"""
Check if any `group in response` is significently different from at least one other `group in response`
"""
function ismodulative(response,gi;alpha=0.05,test=:anova)
    gr = map(i->response[i],gi)
    if test==:ranksum
        h=KruskalWallisTest(gr...)
    else
        h=OneWayANOVATest(gr...)
    end
    pvalue(h) < alpha
end


# ismodulative(response,gi;alpha=0.05) = PyOnewayANOVA.anova_oneway(map(i->response[i],gi),use_var="unequal").pvalue < alpha

"""
Find unique levels for each factor from condition tests, with number of repeat for each level and condition test indices for each repeat.
"""
flin(ctc::DataFrame;sort = true, skipmissing = true) = NamedTuple(f=>condin(ctc[!,[f]];sort,skipmissing) for f in propertynames(ctc))

"""
Find unique conditions of condition tests, with number of repeat for each condition and condition test indices for each repeat.
"""
function condin(ctc::DataFrame;sort = true, skipmissing = true)
    ivn = intersect(propertynames(ctc), (:i,:n))
    if !isempty(ivn)
        error("using factor names: $ivn, i and n are reserved names for condition test indices and number of condition repeat")
    end
    combine(groupby(insertcols(ctc,:i => 1:nrow(ctc)), propertynames(ctc); sort, skipmissing), :i => (x->[x]) => :i, nrow => :n)
end

"Get factor names excluding reserved names"
condfactor(cond)=setdiff(propertynames(cond),(:i,:n,:m,:se,:u))

"Exclude Non-Final version of factor names"
function finalfactor(cond)
    fs = String.(condfactor(cond))
    filter!(f->endswith(f,"_Final") || ("$(f)_Final" ∉ fs),fs)
    Symbol.(fs)
end

"String representation of a condition"
condstring(cond::DataFrameRow;factor=condfactor(cond)) = join(["$f=$(cond[f])" for f in factor],", ")
condstring(cond::DataFrame;factor=condfactor(cond)) = [condstring(c;factor) for c in eachrow(cond)]


"""
Get `Mean` and `SEM` of repeated responses for each condition

1. rs: response of each condition test
2. ci: condition test indices of repeats for each condition
"""
function condresponse(rs,ci)
    crs = [rs[i] for i in ci]
    (m=mean.(crs),se=sem.(crs))
end
function condresponse(rs,cond::DataFrame;u=0,withcond::Bool=true)
    cr = DataFrame(pairs(condresponse(rs,cond.i))...,:u=>u)
    if withcond
        ivn = intersect(propertynames(cond), (:m,:se,:u))
        if !isempty(ivn)
            error("using factor names: $ivn, m, se and u are reserved names for condition response of unit")
        end
        cr = [cr cond[:,condfactor(cond)]]
    end
    cr
end
condresponse(urs::Dict,cond::DataFrame;withcond::Bool=true) = mapreduce(u->condresponse(urs[u],cond;u,withcond),append!,keys(urs))

"""
Get `Mean` and `SEM` of repeated image responses for each condition

1. rs: response of each condition test [height, width, ncondtest]
2. ci: condition test indices of repeats for each condition
"""
function condresponse(rs::AbstractArray{T,3},ci;sfun=nothing,mfun=nothing) where T
    m = Array{Float64}(undef,size(rs)[1:2]...,length(ci))
    se = similar(m)
    for i in eachindex(ci)
        @views r = meanse(rs[:,:,ci[i]];dims=3,sfun,mfun)
        m[:,:,i] = r.m; se[:,:,i] = r.se
    end
    (;m,se)
end

"Get each factor name and its levels as the axes of factor space"
function factoraxis(cond;factor=condfactor(cond))
    fl = flin(cond[!,factor])
    NamedTuple(f=>fl[f][!,f] for f in keys(fl))
end

"Map a variable of conditions into factor space"
function factorspace(cond;fa=factoraxis(cond),col=:i)
    fc = missings(Any,map(length,values(fa)))
    for c in eachrow(cond)
        ci = [findfirst(l->l==c[f],fa[f]) for f in keys(fa)]
        fc[ci...] = c[col]
    end
    (;fc,fa)
end

"""
Get `Mean` and `SEM` of repeated responses for each condition in factor space

1. rs: response of each condition test
2. fi: condition test indices of repeats for each condition in factor space
"""
function factorresponse(rs,fi)
    r = map(i->ismissing(i) ? missing : rs[i],fi)
    m = mean.(r)
    se = sem.(r)
    (;r,m,se)
end
function factorresponse(rs::AbstractMatrix,fi)
    r = map(i->ismissing(i) ? missing : rs[i,:],fi)
    mse = meanse.(r;dims=1)
    (;r,m=first.(mse),se=last.(mse))
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
