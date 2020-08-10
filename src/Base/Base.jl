using LinearAlgebra,Distributions,DataFrames,StatsBase,GLM,LsqFit,Optim,HypothesisTests,Colors,Images
using ImageFiltering,SpecialFunctions,DSP,HCubature,Combinatorics,DataStructures,ANOVA,StatsFuns,Trapz
using ImageSegmentation,ProgressMeter
import Base: vec,range

include("NeuroDataType.jl")
include("CircStats.jl")
include("Spike.jl")
include("LFP.jl")
include("Image.jl")
include("Condition.jl")
include("2P.jl")

vec(x::RGBA)=[x.r,x.g,x.b,x.alpha]
anscombe(x) = 2*sqrt(x+(3/8))

"Check if `response` is significently different from `baseline` by `Wilcoxon Signed Rank Test`"
isresponsive(baseline,response;alpha=0.05) = pvalue(SignedRankTest(baseline,response)) < alpha
"Check if any `sub group of response` is significently different from `baseline` by `Wilcoxon Signed Rank Test`"
isresponsive(baseline,response,gi;alpha=0.05) = any(map(i->isresponsive(baseline[i],response[i],alpha=alpha),gi))
isresponsive(baseline::Vector,response::Matrix;alpha=0.05) = any(isresponsive.(baseline,response,alpha=alpha))

"Check if any factors and their interactions significently modulate response using ANOVA"
function ismodulative(df;alpha=0.05,interact=true)
    xns = filter(i->i!=:Y,propertynames(df))
    categorical!(df,xns)
    if interact
        f = term(:Y) ~ reduce(+,map(i->reduce(&,term.(i)),combinations(xns)))
    else
        f = term(:Y) ~ reduce(+,term.(xns))
    end
    lmr = fit(LinearModel,f,df,contrasts = Dict(x=>EffectsCoding() for x in xns))
    anovatype = length(xns) <= 1 ? 2 : 3
    any(Anova(lmr,anovatype = anovatype).p[1:end-1] .< alpha)
end

"`Gaussian` function"
gaussianf(x;a=1,μ=0,σ=1) = a*exp(-0.5((x-μ)/σ)^2)
function gaussianf(x,y;a=1,μ₁=0,σ₁=1,μ₂=0,σ₂=1,θ=0)
    sinv,cosv = sincos(θ)
    x′ = cosv * x + sinv * y
    y′ = cosv * y - sinv * x
    a*exp(-0.5(((x′-μ₁)/σ₁)^2 + ((y′-μ₂)/σ₂)^2))
end

"""
`von Mises` function [^1]

```math
f(α) =  βe^{κ(cos(n(α - μ)) - 1)}
```

[^1]

Swindale, N.V. (1998). Orientation tuning curves: empirical description and estimation of parameters. Biol Cybern 78, 45–56.

- β: amplitude at μ
- μ: angle of peak
- κ: width parameter
- n: frequency parameter
"""
vmf(α;β=1,μ=0,κ=1,n=1) = β*exp(κ*(cos(n*(α-μ))-1))
"""
`Generalized von Mises` function [^1]

```math
f(α) =  βe^{κ₁(cos(α - μ₁) - 1) + κ₂(cos2(α - μ₂) - 1)}
```

[^1]

Gatto, R., and Jammalamadaka, S.R. (2007). The generalized von Mises distribution. Statistical Methodology 4, 341–353.
"""
gvmf(α;β=1,μ₁=0,κ₁=1,μ₂=0,κ₂=1) = β*exp(κ₁*(cos(α-μ₁)-1) + κ₂*(cos(2(α-μ₂))-1))

"""
`Difference of Gaussians` function
"""
dogf(x;aₑ=2,μₑ=0,σₑ=1,aᵢ=1,μᵢ=0,σᵢ=2) = gaussianf(x,a=aₑ,μ=μₑ,σ=σₑ) - gaussianf(x,a=aᵢ,μ=μᵢ,σ=σᵢ)
function dogf(x,y;aₑ=2,μₑ₁=0,σₑ₁=1,μₑ₂=0,σₑ₂=1,θₑ=0,aᵢ=1,μᵢ₁=0,σᵢ₁=2,μᵢ₂=0,σᵢ₂=2,θᵢ=0)
    sinvₑ,cosvₑ = sincos(θₑ)
    xₑ′ = cosvₑ * x + sinvₑ * y
    yₑ′ = cosvₑ * y - sinvₑ * x
    sinvᵢ,cosvᵢ = sincos(θᵢ)
    xᵢ′ = cosvᵢ * x + sinvᵢ * y
    yᵢ′ = cosvᵢ * y - sinvᵢ * x
    aₑ*exp(-0.5(((xₑ′-μₑ₁)/σₑ₁)^2 + ((yₑ′-μₑ₂)/σₑ₂)^2)) - aᵢ*exp(-0.5(((xᵢ′-μᵢ₁)/σᵢ₁)^2 + ((yᵢ′-μᵢ₂)/σᵢ₂)^2))
end

"""
`sin` grating function

- f: Frequency in cycle/unit_x
- phase: Phase of a cycle in [0, 1] scale
"""
gratingf(x; f=1, phase=0) = sin(2π * (f * x + phase))

"""
`cas` function defined as ``cas(x) = cos(x) + sin(x)``

- f: Frequency in cycle/unit_x
- phase: Phase of a cycle in [0, 1] scale
- isnorm: scale `cas` in [-√2, √2] to [-1, 1]
"""
function cas(x;f=1, phase=0, isnorm::Bool=true)
    r = sum(sincos(2π * (f * x + phase)))
    if isnorm
        r /=sqrt(2)
    end
    return r
end

"""
2D `sin` grating function

- θ: Orientation in radius, 0 is -, increase counter-clock wise
- f: Frequency in cycle/unit_x/y
- phase: Phase of a cycle in [0, 1] scale
"""
function gratingf(x,y; θ=0,f=1,phase=0)
    sinθ,cosθ = sincos(θ)
    y′ = cosθ * y - sinθ * x
    sin(2π * (f * y′ + phase))
end

"""
2D `cas` function defined as ``cas(x+y) = cos(x+y) + sin(x+y)``

- kx: Frequency in cycle/unit_x
- ky: Frequency in cycle/unit_y
- phase: Phase of a cycle in [0, 1] scale
- isnorm: scale `cas` in [-√2, √2] to [-1, 1]
"""
function cas(x,y;kx=1,ky=1, phase=0, isnorm::Bool=true)
    r = sum(sincos(2π * (kx * x + ky * y + phase)))
    isnorm && (r /=sqrt(2))
    return r
end

"""
`cas` phase to `sin` phase, phase is in [0, 1] scale

!!! note
    ``cas(x) = √2 sin(x + π/4)``
"""
cas2sin(phase) = phase + 0.125

"""
2D `cas` to 2D sin `gratingf`
"""
function cas2sin(kx,ky,phase)
    θ = atan(ky,kx) - π/2
    f = sqrt(kx*kx + ky*ky)
    return (θ=θ,f=f,phase=phase + 0.125)
end

"""
`sin` phase to `cas` phase, phase is in [0, 1] scale

!!! note
    ``cas(x) = √2 sin(x + π/4)``
"""
sin2cas(phase) = phase - 0.125

"""
2D sin `gratingf` to 2D `cas`
"""
function sin2cas(θ,f,phase)
    sinθ,cosθ = sincos(θ + π/2)
    kx = cosθ*f
    ky = sinθ*f
    return (kx=kx,ky=ky,phase=phase - 0.125)
end

"`Gabor` function"
gaborf(x;a=1,μ=0,σ=1,f=1,phase=0) = gaussianf(x,a=a,μ=μ,σ=σ)*sin(2π*(f*x+phase))
function gaborf(x,y;a=1,μ₁=0,σ₁=1,μ₂=0,σ₂=1,θ=0,f=1,phase=0)
    sinv,cosv = sincos(θ)
    x′ = cosv * x + sinv * y
    y′ = cosv * y - sinv * x
    a*exp(-0.5(((x′-μ₁)/σ₁)^2 + ((y′-μ₂)/σ₂)^2)) * sin(2π*(f * y′ + phase))
end

"Fit model to data"
function fitmodel(model,x,y)
    alb,aub = abs.(extrema(y))
    ab = max(alb,aub)

    rlt = fun = missing
    if model == :vmn2
        fun = (x,p) -> vmf.(x,β=p[1],μ=p[2],κ=p[3],n=2)
        ofun = (p;x=x,y=y) -> sum((y.-fun(x,p)).^2)

        ub=[1.3ab,   prevfloat(float(2π)),   20]
        lb=[0.3ab,            0,              0]
        p0=[ab,               π,              1]
    elseif model == :gvm
        fun = (x,p) -> gvmf.(x,β=p[1],μ₁=p[2],κ₁=p[3],μ₂=p[4],κ₂=p[5])
        ofun = (p;x=x,y=y) -> sum((y.-fun(x,p)).^2)

        ub=[1.3ab,   prevfloat(float(2π)),   20,    prevfloat(float(2π)),    20]
        lb=[0.3ab,            0,              0,             0,               0]
        p0=[ab,               π,              1,             π,               1]
    end
    if !ismissing(fun)
        ofit = optimize(ofun,lb,ub,p0,SAMIN(rt=0.9),Optim.Options(iterations=200000))
        param=ofit.minimizer

        rlt = (;model,fun,param)
    end
    return rlt
end

function circtuningfeature(mfit;od=π,fn=od==π ? :d : :o)
    x = 0:0.004:2π # 0.004rad = 0.23deg
    circtuningfeature(x,mfit.fun(x,mfit.param),od=od,fn=fn)
end

"""
Properties of Circular Tuning:
    - Prefered Direction/Orientation
    - Selectivity Index
        - version 1: (ResponsePrefered - ResponseOpposing)/ResponsePrefered
        - version 2: (ResponsePrefered - ResponseOpposing)/(ResponsePrefered + ResponseOpposing)
    - Full Width at Half Maximum

1. x: angles in radius
2. y: responses
- od: opposing angle distance, π for DSI, 0.5π for OSI
- fn: factor name
"""
function circtuningfeature(x,y;od=π,fn=od==π ? :d : :o)
    maxi = argmax(y)
    maxr = y[maxi]
    px = x[maxi]
    ox = px+od
    oi = findclosestangle(ox,x)
    or = y[oi]

    si1 = 1-or/maxr
    si2 = (maxr-or)/(maxr+or)

    # minr = minimum(y)
    # hmaxr = minr + (maxr-minr)/2
    (;Symbol(:p,fn)=>rad2deg(mod(px,2od)), Symbol(fn,:si1)=>si1, Symbol(fn,:si2)=>si2)
end

"""
Tuning properties of factor response

1. fl: factor levels
2. fr: factor responses

    HueAngle, Orientation and Direction follow the same convention such that 0 is -/→, then increase counter-clock wise.

    For cases where Orientation and Direction are interlocked(drifting grating):
        - when Orientation is -(0), then Direction is ↑(90)
        - when Direction is →(0), then Orientation is |(-90)
"""
function factorresponsefeature(fl,fr;factor=:Ori,isfit::Bool=true)
    if factor in [:Ori,:Ori_Final]
        θ = deg2rad.(fl)
        d = mean(diff(sort(unique(θ)))) # angle spacing
        # for orientation
        oθ = mod.(θ,π)
        om = circmean(2oθ,fr)
        oo = rad2deg(mod(angle(om),2π)/2)
        ocv = circvar(2oθ,fr,2d)
        # for direction
        dm = circmean(θ.+0.5π,fr)
        od = rad2deg(mod(angle(dm),2π))
        dcv = circvar(θ.+0.5π,fr,d)

        fit = ()
        if isfit
            # fit Generalized von Mises for direction
            try
                mfit = fitmodel(:gvm,θ.+0.5π,fr)
                fit = (circtuningfeature(mfit,od=π,fn=:d)...,gvm=mfit)
            catch
            end
            # fit von Mises for orientation
            try
                mfit = fitmodel(:vmn2,θ,fr)
                fit = (fit...,circtuningfeature(mfit,od=0.5π,fn=:o)...,vmn2=mfit)
            catch
            end
        end

        return (;dm,od,dcv,om,oo,ocv,fit)
    elseif factor == :Dir
        θ = deg2rad.(fl)
        d = mean(diff(sort(unique(θ)))) # angle spacing
        # for orientation
        oθ = mod.(θ.-0.5π,π)
        om = circmean(2oθ,fr)
        oo = rad2deg(mod(angle(om),2π)/2)
        ocv = circvar(2oθ,fr,2d)
        # for direction
        dm = circmean(θ,fr)
        od = rad2deg(mod(angle(dm),2π))
        dcv = circvar(θ,fr,d)
        # fit Generalized von Mises for direction
        fit = ()
        if isfit
            try
                gvmfit = curve_fit((x,p)->gvmf.(x,p...),θ,fr,[1.0,0,1,0,1])
                if gvmfit.converged
                    x = 0:0.004:2π # 0.004rad = 0.23deg
                    y = gvmf.(x,gvmfit.param...)
                    fit = (circtuningstats(x,y,od=π,s=:d)...,gvm=gvmfit)
                end
            catch
            end
            # fit von Mises for orientation
            try
                vmfit = curve_fit((x,p)->vmf.(x,p...,n=2),θ.-0.5π,fr,[1.0,0,1])
                if vmfit.converged
                    x = 0:0.004:2π
                    y = vmf.(x,vmfit.param...,n=2)
                    fit = (fit...,circtuningstats(x,y,od=0.5π,s=:o)...,vm=vmfit)
                end
            catch
            end
        end

        return (dm=dm,od=od,dcv=dcv,om=om,oo=oo,ocv=ocv,fit=fit)
    elseif factor == :SpatialFreq
        osf = 2^(sum(fr.*log2.(fl))/sum(fr)) # weighted average as optimal sf
        fit=()
        if isfit
            # fit difference of gaussians
            try
                dogfit = curve_fit((x,p)->dogf.(x,p...),fl,fr,[1.0,0,1,1,0,1])
                if dogfit.converged
                    fit = (dog=dogfit,)
                end
            catch
            end
        end

        return (osf = osf,fit=fit)
    elseif factor == :ColorID
        # transform colorId to hue angle
        ucid = sort(unique(fl))
        hstep = 2pi/length(ucid)
        ha = map(l->hstep*(findfirst(c->c==l,ucid)-1),fl)
        oh = mod(rad2deg(circmean(ha,fr)),360)
        ohv = circmeanv(ha,fr)
        ohr = circr(ha,fr)
        hcv = circvar(ha,fr)
        hm = circmean(ha,fr)
        oh = mod(rad2deg(angle(hm)),360)
        hcv = circvar(ha,fr,hstep)

        return (hm=hm,oh=oh,hcv=hcv)
    elseif factor == :HueAngle
        θ = deg2rad.(fl)
        d = mean(diff(sort(unique(θ)))) # angle spacing
        # for hue axis
        aθ = mod.(θ,π)
        ham = circmean(2aθ,fr)
        oha = rad2deg(mod(angle(ham),2π)/2)
        hacv = circvar(2aθ,fr,2d)
        # for hue
        hm = circmean(θ,fr)
        oh = rad2deg(mod(angle(hm),2π))
        hcv = circvar(θ,fr,d)
        maxi = argmax(fr)
        maxh = fl[maxi]
        maxr = fr[maxi]

        fit = ()
        if isfit
            # fit Generalized von Mises for hue
            try
                mfit = fitmodel(:gvm,θ,fr)
                fit = (circtuningfeature(mfit,od=π,fn=:h)...,gvm=mfit)
            catch
            end
            # fit von Mises for hue axis
            try
                mfit = fitmodel(:vmn2,θ,fr)
                fit = (fit...,circtuningfeature(mfit,od=0.5π,fn=:ha)...,vmn2=mfit)
            catch
            end
        end

        return (;ham,oha,hacv,hm,oh,hcv,maxh,maxr,fit)
    else
        return ()
    end
end

"""
Spike Triggered Average of Images

1. x: Matrix where each row is one image
2. y: Vector of image response

- norm: normalization factor, default no normalization.
        it could be ``sum(y)`` if y is number of spike or spike rate, then STA would be spiking probability.
- whiten: whiten factor, default no whiten.
        it could be ``(xᵀx)⁻¹`` or inverse of covariance matrix, that decorrelate STA.
"""
function sta(x::AbstractMatrix,y::AbstractVector;norm=nothing,whiten=nothing)
    r = x'*y
    !isnothing(norm) && (r/=norm)
    !isnothing(whiten) && (r=whiten*r)

    # r = x'*x\r
    # r=length(y)*inv(cov(x,dims=1))*r

    return r
end


function psthsts(xs::Vector,binedges::Vector,c;israte::Bool=true,normfun=nothing)
    m,se,x = psth(xs,binedges,israte=israte,normfun=normfun)
    df = DataFrame(x=x,m=m,se=se,c=fill(c,length(x)))
end
function psthsts(xs::Vector,binedges::Vector,cond::DataFrame;israte::Bool=true,normfun=nothing)
    fs = finalfactor(cond)
    vcat([psth(xs[r[:i]],binedges,condstring(r,fs),israte=israte,normfun=normfun) for r in eachrow(cond)]...)
end
function psthsts(xs::Vector,binedges::Vector,ctc::DataFrame,factor;israte::Bool=true,normfun=nothing)
    vf = intersect(names(ctc),factor)
    isempty(vf) && error("No Valid Factor Found.")
    psth(xs,binedges,condin(ctc[:,vf]),israte=israte,normfun=normfun)
end
function psth(ds::DataFrame,binedges::Vector,conds::Vector;normfun=nothing,spike=:spike,isse::Bool=true)
    is,ss = findcond(ds,conds)
    df = psth(map(x->ds[spike][x],is),binedges,ss,normfun=normfun)
    if isse
        df[:ymin] = df[:y]-df[:ysd]./sqrt(df[:n])
        df[:ymax] = df[:y]+df[:ysd]./sqrt(df[:n])
    end
    return df,ss
end
function psthstss(xss::Vector,binedges::Vector,conds;normfun=nothing)
    n = length(xss)
    n!=length(conds) && error("Length of xss and conds don't match.")
    dfs = [psth(xss[i],binedges,conds[i],normfun=normfun) for i=1:n]
    return cat(1,dfs)
end
function spacepsth(unitpsth,unitposition,spacebinedges)
    ys,ns,ws,is = epochspiketrain(unitposition[:,2],spacebinedges)
    x = unitpsth[1][3]
    nbins = length(is)
    um = map(i->i[1],unitpsth)
    use = map(i->i[2],unitpsth)
    spsth = zeros(nbins,length(x))
    for i in 1:nbins
        if length(is[i]) > 0
            spsth[i,:] = mean(hcat(um[is[i]]...),dims=2)
        end
    end
    binwidth = ws[1][2]-ws[1][1]
    bincenters = [ws[i][1]+binwidth/2.0 for i=1:nbins]
    return spsth,x,bincenters,ns
end


"""
Shift(shuffle) corrected, normalized(coincidence/spike), condition/trial-averaged Cross-Correlogram of binary spike trains.
(Bair, W., Zohary, E., and Newsome, W.T. (2001). Correlated Firing in Macaque Visual Area MT: Time Scales and Relationship to Behavior. J. Neurosci. 21, 1676–1697.)
"""
function correlogram(bst1,bst2;lag=nothing,isnorm=true,shiftcorrection=true,condis=nothing)
    if !isnothing(condis)
        cccg=[];x=[]
        for ci in condis
            ccg,x = correlogram(bst1[:,ci],bst2[:,ci],lag=lag,isnorm=isnorm,shiftcorrection=shiftcorrection,condis=nothing)
            push!(cccg,ccg)
        end
        ccg=dropdims(mean(hcat(cccg...),dims=2),dims=2)
        return ccg,x
    end
    n,nepoch = size(bst1)
    lag = floor(Int,isnothing(lag) ? min(n-1, 10*log10(n)) : lag)
    x = -lag:lag;xn=2lag+1
    cc = Array{Float64}(undef,xn,nepoch)
    for k in 1:nepoch
        cc[:,k]=crosscov(bst1[:,k],bst2[:,k],x,demean=false)*n
    end
    ccg = dropdims(mean(cc,dims=2),dims=2)
    if isnorm
        λ1 = mean(mean(bst1,dims=1))
        λ2 = mean(mean(bst2,dims=1))
        gmsr = sqrt(λ1*λ2)
        Θ = n.-abs.(x)
        normfactor = 1 ./ Θ ./ gmsr
    end
    if shiftcorrection
        psth1 = dropdims(mean(bst1,dims=2),dims=2)
        psth2 = dropdims(mean(bst2,dims=2),dims=2)
        s = crosscov(psth1,psth2,x,demean=false)*n
        shiftccg = (nepoch*s .- ccg)/(nepoch-1)
        if isnorm
            ccg .*= normfactor
            shiftccg .*= normfactor
        end
        ccg .-= shiftccg
    elseif isnorm
        ccg .*= normfactor
    end
    ccg,x
end
function circuitestimate(unitbinspike;lag=nothing,maxprojlag=3,minepoch=5,minspike=10,esdfactor=5,isdfactor=3.5,unitid=[],condis=nothing)
    nunit=length(unitbinspike)
    n,nepoch = size(unitbinspike[1])
    lag = floor(Int,isnothing(lag) ? min(n-1, 10*log10(n)) : lag)
    x = -lag:lag;xn=2lag+1

    isenoughspike = (i,j;minepoch=5,minspike=10) -> begin
        vsi = sum(i,dims=1)[:] .>= minspike
        vsj = sum(j,dims=1)[:] .>= minspike
        count(vsi) >= minepoch && count(vsj) >= minepoch
    end

    ccgs=[];ccgis=[];projs=[];eunits=[];iunits=[];projweights=[]
    for (i,j) in combinations(1:nunit,2)
        if isnothing(condis)
            !isenoughspike(unitbinspike[i],unitbinspike[j],minepoch=minepoch,minspike=minspike) && continue
            vcondis=nothing
        else
            vci = map(ci->isenoughspike(unitbinspike[i][:,ci],unitbinspike[j][:,ci],minepoch=minepoch,minspike=minspike),condis)
            all(.!vci) && continue
            vcondis = condis[vci]
        end

        ccg,_ = correlogram(unitbinspike[i],unitbinspike[j],lag=lag,condis=vcondis)
        ps,es,is,pws = projectionfromcorrelogram(ccg,i,j,maxprojlag=maxprojlag,esdfactor=esdfactor,isdfactor=isdfactor)
        if !isempty(ps)
            push!(ccgs,ccg);push!(ccgis,(i,j))
            append!(projs,ps);append!(eunits,es);append!(iunits,is);append!(projweights,pws)
        end
    end
    unique!(eunits);unique!(iunits)
    if length(unitid)==nunit
        map!(t->(unitid[t[1]],unitid[t[2]]),ccgis,ccgis)
        map!(t->(unitid[t[1]],unitid[t[2]]),projs,projs)
        map!(t->unitid[t],eunits,eunits)
        map!(t->unitid[t],iunits,iunits)
    end
    return ccgs,x,ccgis,projs,eunits,iunits,projweights
end
function projectionfromcorrelogram(cc,i,j;maxprojlag=3,minbaselag=maxprojlag+1,esdfactor=5,isdfactor=5)
    midi = Int((length(cc)+1)/2)
    base = vcat(cc[midi+minbaselag:end],cc[1:midi-minbaselag])
    bm,bsd = mean_and_std(base);hl = bm + esdfactor*bsd;ll = bm - isdfactor*bsd
    forwardlags = midi+1:midi+maxprojlag
    backwardlags = midi-maxprojlag:midi-1
    fcc = cc[forwardlags]
    bcc = cc[backwardlags]
    ps=[];ei=[];ii=[];pws=[]

    if ll <= cc[midi] <= hl
        if any(fcc .> hl)
            push!(ps,(i,j));push!(ei,i);push!(pws,(maximum(fcc)-bm)/bsd)
        elseif any(fcc .< ll)
            push!(ps,(i,j));push!(ii,i);push!(pws,(minimum(fcc)-bm)/bsd)
        end
        if any(bcc .> hl)
            push!(ps,(j,i));push!(ei,j);push!(pws,(maximum(bcc)-bm)/bsd)
        elseif any(bcc .< ll)
            push!(ps,(j,i));push!(ii,j);push!(pws,(minimum(bcc)-bm)/bsd)
        end
    end
    return ps,ei,ii,pws
end

function checklayer!(ls::Dict)
    ln=["WM","6","5","5/6","4Cb","4Ca","4C","4B","4A","4A/B","3","2","2/3","1","Out"]
    n = length(ln)
    for i in 1:n-1
        if haskey(ls,ln[i])
            for j in (i+1):n
                if haskey(ls,ln[j])
                    ls[ln[i]][2] = ls[ln[j]][1]
                    break
                end
            end
        end
    end
    return ls
end

"""
Try to locate cell layer.
1. y coordinate of cell postion
2. layer definition
"""
function assignlayer(y,layer)
    l = missing
    for k in keys(layer)
        if layer[k][1] <= y < layer[k][2]
            l=k;break
        end
    end
    return l
end

function checkcircuit(projs,eunits,iunits,projweights)
    ivu = intersect(eunits,iunits)
    veunits = setdiff(eunits,ivu)
    viunits = setdiff(iunits,ivu)
    ivp = map(p->any(i->i in ivu,p),projs)
    vprojs = deleteat!(copy(projs),ivp)
    vprojweights = deleteat!(copy(projweights),ivp)
    return vprojs,veunits,viunits,vprojweights
end
