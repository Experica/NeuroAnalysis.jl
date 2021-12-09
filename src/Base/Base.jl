using LinearAlgebra,Distributions,DataFrames,StatsBase,GLM,LsqFit,Optim,HypothesisTests,Colors,Images,StatsModels,Distances,CategoricalArrays,
ImageFiltering,SpecialFunctions,DSP,HCubature,Combinatorics,DataStructures,ANOVA,StatsFuns,Trapz, ImageSegmentation,ProgressMeter
import Base: vec,range
import StatsBase: predict

include("NeuroDataType.jl")
include("CircStats.jl")
include("Spike.jl")
include("LFP.jl")
include("Image.jl")
include("Condition.jl")
include("2P.jl")

vec(x::RGBA)=[x.r,x.g,x.b,x.alpha]
anscombe(x) = 2*sqrt(x+(3/8))

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

"Fit 1D model to data"
function fitmodel(model,x,y)
    lb,ub = extrema(y)
    bm = (lb+ub)/2
    alb,aub = abs.((lb,ub))
    ab = max(alb,aub)

    xlb,xub = extrema(x)
    xbm = (xlb+xub)/2
    xalb,xaub = abs.((xlb,xub))
    xab = max(xalb,xaub)

    rlt = fun = missing
    if model == :vmn2
        fun = (x,p) -> vmf.(x,β=p[1],μ=p[2],κ=p[3],n=2)
        ofun = (p;x=x,y=y) -> sum((y.-fun(x,p)).^2)

        ub=[1.3ab,   prevfloat(float(2π)),   30]
        lb=[0.3ab,            0,              0]
        p0=[ab,               π,              1]
    elseif model == :gvm
        fun = (x,p) -> gvmf.(x,β=p[1],μ₁=p[2],κ₁=p[3],μ₂=p[4],κ₂=p[5])
        ofun = (p;x=x,y=y) -> sum((y.-fun(x,p)).^2)

        ub=[1.3ab,   prevfloat(float(2π)),   30,    prevfloat(float(2π)),    30]
        lb=[0.3ab,            0,              0,             0,               0]
        p0=[ab,               π,              1,             π,               1]
    elseif model == :dog
        fun = (x,p) -> dogf.(x,aₑ=p[1],μₑ=p[2],σₑ=p[3],aᵢ=p[4],μᵢ=p[5],σᵢ=p[6]) .+ p[7]
        ofun = (p;x=x,y=y) -> sum((y.-fun(x,p)).^2)

        ub=[1.3ab,   10xab,   10xab,    1.3ab,   10xab,   10xab,    ub]
        lb=[0,      -10xab,     0,        0,    -10xab,     0,      lb]
        p0=[ab,        0,       1,       ab,       0,       1,      bm]
    end
    if !ismissing(fun)
        ofit = optimize(ofun,lb,ub,p0,SAMIN(rt=0.9),Optim.Options(iterations=200000))
        param=ofit.minimizer; yy = fun(x,param)

        rlt = (;model,fun,param,goodnessoffit(y,yy,k=length(param))...)
    end
    return rlt
end

"Fit 2D model to image"
function fitmodel2(model,data::Matrix,ppu;w=0.5)
    rpx = (size(data)[1]-1)/2
    radius = rpx/ppu
    x = (mapreduce(i->[i[2] -i[1]],vcat,CartesianIndices(data)) .+ [-(rpx+1) (rpx+1)])/ppu
    y = vec(data)

    # try estimate solution
    roi = peakroi(localcontrast(data,round(Int,w*ppu)))
    alb,aub = abs.(extrema(data[roi.i]))
    ab = max(alb,aub)
    r = roi.radius/ppu
    c = [roi.center[2] - (rpx+1), -roi.center[1] + (rpx+1)]/ppu

    rlt = fun = missing
    if model == :dog
        if aub >= alb
            ai = 5alb
            ae = aub + ai
        else
            ae = 5aub
            ai = alb + ae
        end
        fun = (x,y,p) -> dogf.(x,y,aₑ=p[1],μₑ₁=p[2],σₑ₁=p[3],μₑ₂=p[4],σₑ₂=p[3],θₑ=0,aᵢ=p[5],μᵢ₁=p[2],σᵢ₁=p[6],μᵢ₂=p[4],σᵢ₂=p[6],θᵢ=0)
        ofun = (p;x=x,y=y) -> sum((y.-fun(x[:,1],x[:,2],p)).^2)
        ub=[5ae,    0.5r+c[1],    0.9r,    0.5r+c[2],     5ai,    0.9r]
        lb=[0,     -0.5r+c[1],    0.1r,   -0.5r+c[2],     0,      0.1r]
        p0=[ae,      c[1],        0.3r,     c[2],         ai,     0.3r]
    elseif model == :gabor
        fun = (x,y,p) -> gaborf.(x,y,a=p[1],μ₁=p[2],σ₁=p[3],μ₂=p[4],σ₂=p[5],θ=p[6],f=p[7],phase=p[8])
        ofun = (p;x=x,y=y) -> sum((y.-fun(x[:,1],x[:,2],p)).^2)

        ori,sf = f1orisf(powerspectrum2(data,ppu)...)
        ub=[5ab,   0.5r+c[1],   0.9r,    0.5r+c[2],    0.9r,      prevfloat(float(π)),     12,     prevfloat(1.0)]
        lb=[0,    -0.5r+c[1],   0.1r,   -0.5r+c[2],    0.1r,                0,             0.05,           0]
        p0=[ab,    c[1],        0.3r,     c[2],        0.3r,               ori,            sf,            0.5]
    end
    if !ismissing(fun)
        ofit = optimize(ofun,lb,ub,p0,SAMIN(rt=0.9),Optim.Options(iterations=200000))
        param=ofit.minimizer; yy = fun(x[:,1],x[:,2],param); resid = y .- yy

        rlt = (;model,fun,param,radius,resid,goodnessoffit(y,yy,e=resid,k=length(param))...)
    end
    return rlt
end

predict(fit,x) = fit.fun(x,fit.param)
function predict(fit,x,y;xygrid=true,yflip=false)
    if xygrid
        z = [fit.fun(i,j,fit.param) for j in y, i in x]
        yflip && (z=reverse(z,dims=1))
    else
        z = fit.fun(x,y,fit.param)
    end
    z
end

"""
Goodness of Fit Metrics:

- r: Pearson Correlation Coefficient
- mae: Mean Absolute Error
- rmse: Root Mean Squared Error
- rae: Relative Absolute Error
- rse: Relative Squared Error
- r2: R Squared
- adjr2: Adjusted-R²
- s: Residual Standard Error
- aic: Akaike Information Criterion
- bic: Bayesian Information Criterion

1. y: responses
2. ŷ: model predictions

- n: sample size
- e: errors(y - ŷ)
- k: number of predictors
- df: degree of freedom(n - k - 1)
"""
function goodnessoffit(y,ŷ;n = length(y),e = y .- ŷ,k=missing,df = n-k-1)
    r = cor(y,ŷ)
    ae = abs.(e)
    e2 = e.^2
    ssᵣ = sum(e2)
    mae = mean(ae)
    rmse = sqrt(mean(e2))
    ydm = y .- mean(y)
    ssₜ = sum(ydm.^2)
    rae = sum(ae)/sum(abs.(ydm))
    fvu = ssᵣ/ssₜ
    rse = sqrt(fvu)
    r2 = 1 - fvu
    s = ssᵣ/df
    s = s < 0 ? missing : sqrt(s)
    adjr2 = 1 - fvu*(n-1)/df
    aic = n*log(ssᵣ) + 2k
    bic = n*log(ssᵣ/n) + k*log(n)
    (;r,mae,rmse,rae,rse,r2,adjr2,s,aic,bic)
end

function searchclosest(v,vs;start::Integer=1,step::Integer=1,circ=false)
    n=length(vs);ssign = sign(vs[start]-v)
    i = start
    for _ in 1:n
        if !circ
            i<1 && return -Inf
            n<i && return Inf
        end
        sign(vs[i]-v) != ssign && return i
        i += step
        if circ
            i<1 && (i+=n)
            n<i && (i-=n)
        end
    end
    Inf
end

"""
left and right half width of `v` relative to `y[start]`, -Inf/Inf when no `v` is found.

- start: index of `y` which is the center of the width
- v: value on which width is cutoff
- circ: whether `y` is defined on circular domain and width can wrap around
- x: domain of `y`, return width when `x` provided, otherwise return cutoff indices
"""
function halfwidth(y;start=argmax(y),v=y[start]/2,circ=false,x=nothing)
    li = searchclosest(v,y;start,step=-1,circ)
    ri = searchclosest(v,y;start,step=1,circ)
    if isnothing(x)
        return li,ri
    else
        if circ
            lw = isinf(li) ? li : li<=start ? abs(x[start]-x[li]) : abs(x[start]-x[1])+abs(x[end]-x[li])
            rw = isinf(ri) ? ri : ri>=start ? abs(x[ri]-x[start]) : abs(x[ri]-x[1])+abs(x[end]-x[start])
        else
            lw = isinf(li) ? li : abs(x[start]-x[li])
            rw = isinf(ri) ? ri : abs(x[ri]-x[start])
        end
        return lw,rw
    end
end

function circtuningfeature(mfit;od=π,fn=od==π ? :d : :o)
    x = 0:0.002:2π # 0.002rad ≈ 0.11deg
    circtuningfeature(x,predict(mfit,x),od=od,fn=fn)
end

"""
Properties of Circular Tuning

    - Prefered Direction/Orientation
    - Selectivity Index
        - version 1: (ResponsePrefered - ResponseOpposing)/ResponsePrefered
        - version 2: (ResponsePrefered - ResponseOpposing)/(ResponsePrefered + ResponseOpposing)
    - Half Width at Half Peak-to-Trough

1. x: angles in radius
2. y: responses
- od: opposing angle distance, π for DSI, 0.5π for OSI
- fn: factor name
"""
function circtuningfeature(x,y;od=π,fn=od==π ? :d : :o)
    maxi = argmax(y)
    mini = argmin(y)
    maxr = y[maxi]
    minr = y[mini]
    px = x[maxi]
    ox = px+od
    oi = findclosestangle(ox,x)
    or = y[oi]

    si1 = 1-or/maxr
    si2 = (maxr-or)/(maxr+or)
    hw = halfwidth(y,start=maxi,v=(maxr+minr)/2,circ=true,x=x)

    (;Symbol(:p,fn)=>rad2deg(mod(px,2od)), Symbol(fn,:hw)=>rad2deg.(hw), Symbol(fn,:si1)=>si1, Symbol(fn,:si2)=>si2)
end

function sftuningfeature(mfit)
    x = 0:0.002:10
    sftuningfeature(x,predict(mfit,x))
end

"""
Properties of Spatial Frequency Tuning

    - Prefered Spatial Frequency
    - Half Width at Half Peak-to-Trough
    - Freq Passing Type {A:All Pass, H:High Pass, L:Low Pass, B:Band Pass}
    - Bandwidth ``log2(H_cut/L_cut)``
    - Passwidth at Half Peak-to-Trough constrained by low freq lim and high freq lim

1. x: sf in cycle/degree
2. y: responses
"""
function sftuningfeature(x,y;low=minimum(x),high=maximum(x))
    maxi = argmax(y)
    mini = argmin(y)
    maxr = y[maxi]
    minr = y[mini]
    px = x[maxi]

    hw = halfwidth(y,start=maxi,v=(maxr+minr)/2,circ=false,x=x)
    pt = all(isinf.(hw)) ? 'A' : isinf(hw[1]) ? 'L' : isinf(hw[2]) ? 'H' : 'B'
    bw = log2((px+hw[2])/(px-hw[1]))
    pw = pt == 'A' ? high-low : pt == 'L' ? px-low+hw[2] : pt == 'H' ? high-px+hw[1] : sum(hw)

    (;psf=px,sfhw=hw,sftype=pt,sfbw=bw,sfpw=pw)
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
                mfit = fitmodel(:gvm,θ,fr)
                fit = (circtuningfeature(mfit,od=π,fn=:d)...,gvm=mfit)
            catch
            end
            # fit von Mises for orientation
            try
                mfit = fitmodel(:vmn2, θ.-0.5π,fr)
                fit =(fit...,circtuningfeature(mfit, od=0.5π,fn=:o)...,vmn2=mfit)
            catch
            end
        end

        return (;dm,od,dcv,om,oo,ocv,fit)
    elseif factor == :SpatialFreq
        osf = 2^(sum(fr.*log2.(fl))/sum(fr)) # weighted average as optimal sf
        maxi = argmax(fr)
        maxsf = fl[maxi]
        maxr = fr[maxi]

        fit = ()
        if isfit
            # fit difference of gaussians for SpatialFreq
            try
                mfit = fitmodel(:dog,fl,fr)
                fit = (sftuningfeature(mfit)...,dog=mfit)
            catch
            end
        end

        return (;osf,maxsf,maxr,fit)
    elseif factor == :ColorID
        # transform colorId to hue angle
        ucid = sort(unique(fl))
        hstep = 2pi/length(ucid)
        ha = map(l->hstep*(findfirst(c->c==l,ucid)-1),fl)
        # for hue direction
        hm = circmean(ha,fr)
        oh = mod(rad2deg(angle(hm)),360)
        hcv = circvar(ha,fr,hstep)
        maxi = argmax(fr)
        maxh = rad2deg(ha[maxi])
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
    isnothing(norm) || (r/=norm)
    isnothing(whiten) || (r=whiten*r)

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
function unitdensity(pos;w=ones(length(pos)),lim=extrema(pos),bw=0.01(lim[2]-lim[1]),step=bw/2,r=nothing,wfun=sum)
    hbw = bw/2
    y = lim[1]:step:lim[2]
    n = [wfun(w[i-hbw .<=pos.< i+hbw]) for i in y]
    if !isnothing(r)
        n = n/(bw*π*r^2)
    end
    return (;n,y)
end
function spacepsth(unitpsth,unitposition;w=ones(size(unitposition,1)),lim=extrema(unitposition),bw=0.01(lim[2]-lim[1]),step=bw/2)
    hbw = bw/2
    x = unitpsth[1].x
    y = lim[1]:step:lim[2]

    n = zeros(length(y))
    psth = zeros(length(y),length(x))
    for i in eachindex(y)
        @views j = y[i]-hbw .<=unitposition[:,2].< y[i]+hbw
        n[i] = sum(w[j])
        if n[i] > 0
            @views psth[i,:] = mapreduce(u->u.m,.+,unitpsth[j])/n[i]
        end
    end
    return (;psth,x,y,n)
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
    base = [cc[midi+minbaselag:end];cc[1:midi-minbaselag]]
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

"Check Layer Boundaries"
function checklayer!(ls::Dict;ln=["WM","6","5","56","4Cb","4Ca","4C","4B","4A","4AB","3","2","23","1","Out"])
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
function assignlayer(y,layer::Dict)
    l = missing
    for k in keys(layer)
        if layer[k][1] <= y < layer[k][2]
            l=k;break
        end
    end
    return l
end

"Check circuit consistency and remove duplicates"
function checkcircuit(projs,eunits,iunits,projweights)
    ivu = intersect(eunits,iunits)
    veunits = setdiff(eunits,ivu)
    viunits = setdiff(iunits,ivu)
    ivp = map(p->p[1] in ivu,projs)
    vprojs = projs[.!ivp]
    vprojweights = projweights[.!ivp]
    ui = indexin(unique(vprojs),vprojs)
    return vprojs[ui],veunits,viunits,vprojweights[ui]
end
