using LinearAlgebra,Distributions,DataFrames,StatsBase,GLM,LsqFit,Optim,BlackBoxOptim,HypothesisTests,Colors,Images,StatsModels,CategoricalArrays,
ImageFiltering,SpecialFunctions,DSP,HCubature,Combinatorics,DataStructures,StatsFuns,Trapz,CircStats,ImageSegmentation,ProgressMeter,#PyCall,
Dierckx,BandedMatrices,OffsetArrays#,ANOVA
import Base: vec,range
import StatsBase: predict

include("NeuroDataType.jl")
include("Stats.jl")
include("Spike.jl")
include("Jitter.jl")
include("LFP.jl")
include("Image.jl")
include("Condition.jl")
include("2P.jl")

vec(x::RGBA)=[x.r,x.g,x.b,x.alpha]
anscombe(x) = 2*sqrt(x+(3/8))

"`Ellipse` function"
function ellipsef(α;a=2,b=1,θ=π/4,μ₁=0,μ₂=0)
    sinθ,cosθ = sincos(θ)
    sinα,cosα = sincos(α)
    x′ = a*cosα
    y′ = b*sinα
    cosθ*x′ - sinθ*y′ + μ₁, sinθ*x′ + cosθ*y′ + μ₂
end
"`Gaussian` function"
gaussianf(x;a=1,μ=0,σ=1) = a*exp(-0.5((x-μ)/σ)^2)
"2D `Gaussian` function"
function gaussianf(x,y;a=1,μ₁=0,σ₁=1,μ₂=0,σ₂=1,θ=0)
    sinθ,cosθ = sincos(θ)
    x₀ = x-μ₁
    y₀ = y-μ₂
    x′ = cosθ * x₀ + sinθ * y₀
    y′ = cosθ * y₀ - sinθ * x₀
    a*exp(-0.5((x′/σ₁)^2 + (y′/σ₂)^2))
end
"`Gaussian` contour"
gaussiancontour(α;fσ=2.5,μ₁=0,σ₁=1,μ₂=0,σ₂=1,θ=0) = ellipsef(α;a=fσ*σ₁,b=fσ*σ₂,θ,μ₁,μ₂)
"2D `Gaussian` envelope"
function gaussianenvelope(x,y;fσ=2.5,μ₁=0,σ₁=1,μ₂=0,σ₂=1,θ=0)
    sinθ,cosθ = sincos(θ)
    x₀ = x-μ₁
    y₀ = y-μ₂
    x′ = cosθ * x₀ + sinθ * y₀
    y′ = cosθ * y₀ - sinθ * x₀
    (x′/(fσ*σ₁))^2 + (y′/(fσ*σ₂))^2 - 1 # ellipse function around tail of gaussian
end
"2D `Gaussian` envelope binary mask"
gaussianenvelopemask(x,y;fσ=2.5,μ₁=0,σ₁=1,μ₂=0,σ₂=1,θ=0)=gaussianenvelope(x,y;fσ,μ₁,σ₁,μ₂,σ₂,θ) <= 0 ? true : false

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
gvmf(α;β=1,μ₁=0,κ₁=1,μ₂=0,κ₂=1) = β*exp(κ₁*cos(α-μ₁) + κ₂*cos(2(α-μ₂)))

"`Difference of Gaussians` function"
dogf(x;aₑ=2,μₑ=0,σₑ=1,aᵢ=1,μᵢ=0,σᵢ=2) = aₑ*exp(-0.5((x-μₑ)/σₑ)^2) - aᵢ*exp(-0.5((x-μᵢ)/σᵢ)^2)
"2D `Difference of Gaussians` function"
function dogf(x,y;aₑ=2,μₑ₁=0,σₑ₁=1,μₑ₂=0,σₑ₂=1,θₑ=0,aᵢ=1,μᵢ₁=0,σᵢ₁=2,μᵢ₂=0,σᵢ₂=2,θᵢ=0)
    sinθₑ,cosθₑ = sincos(θₑ)
    xₑ₀ = x-μₑ₁
    yₑ₀ = y-μₑ₂
    xₑ′ = cosθₑ * xₑ₀ + sinθₑ * yₑ₀
    yₑ′ = cosθₑ * yₑ₀ - sinθₑ * xₑ₀
    sinθᵢ,cosθᵢ = sincos(θᵢ)
    xᵢ₀ = x-μᵢ₁
    yᵢ₀ = y-μᵢ₂
    xᵢ′ = cosθᵢ * xᵢ₀ + sinθᵢ * yᵢ₀
    yᵢ′ = cosθᵢ * yᵢ₀ - sinθᵢ * xᵢ₀
    aₑ*exp(-0.5((xₑ′/σₑ₁)^2 + (yₑ′/σₑ₂)^2)) - aᵢ*exp(-0.5((xᵢ′/σᵢ₁)^2 + (yᵢ′/σᵢ₂)^2))
end

"""
`sin` grating function

- μ: x offset
- f: Frequency in cycle/unit_x
- phase: Phase of a cycle in [0, 1] scale
"""
gratingf(x;μ=0, f=1, phase=0) = sin(2π * (f * (x-μ) + phase))

"""
2D `sin` grating function

- μ₁: x offset
- μ₂: y offset
- θ: Orientation in radius, 0 is -, increase counter-clock wise
- f: Frequency in cycle/unit_distance orthogonal to orientation
- phase: Phase of a cycle in [0, 1] scale
"""
function gratingf(x,y;μ₁=0,μ₂=0,θ=0,f=1,phase=0)
    sinθ,cosθ = sincos(θ)
    y′ = cosθ * (y-μ₂) - sinθ * (x-μ₁)
    sin(2π * (f * y′ + phase))
end

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
gaborf(x;a=1,μ=0,σ=1,f=1,phase=0) = a*exp(-0.5((x-μ)/σ)^2)*sin(2π*(f*(x-μ)+phase))
"2D `Gabor` function"
function gaborf(x,y;a=1,μ₁=0,σ₁=1,μ₂=0,σ₂=1,θ=0,f=1,phase=0)
    sinθ,cosθ = sincos(θ)
    x₀ = x-μ₁
    y₀ = y-μ₂
    x′ = cosθ * x₀ + sinθ * y₀
    y′ = cosθ * y₀ - sinθ * x₀
    a*exp(-0.5((x′/σ₁)^2 + (y′/σ₂)^2)) * sin(2π*(f * y′ + phase))
end

"`gabor` contour"
gaborcontour(α;fσ=2.5,μ₁=0,σ₁=1,μ₂=0,σ₂=1,θ=0) = ellipsef(α;a=fσ*σ₁,b=fσ*σ₂,θ,μ₁,μ₂)
"2D `gabor` envelope"
gaborenvelope(x,y;fσ=2.5,μ₁=0,σ₁=1,μ₂=0,σ₂=1,θ=0) = gaussianenvelope(x,y;fσ,μ₁,σ₁,μ₂,σ₂,θ)
"2D `gabor` envelope binary mask"
gaborenvelopemask(x,y;fσ=2.5,μ₁=0,σ₁=1,μ₂=0,σ₂=1,θ=0) = gaussianenvelopemask(x,y;fσ,μ₁,σ₁,μ₂,σ₂,θ)

"concentric circular dog contour"
function dogcontour(α;fσ=2.5,μ₁=0,σₑ₁=1,μ₂=0,rσᵢₑ=2)
    a=fσ*max(σₑ₁,rσᵢₑ*σₑ₁)
    ellipsef(α;a,b=a,θ=0,μ₁,μ₂)
end
"2D concentric circular dog envelope"
function dogenvelope(x,y;fσ=2.5,μ₁=0,σₑ₁=1,μ₂=0,rσᵢₑ=2)
    σ₁ = max(σₑ₁,rσᵢₑ*σₑ₁)
    gaussianenvelope(x,y;fσ,μ₁,σ₁,μ₂,σ₂=σ₁,θ=0)
end
"2D concentric circular dog envelope binary mask"
function dogenvelopemask(x,y;fσ=2.5,μ₁=0,σₑ₁=1,μ₂=0,rσᵢₑ=2)
    σ₁ = max(σₑ₁,rσᵢₑ*σₑ₁)
    gaussianenvelopemask(x,y;fσ,μ₁,σ₁,μ₂,σ₂=σ₁,θ=0)
end

"concentric orientated elliptical dog contour"
function edogcontour(α;fσ=2.5,μ₁=0,σₑ₁=1,rσ₂₁=1,μ₂=0,rσᵢₑ=2,θ=0)
    σ₁ = max(σₑ₁,rσᵢₑ*σₑ₁)
    σ₂ = rσ₂₁*σ₁
    ellipsef(α;a=fσ*σ₁,b=fσ*σ₂,θ,μ₁,μ₂)
end
"2D concentric orientated elliptical dog envelope"
function edogenvelope(x,y;fσ=2.5,μ₁=0,σₑ₁=1,rσ₂₁=1,μ₂=0,rσᵢₑ=2,θ=0)
    σ₁ = max(σₑ₁,rσᵢₑ*σₑ₁)
    σ₂ = rσ₂₁*σ₁
    gaussianenvelope(x,y;fσ,μ₁,σ₁,μ₂,σ₂,θ)
end
"2D concentric orientated elliptical dog envelope binary mask"
function edogenvelopemask(x,y;fσ=2.5,μ₁=0,σₑ₁=1,rσ₂₁=1,μ₂=0,rσᵢₑ=2,θ=0)
    σ₁ = max(σₑ₁,rσᵢₑ*σₑ₁)
    σ₂ = rσ₂₁*σ₁
    gaussianenvelopemask(x,y;fσ,μ₁,σ₁,μ₂,σ₂,θ)
end

"Fit 1D model to data"
function fitmodel(model,x,y)
    lb,ub = extrema(y)
    bm = (lb+ub)/2
    br = (ub-lb)/2
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

        ub=[1.8ab,   prevfloat(float(π)),   200]
        lb=[0.2ab,            0,              0]
        p0=[ab,               0,              1]
    elseif model == :gvm
        # fun = (x,p) -> gvmf.(x,β=p[1],μ₁=p[2],κ₁=p[3],μ₂=p[4],κ₂=p[5])
        fun = gvmff
        ofun = (p;x=x,y=y) -> sum((y.-fun(x,p)).^2)

        ub=[1.8ab,   prevfloat(float(2π)),   40,    prevfloat(float(π)),     40]
        lb=[nextfloat(0.0),   0,              0,             0,               0]
        p0=[ab,               0,              1,             0,               1]
    elseif model == :dog
        fun = (x,p) -> dogf.(x,aₑ=p[1],μₑ=p[2],σₑ=p[3],aᵢ=p[4],μᵢ=p[5],σᵢ=p[6])
        ofun = (p;x=x,y=y) -> sum((y.-fun(x,p)).^2)

        ub=[1.8ab,   10xab,            10xab,    1.8ab,     10xab,              10xab]
        lb=[0,      -10xab,   nextfloat(0.0),        0,    -10xab,     nextfloat(0.0)]
        p0=[ab,        0,                  1,       ab,         0,                  1]
    elseif model == :sfdog
        # fun = (x,p) -> dogf.(x,aₑ=p[1],μₑ=p[2],σₑ=p[3],aᵢ=p[4],μᵢ=p[5],σᵢ=p[6]) .+ p[7]
        fun = sfdogff
        ofun = (p;x=x,y=y) -> sum((y.-fun(x,p)).^2)

        ub=[1.5ab,    10,             10,          1.5ab,       10,                10,         bm+br/3]
        lb=[0,         0,      nextfloat(0.0),       0,          0,        nextfloat(0.0),        0]
        p0=[ab,     0.5xab,         0.5xab,          0,        0.5xab,           0.5xab,       bm-br/3]
    elseif model == :sfgaussian
        # fun = (x,p) -> gaussianf.(log2.(x),a=p[1],μ=p[2],σ=p[3]) .+ p[4]
        fun = sfgaussianff
        ofun = (p;x=x,y=y) -> sum((y.-fun(x,p)).^2)

        # limit gaussian center ~[-8 10], sigma ~[1.4 8]
        ub=[1.5ab,     3.3,          3,       bm]
        lb=[0,          -3,        0.5,        0]
        p0=[ab,          0,          1,        0]
    end
    if !ismissing(fun)
        # ofit = optimize(ofun,lb,ub,p0,SAMIN(rt=0.92),Optim.Options(iterations=220000))
        # param=ofit.minimizer
        ofit = bboptimize(ofun,p0;SearchRange=collect(zip(lb,ub)),Method=:adaptive_de_rand_1_bin_radiuslimited,MaxSteps=100000)
        param = best_candidate(ofit)

        rlt = (;model,fun,param, goodnessoffit(y,fun(x,param),k=length(param))...)
    end
    return rlt
end
gvmff(x,p) = gvmf.(x,β=p[1],μ₁=p[2],κ₁=p[3],μ₂=p[4],κ₂=p[5])
sfdogff(x,p) = dogf.(x,aₑ=p[1],μₑ=p[2],σₑ=p[3],aᵢ=p[4],μᵢ=p[5],σᵢ=p[6]) .+ p[7]
sfgaussianff(x,p) = gaussianf.(log2.(x),a=p[1],μ=p[2],σ=p[3]) .+ p[4]


"Fit 2D model to image"
function fitmodel2(model,data::AbstractMatrix,ppu;w=0.5,fun=rms)
    rspx = (size(data).-1)./2 # data should have odd pixels
    radii = rspx./ppu
    # standard 2D coordinates(rightward:x, upward:y) vector in unit
    x = (vcat(map(i->[i[2] -i[1]],CartesianIndices(data))...) .- [(rspx[2]+1) -(rspx[1]+1)]) / ppu
    y = vec(data)

    # prepare bounds and inital param for model solution
    roi = peakroi(localcontrast(data,round(Int,w*ppu);fun)) # locate meaningful data region
    alb,aub = abs.(extrema(data[roi.i]))
    ab = max(alb,aub)
    r = roi.radius / ppu
    c = [roi.center[2] - (rspx[2]+1), -roi.center[1] + (rspx[1]+1)] / ppu

    rlt = fun = missing
    if model == :dog # concentric circular dog
        fun = dogff
        mfun = dogfmf
        cfun = dogfcf
        ofun = (p;x=x,y=y) -> @views sum((y.-fun(x[:,1],x[:,2],p)).^2)

        if aub >= alb
            ai = 5alb
            ae = aub + ai
        else
            ae = 5aub
            ai = alb + ae
        end
        ub=[5ae,    0.5r+c[1],    0.9r,    0.5r+c[2],     5ai,    4]
        lb=[0,     -0.5r+c[1],    0.1r,   -0.5r+c[2],     0,      0.25]
        p0=[ae,      c[1],        0.3r,     c[2],         ai,     1]
    elseif model == :edog # concentric orientated elliptical dog
        fun = edogff
        mfun = edogfmf
        cfun = edogfcf
        ofun = (p;x=x,y=y) -> @views sum((y.-fun(x[:,1],x[:,2],p)).^2)

        if aub >= alb
            ai = 5alb
            ae = aub + ai
        else
            ae = 5aub
            ai = alb + ae
        end
        ub=[5ae,    0.5r+c[1],    0.9r,    0.5r+c[2],       1,       prevfloat(float(π)),    5ai,     4]
        lb=[0,     -0.5r+c[1],    0.1r,   -0.5r+c[2],      0.5,              0,               0,      0.25]
        p0=[ae,      c[1],        0.3r,     c[2],           1,               0,               ai,     1]
    elseif model == :gabor
        fun = gaborff
        mfun = gaborfmf
        cfun = gaborfcf
        ofun = (p;x=x,y=y) -> @views sum((y.-fun(x[:,1],x[:,2],p)).^2)

        ori,sf = f1orisf(powerspectrum2(data,ppu)...)
        ub=[5ab,   0.5r+c[1],   0.9r,    0.5r+c[2],    6.0,     prevfloat(float(π)),     12,            prevfloat(1.0)]
        lb=[0,    -0.5r+c[1],   0.1r,   -0.5r+c[2],    0.2,             0,              0.05,                  0]
        p0=[ab,    c[1],        0.3r,     c[2],          1,            ori,     clamp(sf,lb[7],ub[7]),       0.5]
    end
    if !ismissing(fun)
        # ofit = optimize(ofun,lb,ub,p0,SAMIN(rt=0.92),Optim.Options(iterations=220000))
        # param = ofit.minimizer
        ofit = bboptimize(ofun,p0;SearchRange=collect(zip(lb,ub)),Method=:adaptive_de_rand_1_bin_radiuslimited,MaxSteps=200000)
        param = best_candidate(ofit)

        @views rlt = (;model,fun,mfun,cfun,param,radii,ppu, goodnessoffit(y,fun(x[:,1],x[:,2],param),k=length(param))...)
    end
    return rlt
end

dogff(x,y,p) = dogf.(x,y,aₑ=p[1],μₑ₁=p[2],σₑ₁=p[3],μₑ₂=p[4],σₑ₂=p[3],θₑ=0,aᵢ=p[5],μᵢ₁=p[2],σᵢ₁=p[6]*p[3],μᵢ₂=p[4],σᵢ₂=p[6]*p[3],θᵢ=0)
dogfmf(x,y,p) = dogenvelopemask.(x,y;fσ=2.5,μ₁=p[2],μ₂=p[4],σₑ₁=p[3],rσᵢₑ=p[6])
dogfcf(x,p) = dogcontour.(x;fσ=2.5,μ₁=p[2],μ₂=p[4],σₑ₁=p[3],rσᵢₑ=p[6])

edogff(x,y,p) = dogf.(x,y,aₑ=p[1],μₑ₁=p[2],σₑ₁=p[3],μₑ₂=p[4],σₑ₂=p[5]*p[3],θₑ=p[6],aᵢ=p[7],μᵢ₁=p[2],σᵢ₁=p[8]*p[3],μᵢ₂=p[4],σᵢ₂=p[5]*p[8]*p[3],θᵢ=p[6])
edogfmf(x,y,p) = edogenvelopemask.(x,y;fσ=2.5,μ₁=p[2],μ₂=p[4],σₑ₁=p[3],rσ₂₁=p[5],rσᵢₑ=p[8],θ=p[6])
edogfcf(x,p) = edogcontour.(x;fσ=2.5,μ₁=p[2],μ₂=p[4],σₑ₁=p[3],rσ₂₁=p[5],rσᵢₑ=p[8],θ=p[6])

gaborff(x,y,p) = gaborf.(x,y;a=p[1],μ₁=p[2],σ₁=p[3],μ₂=p[4],σ₂=p[5]*p[3],θ=p[6],f=p[7],phase=p[8])
gaborfmf(x,y,p) = gaborenvelopemask.(x,y;fσ=2.5,μ₁=p[2],σ₁=p[3],μ₂=p[4],σ₂=p[5]*p[3],θ=p[6])
gaborfcf(x,p) = gaborcontour.(x;fσ=2.5,μ₁=p[2],σ₁=p[3],μ₂=p[4],σ₂=p[5]*p[3],θ=p[6])


predict(fit,x) = fit.fun(x,fit.param)
function predict(fit,x,y;xygrid=true,yflip=false)
    if xygrid
        z = [fit.fun(i,j,fit.param) for j in y, i in x]
        yflip && reverse!(z,dims=1)
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
- k: number of model parameters
- df: degree of freedom(n - k)
"""
function goodnessoffit(y,ŷ;n = length(y),e = y .- ŷ,k = missing,df = n - k,corrected=false)
    r = cor(y,ŷ)
    ae = abs.(e)
    e2 = e.^2
    ssᵣ = sum(e2)
    mae = mean(ae)
    mse = ssᵣ/n
    rmse = sqrt(mse)
    yₘ = y .- mean(y)
    ssₜ = sum(yₘ.^2)
    rae = sum(ae)/sum(abs.(yₘ))
    fvu = ssᵣ/ssₜ
    rse = sqrt(fvu)
    r2 = 1 - fvu
    varᵣ = ssᵣ/df
    s = varᵣ < 0 ? missing : sqrt(varᵣ)
    adjr2 = 1 - fvu*(n-1)/df
    # under iid normal error distribution
    ll = corrected ? n*log(varᵣ) : n*log(mse)
    aic = ll + 2k
    bic = ll + k*log(n)
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

circtuningfeature(mfit;od=[π,0.5π],fn=:a,x = 0:0.002:2π) = circtuningfeature(x,predict(mfit,x);od,fn) # 0.002rad ≈ 0.11deg
"""
Properties of Circular Tuning

    - Prefered Angle with Peak Response
    - Selectivity Index
        - version 1: (ResponsePeak - ResponseOpposing)/ResponsePeak
        - version 2: (ResponsePeak - ResponseOpposing)/(ResponsePeak + ResponseOpposing)
    - Half Width at Half Peak-to-Trough

1. x: angles in radius
2. y: responses
- od: opposing angle distance to prefered angle, e.g. π for DSI, 0.5π for OSI
- fn: factor name
"""
function circtuningfeature(x,y;od=[π,0.5π],fn=:a)
    maxr,maxi = findmax(y)
    minr,mini = findmin(y)
    maxx = x[maxi]
    ox = maxx.+od
    _,oi = findclosestangle(x,ox)
    or = y[oi]

    si1 = 1 .- or./maxr
    si2 = (maxr .- or)./(maxr .+ or)
    hw = halfwidth(y,start=maxi,v=(maxr+minr)/2,circ=true,x=x)

    (;Symbol(:p,fn)=>rad2deg(mod2pi(maxx)),Symbol(fn,:hw)=>rad2deg.(hw),Symbol(fn,:si1)=>si1,Symbol(fn,:si2)=>si2,Symbol(fn,:od)=>od)
end

sftuningfeature(mfit;x = 0:0.003:10) = sftuningfeature(x,predict(mfit,x))
"""
Properties of Spatial Frequency Tuning

    - Prefered Spatial Frequency with Peak Response
    - Half Width at Half Peak-to-Trough
    - Freq Passing Type {A:All Pass, H:High Pass, L:Low Pass, B:Band Pass}
    - Bandwidth ``log2(H_cut/L_cut)``
    - Passwidth at Half Peak-to-Trough constrained by low freq lim and high freq lim

1. x: sf in cycle/degree
2. y: responses
"""
function sftuningfeature(x,y;low=minimum(x),high=maximum(x))
    maxr,maxi = findmax(y)
    minr,mini = findmin(y)
    maxx = x[maxi]

    hw = halfwidth(y,start=maxi,v=(maxr+minr)/2,circ=false,x=x)
    pt = all(isinf.(hw)) ? 'A' : isinf(hw[1]) ? 'L' : isinf(hw[2]) ? 'H' : 'B'
    bw = log2((maxx+hw[2])/(maxx-hw[1]))
    pw = pt == 'A' ? high-low : pt == 'L' ? maxx-low+hw[2] : pt == 'H' ? high-maxx+hw[1] : sum(hw)

    (;psf=maxx,sfhw=hw,sftype=pt,sfbw=bw,sfpw=pw)
end

"""
Tuning properties of factor response

1. fl: factor levels
2. fr: factor responses for each level

    Angle, Orientation and Direction follow the same convention such that 0 is -/→, then increase counter-clock wise.

    For cases where Orientation and Direction are interlocked(drifting grating):
        - when Orientation is -(0), then Direction is ↑(90)
        - when Direction is →(0), then Orientation is |(-90)
"""
function factorresponsefeature(fl,fr;fm=mean.(fr),factor=:Ori,isfit::Bool=true)
    i = .!ismissing.(fr)
    fl = fl[i];fr=fr[i];fm=fm[i]
    # ls = mapreduce((l,r)->fill(l,length(r)),append!,fl,fr)
    # rs = mapreduce(deepcopy,append!,fr)

    if factor in [:Ori,:Ori_Final]
        # θ = mod2pi.(deg2rad.(ls))
        α = mod2pi.(deg2rad.(fl))
        d = mean(diff(sort(unique(α)))) # angle spacing
        # for orientation
        oα = mod.(α,π)
        ocv, = circ_var(2oα,w=fm,d=2d)
        ocm, = circ_mean(2oα,w=fm)
        ocm = rad2deg(mod2pi(ocm)/2)
        oup, = circ_otest(2oα,w=fm) # Omnibus test for orientation non-uniformity
        # for direction
        dcv, = circ_var(α;w=fm,d)
        dcm, = circ_mean(α,w=fm)
        dcm = rad2deg(mod2pi(dcm+0.5π))
        dup, = circ_otest(α,w=fm) # Omnibus test for direction non-uniformity
        maxr,maxi = findmax(fm)
        maxl = fl[maxi]

        fit = ()
        if isfit
            try
                mfit = fitmodel(:gvm,α,fm) # fit Generalized von Mises
                fit = (;circtuningfeature(mfit,od=[π,0.5π],fn=:o)...,mfit)
            catch
                display.(stacktrace(catch_backtrace()))
            end
        end

        return (;oup,ocv,ocm,dup,dcv,dcm,max=maxl=>maxr,fit)
    elseif factor == :Dir
        θ = deg2rad.(fl)
        d = mean(diff(sort(unique(θ)))) # angle spacing
        # for orientation
        oθ = mod.(θ.-0.5π,π)
        om, = circ_mean(2oθ,w=fr)
        oo = rad2deg(mod(angle(om),2π)/2)
        ocv = circ_var(2oθ,w=fr,d=2d)
        # for direction
        dm, = circ_mean(θ,w=fr)
        od = rad2deg(mod(angle(dm),2π))
        dcv = circ_var(θ;w=fr,d)
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
        up = PyOnewayANOVA.anova_oneway(fr,use_var="unequal").pvalue
        msf = 2^(sum(fm.*log2.(fl))/sum(fm)) # weighted average
        maxr,maxi = findmax(fm)
        maxl = fl[maxi]

        fit = ()
        if isfit
            try
                # mfit = fitmodel(:sfdog,fl,fm) # fit Difference of Gaussians
                mfit = fitmodel(:sfgaussian,fl,fm) # fit Gaussian of logarithm sf
                fit = (;sftuningfeature(mfit,x = range(extrema(fl)...,step=0.003))...,mfit)
            catch
                display.(stacktrace(catch_backtrace()))
            end
        end

        return (;up,msf,max=maxl=>maxr,fit)
    elseif factor == :ColorID
        # transform colorId to hue angle
        ucid = sort(unique(fl))
        hstep = 2pi/length(ucid)
        ha = map(l->hstep*(findfirst(c->c==l,ucid)-1),fl)
        # for hue direction
        hm, = circ_mean(ha,w=fr)
        oh = mod(rad2deg(angle(hm)),360)
        hcv = circ_var(ha,w=fr,d=hstep)
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
    elseif factor in [:HueAngle,:Angle]
        # θ = mod2pi.(deg2rad.(ls))
        α = mod2pi.(deg2rad.(fl))
        d = mean(diff(sort(unique(α)))) # angle spacing
        # for axis
        aα = mod.(α,π)
        acv, = circ_var(2aα,w=fm,d=2d)
        acm, = circ_mean(2aα,w=fm)
        acm = rad2deg(mod2pi(acm)/2)
        aup, = circ_otest(2aα,w=fm) # Omnibus test for axis non-uniformity
        # for angle
        cv, = circ_var(α;w=fm,d)
        cm, = circ_mean(α,w=fm)
        cm = rad2deg(mod2pi(cm))
        up, = circ_otest(α,w=fm) # Omnibus test for angle non-uniformity
        maxr,maxi = findmax(fm)
        maxl = fl[maxi]

        fit = ()
        if isfit
            try
                mfit = fitmodel(:gvm,α,fm) # fit Generalized von Mises
                fit = (;circtuningfeature(mfit,od=[π,0.5π],fn=:a)...,mfit)
            catch
                display.(stacktrace(catch_backtrace()))
            end
        end

        return (;aup,acv,acm,up,cv,cm,max=maxl=>maxr,fit)
    else
        return ()
    end
end

"""
Spike Triggered Average of Images based on linear model of response,

```math
yᵢ = KᵀXᵢ + ϵᵢ, ϵᵢ ~ N(0, σ²)
```

then the maximum likelihood estimator and least square estimator have the same analytical form.

```math
K̂ₘₗ = argmax[P(Y|X,K)] = (XᵀX)⁻¹XᵀY = K̂ₗₛ
```

1. x: Matrix where each row is one image
2. y: Vector of image response

- norm: normalization factor, default no normalization.
        it could be ``sum(y)`` if y is number of spike or spike rate, then STA would be spiking probability.
- whiten: whiten factor, default no whiten.
        it could be ``(XᵀX)⁻¹`` or inverse of covariance matrix, that decorrelate STA.
"""
function sta(x::AbstractMatrix,y::AbstractVector;norm=nothing,whiten=nothing)
    r = x'*y
    isnothing(norm) || (r=r/norm)
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
function unitdensity(pos;w=ones(length(pos)),spacerange=extrema(pos),bw=0.01(last(spacerange)-first(spacerange)),
                    step=bw/2,r=nothing,wfun=sum,s=nothing)
    hbw = bw/2
    y = first(spacerange):step:last(spacerange)
    n = [wfun(w[i-hbw .<=pos.< i+hbw]) for i in y]
    if !isnothing(r)
        n = n/(bw*π*r^2)
    end
    i = isnan.(n) .| isinf.(n)
    if any(i)
        n = Spline1D(y[.!i],n[.!i],k=3,bc="extrapolate")(y)
    end
    if !isnothing(s)
        g = gaussianf.(-5:5,σ=s)
        g ./= sum(g)
        n = imfilter(n,centered(g),ImageFiltering.Fill(0))
    end
    return (;n,y)
end
function spacepsth(unitpsth,unitposition;dims=1,spacerange=extrema(unitposition[:,dims]),
                    bw=0.01(last(spacerange)-first(spacerange)),step=bw/2)
    hbw = bw/2
    y = first(spacerange):step:last(spacerange)

    n = zeros(length(y))
    psth = zeros(length(y),length(unitpsth[1]))
    for i in eachindex(y)
        @views j = y[i]-hbw .<=unitposition[:,dims].< y[i]+hbw
        n[i] = count(j)
        if n[i] > 0
            @views psth[i,:] = reduce(.+,unitpsth[j])/n[i]
        end
    end
    return (;psth,y,n)
end


"""
Normalized(coincidence/spike), condition and trial averaged Cross-Correlogram of binary spike trains (bins x trials).

- Correction:
    - (shuffle=true), "all-way" nonsimultaneous trials shuffle

        (Bair, W., Zohary, E., and Newsome, W.T. (2001). Correlated Firing in Macaque Visual Area MT: Time Scales and Relationship to Behavior. J. Neurosci. 21, 1676-1697.)

    - (shufflejitter=true, l=25), shuffle jittered spikes across trials in consecutive intervals of length l

        (Smith, Matthew A., and Adam Kohn (2008). Spatial and Temporal Scales of Neuronal Correlation in Primary Visual Cortex. J. Neurosci. 28.48 : 12591-12603.)
"""
function correlogram(bst1,bst2;lag=nothing,norm=true,correction=(;shuffle=true),condis=nothing)
    if !isnothing(condis)
        cccg=[];τ=[]
        @views for ci in condis
            ccg,τ = correlogram(bst1[:,ci],bst2[:,ci];lag,norm,correction,condis=nothing)
            push!(cccg,ccg)
        end
        ccg = reduce(.+,cccg)/length(cccg)
        return ccg,τ
    end

    nbin,ntrial = size(bst1)
    lag = floor(Int,isnothing(lag) ? min(nbin-1, 10*log10(nbin)) : lag)
    τ = -lag:lag
    cc = Array{Float64}(undef,length(τ),ntrial)
    @views for i in 1:ntrial
        cc[:,i]=crosscov(bst1[:,i],bst2[:,i],τ,demean=false)*nbin
    end
    ccg = dropdims(mean(cc,dims=2),dims=2)

    if haskey(correction,:shuffle) && correction.shuffle
        psth1 = dropdims(mean(bst1,dims=2),dims=2)
        psth2 = dropdims(mean(bst2,dims=2),dims=2)
        s = crosscov(psth1,psth2,τ,demean=false)*nbin
        shiftccg = (ntrial*s .- ccg)/(ntrial-1)
        ccg .-= shiftccg
    elseif haskey(correction,:shufflejitter) && correction.l > 0
        jbst1 = shufflejitter(bst1;l=correction.l)
        jbst2 = shufflejitter(bst2;l=correction.l)
        jitterccg,_ = correlogram(jbst1,jbst2;lag,norm=false,correction=(;),condis=nothing)
        ccg .-= jitterccg
    end

    if norm
        λ1 = mean(mean(bst1,dims=1))
        λ2 = mean(mean(bst2,dims=1))
        gmsr = sqrt(λ1*λ2)
        Θ = nbin.-abs.(τ)
        normfactor = 1 ./ Θ ./ gmsr
        ccg .*= normfactor
    end
    ccg,τ
end
function circuitestimate(unitbinspike;lag=nothing,correction=(;shuffle=true),maxprojlag=5,minbaselag=10,minepoch=4,minspike=4,esdfactor=5,isdfactor=3.5,nosync=true,unitid=[],condis=nothing)
    nunit=length(unitbinspike)
    n,nepoch = size(unitbinspike[1])
    lag = floor(Int,isnothing(lag) ? min(n-1, 10*log10(n)) : lag)
    x = -lag:lag;xn=2lag+1

    isenough = (i,j;minepoch=4,minspike=4) -> begin
        vei = sum(i,dims=1) .>= minspike
        vej = sum(j,dims=1) .>= minspike
        count(vei) >= minepoch && count(vej) >= minepoch
    end

    ccgs=[];ccgis=[];projs=[];eunits=[];iunits=[];projlags=[];projweights=[]
    for (i,j) in combinations(1:nunit,2)
        if isnothing(condis)
            isenough(unitbinspike[i],unitbinspike[j];minepoch,minspike) || continue
            vcondis=nothing
        else
            @views vci = map(ci->isenough(unitbinspike[i][:,ci],unitbinspike[j][:,ci];minepoch,minspike),condis)
            all(.!vci) && continue
            vcondis = condis[vci]
        end

        ccg,_ = correlogram(unitbinspike[i],unitbinspike[j];lag,correction,condis=vcondis)
        ps,es,is,pls,pws = projectionfromcorrelogram(ccg,i,j;maxprojlag,minbaselag,esdfactor,isdfactor,nosync)
        if !isempty(ps)
            push!(ccgs,ccg);push!(ccgis,(i,j))
            append!(projs,ps);append!(eunits,es);append!(iunits,is);append!(projlags,pls);append!(projweights,pws)
        end
    end
    unique!(eunits);unique!(iunits)
    if length(unitid)==nunit
        map!(t->(unitid[t[1]],unitid[t[2]]),ccgis,ccgis)
        map!(t->(unitid[t[1]],unitid[t[2]]),projs,projs)
        map!(t->unitid[t],eunits,eunits)
        map!(t->unitid[t],iunits,iunits)
    end
    return (;ccgs,x,ccgis,projs,eunits,iunits,projlags,projweights)
end
"excitatory and inhibitory projections between two spiking units from cross-correlogram"
function projectionfromcorrelogram(cc,i,j;maxprojlag=5,minbaselag=10,esdfactor=5,isdfactor=5,nosync=true)
    midi = ceil(Int,length(cc)/2)
    base = [cc[midi+minbaselag:end];cc[1:midi-minbaselag]]
    bm,bsd = mean_and_std(base);hl = bm + esdfactor*bsd;ll = bm - isdfactor*bsd
    forwardlags = midi.+(1:maxprojlag)
    backwardlags = midi.-(1:maxprojlag)
    @views fcc = cc[forwardlags]
    @views bcc = cc[backwardlags]
    ps=[];ei=[];ii=[];pls=[];pws=[]

    if nosync && !(ll <= cc[midi] <= hl)
        return (;ps,ei,ii,pls,pws)
    end
    if any(fcc .> hl)
        push!(ps,(i,j));push!(ei,i)
        mi = argmax(fcc)
        push!(pls,mi+1)
        push!(pws,(fcc[mi]-bm)/bsd)
    elseif any(fcc .< ll)
        push!(ps,(i,j));push!(ii,i)
        mi = argmin(fcc)
        push!(pls,mi+1)
        push!(pws,(fcc[mi]-bm)/bsd)
    end
    if any(bcc .> hl)
        push!(ps,(j,i));push!(ei,j)
        mi = argmax(bcc)
        push!(pls,mi+1)
        push!(pws,(bcc[mi]-bm)/bsd)
    elseif any(bcc .< ll)
        push!(ps,(j,i));push!(ii,j)
        mi = argmin(bcc)
        push!(pls,mi+1)
        push!(pws,(bcc[mi]-bm)/bsd)
    end
    return (;ps,ei,ii,pls,pws)
end

"Check circuit consistency, remove conflicts and duplicates"
function checkcircuit(projs,eunits,iunits,projlags,projweights;debug=false)
    ivu = intersect(eunits,iunits)
    if debug && !isempty(ivu)
        @info "E&I Unit: $ivu"
    end
    veunits = setdiff(eunits,ivu)
    viunits = setdiff(iunits,ivu)
    ivp = map(p->p[1] in ivu,projs)
    vprojs = projs[.!ivp]
    vprojlags = projlags[.!ivp]
    vprojweights = projweights[.!ivp]
    ui = indexin(unique(vprojs),vprojs)
    return vprojs[ui],veunits,viunits,vprojlags[ui],vprojweights[ui]
end

"""
Check layer boundaries, make sure no gap and overlap between layers.
Layers are reversely ordered in layer names in `ln`, and start boundary is set to the same as the stop boundary of previous layer.
"""
function checklayer!(ls::Dict;ln=["1", "2", "3A", "2/3A", "3B", "3", "2/3", "4A", "4B", "4A/4B", "4Cα", "4Cβ", "4C", "4", "5A", "5B", "5", "6A", "6B", "6", "5/6", "WM"])
    n = length(ln)
    for i in 1:(n-1)
        if haskey(ls,ln[i])
            for j in (i+1):n
                if haskey(ls,ln[j])
                    ls[ln[i]][begin] = ls[ln[j]][end]
                    break
                end
            end
        end
    end
    return ls
end

"""
Try to locate cell's layer.

1. coordinate of cell on the axis perpendicular to layers
2. layers definition
"""
function assignlayer(y,layer::Dict)
    l = missing
    for k in keys(layer)
        if layer[k][begin] <= y < layer[k][end]
            l=k;break
        end
    end
    return l
end
