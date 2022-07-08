function huecolors(n::Integer=100;alpha=0.8,saturation=1,brightness=1,precolors=[HSVA(0,1,0,alpha)],sufcolors=[])
    hc = [HSVA(360(i-1)/n,saturation,brightness,alpha) for i in 1:n]
    prepend!(hc,precolors)
    append!(hc,sufcolors)
    hc
end

"Generates linearly interpolated `ColorGradient`"
function cgrad(start::T,stop::T;length=100) where T<:Colorant
    cgrad(range(start,stop,length=length))
end

"Generates linearly interpolated `ColorGradient` between key colors."
function cgrad(x::T...;length=100) where T<:Colorant
    cgrad(range(x...,length=length))
end

"""
Generates linearly interpolated colors between key colors inclusively.

return an `Array` of colors.
"""
function range(x::T...;length=100) where T<:Colorant
    nk = Base.length(x)
    sn = round(Int,length/(nk-1))
    segs = map(i->range(x[i],x[i+1],length=sn),1:nk-1)
    cs = push!(mapreduce(i->i[1:end-1],append!,segs),segs[end][end])
end

"Predefined ColorMaps"
ColorMaps=Dict()

function plotcolormap(cm;title="",xlabel="",ylabel="",markersize=12,shape=:circle,frame=:none,theta=45)
    plot()
    plotcolormap!(cm;title,xlabel,ylabel,markersize,shape,frame,theta)
end
function plotcolormap!(cm;title="",xlabel="",ylabel="",markersize=12,shape=:circle,frame=:none,theta=45)
    if shape == :circle
        yx = sincos.(2π*range(0,1,length=length(cm)))
        x = map(i->i[2],yx)
        y = map(i->i[1],yx)
        marker = :circle
        xlims=ylims=[-1.2,1.2]
        color=cm
    elseif shape == :hline
        x = range(-1,1,length=length(cm))
        y = zeros(length(cm))
        marker = :rect
        ylims=[-0.5,0.5]
        xlims = [-1.2,1.2]
        color=cm
    elseif shape == :vline
        y = range(-1,1,length=length(cm))
        x = zeros(length(cm))
        marker = :rect
        xlims=[-0.5,0.5]
        ylims = [-1.2,1.2]
        color=cm
    elseif shape == :diag
        s,c = sincosd(theta)
        x = range(-1,1,length=length(cm))
        y = s*x
        x = c*x
        marker = :circle
        xlims=ylims=[-1.2,1.2]
        color=cm
    elseif shape == :sin
        x = range(-1,1,length=length(cm))
        y = sin.(π*(x.+1))
        color=cgrad(cm)[(y.+1)/2]
        xlims=ylims=[-1.2,1.2]
        marker = :circle
    end
    scatter!(x,y;aspectratio=:equal,color,markersize,frame,marker,xlims,ylims,msw=0,leg=false,xlabel,ylabel,title)
end
