huecolors(n::Int;alpha=0.8,saturation=1,brightness=1)=[HSVA(((i-1)/n)*360,saturation,brightness,alpha) for i=1:n]

# function minmaxcolormap(cname,min,max;isreverse=false)
#     c=colormap(cname,101,mid=max/(abs(min)+abs(max)))
#     if isreverse
#         c=reverse(c)
#     end
#     ColorGradient(c,0:0.01:1)
# end
# function minmaxcolorgradient(minc,maxc;n=100)
#     d = maxc-minc
#     r = range(0,1,length=n)
#     ColorGradient(map(i->minc+i*d,r),r)
# end
# function mapcolor(data,cg::ColorGradient)
#     minv,maxv = extrema(data)
#     r=maxv-minv
#     map(i->RGBA(cg[(i-minv)/r]),data)
# end

function unitcolors(uids=[];n=5,alpha=0.8,saturation=1,brightness=1)
    uc = huecolors(n,alpha=alpha,saturation=saturation,brightness=brightness)
    insert!(uc,1,HSVA(0,0,0,alpha))
    if !isempty(uids)
        uc=uc[sort(unique(uids)).+1]
    end
    return uc
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

function plotcolormap(cm;title="",xlabel="",ylabel="",markersize=12,shape=:circle)
    if shape == :circle
        yx = sincos.(2π*range(0,1,length=length(cm)))
        x = map(i->i[2],yx)
        y = map(i->i[1],yx)
        marker = :circle
        ylims=[-1.2,1.2]
    else
        x = range(0,1,length=length(cm))
        y = ones(length(cm))
        marker = :rect
        ylims=[0.5,1.5]
    end
    scatter(x,y,aspectratio=:equal,color=cm,markersize=markersize,marker=marker,ylims=ylims,
    markerstrokewidth=0,legend=false,xlabel=xlabel,ylabel=ylabel,title=title)
end
