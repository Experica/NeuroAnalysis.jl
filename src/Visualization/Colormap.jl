huecolors(n::Int;alpha=0.8,saturation=1,brightness=1)=[HSVA(((i-1)/n)*360,saturation,brightness,alpha) for i=1:n]

function minmaxcolormap(cname,min,max;isreverse=false)
    c=colormap(cname,101,mid=max/(abs(min)+abs(max)))
    if isreverse
        c=reverse(c)
    end
    ColorGradient(c,0:0.01:1)
end
function minmaxcolorgradient(minc,maxc;n=100)
    d = maxc-minc
    r = range(0,1,length=n)
    ColorGradient(map(i->minc+i*d,r),r)
end
function mapcolor(data,cg::ColorGradient)
    minv,maxv = extrema(data)
    r=maxv-minv
    map(i->RGBA(cg[(i-minv)/r]),data)
end

function unitcolors(uids=[];n=5,alpha=0.8,saturation=1,brightness=1)
    uc = huecolors(n,alpha=alpha,saturation=saturation,brightness=brightness)
    insert!(uc,1,HSVA(0,0,0,alpha))
    if !isempty(uids)
        uc=uc[sort(unique(uids)).+1]
    end
    return uc
end

function cgrad(start::T,stop::T;length=100) where T<:Colorant
    cgrad(range(start,stop,length=length))
end

function range(x::T...;length=100) where T<:Colorant
    n = Base.length(x)
    nps = round(Int,(length+n)/(n-1))
    segs = map(i->range(x[i],x[i+1],length=nps),1:n-1)
    cs = push!(mapreduce(i->i[1:end-1],append!,segs),segs[end][end])
end

function plotcolormap(cm::Dict;title="",xlabel="",ylabel="",markersize=12)
    yx = sincos.(2π*cm["values"])
    scatter(map(i->i[2],yx),map(i->i[1],yx),aspectratio=:equal,color=cm["colors"],markersize=markersize,marker=:circle,
    markerstrokewidth=0,legend=false,xlabel=xlabel,ylabel=ylabel,title=title)
end
