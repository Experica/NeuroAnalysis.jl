using Gadfly,Plots,StatPlots,Rsvg

export factorunit,huecolors,unitcolors,plotspiketrain,plotpsth,plotcondresponse,plotsta,savefig

plotlyjs()

factorunit(fs::Vector{Symbol};timeunit=SecondPerUnit)=join(factorunit.(fs,timeunit=timeunit),", ")
function factorunit(f::Symbol;timeunit=SecondPerUnit)
    fu=String(f)
    if occursin("Ori",fu)
        fu="$fu (deg)"
    elseif fu=="Diameter"
        fu="$fu (deg)"
    elseif fu=="SpatialFreq"
        fu="$fu (cycle/deg)"
    elseif fu=="TemporalFreq"
        fu = "$fu (cycle/sec)"
    elseif fu=="Time"
        if timeunit ==1
            fu="$fu (s)"
        elseif timeunit == 0.001
            fu="$fu (ms)"
        end
    elseif fu=="Response"
        fu="$fu (spike/s)"
    end
    return fu
end

huecolors(n::Int;alpha=0.8,saturation=1,brightness=1)=[HSVA(((i-1)/n)*360,saturation,brightness,alpha) for i=1:n]

function unitcolors(uids=[];n=5,alpha=0.8,saturation=1,brightness=1)
    uc = huecolors(n,alpha=alpha,saturation=saturation,brightness=brightness)
    insert!(uc,1,HSVA(0,0,0,alpha))
    if !isempty(uids)
        uc=uc[sort(unique(uids)).+1]
    end
    return uc
end

function plotspiketrain(x,y;group::Vector=[],timeline=[],colors=unitcolors(),title="")
    if isempty(group)
        scatter(x,y,label="SpikeTrain",markershape=:vline,markersize=1,markerstrokecolor=RGBA(0.0,0.1,0.2,0.8),markerstrokewidth = 1)
    else
        scatter(x,y,group=group,markershape=:vline,markersize=1,markerstrokewidth = 1,markerstrokecolor=reshape(colors,1,:))
    end
    vline!(timeline,line=(:grey),label="TimeLine",grid=false,xaxis=(factorunit(:Time)),yaxis=("Trial"),title=(title))
end
function plotspiketrain(sts::RVVector;uids::RVVector=RealVector[],sortvalues=[],timeline=[0],colors=unitcolors(),title="")
    if isempty(uids)
        g=uids;uc=colors
    else
        fuids = flatrvv(uids,sortvalues)[1]
        g=map(i->"U$i",fuids);uc=colors[sort(unique(fuids)).+1]
    end
    plotspiketrain(flatrvv(sts,sortvalues)[1:2]...,group=g,timeline=timeline,colors=uc,title=title)
end

function plotspiketrain1(x::Vector,y::Vector,c::Vector=[];xmin=minimum(x)-10,xmax=maximum(x)+10,xgroup::Vector=[],
    ymin=minimum(y)-1,ymax=maximum(y)+1,timemark=[0],theme=Theme(),
    colorkey="",colorfun=Scale.lab_gradient(colorant"white",colorant"red"),colorminv=[],colormaxv=[])
    xl="Time (ms)";yl="Trial"
    if isempty(c)
        if isempty(xgroup)
            plot(x=x,y=y,xintercept=timemark,theme,Geom.point,Geom.vline(color="gray",size=1pt),
            Coord.Cartesian(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),Guide.xlabel(xl),Guide.ylabel(yl))
        else
            plot(x=x,y=y,xgroup=xgroup,xintercept=fill(timemark[1],length(x)),theme,
            Geom.subplot_grid(Geom.point,Geom.vline(color="gray",size=1pt),free_x_axis=true),
            Guide.xlabel(xl),Guide.ylabel(yl))
        end
    else
        yl="$yl Sorted"
        if isempty(colorminv);colorminv=minimum(c);end
        if isempty(colormaxv);colormaxv=maximum(c);end
        if isempty(xgroup)
            plot(x=x,y=y,color=c,xintercept=timemark,theme,Geom.point,Geom.vline(color="gray",size=1pt),
            Coord.Cartesian(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),Guide.xlabel(xl),Guide.ylabel(yl),Guide.colorkey(colorkey),
            Scale.ContinuousColorScale(colorfun,minvalue=colorminv,maxvalue=colormaxv))
        else
            plot(x=x,y=y,color=c,xgroup=xgroup,xintercept=fill(timemark[1],length(x)),theme,
            Geom.subplot_grid(Geom.point,Geom.vline(color="gray",size=1pt),free_x_axis=true),
            Guide.xlabel(xl),Guide.ylabel(yl),Guide.colorkey(colorkey),
            Scale.ContinuousColorScale(colorfun,minvalue=colorminv,maxvalue=colormaxv))
        end
    end
end
plotspiketrain1(rvs::RVVector;sortvar=[],xgroup::Vector=[],timemark=[0],theme=Theme(),colorkey="",colorfun=Scale.lab_gradient(colorant"white",colorant"red"),colorminv=[],colormaxv=[]) = plotspiketrain1(flatrvs(rvs,sortvar)...,xgroup=xgroup,timemark=timemark,theme=theme,colorkey=colorkey,colorfun=colorfun,colorminv=colorminv,colormaxv=colormaxv)

function plotpsth1(rvs::RVVector,binedges::RealVector;theme=Theme(),timeline=[0],title="")
    m,sd,n,x = psth(rvs,binedges)
    xl = "Time (ms)";yl = "Response (spike/s)"
    if n>1
        Gadfly.plot(y=m,x=x,ymin=m-sd/sqrt(n),ymax=m+sd/sqrt(n),xintercept=timeline,theme,Geom.line,Geom.ribbon,
        Geom.vline(color="gray",size=1pt),Coord.Cartesian(xmin=binedges[1],xmax=binedges[end],ymin=0),
        Guide.xlabel(xl),Guide.ylabel(yl),Guide.title(title))
    else
        Gadfly.plot(y=m,x=x,xintercept=timeline,theme,Geom.line,Geom.vline(color="gray",size=1pt),
        Coord.Cartesian(xmin=binedges[1],xmax=binedges[end],ymin=0),Guide.xlabel(xl),Guide.ylabel(yl),Guide.title(title))
    end
end
function plotpsth1(rvs::RVVector,binedges::RealVector,rvsidx,condstr;theme=Theme(),timeline=[0],title="")
    df = psth(rvs,binedges,rvsidx,condstr)
    df[:ymin] = df[:y]-df[:ysd]./sqrt(df[:n])
    df[:ymax] = df[:y]+df[:ysd]./sqrt(df[:n])
    xl = "Time (ms)";yl = "Response (spike/s)"
    Gadfly.plot(df,y=:y,x=:x,ymin=:ymin,ymax=:ymax,color=:condition,xintercept=timeline,theme,Geom.line,Geom.ribbon,
    Geom.vline(color="gray",size=1pt),Coord.Cartesian(xmin=binedges[1],xmax=binedges[end],ymin=0),
    Guide.xlabel(xl),Guide.ylabel(yl),Guide.title(title))
end

function plotpsth1(ds::DataFrame,binedges::RealVector,conds::Vector{Vector{Any}};spike=:spike,theme=Theme(),timeline=[0])
    df,ss = psth(ds,binedges,conds,spike=spike)
    xl = "Time (ms)";yl = "Response (spike/s)"
    Gadfly.plot(df,x=:x,y=:y,ymin=:ymin,ymax=:ymax,color=:condition,xintercept=timeline,
    theme,Geom.line,Geom.ribbon,Geom.vline(color="gray",size=1pt),
    Coord.Cartesian(xmin=binedges[1],xmax=binedges[end],ymin=0),Guide.xlabel(xl),Guide.ylabel(yl))
end

plotcondresponse(rs,ctc,factor;u=0,style=:path,title="",projection=[],linewidth=:auto,legend=:best)=plotcondresponse(Dict(u=>rs),ctc,factor,style=style,title=title,projection=projection,linewidth=linewidth,legend=legend)
function plotcondresponse(urs::Dict,ctc::DataFrame,factor;colors=unitcolors(collect(keys(urs))),style=:path,projection=[],title="",linewidth=:auto,legend=:best)
    mseuc = condresponse(urs,ctc,factor)
    plotcondresponse(mseuc,colors=colors,style=style,title=title,projection=projection,linewidth=linewidth,legend=legend)
end
function plotcondresponse(urs::Dict,cond::DataFrame;colors=unitcolors(collect(keys(urs))),style=:path,projection=[],title="",linewidth=:auto,legend=:best)
    mseuc = condresponse(urs,cond)
    plotcondresponse(mseuc,colors=colors,style=style,title=title,projection=projection,linewidth=linewidth,legend=legend)
end
function plotcondresponse(mseuc::DataFrame;colors=unitcolors(unique(mseuc[:u])),style=:path,projection=[],title="",linewidth=:auto,legend=:best)
    us = sort(unique(mseuc[:u]))
    factor=setdiff(names(mseuc),[:m,:se,:u])
    nfactor=length(factor)
    if nfactor==1
        factor=factor[1]
        if typeof(mseuc[factor][1]) <: Array
            map!(string,mseuc[factor],mseuc[factor])
            style=:bar
        end
    else
        mseuc[:Condition]=condstring(mseuc[:,factor])
        factor=:Condition
        style=:bar
    end
    if projection==:polar
        c0 = mseuc[mseuc[factor].==0,:]
        c0[factor]=360
        mseuc = [mseuc;c0]
        mseuc[factor]=deg2rad.(mseuc[factor])
    end
    sort!(mseuc,factor)
    if projection==:polar
        @df mseuc Plots.plot(cols(factor),:m,group=:u,line=style,markerstrokecolor=:auto,color=reshape(colors,1,:),label=reshape(["U$k" for k in us],1,:),
        grid=false,projection=projection,legend=legend,xaxis=(factorunit(factor)),yaxis=(factorunit(:Response)),title=(title),linewidth=linewidth)
    else
        @df mseuc Plots.plot(cols(factor),:m,yerror=:se,group=:u,line=style,markerstrokecolor=:auto,color=reshape(colors,1,:),label=reshape(["U$k" for k in us],1,:),
        grid=false,projection=projection,legend=legend,xaxis=(factorunit(factor)),yaxis=(factorunit(:Response)),title=(title),linewidth=linewidth)
    end
end

plotpsth(rvv::RVVector,binedges::RealVector;timeline=[0],colors=[:auto],title="")=plotpsth(rvv,binedges,DataFrame(Factor="Value",i=[1:length(rvv)]),timeline=timeline,colors=colors,title=title)
function plotpsth(rvv::RVVector,binedges::RealVector,ctc::DataFrame,factor;timeline=[0],colors=nothing,title="")
    msexc = psth(rvv,binedges,ctc,factor)
    plotpsth(msexc,timeline=timeline,colors=colors==nothing ? huecolors(length(levels(msexc[:c]))) : colors,title=title)
end
function plotpsth(rvv::RVVector,binedges::RealVector,cond::DataFrame;timeline=[0],colors=huecolors(nrow(cond)),title="")
    msexc = psth(rvv,binedges,cond)
    plotpsth(msexc,timeline=timeLine,colors=colors,title=title)
end
function plotpsth(msexc::DataFrame;timeline=[0],colors=[:auto],title="")
    @df msexc Plots.plot(:x,:m,ribbon=:se,group=:c,fillalpha=0.2,color=reshape(colors,1,:))
    vline!(timeline,line=(:grey),label="TimeLine",grid=false,xaxis=(factorunit(:Time)),yaxis=(factorunit(:Response)),title=(title))
end

function plotsta(α;delay=nothing,decor=false,savedir=nothing)
    ds = delay==nothing ? "" : "_$(delay)ms"
    t = (decor ? "d" : "") * "STA$ds"
    p=Plots.plot(α,seriestype=:heatmap,color=:fire,ratio=:equal,yflip=true,leg=false,framestyle=:none,title=t)
    if savedir!=nothing
        !isdir(savedir) && mkpath(savedir)
        png(p,joinpath(savedir,t))
    end
    p
end

function savefig(fig,filename::AbstractString;path::AbstractString="",format::AbstractString="svg")
    f = joinpath(path,"$filename.$format")
    if !ispath(path)
        mkpath(path)
    end
    if format=="svg"
        format = "$format+xml"
    end
    open(f, "w") do io
        writemime(io,"image/$format",fig)
    end
end

function savefig(fig::Gadfly.Plot,filename::AbstractString;path::AbstractString="",format::AbstractString="svg",width=22cm,height=13cm,dpi=300)
    f = joinpath(path,"$filename.$format")
    if !ispath(path)
        mkpath(path)
    end
    if format=="svg"
        format = SVG(f,width,height)
    elseif format=="png"
        format = PNG(f,width,height,dpi=dpi)
    elseif format=="pdf"
        format = PDF(f,width,height,dpi=dpi)
    end
    draw(format,fig)
end
