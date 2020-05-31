using Plots,StatsPlots,VegaLite
import Plots: cgrad

include("Colormap.jl")

factorunit(fs::Vector{Symbol};timeunit=SecondPerUnit)=join(factorunit.(fs,timeunit=timeunit),", ")
function factorunit(f::Symbol;timeunit=SecondPerUnit)
    fu=String(f)
    if any(occursin.(["Ori","Angle"],fu))
        fu="$fu (deg)"
    elseif fu=="dir"
        fu="Direction (deg)"
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
    elseif fu=="ResponseF"
        fu="Response (% \\DeltaF / F)"
    end
    return fu
end

"scatter plot of spike trains"
function plotspiketrain(x,y;group::Vector=[],timeline=[0],colors=unitcolors(),title="",size=(800,550))
    nt = isempty(x) ? 0 : maximum(y)
    s = min(size[2]/nt,1)
    if isempty(group)
        scatter(x,y,label="SpikeTrain",markershape=:vline,size=size,markersize=s,markerstrokewidth = s,markerstrokecolor=RGBA(0.1,0.1,0.3,0.8),legend=false)
    else
        scatter(x,y,group=group,markershape=:vline,size=size,markersize=s,markerstrokewidth = s,markerstrokecolor=reshape(colors,1,:))
    end
    vline!(timeline,line=(:grey),label="TimeLine",grid=false,xaxis=(factorunit(:Time)),yaxis=("Trial"),title=(title),legend=false)
end
function plotspiketrain(sts::Vector;uids::Vector=[],sortvalues=[],timeline=[0],colors=unitcolors(),title="",size=(800,550))
    if isempty(uids)
        g=uids;uc=colors
    else
        fuids = flatspiketrains(uids,sortvalues)[1]
        g=map(i->"U$i",fuids);uc=colors[sort(unique(fuids)).+1]
    end
    plotspiketrain(flatspiketrains(sts,sortvalues)[1:2]...,group=g,timeline=timeline,colors=uc,title=title,size=size)
end

# function plotspiketrain1(x::Vector,y::Vector,c::Vector=[];xmin=minimum(x)-10,xmax=maximum(x)+10,xgroup::Vector=[],
#     ymin=minimum(y)-1,ymax=maximum(y)+1,timemark=[0],theme=Theme(),
#     colorkey="",colorfun=Scale.lab_gradient(colorant"white",colorant"red"),colorminv=[],colormaxv=[])
#     xl="Time (ms)";yl="Trial"
#     if isempty(c)
#         if isempty(xgroup)
#             plot(x=x,y=y,xintercept=timemark,theme,Geom.point,Geom.vline(color="gray",size=1pt),
#             Coord.Cartesian(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),Guide.xlabel(xl),Guide.ylabel(yl))
#         else
#             plot(x=x,y=y,xgroup=xgroup,xintercept=fill(timemark[1],length(x)),theme,
#             Geom.subplot_grid(Geom.point,Geom.vline(color="gray",size=1pt),free_x_axis=true),
#             Guide.xlabel(xl),Guide.ylabel(yl))
#         end
#     else
#         yl="$yl Sorted"
#         if isempty(colorminv);colorminv=minimum(c);end
#         if isempty(colormaxv);colormaxv=maximum(c);end
#         if isempty(xgroup)
#             plot(x=x,y=y,color=c,xintercept=timemark,theme,Geom.point,Geom.vline(color="gray",size=1pt),
#             Coord.Cartesian(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),Guide.xlabel(xl),Guide.ylabel(yl),Guide.colorkey(colorkey),
#             Scale.ContinuousColorScale(colorfun,minvalue=colorminv,maxvalue=colormaxv))
#         else
#             plot(x=x,y=y,color=c,xgroup=xgroup,xintercept=fill(timemark[1],length(x)),theme,
#             Geom.subplot_grid(Geom.point,Geom.vline(color="gray",size=1pt),free_x_axis=true),
#             Guide.xlabel(xl),Guide.ylabel(yl),Guide.colorkey(colorkey),
#             Scale.ContinuousColorScale(colorfun,minvalue=colorminv,maxvalue=colormaxv))
#         end
#     end
# end
# plotspiketrain1(rvs::RVVector;sortvar=[],xgroup::Vector=[],timemark=[0],theme=Theme(),colorkey="",colorfun=Scale.lab_gradient(colorant"white",colorant"red"),colorminv=[],colormaxv=[]) = plotspiketrain1(flatrvs(rvs,sortvar)...,xgroup=xgroup,timemark=timemark,theme=theme,colorkey=colorkey,colorfun=colorfun,colorminv=colorminv,colormaxv=colormaxv)


plotcondresponse(rs,ctc;factors=names(ctc),u=0,style=:path,title="",projection=[],linewidth=:auto,legend=:best,responseline=[])=plotcondresponse(Dict(u=>rs),ctc,factors,style=style,title=title,projection=projection,linewidth=linewidth,legend=legend,responseline=responseline)
function plotcondresponse(urs::Dict,ctc::DataFrame,factors;colors=unitcolors(collect(keys(urs))),style=:path,projection=[],title="",linewidth=:auto,legend=:best,responseline=[])
    mseuc = condresponse(urs,ctc,factors)
    plotcondresponse(mseuc,colors=colors,style=style,title=title,projection=projection,linewidth=linewidth,legend=legend,responseline=responseline)
end
function plotcondresponse(urs::Dict,cond::DataFrame;colors=unitcolors(collect(keys(urs))),style=:path,projection=[],title="",linewidth=:auto,legend=:best,responseline=[])
    mseuc = condresponse(urs,cond)
    plotcondresponse(mseuc,colors=colors,style=style,title=title,projection=projection,linewidth=linewidth,legend=legend,responseline=responseline)
end
function plotcondresponse(mseuc::DataFrame;colors=unitcolors(unique(mseuc[:ug])),style=:path,projection=[],title="",linewidth=:auto,legend=:best,responseline=[],responsetype=:Response)
    ugs = sort(unique(mseuc[:,[:u,:ug]]))
    factors=setdiff(names(mseuc),[:m,:se,:u,:ug])
    nfactor=length(factors)
    if nfactor==1
        factor=factors[1]
        if typeof(mseuc[!,factor][1]) <: Array
            map!(string,mseuc[factor],mseuc[factor])
            style=:bar
        end
    elseif nfactor==2
        fm,fse,fa = factorresponse(mseuc)
        clim=maximum(skipmissing(fm))
        yfactor,xfactor = collect(keys(fa))
        y,x = collect(values(fa))
    else
        mseuc[:Condition]=condstring(mseuc[:,factors])
        factor=:Condition
        style=:bar
    end
    if nfactor==2
        x=float.(x)
        y=float.(y)
        heatmap(x,y,fm,color=:fire,title=title,legend=legend,xaxis=(factorunit(xfactor)),yaxis=(factorunit(yfactor)),colorbar_title=factorunit(responsetype),clims=(0,clim))
    else
        if projection==:polar
            c0 = mseuc[mseuc[!,factor].==0,:]
            c0[:,factor].=360
            mseuc = [mseuc;c0]
            mseuc[!,factor]=deg2rad.(mseuc[!,factor])
        end
        sort!(mseuc,factor)
        if projection==:polar
            p = @df mseuc Plots.plot(cols(factor),:m,yerror=:se,group=:ug,line=style,markerstrokecolor=:auto,color=reshape(colors,1,:),label=reshape(["$(k.ug)$(k.u)" for k in eachrow(ugs)],1,:),
            grid=false,projection=projection,legend=legend,xaxis=(factorunit(factor)),yaxis=(factorunit(responsetype)),title=(title),linewidth=linewidth)
        else
            p = @df mseuc plot(cols(factor),:m,yerror=:se,group=:ug,line=style,markerstrokecolor=:auto,color=reshape(colors,1,:),label=reshape(["$(k.ug)$(k.u)" for k in eachrow(ugs)],1,:),
            grid=false,projection=projection,legend=legend,xaxis=(factorunit(factor)),yaxis=(factorunit(responsetype)),title=(title),linewidth=linewidth)
        end
        if !isempty(responseline)
            for i in responseline
                hline!(p,[i[1]],ribbon=[i[2]],color=colors,legend=false)
            end
        end
        p
    end
end

plotpsth(rvv::RVVector,binedges::RealVector;timeline=[0],colors=[:auto],title="")=plotpsth(rvv,binedges,DataFrame(Factor="Value",i=[1:length(rvv)]),timeline=timeline,colors=colors,title=title)
function plotpsth(rvv::RVVector,binedges::RealVector,ctc::DataFrame,factor;timeline=[0],colors=nothing,title="")
    msexc = psth(rvv,binedges,ctc,factor)
    plotpsth(msexc,timeline=timeline,colors=colors==nothing ? huecolors(length(levels(msexc[:c]))) : colors,title=title)
end
function plotpsth(rvv::RVVector,binedges::RealVector,cond::DataFrame;timeline=[0],colors=huecolors(nrow(cond)),title="")
    msexc = psth(rvv,binedges,cond)
    plotpsth(msexc,timeline=timeline,colors=colors,title=title)
end
function plotpsth(msexc::DataFrame;timeline=[0],colors=[:auto],title="")
    @df msexc Plots.plot(:x,:m,ribbon=:se,group=:c,fillalpha=0.2,color=reshape(colors,1,:))
    vline!(timeline,line=(:grey),label="TimeLine",grid=false,xaxis=(factorunit(:Time)),yaxis=(factorunit(responsetype)),title=(title))
end
function plotpsth(data::RealMatrix,x,y;color=:Reds,timeline=[0],hlines=[],layer=nothing,n=[])
    xms = x*SecondPerUnit*1000
    p=heatmap(xms,y,data,color=color,colorbar_title="Spike/Sec",xlabel="Time (ms)",ylabel="Depth (um)")
    vline!(p,timeline,color=:gray,label="TimeLine")
    if !isempty(n)
        pn = n./maximum(n) .* maximum(x) .* 0.2 .+ minimum(x)
        plot!(p,pn,y,label="Number of Units",color=:seagreen)
    end
    if !isnothing(layer)
        lx = minimum(xms)+5
        hline!(p,[layer[k][1] for k in keys(layer)],linestyle=:dash,annotations=[(lx,layer[k][1],text(k,6,:gray20,:bottom)) for k in keys(layer)],linecolor=:gray30,legend=false)
    end
    return p
end

function plotsta(ps;sizepx=size(ps),sizedeg=nothing,ppd=45,index=nothing,filter=Kernel.gaussian(1),title="",color="redblue",r=[extrema(ps)...],bg="white")
    nd = ndims(ps)
    if nd==1
        if isnothing(index)
            sta = reshape(ps,sizepx)
        else
            sta = fill(mean(ps),sizepx)
            sta[index] = ps
        end
    elseif nd==2
        sta = ps
    end
    if !isnothing(filter)
        sta = imfilter(sta,filter)
    end
    if !isnothing(sizedeg)
        ppd = first(sizepx./sizedeg)
    end

    x = vec([(j-1)/ppd for i in 1:sizepx[1], j in 1:sizepx[2]])
    y = vec([(i-1)/ppd for i in 1:sizepx[1], j in 1:sizepx[2]])
    z = vec(sta)
    if !isnothing(index)
        x=x[index];y=y[index];z=z[index]
    end
    xlim=[extrema(x)...];ylim=[extrema(y)...];psize=2;width=psize*length(unique(x));height=psize*length(unique(y))

    DataFrame(z=z,x=x,y=y) |> @vlplot(mark={:rect,size=psize,strokeWidth=0},width=width,height=height,background=bg,
    x={"x:o",title="X (deg)",axis={values=xlim,format=".1",labelAngle=0}},
    y={"y:o",title="Y (deg)",axis={values=ylim,format=".1"}},
    color={"z:q",title="",scale={domain=r}},title={text=title},
    config={range={heatmap={scheme=color,extent=[1,0]}},
    view={strokeWidth=0,fill=bg}})

    # x=[(i-1)/ppd for i in 1:imagesize[1]]
    #
    # color=:coolwarm
    # if true
    #     cg = cgrad(color)
    #     minv,maxv = extrema(sta)
    #     vr=maxv-minv
    #     sta = map(i->coloralpha(RGB(cg[(i-minv)/vr]),0),sta)
    #     sta[index] = coloralpha.(sta[index],1)
    # end
    # Plots.heatmap(x,x,sta,color=:coolwarm,ratio=:equal,yflip=true,leg=true,framestyle=:grid,title=title,
    # xlabel="Position_X (Deg)",ylabel="Position_Y (Deg)",xtick=xlim,ytick=[])
end

function plotanalog(data;x=nothing,y=nothing,fs=0,xext=0,timeline=[0],xlabel="Time",xunit=:ms,cunit=:v,plottype=:heatmap,ystep=20,color=:coolwarm,layer=nothing)
    nd=ndims(data)
    if nd==1
        x=1:length(y)
        if fs>0
            x = x./fs.-xext
            if timeunit==:ms
                x*=1000
            end
        end
        ylim=maximum(abs.(y))
        p=plot(x,y,ylims=(-ylim,ylim))
    elseif nd==2
        if isnothing(x)
            x=1:size(data,2)
            if fs>0
                x = x./fs.-xext
                if xunit==:ms
                    x*=1000
                end
            end
        end
        df=1
        if cunit==:v
            lim = maximum(abs.(data))
            clim = (-lim,lim)
        elseif cunit==:uv
            df = 1e6
            lim = maximum(abs.(data))*df
            clim = (-lim,lim)
        elseif cunit == :fr
            lim = maximum(data)
            clim = (0,lim)
        else
            clim = :auto
        end
        if plottype==:heatmap
            if isnothing(y)
                y = (1:size(data,1))*ystep
            end
            p=heatmap(x,y,data.*df,color=color,clims=clim,xlabel="$xlabel ($xunit)")
        else
            p=plot(x,data'.*df,legend=false,color_palette=color,grid=false,ylims=clim,xlabel="$xlabel ($xunit)")
        end
    end
    !isempty(timeline) && vline!(p,timeline,line=(:grey),label="TimeLine")
    if !isnothing(layer)
        lx = minimum(x)+5
        hline!(p,[layer[k][1] for k in keys(layer)],linestyle=:dash,annotations=[(lx,layer[k][1],text(k,6,:gray20,:bottom)) for k in keys(layer)],linecolor=:gray30,legend=false)
    end
    return p
end

plotunitposition(spike::Dict;layer=nothing,color=nothing,alpha=0.4,title="",markersize=5,unitidsize=3) = plotunitposition(spike["unitposition"],unitgood=spike["unitgood"],chposition=spike["chposition"],unitid=spike["unitid"],layer=layer,color=color,alpha=alpha,title=title,markersize=markersize,unitidsize=unitidsize)
function plotunitposition(unitposition;unitgood=[],chposition=[],unitid=[],layer=nothing,color=nothing,alpha=0.4,title="",markersize=5,unitidsize=3)
    nunit = size(unitposition,1);ngoodunit = isempty(unitgood) ? nunit : count(unitgood);us = "$ngoodunit/$nunit"
    xlim = isempty(chposition) ? (minimum(unitposition[:,1])-4,maximum(unitposition[:,1])+2) : (minimum(chposition[:,1])-5,maximum(chposition[:,1])+5)
    p = plot(legend=:topright,xlabel="Position_X (um)",ylabel="Position_Y (um)",xlims=xlim)
    if !isempty(chposition)
        scatter!(p,chposition[:,1],chposition[:,2],markershape=:rect,markerstrokewidth=0,markersize=2,color=:grey60,label="Electrode")
    end
    if isnothing(color)
        if !isempty(unitgood)
            color = map(i->i ? :darkgreen : :gray30,unitgood)
        else
            color = :gray30
        end
        if !isnothing(alpha)
            color = coloralpha.(parse.(RGB,color),alpha)
        end
    end
    if !isnothing(layer)
        lx = xlim[1]+2
        hline!(p,[layer[k][1] for k in keys(layer)],linestyle=:dash,annotations=[(lx,layer[k][1],text(k,5,:gray20,:bottom)) for k in keys(layer)],linecolor=:gray30,legend=false)
    end
    if !isempty(unitid)
        scatter!(p,unitposition[:,1],unitposition[:,2],label=us,color=color,markerstrokewidth=0,markersize=markersize,series_annotations=text.(unitid,unitidsize,:gray10,:center),title=title)
    else
        scatter!(p,unitposition[:,1],unitposition[:,2],label=us,color=color,markerstrokewidth=0,markersize=markersize,title=title)
    end
    return p
end

function plotunitpositionproperty(unitposition;ori=nothing,os=nothing,dir=nothing,ds=nothing,sf=nothing,width=500,height=400,title="",layer=nothing)
    df = DataFrame(x=unitposition[:,1],y=unitposition[:,2],m=Any[:circle for _ in 1:size(unitposition,1)],a=0.0,sw=0.0,c=0.0,s=10.0)
    l = DataFrame()
    xlim = [minimum(df[:,:x])-4,maximum(df[:,:x])+2]
    ylim = [minimum(df[:,:y])-100,maximum(df[:,:y])+100]
    if !isnothing(sf)
        df[!,:c]=sf
    end
    if !isnothing(ori)
        df[!,:m] .= :stroke
        df[!,:a] = -ori
        df[!,:sw] .= 1
        if !isnothing(os)
            df[!,:s]=os
        else
            df[!,:s].=160
        end
    end
    if !isnothing(dir)
        t = copy(df)
        arrowpath = "M 0 -0.1 H 1 V -0.3 L 1.6 0 L 1 0.3 V 0.1 H 0 Z"
        t[!,:m] .= arrowpath
        t[!,:a] = -dir
        t[!,:sw] .= 0
        if !isnothing(ds)
            t[!,:s]=ds
        else
            t[!,:s].=200
        end
        df = [df;t]
    end
    if !isnothing(layer)
        l[!,:x] = fill(xlim[1],length(layer))
        l[!,:x2] .= xlim[2]
        l[!,:y] = [v[1] for v in values(layer)]
        l[!,:y2] = l[!,:y]
        l[!,:l] = collect(keys(layer))
        ylim = [extrema([ylim;l[:,:y]])...].+[-50,100]
    end
    @vgplot(height=height,width=width,padding=5,data=[:df=>df,:l=>l],
    marks=[
    {
    type="rule",
    from={data="l"},
    encode={
        update={
        x={field="x",scale="x"},
        y={field="y",scale="y"},
        x2={field="x2",scale="x"},
        y2={field="y2",scale="y"},
        strokeDash={value=[4,2]},
        strokeWidth={value=0.5},
        stroke={value="dimgray"}
        }}
    },
    {
    type="text",
    from={data="l"},
    encode={
        update={
        x={field="x",scale="x"},
        y={field="y",scale="y"},
        text={field="l"},
        align={value="left"},
        baseline={value="bottom"},
        fontSize={value=7},
        dx={value=15}
        }}
    },
    {
    type="symbol",
    from={data="df"},
    encode={
        update={
        x={field="x",scale="x"},
        y={field="y",scale="y"},
        shape={field="m"},
        angle={field="a"},
        size={field="s",scale="s"},
        strokeWidth={field="sw"},
        stroke={field="c",scale="c"},
        fill={field="c",scale="c"}
        }}
    }
    ],
    scales=[
    {
        name="x",
        nice=false,
        zero=false,
        range="width",
        domain=xlim,
        type="linear",
        round=true
    },
    {
        name="y",
        nice=false,
        zero=false,
        range="height",
        domain=ylim,
        type="linear",
        round=true
    },
    {
        name="c",
        nice=false,
        zero=false,
        round=false,
        type="linear",
        range={scheme = "plasma",extent=[0.8,0.2]},
        domain={data="df",field="c"}
    },
    {
        name="s",
        nice=false,
        zero=false,
        domain={data="df",field="s"},
        type="linear",
        round=true,
        range=[10, 200]
    }
    ],
    axes=[
    {
        domain=true,
        tickCount=5,
        grid=false,
        title="Position_X (μm)",
        scale="x",
        orient="bottom"
    },
    {
        domain=true,
        tickCount=5,
        grid=false,
        title="Position_Y (μm)",
        scale="y",
        orient="left"
    }
    ],
    title={
    text = title
    },
    legends=[
    {
    type="gradient",
    fill="c",
    title="SF"
    },
    {
    type="symbol",
    symbolType="stroke",
    size="s",
    title="1-CV"
    }
    ])
end

function plotunitpositionimage(unitposition,image;width=800,height=600,markersize=20,title="",layer=nothing)
    tempimagedir = joinpath(pwd(),"tempimage");nu=size(unitposition,1)
    isdir(tempimagedir) || mkpath(tempimagedir)
    urls = map(i->"tempimage/$(i).png",1:nu)
    foreach(i->save(urls[i],image[i]),1:nu)

    df = DataFrame(x=unitposition[:,1],y=unitposition[:,2],m=urls,s=markersize)
    l = DataFrame()
    xlim = [minimum(df[:,:x])-4,maximum(df[:,:x])+2]
    ylim = [minimum(df[:,:y])-100,maximum(df[:,:y])+100]
    if !isnothing(layer)
        l[!,:x] = fill(xlim[1],length(layer))
        l[!,:x2] .= xlim[2]
        l[!,:y] = [v[1] for v in values(layer)]
        l[!,:y2] = l[!,:y]
        l[!,:l] = collect(keys(layer))
        ylim = [extrema([ylim;l[:,:y]])...].+[-50,100]
    end
    @vgplot(height=height,width=width,padding=5,data=[:df=>df,:l=>l],
    marks=[
    {
    type="rule",
    from={data="l"},
    encode={
        update={
        x={field="x",scale="x"},
        y={field="y",scale="y"},
        x2={field="x2",scale="x"},
        y2={field="y2",scale="y"},
        strokeDash={value=[4,2]},
        strokeWidth={value=0.5},
        stroke={value="dimgray"}
        }}
    },
    {
    type="text",
    from={data="l"},
    encode={
        update={
        x={field="x",scale="x"},
        y={field="y",scale="y"},
        text={field="l"},
        align={value="left"},
        baseline={value="bottom"},
        fontSize={value=7},
        dx={value=15}
        }}
    },
    {
    type="image",
    from={data="df"},
    encode={
        update={
        x={field="x",scale="x"},
        y={field="y",scale="y"},
        url={field="m"},
        width={field="s"},
        height={field="s"},
        align={value="center"},
        baseline={value="middle"}
        }}
    }
    ],
    scales=[
    {
        name="x",
        nice=false,
        zero=false,
        range="width",
        domain=xlim,
        type="linear",
        round=true
    },
    {
        name="y",
        nice=false,
        zero=false,
        range="height",
        domain=ylim,
        type="linear",
        round=true
    }
    ],
    axes=[
    {
        domain=true,
        tickCount=5,
        grid=false,
        title="Position_X (μm)",
        scale="x",
        orient="bottom"
    },
    {
        domain=true,
        tickCount=5,
        grid=false,
        title="Position_Y (μm)",
        scale="y",
        orient="left"
    }
    ],
    title={
    text = title
    }
    )
end

function plotunitlayerimage(unitlayer,image;width=800,height=600,markersize=40,title="",unitid=nothing)
    tempimagedir = joinpath(pwd(),"tempimage");nu=length(unitlayer)
    isdir(tempimagedir) || mkpath(tempimagedir)
    urls = map(i->"tempimage/$(i).png",1:nu)
    foreach(i->save(urls[i],image[i]),1:nu)

    x = zeros(nu)
    for r in eachrow(condin(DataFrame(l=unitlayer)))
        x[r.i] = 1:r.n
    end
    df = DataFrame(x=x,y=unitlayer,m=urls)
    if isnothing(unitid)
        df[:u] = ""
    else
        df[:u] = unitid
    end

    @vlplot(height=height,width=width,padding=5,title=title,data=df,
    y={"y:o",axis={title="Layer"}},
    x={"x",axis={grid=false,title=nothing,ticks=false,domain=false,labels=false}}) +
    @vlplot(mark={:image,width=markersize,height=markersize,align=:center,baseline=:middle}, url=:m) +
    @vlplot(mark={:text,align=:center,baseline=:bottom,dy=30,fontSize=7},text=:u)
end

function plotcircuit(unitposition,unitid,projs;unitgood=[],eunits=[],iunits=[],projweights=[],layer=nothing,showuid=true,showmode=:none)
    vpi = indexin(unique(projs),projs)
    projs=projs[vpi];np=length(projs)
    if !isempty(projweights)
        projweights=projweights[vpi]
        t = abs.(projweights)
        pcolor = coloralpha.(RGB(0.5,0.5,0.5),t./maximum(t))
    else
        pcolor=:gray50
    end
    nu = length(unitid)
    nsu = isempty(unitgood) ? nu : count(unitgood)
    uidi = Dict(i=>findfirst(i.==unitid) for i in unitid)
    nsupair=binomial(nsu,2);ss = "$nsu: $np/$nsupair($(round(np/nsupair*100,digits=3))%)"
    xlim = (minimum(unitposition[:,1])-5,maximum(unitposition[:,1])+5)
    ylim = (minimum(unitposition[:,2]),maximum(unitposition[:,2]))
    p = plot(legend=:topright,xlabel="Position_X (um)",ylabel="Position_Y (um)",xlims=xlim,grid=false)

    for i in projs
        t=hcat(unitposition[uidi[i[1]],:],unitposition[uidi[i[2]],:])
        plot!(p,t[1,:],t[2,:],linewidth=0.3,color=pcolor,arrow=arrow(:closed,:head,0.45,0.12))
    end

    color = fill(:gray30,nu)
    color[unitgood] .= :darkgreen
    if !isempty(eunits)
        color[map(i->uidi[i],eunits)] .= :darkred
    end
    if !isempty(iunits)
        color[map(i->uidi[i],iunits)] .= :darkblue
    end
    if showmode == :su
        showuidi = unitgood
    elseif showmode == :circuit
        cuid=[]
        foreach(p->append!(cuid,[p...]),projs)
        showuidi=indexin(unique(cuid),unitid)
    else
        showuidi = 1:nu
    end
    if showuid
        scatter!(p,unitposition[showuidi,1],unitposition[showuidi,2],label=ss,color=color[showuidi],alpha=0.4,markerstrokewidth=0,markersize=6,series_annotations=text.(unitid[showuidi],3,:gray10,:center),legend=false)
    else
        scatter!(p,unitposition[showuidi,1],unitposition[showuidi,2],label=ss,color=color[showuidi],alpha=0.4,markerstrokewidth=0,markersize=5,legend=false)
    end
    if !isnothing(layer)
        hline!(p,[layer[k][1] for k in keys(layer)],linestyle=:dash,annotations=[(xlim[1]+2,layer[k][1],text(k,5,:gray20,:bottom)) for k in keys(layer)],linecolor=:gray30,legend=false)
    end
    annotate!(p,[(xlim[2]-6,ylim[2]-100,text(ss,6,:gray20,:bottom))])
    return p
end

function plothartleyspace(space,nk,dk)
    n = 2nk+1
    p=plot(layout=(n,n),size=(120n,120n),leg=false,clims=(-1,1),frame=:none,aspect_ratio=:equal,yflip=true)
    x=y=0:0.01:1
    for (kx,ky,phase) in space
        r = -Int(ky/dk)+nk+1; c = Int(kx/dk)+nk+1
        cg = [cas(i,j,kx=kx,ky=ky,phase=phase) for j in reverse(y),i in x]
        heatmap!(p,subplot=c+n*(r-1),cg,color=:grays)
    end
    p
end

# function savefig(fig,filename::AbstractString;path::AbstractString="",format::AbstractString="svg")
#     f = joinpath(path,"$filename.$format")
#     if !ispath(path)
#         mkpath(path)
#     end
#     if format=="svg"
#         format = "$format+xml"
#     end
#     open(f, "w") do io
#         writemime(io,"image/$format",fig)
#     end
# end
#
# function savefig(fig::Gadfly.Plot,filename::AbstractString;path::AbstractString="",format::AbstractString="svg",width=22cm,height=13cm,dpi=300)
#     f = joinpath(path,"$filename.$format")
#     if !ispath(path)
#         mkpath(path)
#     end
#     if format=="svg"
#         format = SVG(f,width,height)
#     elseif format=="png"
#         format = PNG(f,width,height,dpi=dpi)
#     elseif format=="pdf"
#         format = PDF(f,width,height,dpi=dpi)
#     end
#     draw(format,fig)
# end
