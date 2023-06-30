using Plots,StatsPlots,VegaLite
import Plots: cgrad

include("Color.jl")

iscircfactor(f::Symbol) = f in [:Ori,:Dir,:Ori_Final,:Dir_Final,:HueAngle,:Angle]
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
factorunitstring(f,u) = factorunitstring(string(f),string(u))
function factorunitstring(f::String,u::String)
    isempty(f) && return f
    isempty(u) && return f
    "$f ($u)"
end

"scatter plot of spikes"
function plotspiketrain(x,y;group=[],timeline=[],timespan=[],color=RGBA(0,0,0,0.7),title="",size=(800,600),spancolor=RGBA(0.5,0.5,0.5,0.3),
                        linewidth=:auto,linecolor=:gray,linestyle=:solid,markershape=:vline,markersize=:auto,markerstrokewidth=1,frame=:axes)
    if markersize==:auto
        nt = isempty(x) ? 0 : maximum(y)
        markersize =0.5*size[2]/(nt+10)
    end

    p = plot(;frame,size,leg=false,title,grid=false,tickdir=:out)
    isempty(timespan) || vspan!(p,timespan;color=spancolor,label="TimeSpan")
    isempty(timeline) || vline!(p,timeline;linewidth,linecolor,linestyle,label="TimeLine")
    if isempty(group)
        scatter!(p,x,y;markershape,markersize,markerstrokewidth,color,xlabel=factorunit(:Time),ylabel="Trial",label="SpikeTrain")
    else
        color=huecolors(length(unique(group)))'
        scatter!(p,x,y;group,markershape,markersize,markerstrokewidth,color,xlabel=factorunit(:Time),ylabel="Trial",label="SpikeTrain")
    end
end
"scatter plot of spike trains"
function plotspiketrain(sts::AbstractVecOrMat;uids=[],trialorder=[],timeline=[0],timespan=[],color=RGBA(0,0,0,0.7),title="",size=(800,600),spancolor=RGBA(0.5,0.5,0.5,0.3),
                        linewidth=:auto,linecolor=:gray,linestyle=:solid,markershape=:vline,markersize=:auto,markerstrokewidth=1,frame=:axes)
    if isempty(uids)
        group = uids
    else
        group = flatspiketrains(uids;trialorder)[1]
    end
    plotspiketrain(flatspiketrains(sts;trialorder)[1:2]...;group,timeline,timespan,color,title,size,spancolor,
                    linewidth,linecolor,linestyle,markershape,markersize,markerstrokewidth,frame)
end

"Plot `Mean` and `SEM` of repeated responses for each condition"
function plotcondresponse(mseuc::DataFrame;group=:u,color=:auto,style=:path,projection=:none,title="",grid=false,
                            linestyle=:solid,linewidth=:auto,legend=:best,response=[],responsetype=:Response)
    ug = unique(mseuc.u)
    nug = length(ug)
    if nug>1
        # sort to the group order plot uses
        gi = sortperm(ug)
        color isa Vector && (color=permutedims(color[gi]))
        linewidth isa Vector && (linewidth=permutedims(linewidth[gi]))
        linestyle isa Vector && (linestyle=permutedims(linestyle[gi]))
    end
    factor = condfactor(mseuc)
    nfactor = length(factor)
    if nfactor==1
        factor=factor[1]
        # if typeof(mseuc[!,factor][1]) <: Array
        #     map!(string,mseuc[factor],mseuc[factor])
        #     style=:bar
        # end
    elseif nfactor==2
        p = nothing
        if nug==1
            fm,fa = factorspace(mseuc;col=:m)
            ab = skipmissing(fm) |> extrema .|> abs |> maximum
            lb = skipmissing(fm) |> minimum
            yfactor,xfactor = collect(keys(fa))
            y,x = collect(values(fa))
            clims = lb < 0 ? (-ab,ab) : (0,ab)
            p = heatmap(x,y,fm;color,title,legend,xlabel=factorunit(xfactor),ylabel=factorunit(yfactor),tickdir=:out,
                        colorbar_title=factorunit(responsetype),clims)
        else
            @info "Plot >1 unit group response is not implemented for 2D heatmap"
        end
        return p
    else
        mseuc.Condition = condstring(mseuc[!,factor])
        factor=:Condition
        style=:bar
    end

    if projection==:polar
        # close curve
        c0 = mseuc[mseuc[!,factor].==0,:]
        c0[!,factor].=360
        mseuc = [mseuc;c0]
        mseuc[!,factor] = deg2rad.(mseuc[!,factor])
        be = backend_name()
        if be == :gr
            legend=(0.86,0.96)
        elseif be == :pythonplot
            legend=(0.91,0.81)
        end
        sort!(mseuc,factor)
        if be == :gr
            p = @df mseuc plot(cols(factor),:m;yerror=:se,group=cols(group),line=style,markerstrokecolor=color,color,
            grid,projection,legend,xlabel=factorunit(factor),ylabel=factorunit(responsetype),title,linestyle,linewidth)
        elseif be == :pythonplot
            p = @df mseuc plot(cols(factor),:m;ribbon=:se,group=cols(group),line=style,markerstrokecolor=color,color,
            grid,projection,legend,xticks=0:0.25π:1.75π,xformatter=x->round(Int,rad2deg(x)),title,linestyle,linewidth)
        end
    else
        sort!(mseuc,factor)
        p = @df mseuc plot(cols(factor),:m;ribbon=:se,group=cols(group),line=style,markerstrokecolor=color,color,
        grid,projection,legend,xlabel=factorunit(factor),ylabel=factorunit(responsetype),title,linestyle,linewidth)
    end
    if !isempty(response)
        for i in response
            hline!(p,[i[1]],ribbon=[i[2]],color=colors,legend=false)
        end
    end
    p
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

"Plot Analog Signals"
function plotanalog(data;x=nothing,y=nothing,fs=0,xext=0,timeline=[],xlabel="Time",ylabel="Depth",clims=nothing,
    xunit=:ms,yunit=:μm,cunit=:v,plottype=:heatmap,hy=0,color=:coolwarm,layer=nothing,n=nothing,aspectratio=:auto)
    nd=ndims(data);ds=1
    if nd==1
        if isnothing(x)
            x=1:length(data)
            if fs>0
                x = x./fs.-xext
                xunit==:ms && (x*=1000)
            end
        end
        if isnothing(clims)
            if cunit==:v
                lim = maximum(abs.(data))
                clims = (-lim,lim)
            elseif cunit==:uv
                ds = 1e6
                lim = maximum(abs.(data))*ds
                clims = (-lim,lim)
            elseif cunit == :fr
                lim = maximum(data)
                clims = (0,lim)
            else
                clims = :auto
            end
        end
        xlabel = factorunitstring(xlabel,xunit)
        ylabel = factorunitstring(ylabel,yunit)
        p=plot(x,data*ds;legend=false,color_palette=color,grid=false,ylims=clims,xlabel,ylabel,aspectratio)
    elseif nd==2
        if isnothing(x)
            x=1:size(data,2)
            if fs>0
                x = x./fs.-xext
                xunit==:ms && (x*=1000)
            end
        end
        if isnothing(clims)
            if cunit==:v
                lim = maximum(abs.(data))
                clims = (-lim,lim)
            elseif cunit==:uv
                ds = 1e6
                lim = maximum(abs.(data))*ds
                clims = (-lim,lim)
            elseif cunit == :fr
                lim = maximum(data)
                clims = (0,lim)
            else
                clims = :auto
            end
        end
        xlabel = factorunitstring(xlabel,xunit)
        ylabel = factorunitstring(ylabel,yunit)
        if plottype==:heatmap
            if isnothing(y)
                y = (1:size(data,1))
                hy>0 && (y*=hy)
            end
            p=heatmap(x,y,data*ds;color=color,clims=clims,xlabel,ylabel,aspectratio,xlim=extrema(x),ylim=extrema(y))
            if !isnothing(n)
                pn = clampscale(n) .* (maximum(x)-minimum(x)) .* 0.2 .+ minimum(x)
                plot!(p,pn,y,color=:seagreen,leg=false,xlim=extrema(x),ylim=extrema(y))
            end
        else
            p=plot(x,data'*ds;legend=false,color_palette=color,grid=false,ylims=clims,xlabel,ylabel,aspectratio)
        end
    end
    isempty(timeline) || vline!(p,timeline,line=(:grey),label="TimeLine",leg=false)
    if !isnothing(layer)
        anno = [(minimum(x)+3,mean(layer[k]),text(k,6,:gray10,:center,:left)) for k in keys(layer)]
        hline!(p,[l[2] for l in values(layer)],linecolor=:gray25,legend=false,lw=0.5,annotations=anno)
    end
    return p
end

plotunitposition(spike::Dict;layer=nothing,color=nothing,alpha=0.4,title="",markersize=5,unitidsize=3,size=(600,450)) = plotunitposition(spike["unitposition"];unitgood=spike["unitgood"],chposition=spike["chposition"],unitid=spike["unitid"],layer,color,alpha,title,markersize,unitidsize,size)
function plotunitposition(unitposition;unitgood=[],chposition=[],unitid=[],layer=nothing,color=nothing,alpha=0.4,title="",markersize=5,unitidsize=3,size=(600,450))
    nunit = Base.size(unitposition,1);ngoodunit = isempty(unitgood) ? nunit : count(unitgood);ustr = "$ngoodunit / $nunit"
    xlims = isempty(chposition) ? (minimum(unitposition[:,1])-4,maximum(unitposition[:,1])+2) : (minimum(chposition[:,1])-5,maximum(chposition[:,1])+5)
    p = plot(;legend=:topright,xlabel="X (μm)",ylabel="Y (μm)",xlims,size,tickdir=:out,leftmargin=4Plots.mm)
    if !isempty(chposition)
        scatter!(p,chposition[:,1],chposition[:,2],markershape=:rect,markerstrokewidth=0,markersize=1.5,color=:grey70,label="Electrode")
    end
    if isnothing(color)
        if isempty(unitgood)
            color = :gray30
        else
            color = map(i->i ? :darkgreen : :gray30,unitgood)
        end
    end
    if !isnothing(layer)
        xmin,xmax = xlims
        ann = [(xmin+0.02(xmax-xmin),mean(layer[k]),text(k,6,:gray10,:left,:vcenter)) for k in keys(layer)]
        hline!(p,[l[2] for l in values(layer)];linecolor=:gray25,leg=false,lw=0.5,ann)
    end
    if isempty(unitid)
        series_annotations=[]
    else
        series_annotations=text.(unitid,unitidsize,:gray10,:center)
    end
    scatter!(p,unitposition[:,1],unitposition[:,2];label=ustr,color,alpha,msw=0,markersize,series_annotations,title)
    return p
end

# function plotunitpositionproperty(unitposition;ori=nothing,os=nothing,dir=nothing,ds=nothing,sf=nothing,width=500,height=400,title="",layer=nothing)
#     df = DataFrame(x=unitposition[:,1],y=unitposition[:,2],m=Any[:circle for _ in 1:size(unitposition,1)],a=0.0,sw=0.0,c=0.0,s=10.0)
#     l = DataFrame()
#     xlim = [minimum(df[:,:x])-4,maximum(df[:,:x])+2]
#     ylim = [minimum(df[:,:y])-100,maximum(df[:,:y])+100]
#     if !isnothing(sf)
#         df[!,:c]=sf
#     end
#     if !isnothing(ori)
#         df[!,:m] .= :stroke
#         df[!,:a] = -ori
#         df[!,:sw] .= 1
#         if !isnothing(os)
#             df[!,:s]=os
#         else
#             df[!,:s].=160
#         end
#     end
#     if !isnothing(dir)
#         t = copy(df)
#         arrowpath = "M 0 -0.1 H 1 V -0.3 L 1.6 0 L 1 0.3 V 0.1 H 0 Z"
#         t[!,:m] .= arrowpath
#         t[!,:a] = -dir
#         t[!,:sw] .= 0
#         if !isnothing(ds)
#             t[!,:s]=ds
#         else
#             t[!,:s].=200
#         end
#         df = [df;t]
#     end
#     if !isnothing(layer)
#         l[!,:x] = fill(xlim[1],length(layer))
#         l[!,:x2] .= xlim[2]
#         l[!,:y] = [v[1] for v in values(layer)]
#         l[!,:y2] = l[!,:y]
#         l[!,:l] = collect(keys(layer))
#         ylim = [extrema([ylim;l[:,:y]])...].+[-50,100]
#     end
#     @vgplot(height=height,width=width,padding=5,data=[:df=>df,:l=>l],
#     marks=[
#     {
#     type="rule",
#     from={data="l"},
#     encode={
#         update={
#         x={field="x",scale="x"},
#         y={field="y",scale="y"},
#         x2={field="x2",scale="x"},
#         y2={field="y2",scale="y"},
#         strokeDash={value=[4,2]},
#         strokeWidth={value=0.5},
#         stroke={value="dimgray"}
#         }}
#     },
#     {
#     type="text",
#     from={data="l"},
#     encode={
#         update={
#         x={field="x",scale="x"},
#         y={field="y",scale="y"},
#         text={field="l"},
#         align={value="left"},
#         baseline={value="bottom"},
#         fontSize={value=7},
#         dx={value=15}
#         }}
#     },
#     {
#     type="symbol",
#     from={data="df"},
#     encode={
#         update={
#         x={field="x",scale="x"},
#         y={field="y",scale="y"},
#         shape={field="m"},
#         angle={field="a"},
#         size={field="s",scale="s"},
#         strokeWidth={field="sw"},
#         stroke={field="c",scale="c"},
#         fill={field="c",scale="c"}
#         }}
#     }
#     ],
#     scales=[
#     {
#         name="x",
#         nice=false,
#         zero=false,
#         range="width",
#         domain=xlim,
#         type="linear",
#         round=true
#     },
#     {
#         name="y",
#         nice=false,
#         zero=false,
#         range="height",
#         domain=ylim,
#         type="linear",
#         round=true
#     },
#     {
#         name="c",
#         nice=false,
#         zero=false,
#         round=false,
#         type="linear",
#         range={scheme = "plasma",extent=[0.8,0.2]},
#         domain={data="df",field="c"}
#     },
#     {
#         name="s",
#         nice=false,
#         zero=false,
#         domain={data="df",field="s"},
#         type="linear",
#         round=true,
#         range=[10, 200]
#     }
#     ],
#     axes=[
#     {
#         domain=true,
#         tickCount=5,
#         grid=false,
#         title="Position_X (μm)",
#         scale="x",
#         orient="bottom"
#     },
#     {
#         domain=true,
#         tickCount=5,
#         grid=false,
#         title="Position_Y (μm)",
#         scale="y",
#         orient="left"
#     }
#     ],
#     title={
#     text = title
#     },
#     legends=[
#     {
#     type="gradient",
#     fill="c",
#     title="SF"
#     },
#     {
#     type="symbol",
#     symbolType="stroke",
#     size="s",
#     title="1-CV"
#     }
#     ])
# end

# function plotunitpositionimage(unitposition,image;width=800,height=600,markersize=20,title="",layer=nothing)
#     tempimagedir = "tempimage";nu=size(unitposition,1)
#     isdir(tempimagedir) || mkpath(tempimagedir)
#     urls = map(i->"tempimage/$(i).png",1:nu)
#     foreach(i->save(urls[i],image[i]),1:nu)

#     df = DataFrame(x=unitposition[:,1],y=unitposition[:,2],m=urls,s=markersize)
#     l = DataFrame()
#     # xlim = [minimum(df[:,:x])-4,maximum(df[:,:x])+2]
#     # ylim = [minimum(df[:,:y])-100,maximum(df[:,:y])+100]
#     xlim = [minimum(df[:,:x])-0.2,maximum(df[:,:x])+0.2]
#     ylim = [minimum(df[:,:y])-0.2,maximum(df[:,:y])+0.2]
#     # xlim = [minimum(df[:,:x])-2,maximum(df[:,:x])+2]
#     # ylim = [minimum(df[:,:y])-2,maximum(df[:,:y])+2]
#     if !isnothing(layer)
#         l[!,:x] = fill(xlim[1],length(layer))
#         l[!,:x2] .= xlim[2]
#         l[!,:y] = [v[1] for v in values(layer)]
#         l[!,:y2] = l[!,:y]
#         l[!,:l] = collect(keys(layer))
#         ylim = [extrema([ylim;l[:,:y]])...].+[-50,100]
#     end
#     @vgplot(height=height,width=width,padding=5,data=[:df=>df,:l=>l],
#     marks=[
#     {
#     type="rule",
#     from={data="l"},
#     encode={
#         update={
#         x={field="x",scale="x"},
#         y={field="y",scale="y"},
#         x2={field="x2",scale="x"},
#         y2={field="y2",scale="y"},
#         strokeDash={value=[4,2]},
#         strokeWidth={value=0.5},
#         stroke={value="dimgray"}
#         }}
#     },
#     {
#     type="text",
#     from={data="l"},
#     encode={
#         update={
#         x={field="x",scale="x"},
#         y={field="y",scale="y"},
#         text={field="l"},
#         align={value="left"},
#         baseline={value="bottom"},
#         fontSize={value=7},
#         dx={value=15}
#         }}
#     },
#     {
#     type="image",
#     from={data="df"},
#     encode={
#         update={
#         x={field="x",scale="x"},
#         y={field="y",scale="y"},
#         url={field="m"},
#         width={field="s"},
#         height={field="s"},
#         align={value="center"},
#         baseline={value="middle"}
#         }}
#     }
#     ],
#     scales=[
#     {
#         name="x",
#         nice=false,
#         zero=false,
#         range="width",
#         domain=xlim,
#         type="linear",
#         round=true
#     },
#     {
#         name="y",
#         nice=false,
#         zero=false,
#         range="height",
#         domain=ylim,
#         type="linear",
#         round=true
#     }
#     ],
#     axes=[
#     {
#         domain=false,
#         tickCount=0,
#         grid=false,
#         # title="Position_X (μm)",
#         title="",
#         scale="x",
#         orient="bottom"
#     },
#     {
#         domain=false,
#         tickCount=0,
#         grid=false,
#         # title="Position_Y (μm)",
#         title="",
#         scale="y",
#         orient="left"
#     }
#     ],
#     title={
#     text = title
#     }
#     )
# end

# function plotunitlayerimage(unitlayer,image;width=800,height=600,markersize=40,title="",unitid=nothing)
#     tempimagedir = joinpath(pwd(),"tempimage");nu=length(unitlayer)
#     isdir(tempimagedir) || mkpath(tempimagedir)
#     urls = map(i->"tempimage/$(i).png",1:nu)
#     foreach(i->save(urls[i],image[i]),1:nu)

#     x = zeros(nu)
#     for r in eachrow(condin(DataFrame(l=unitlayer)))
#         x[r.i] = 1:r.n
#     end
#     df = DataFrame(x=x,y=unitlayer,m=urls)
#     if isnothing(unitid)
#         df[:u] = ""
#     else
#         df[:u] = unitid
#     end

#     @vlplot(height=height,width=width,padding=5,title=title,data=df,
#     y={"y:o",axis={title="Layer"}},
#     x={"x",axis={grid=false,title=nothing,ticks=false,domain=false,labels=false}}) +
#     @vlplot(mark={:image,width=markersize,height=markersize,align=:center,baseline=:middle}, url=:m) +
#     @vlplot(mark={:text,align=:center,baseline=:bottom,dy=30,fontSize=7},text=:u)
# end

function plotcircuit(unitposition,unitid,projs;unitgood=[],projtypes=[],projweights=[],layer=nothing,
    showuid=true,alpha=0.4,title="",markersize=6,unitidsize=3,size=(600,800),showunit=:all,plw=0.4,pcolor=:gray50)
    np=length(projs);nu = length(unitid)
    nsu = isempty(unitgood) ? nu : count(unitgood)
    uidi = Dict(unitid[i]=>i for i in eachindex(unitid))
    nsupair=binomial(nsu,2);annstr = "$nsu/$nu, $np/$nsupair($(round(np/nsupair*100,digits=3))%)"
    p = plot(;legend=false,xlabel="X (μm)",ylabel="Y (μm)",grid=false,size,title,tickdir=:out,leftmargin=4Plots.mm)

    pii = map(p->uidi[p[1]],projs);pji = map(p->uidi[p[2]],projs)
    px = [unitposition[pii,1]';unitposition[pji,1]']
    py = [unitposition[pii,2]';unitposition[pji,2]']
    if !isempty(projweights)
        w = clampscale(abs.(projweights),3) .+ 0.4
        plw = (plw.*w)'
    end
    plot!(p,px,py,lw=plw,color=pcolor,arrow=arrow(:closed,:head,0.46,0.14))

    color = fill(:gray30,nu)
    color[unitgood] .= :darkgreen
    if !isempty(projtypes)
        color[pii] = replace(projtypes,true=>:darkred,false=>:darkblue)
    end

    if showunit == :su
        showui = unitgood
    elseif showunit == :circuit
        cuid = mapreduce(p->[p...],append!,projs)
        showui=map(i->uidi[i],unique!(cuid))
    elseif showunit == :all
        showui = 1:nu
    end

    showunitposx = unitposition[showui,1]
    showunitposy = unitposition[showui,2]
    xmin,xmax = extrema(showunitposx);xr=xmax-xmin
    xmin -= 0.08xr;xmax += 0.02xr
    annx = xmin+0.02xr
    xlims=(xmin,xmax)
    ymin,ymax = extrema(showunitposy)
    if !isnothing(layer)
        ann = [(annx,mean(layer[k]),text(k,6,:gray10,:left,:vcenter)) for k in keys(layer)]
        ls = [l[2] for l in values(layer)]
        hline!(p,ls;linecolor=:gray25,lw=0.5,ann)
        ymin = min(ymin,minimum(ls))
        ymax = max(ymax,maximum(ls))
    end
    yr = ymax-ymin
    ymin -= 0.05yr;ymax += 0.02yr
    anny = ymin+0.02yr
    ylims=(ymin,ymax)

    if showuid
        series_annotations=text.(unitid[showui],unitidsize,:gray10,:center)
    else
        series_annotations=[]
    end
    scatter!(p,showunitposx,showunitposy;color=color[showui],alpha,msw=0,markersize,series_annotations,xlims,ylims)
    annotate!(p,[(annx,anny,text(annstr,9,:gray10,:left,:vcenter))])
    p
end

"plot hartley subspace"
function plothartleysubspace(ps,nk,dk;color=:grays)
    n = 2nk+1
    p=plot(layout=(n,n),size=(120n,120n),leg=false,clims=(-1,1),frame=:none,aspect_ratio=:equal)
    xy=0:0.01:1
    for (kx,ky,phase) in ps
        r = -round(Int,ky/dk)+nk+1; c = round(Int,kx/dk)+nk+1
        cg = [cas(i,j,kx=kx,ky=ky,phase=phase) for j in xy,i in xy]
        heatmap!(p,subplot=c+n*(r-1),cg,color=color)
    end
    p
end
