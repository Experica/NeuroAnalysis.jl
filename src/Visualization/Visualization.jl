export factorunit,unitcolor,plotspiketrain,plotpsth,plottuning,savefig

using Gadfly,Plots

function factorunit(f)
    fu=""
    if f=="Ori"
        fu="Deg"
      elseif f=="SpatialFreq"
        fu="Cycle/Deg"
      elseif f=="TemporalFreq"
      fu = "Cycle/Sec"
    end
    return fu
end

function unitcolor(;n=5,alpha=0.8,saturation=1,brightness=1)
    [HSVA(0,0,0,alpha) [HSVA(((i-1)/n)*360,saturation,brightness,alpha) for j=1:1, i=1:n]]
end
function plotspiketrain(x,y,g::Vector=[];timeline=[0],colorseq=unitcolor(),title=title)
    xl="Time (ms)";yl="Trial"
    if isempty(g)
    scatter(x,y,label="SpikeTrain",markershape=:vline,markersize=1,markerstrokecolor=RGBA(0.0,0.1,0.2,0.8),markerstrokewidth = 1)
    else
       scatter(x,y,group=g,markershape=:vline,markersize=1,markerstrokewidth = 1,
        markerstrokecolor=colorseq)
    end
        vline!(timeline,line=(:gray),label="TimeLine",xaxis=(xl),yaxis=(yl),title=(title))
end
plotspiketrain(sts::RVVector;sv=[],timeline=[0],colorseq=unitcolor(),title="")=plotspiketrain(flatrvs(sts,sv)...,timeline=timeline,colorseq=colorseq,title=title)
plotspiketrain(sts::RVVector,uids::RVVector;sv=[],timeline=[0],colorseq=unitcolor(),title="")=plotspiketrain(flatrvs(sts,sv)[1:2]...,map(i->"U$i",flatrvs(uids,sv)[1]),timeline=timeline,colorseq=colorseq,title=title)

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

function plotpsth(rvs::RVVector,binedges::RealVector;theme=Theme(),timeline=[0],title="")
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
function plotpsth(rvs::RVVector,binedges::RealVector,rvsidx,condstr;theme=Theme(),timeline=[0],title="")
    df = psth(rvs,binedges,rvsidx,condstr)
    df[:ymin] = df[:y]-df[:ysd]./sqrt(df[:n])
    df[:ymax] = df[:y]+df[:ysd]./sqrt(df[:n])
  xl = "Time (ms)";yl = "Response (spike/s)"
        Gadfly.plot(df,y=:y,x=:x,ymin=:ymin,ymax=:ymax,color=:condition,xintercept=timeline,theme,Geom.line,Geom.ribbon,
    Geom.vline(color="gray",size=1pt),Coord.Cartesian(xmin=binedges[1],xmax=binedges[end],ymin=0),
    Guide.xlabel(xl),Guide.ylabel(yl),Guide.title(title))
end

function plotpsth(ds::DataFrame,binedges::RealVector,conds::Vector{Vector{Any}};spike=:spike,theme=Theme(),timeline=[0])
  df,ss = psth(ds,binedges,conds,spike=spike)
  xl = "Time (ms)";yl = "Response (spike/s)"
  Gadfly.plot(df,x=:x,y=:y,ymin=:ymin,ymax=:ymax,color=:condition,xintercept=timeline,
       theme,Geom.line,Geom.ribbon,Geom.vline(color="gray",size=1pt),
       Coord.Cartesian(xmin=binedges[1],xmax=binedges[end],ymin=0),Guide.xlabel(xl),Guide.ylabel(yl))
end

function plottuning(rs,ridx,conds;title="")
    m,sd,n=condmean(rs,ridx,conds)
    nf = length(conds[1])
    if nf==1
        f=conds[1][1][1]
        xl="$f ($(factorunit(f)))"
        yl="Response (spike/s)"
        x=map(i->i[1][2],conds)
        Plots.plot(x,m,yerror=sd./sqrt(n),xaxis=(xl),yaxis=(yl),title=(title))
    end
end
function plottuning(rs,us,ridx,conds;title="")
    m,sd,n=condmean(rs,us,ridx,conds)
    nf = length(conds[1])
    if nf==1
        f=conds[1][1][1]
        xl="$f ($(factorunit(f)))"
        yl="Response (spike/s)"
        x=map(i->i[1][2],conds)
        Plots.plot(x,m,yerror=sd./sqrt(n),xaxis=(xl),yaxis=(yl),title=(title))
    end
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
