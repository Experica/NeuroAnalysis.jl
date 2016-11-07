export unitcolor,plotspiketrain,plotpsth,savefig

using Gadfly,Plots

function unitcolor(;n=5,alpha=0.8,saturation=1,brightness=1)
    [HSVA(0,0,0,alpha) [HSVA(((i-1)/n)*360,saturation,brightness,alpha) for j=1:1, i=1:n]]
end
function plotspiketrain(x,y,g::Vector=[];timeline=[0],colorseq=unitcolor())
    xl="Time (ms)";yl="Trial"
    if isempty(g)
    scatter(x,y,label="SpikeTrain",markershape=:vline,markersize=1,markerstrokecolor=RGBA(0.0,0.1,0.2,0.8),markerstrokewidth = 1)
    else
       scatter(x,y,group=g,markershape=:vline,markersize=1,markerstrokewidth = 1,
        markerstrokecolor=colorseq)
    end
        vline!(timeline,line=(:gray),label="TimeLine",xaxis=(xl),yaxis=(yl))
end
plotspiketrain(sts::RVVector;sv=[],timeline=[0],colorseq=unitcolor())=plotspiketrain(flatrvs(sts,sv)...,timeline=timeline,colorseq=colorseq)
plotspiketrain(sts::RVVector,uids::RVVector;sv=[],timeline=[0],colorseq=unitcolor())=plotspiketrain(flatrvs(sts,sv)[1:2]...,map(i->"U$i",flatrvs(uids,sv)[1]),timeline=timeline,colorseq=colorseq)

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

function plotpsth(rvs::RVVector,binedges::RealVector;theme=Theme(),timeline=[0])
  m,sd,n,x = psth(rvs,binedges)
  if n>1
    Gadfly.plot(y=m,x=x,ymin=m-sd/sqrt(n),ymax=m+sd/sqrt(n),xintercept=timeline,theme,Geom.line,Geom.ribbon,
    Geom.vline(color="gray",size=1pt),Coord.Cartesian(xmin=binedges[1],xmax=binedges[end],ymin=0),Guide.xlabel("Time (ms)"),Guide.ylabel("Response (spike/s)"))
  else
    Gadfly.plot(y=m,x=x,xintercept=timeline,theme,Geom.line,Geom.vline(color="gray",size=1pt),
    Coord.Cartesian(xmin=binedges[1],xmax=binedges[end],ymin=0),Guide.xlabel("Time (ms)"),Guide.ylabel("Response (spike/s)"))
  end
end

function plotpsth(ds::DataFrame,binedges::RealVector,conds::Vector{Vector{Any}};spike=:spike,theme=Theme(),timemark=[0])
  df,ss = psth(ds,binedges,conds,spike=spike)
  Gadfly.plot(df,x=:x,y=:y,ymin=:ymin,ymax=:ymax,color=:condition,xintercept=timemark,
       theme,Geom.line,Geom.ribbon,Geom.vline(color="gray",size=1pt),
       Coord.Cartesian(xmin=binedges[1],xmax=binedges[end],ymin=0),Guide.xlabel("Time (ms)"),Guide.ylabel("Response (spike/s)"))
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
