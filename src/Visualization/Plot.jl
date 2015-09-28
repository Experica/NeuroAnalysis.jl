using Gadfly,DataFrames
export plotspiketrain,plotpsth


<<<<<<< HEAD
function plotspiketrain(x::Vector,y::Vector,c::Vector=[];xgroup::Vector=[],timemark=[0],theme=Theme(),colorkey="",colorfun=Scale.lab_gradient(color("white"),color("red")),colorminv=[],colormaxv=[])
  xl="Time (ms)";yl="Trial"
  if isempty(c)
    if isempty(xgroup)
      plot(x=x,y=y,xintercept=timemark,theme,Geom.point,Geom.vline(color="gray",size=1pt),
           Guide.xlabel(xl),Guide.ylabel(yl))
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
           Guide.xlabel(xl),Guide.ylabel(yl),Guide.colorkey(colorkey),
           Scale.ContinuousColorScale(colorfun,minvalue=colorminv,maxvalue=colormaxv))
    else
      plot(x=x,y=y,color=c,xgroup=xgroup,xintercept=fill(timemark[1],length(x)),theme,
           Geom.subplot_grid(Geom.point,Geom.vline(color="gray",size=1pt),free_x_axis=true),
           Guide.xlabel(xl),Guide.ylabel(yl),Guide.colorkey(colorkey),
           Scale.ContinuousColorScale(colorfun,minvalue=colorminv,maxvalue=colormaxv))
    end
  end
end

plotspiketrain(rvs::RVVector;sortvar=[],xgroup::Vector=[],timemark=[0],theme=Theme(),colorkey="",colorfun=Scale.lab_gradient(color("white"),color("red")),colorminv=[],colormaxv=[]) = plotspiketrain(flatrv(rvs,sortvar)...,xgroup=xgroup,timemark=timemark,theme=theme,colorkey=colorkey,colorfun=colorfun,colorminv=colorminv,colormaxv=colormaxv)

function plotpsth(rvs::RVVector,binedges::RealVector;theme=Theme())
  m,sd,n,x = psth(rvs,binedges)
  plot(y=m,x=x,theme,Geom.line,Guide.xlabel("Time (ms)"),Guide.ylabel("Response (spike/s)"))
end

=======
function plotspiketrain(rvs::RVVector;timemark=[],theme=Theme(),sorttrial::Bool=false,sortparam=[],sortparamname="")
  rvn = length(rvs)
  colorsort=false
  yl = "Trial"
  if length(sortparam) > 0
    if length(sortparam) == rvn
      colorsort = true
    else
      warn("Lengths of rvs and sortparam do not match.")
    end
  end
  if sorttrial && colorsort
    rvs=rvs[sortperm(sortparam)]
    sortparam=sort(sortparam)
    yl = "$yl Sorted"
  end
  x=[];y=[];s=[]
  for i in 1:rvn
    rv = rvs[i]
    if isempty(rv);continue;end
    n=length(rv)
    x = [x,rv]
    y = [y,ones(n)*i]
    if colorsort
      s=[s,fill(sortparam[i],n)]
    end
  end
  if colorsort
    plot(x=x,y=y,color=s,xintercept=timemark,theme,Geom.point,Geom.vline(color="gray",size=1pt),
         Guide.xlabel("Time (ms)"),Guide.ylabel(yl),Guide.colorkey(sortparamname))
  else
    plot(x=x,y=y,xintercept=timemark,theme,Geom.point,Geom.vline(color="gray",size=1pt),
         Guide.xlabel("Time (ms)"),Guide.ylabel(yl))
  end
end

function plotpsth(rvs::RVVector,binedges::RealVector;theme=Theme(),timemark=[0])
  m,sd,n,x = psth(rvs,binedges)
  plot(y=m,x=x,ymin=m-sd/sqrt(n),ymax=m+sd/sqrt(n),xintercept=timemark,theme,
    Geom.line,Geom.ribbon,Geom.vline(color="gray",size=1pt),Coord.Cartesian(xmin=binedges[1],xmax=binedges[end],ymin=0),
    Guide.xlabel("Time (ms)"),Guide.ylabel("Response (spike/s)"))
end

>>>>>>> 107f77da1a310ca40b52df72ae4d71ec825ec5d1
function plotpsth(ds::DataFrame,binedges::RealVector,conds::Vector{Vector{Any}};theme=Theme(),timemark=[0])
  df,ss = psth(ds,binedges,conds)
  plot(df,x=:x,y=:y,ymin=:ymin,ymax=:ymax,color=:condition,xintercept=timemark,
       theme,Geom.line,Geom.ribbon,Geom.vline(color="gray",size=1pt),
       Coord.Cartesian(xmin=binedges[1],xmax=binedges[end],ymin=0),Guide.xlabel("Time (ms)"),Guide.ylabel("Response (spike/s)"))
end
