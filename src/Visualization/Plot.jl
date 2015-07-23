using Gadfly,DataFrames
export plotspiketrain,plotpsth


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

function plotpsth(ds::DataFrame,binedges::RealVector,conds::Vector{Vector{Any}};theme=Theme(),timemark=[0])
  df,ss = psth(ds,binedges,conds)
  plot(df,x=:x,y=:y,ymin=:ymin,ymax=:ymax,color=:condition,xintercept=timemark,
       theme,Geom.line,Geom.ribbon,Geom.vline(color="gray",size=1pt),
       Coord.Cartesian(xmin=binedges[1],xmax=binedges[end],ymin=0),Guide.xlabel("Time (ms)"),Guide.ylabel("Response (spike/s)"))
end
