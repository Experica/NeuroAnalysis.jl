using Gadfly
export plotspiketrain,plotpsth


function plotspiketrain(sts::TPsVector;timemark=[],theme=Theme())
  x=[];y=[]
  xmin=0;xmax=0
  for i in 1:length(sts)
    cst = sts[i]
    xmin = minimum([xmin,minimum(cst)])
    xmax = maximum([xmax,maximum(cst)])
    x = [x,cst]
    y = [y,ones(length(cst))*i]
  end
  plot(x=x,y=y,xintercept=timemark,theme,Geom.point,Geom.vline(color="gray",size=1pt),
       Coord.Cartesian(xmin=xmin,xmax=xmax),Guide.xlabel("Time (ms)"),Guide.ylabel("Trials"))
end

function plotpsth(sts::TPsVector,binedges::TimePoints;theme=Theme())
  vhist,wins,vsubs,vsis = histtps(sts,binedges)
  hm,x = histmatrix(vhist,wins)
  y = mean(hm,1)
  plot(y=y,x=x,theme,Geom.bar,Guide.xlabel("Time (ms)"),Guide.ylabel("Counts"))
end
