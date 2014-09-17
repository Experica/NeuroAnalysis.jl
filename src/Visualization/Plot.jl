using Gadfly
export plotspiketrain,plotpsth


function plotspiketrain(sts::TPsVector;theme=Theme())
  stn = length(sts)
  ls = []
  for i in 1:stn
    cst = sts[i]
    cstn = length(cst)
    y = ones(cstn)*i
    ls = [ls,layer(x=cst,y=y,theme,Geom.point)]
  end
  plot(ls,Guide.xlabel("Time (ms)"),Guide.ylabel("Trials"))
end

function plotpsth(sts::TPsVector,binedges::TimePoints;theme=Theme())
  vhist,wins,vsubs,vsis = histtps(sts,binedges)
  hm,x = histmatrix(vhist,wins)
  y = mean(hm,1)
  plot(y=y,x=x,theme,Geom.bar,Guide.xlabel("Time (ms)"),Guide.ylabel("Counts"))
end
