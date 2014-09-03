using Gadfly
export plotspiketrain,plotpsth


function plotspiketrain(sts::TPsVector;theme=Theme())
  stsn = length(sts)
  ls = []
  for i=1:stsn
    csubst = sts[i]
    csubstn = length(csubst)
    y = ones(csubstn)*i
    ls = [ls,layer(x=csubst,y=y,theme,Geom.point)]
  end
  plot(ls,Guide.xlabel("Time (ms)"),Guide.ylabel("Trials"))
end

function plotpsth(psths::AbstractVector{AbstractVector{Integer}},bins::TPsVector;theme=Theme())
  parray,x = histarray(psths,bins)
  y = mean(parray,1)
  plot(y=y,x=x,theme,Geom.bar,Guide.xlabel("Time (ms)"),Guide.ylabel("Counts"))
end
