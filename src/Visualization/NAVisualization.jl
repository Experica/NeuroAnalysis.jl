module NAVisualization

using NeuroAnalysis.NACore, NeuroAnalysis.NABase

include("Plot.jl")

export savefig

function savefig(fig,fname::String;path::String="",format::String="svg")
  fn = joinpath(path,"$fname.$format")
  if !ispath(path)
    mkpath(path)
  end
  if format=="svg"
    format = "$format+xml"
  end
  open(fn, "w") do io
    writemime(io,"image/$format",fig)
  end
end

function savefig(fig::Plot,fname::String;path::String="",format::String="svg",width=22cm,height=13cm)
  fn = joinpath(path,"$fname.$format")
  if !ispath(path)
    mkpath(path)
  end
  if format=="svg"
    format = SVG(fn,width,height)
  elseif format=="png"
    format = PNG(fn,width,height)
  elseif format=="pdf"
    format = PDF(fn,width,height)
  end
  draw(format,fig)
end

end # module
