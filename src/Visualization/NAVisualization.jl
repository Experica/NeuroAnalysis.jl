module NAVisualization
using NeuroAnalysis.NACore, NeuroAnalysis.NABase

include("Plot.jl")

export savefig

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

function savefig(fig::Plot,filename::AbstractString;path::AbstractString="",format::AbstractString="svg",width=22cm,height=13cm,dpi=300)
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

end # module
