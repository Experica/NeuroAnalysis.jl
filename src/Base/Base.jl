module Base

using NeuroAnalysis.Core

include("NeuroDataQuery.jl")
include("NeuroDataPrepare.jl")

export anscombe

function anscombe(x)
  2*sqrt(x+(3/8))
end

end # module
