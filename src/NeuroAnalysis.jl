module NeuroAnalysis

include("Base/Base.jl")
include("Visualization/Visualization.jl")
include("IO/IO.jl")

# export all symbols
for n in names(@__MODULE__, all=true)
    if Base.isidentifier(n) && n ∉ (nameof(@__MODULE__), :eval, :include)
        @eval export $n
    end
end

function __init__()
    f = joinpath(@__DIR__,"Visualization","Colors","colormaps.jld2")
    global ColorMaps = load(f)
    # global PyOnewayANOVA = pyimport("statsmodels.stats.oneway")
end

end
