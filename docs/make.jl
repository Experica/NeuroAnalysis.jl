using Documenter,NeuroAnalysis

uname = "Experica"
pname = "NeuroAnalysis"
makedocs(sitename="$pname.jl",modules=[NeuroAnalysis],
    pages=[])

deploydocs(repo = "github.com/$uname/$pname.jl.git")
