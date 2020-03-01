using Documenter,NeuroAnalysis

username = "Experica"
pkgname = "NeuroAnalysis"

makedocs(sitename="$pkgname.jl",modules=[NeuroAnalysis],
    pages=[])

deploydocs(repo = "github.com/$username/$pkgname.jl.git")
