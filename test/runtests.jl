using NeuroAnalysis, Test, BenchmarkTools, DataFrames, Plots

@testset "NeuroAnalysis" begin

    include("spiketest.jl")
    include("imagetest.jl")
    include("conditiontest.jl")


## 1D grating
# @manipulate for f in 0:0.1:1, p in 0:0.001:1
#     plot([x->gratingf(x,f=f,phase=p),x->cas(x,f=f,phase=sin2cas(p))],-2,2,labels=["sin" "cas"])
# end
## 2D grating
# @manipulate for θ in 0:0.01:π, f in 0:0.05:1, sinp in 0:0.001:1, kx in -1:0.05:1,ky in -1:0.05:1, casp in 0:0.001:1
#     x=y=-2:0.05:2
#     sg = [gratingf(i,j,θ=θ,f=f,phase=sinp) for j in reverse(y),i in x]
#     scg = [cas(i,j;sin2cas(θ,f,sinp)...) for j in reverse(y),i in x]
#     cg = [cas(i,j,kx=kx,ky=ky,phase=casp) for j in reverse(y),i in x]
#     csg = [gratingf(i,j;cas2sin(kx,ky,casp)...) for j in reverse(y),i in x]
#
#     p = plot(layout=(2,3),leg=false,size=(600,400),clims=(-1,1),frame=:none,aspect_ratio=:equal,yflip=true)
#     heatmap!(p,subplot=1,sg,color=:grays,title="sin")
#     heatmap!(p,subplot=2,scg,color=:grays,title="sin2cas")
#     heatmap!(p,subplot=3,sg.-scg,color=:grays,title="sin - sin2cas")
#     heatmap!(p,subplot=4,cg,color=:grays,title="cas")
#     heatmap!(p,subplot=5,csg,color=:grays,title="cas2sin")
#     heatmap!(p,subplot=6,cg.-csg,color=:grays,title="cas - cas2sin")
#     p
# end


## colormap
cgrad(RGB(1,0.0,0),RGB(0,1.0,0)).colors
cgrad(RGB(1,0.0,0),RGB(0,1.0,0),RGB(0,0.0,1)).colors

end
