using Test, NeuroAnalysis, BenchmarkTools, DataFrames, Plots, FileIO

@testset "NeuroAnalysis" begin

    include("spiketest.jl")
    include("imagetest.jl")
    include("conditiontest.jl")
    include("regressiontest.jl")
    include("visualizationtest.jl")
    



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

## functions
plot(vmf,-π,π)
plot(gvmf,-π,π)
plot(gratingf,-3,3)
plot(gaussianf,-3,3)
plot(gaborf,-3,3)
plot(dogf,-3,3)

funnames = ["grating","gaussian","gabor","dog"]
plot([gratingf,gaussianf,gaborf,dogf],-3,3,labels=permutedims(funnames),lw=[2 2 3 3],xlabel="y′")


x=y=-3:0.05:3;z=[]
push!(z,[gratingf(i,j,θ=0.25π,f=0.5) for j in y,i in x])
push!(z,[gaussianf(i,j,σ₁=0.8,σ₂=0.5,θ=0.25π) for j in y,i in x])
push!(z,[gaborf(i,j,σ₁=0.8,σ₂=0.5,θ=0.25π,f=0.5) for j in y,i in x])
push!(z,[dogf(i,j,σₑ₁=0.8,σₑ₂=0.5,σᵢ₁=1,σᵢ₂=0.7,θₑ=0.25π,θᵢ=0.25π) for j in y,i in x])

p = plot(layout=(2,2),legend=false,size=(600,600))
foreach(i->heatmap!(p,z[i],subplot=i,aspect_ratio=:equal,frame=:none,color=:coolwarm,clims=(-1,1),title=funnames[i]),1:4)
p


end
