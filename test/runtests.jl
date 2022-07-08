using Test, NeuroAnalysis, BenchmarkTools, DataFrames, Plots, FileIO, FFTW

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
push!(z,[gratingf(i,j,μ₁=0,μ₂=0.5,θ=0.25π,f=0.5) for j in y,i in x])
push!(z,[gaussianf(i,j,μ₁=0,μ₂=0.5,σ₁=0.8,σ₂=0.5,θ=0.25π) for j in y,i in x])
push!(z,[gaborf(i,j,μ₁=0,μ₂=0.5,σ₁=1,σ₂=0.5,θ=0.25π,f=0.5) for j in y,i in x])
push!(z,[dogf(i,j,aₑ=1,aᵢ=2,μₑ₁=0,μₑ₂=0.5,μᵢ₁=0,μᵢ₂=0.5,σₑ₁=1,σₑ₂=0.6,σᵢ₁=0.7,σᵢ₂=0.42,θₑ=0.25π,θᵢ=0.25π) for j in y,i in x])

p = plot(layout=(2,2),legend=false,size=(600,600))
foreach(i->heatmap!(p[i],z[i],aspect_ratio=:equal,frame=:none,color=:coolwarm,clims=(-1,1),title=funnames[i]),1:4)
p

z[2] = [gaussianenvelopemask(i,j;fσ=2.5,μ₁=0,σ₁=0.8,μ₂=0.5,σ₂=0.5,θ=0.25π) for j in y,i in x]
z[3] = [gaborenvelopemask(i,j;fσ=2.5,μ₁=0,σ₁=1,μ₂=0.5,σ₂=0.5,θ=0.25π) for j in y,i in x]
z[4] = [edogenvelopemask(i,j;fσ=2.5,μ₁=0,σₑ₁=1,rσ₂₁=0.6,μ₂=0.5,rσᵢₑ=0.7,θ=0.25π) for j in y,i in x]

foreach(i->heatmap!(p[i],z[i],aspect_ratio=:equal,frame=:none,alpha=0.1,color=:coolwarm,clims=(-1,1),title=funnames[i]),2:4)
p

t = 0:0.02:2π
z[2] = [gaussiancontour(i;fσ=2.5,μ₁=0,σ₁=0.8,μ₂=0.5,σ₂=0.5,θ=0.25π) for i in t]
z[3] = [gaborcontour(i;fσ=2.5,μ₁=0,σ₁=1,μ₂=0.5,σ₂=0.5,θ=0.25π) for i in t]
z[4] = [edogcontour(i;fσ=2.5,μ₁=0,σₑ₁=1,rσ₂₁=0.6,μ₂=0.5,rσᵢₑ=0.7,θ=0.25π) for i in t]

p = plot(layout=(2,2),legend=false,size=(600,600),frame=:none)
foreach(i->plot!(p[i],z[i],aspect_ratio=:equal,frame=:box,xlims=(-3,3),ylims=(-3,3),title=funnames[i]),2:4)
p

# dog pattern
x = -6:0.01:6
p=plot(layout=(3,3),leg=false,frame=:origin,yticks=[],xticks=[],link=:all,size=(600,500),leftmargin=4mm)
plot!(p[1,1],x->dogf(x,aₑ=2,σₑ=2,aᵢ=1,σᵢ=1),x,lw=2,ylabel="Aₑ/Aᵢ = 2")
plot!(p[1,2],x->dogf(x,aₑ=2,σₑ=1,aᵢ=1,σᵢ=1),x,lw=2)
plot!(p[1,3],x->dogf(x,aₑ=2,σₑ=1,aᵢ=1,σᵢ=2),x,lw=2)

plot!(p[2,1],x->dogf(x,aₑ=1,σₑ=2,aᵢ=1,σᵢ=1),x,lw=2,ylabel="Aₑ/Aᵢ = 1")
plot!(p[2,2],x->dogf(x,aₑ=1,σₑ=1,aᵢ=1,σᵢ=1),x,lw=2)
plot!(p[2,3],x->dogf(x,aₑ=1,σₑ=1,aᵢ=1,σᵢ=2),x,lw=2)

plot!(p[3,1],x->dogf(x,aₑ=1,σₑ=2,aᵢ=2,σᵢ=1),x,lw=2,xlabel="σᵢ/σₑ = 0.5",ylabel="Aₑ/Aᵢ = 0.5")
plot!(p[3,2],x->dogf(x,aₑ=1,σₑ=1,aᵢ=2,σᵢ=1),x,lw=2,xlabel="σᵢ/σₑ = 1")
plot!(p[3,3],x->dogf(x,aₑ=1,σₑ=1,aᵢ=2,σᵢ=2),x,lw=2,xlabel="σᵢ/σₑ = 2")
p

# gabor pattern
x = -3:0.01:3
p=plot(layout=(3,4),leg=false,frame=:origin,yticks=[],xticks=[],link=:all,size=(600,500),leftmargin=4mm)
plot!(p[1,1],x->gaborf(x,f=1.5/5,phase=0),x,lw=2,ylabel="cyc = 1.5")
plot!(p[1,2],x->gaborf(x,f=1.5/5,phase=0.25),x,lw=2)
plot!(p[1,3],x->gaborf(x,f=1.5/5,phase=0.5),x,lw=2)
plot!(p[1,4],x->gaborf(x,f=1.5/5,phase=0.75),x,lw=2)

plot!(p[2,1],x->gaborf(x,f=1/5,phase=0),x,lw=2,ylabel="cyc = 1")
plot!(p[2,2],x->gaborf(x,f=1/5,phase=0.25),x,lw=2)
plot!(p[2,3],x->gaborf(x,f=1/5,phase=0.5),x,lw=2)
plot!(p[2,4],x->gaborf(x,f=1/5,phase=0.75),x,lw=2)

plot!(p[3,1],x->gaborf(x,f=0.5/5,phase=0),x,lw=2,xlabel="p = 0",ylabel="cyc = 0.5")
plot!(p[3,2],x->gaborf(x,f=0.5/5,phase=0.25),x,lw=2,xlabel="p = 0.25")
plot!(p[3,3],x->gaborf(x,f=0.5/5,phase=0.5),x,lw=2,xlabel="p = 0.5")
plot!(p[3,4],x->gaborf(x,f=0.5/5,phase=0.75),x,lw=2,xlabel="p = 0.75")
p

# dog opponency
x = -10:0.01:10
plot([x->dogf(x,aₑ=2,σₑ=1,aᵢ=1,σᵢ=2),
    x->dogf(x,aₑ=4,σₑ=1,aᵢ=1,σᵢ=2),
    x->dogf(x,aₑ=2,σₑ=1,aᵢ=1,σᵢ=4)],x,lw=2,frame=:origin,fill=0,alpha=0.2,
    color=[:green :red :blue],
    label=["ρa=2, ρ=2" "ρa=4, ρ=2" "ρa=2, ρ=4"])


# dft
fs = 20
f = 2
t = 0:1/fs:10
x = cos.(2π*(f*t))
plot(t,x,xlabel="Second")
Fs = dft(x,fs,0,f)
FF = rfft(x)
plot(abs.(FF))

i = round(Int,f*length(x)/fs) + 1
@test FF[i] ≈ Fs[2]
@test FF[1] ≈ Fs[1]

end
