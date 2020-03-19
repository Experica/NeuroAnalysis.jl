using NeuroAnalysis, Test, BenchmarkTools

@testset "Spike" begin
## subrvr and subrvr_ono
spike = sort(rand(1:100000,1000))
on = collect(10:100:95000)
off = collect(90:100:95080)
n1 = subrvr(spike,on,off,israte=false)
n2 = subrvr_ono(spike,on,off,israte=false)
@test n1==n2
# @btime n1 = subrvr($(spike),$(on),$(off),israte=false)
# @btime n2 = subrvr_ono($(spike),$(on),$(off),israte=false)

end

@testset "Image" begin
## image and freqimage
ppd = 30;ori=0.5π;sf=2
img = grating(θ=ori,sf=sf,ppd=ppd)
# heatmap(img,yflip=true)
ps,f1,f2 = powerspectrum(img,ppd,freqrange=[-6,6])
# heatmap(f1,f2,ps)
eori,esf = freqimagestats(ps,f1,f2)
@test sf≈esf
@test ori≈eori
## hartley space
# @manipulate for kb in 0:4, ke in 1:5, dk in 0.4:0.2:1, p in 0:0.005:1, shape in [:square, :circle]
#     plothartleyspace(hartleysubspace(kbegin=kb,kend=ke,dk=dk,phase=p,shape=shape),floor(Int,ke/dk),dk)
# end

end

@testset "Function" begin
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

end
