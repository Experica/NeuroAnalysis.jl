# hartley space
# @manipulate for kb in 0:4, ke in 1:5, dk in 0.4:0.2:1, p in 0:0.005:1, shape in [:square, :circle]
#     plothartleysubspace(hartleysubspace(kbegin=kb,kend=ke,dk=dk,phase=p,shape=shape),floor(Int,ke/dk),dk,color=:coolwarm)
# end

hs = hartleysubspace(kbegin=0.2,kend=6.6,dk=0.2,addhalfcycle=true)
hgs = map(i -> begin
          ss = cas2sin(i...)
          grating(θ = ss.θ, sf = ss.f, phase = ss.phase, size = (5, 5), ppd = 30)
          end, hs)
# mean of gratings[0, 1] of a hartley subspace should be a uniform gray
minv,maxv = extrema(reduce((i,j)->i.+j,hgs)/length(hs))
@test minv≈maxv≈0.5

hs = hartleysubspace(kbegin=0.2,kend=6.6,dk=0.2,addhalfcycle=true)
hgs = map(i -> begin
          ss = cas2sin(i...)
          grating(θ = ss.θ, sf = ss.f, phase = ss.phase, size = (2, 2), ppd = 10)
          end, hs)
Xᵀ = mapreduce(vec,hcat,hgs)
X = Xᵀ'
a=Xᵀ*X
b=inv(Xᵀ*X)


ggs = [randn(30,30) for _ in 1:10000]
Xᵀ = mapreduce(vec,hcat,ggs)
X = Xᵀ'
a=Xᵀ*X
b=inv(Xᵀ*X)



heatmap(a,yflip=true,aspect_ratio=1,frame=:none,color=:fire)

## grating
ori=0.5π;sf=2;ppd=30
img = grating(θ=ori,sf=sf,ppd=ppd)
img .+= rand(size(img)...)/3
heatmap(img,yflip=true,aspect_ratio=:equal,frame=:none)
ps,f1,f2 = powerspectrum2(img,ppd,freqrange=[-6,6])
heatmap(f2,f1,ps,aspect_ratio=:equal,frame=:grid)
eori,esf = f1orisf(ps,f1,f2)
@test sf≈esf
@test ori≈eori


## image resize
img = grating(θ=0.25π,sf=6,ppd=60,size=(3,3))
heatmap(img,yflip=true,aspect_ratio=:equal,frame=:none)

imgra = imresize_antialiasing(img,(64,64))
heatmap(imgra,yflip=true,aspect_ratio=:equal,frame=:none)

## image mask and blend
@test alphablend(Gray(0.7),Gray(0.5)) == Gray(0.7)
@test alphablend(GrayA(0.7,0),Gray(0.5)) == Gray(0.5)
@test alphablend(RGBA(0.5,1,1,1),RGB(0.6,0.3,0.2)) == RGB(0.5,1,1)

@test alphablend(10,5,1) == 10
@test alphablend(fill(10,10,10),fill(5,10,10),0.5) == fill(7.5,10,10)

alphablend.(alphamask(Gray.(img)).y,Gray(0.5))
alphablend.(alphamask(Gray.(img),masktype="Gaussian",sigma=0.4).y,Gray(0.5))
alphablend.(alphamask(Gray.(img),masktype="DiskFade",sigma=15,radius=0.45).y,Gray(0.5))
