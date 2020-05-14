## hartley space
# @manipulate for kb in 0:4, ke in 1:5, dk in 0.4:0.2:1, p in 0:0.005:1, shape in [:square, :circle]
#     plothartleyspace(hartleysubspace(kbegin=kb,kend=ke,dk=dk,phase=p,shape=shape),floor(Int,ke/dk),dk)
# end

hs = hartleysubspace(kbegin=0.2,kend=6.6,dk=0.2,addhalfcycle=true)
mhg = mapreduce(i -> begin
                ss = cas2sin(i...)
                grating(θ = ss.θ, sf = ss.f, phase = ss.phase, size = (5, 5), ppd = 30)
                end, (i, j) -> i .+ j, hs) / length(hs)
# mean of a hartley subspace gratings should be a uniform gray
minv,maxv = extrema(mhg)
@test minv≈maxv≈0.5

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
