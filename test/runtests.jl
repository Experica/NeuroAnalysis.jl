using NeuroAnalysis, Test, BenchmarkTools

@testset begin
## subrvr and subrvr_ono
spike = sort(rand(1:100000,1000))
on = collect(10:100:95000)
off = collect(90:100:95080)
n1 = subrvr(spike,on,off,israte=false)
n2 = subrvr_ono(spike,on,off,israte=false)
@test n1==n2
# @btime n1 = subrvr($(spike),$(on),$(off),israte=false)
# @btime n2 = subrvr_ono($(spike),$(on),$(off),israte=false)
## image and freqimage
ppd = 30;ori=0.5π;sf=2
img = grating(θ=ori,sf=sf,ppd=ppd)
# heatmap(img,yflip=true)
ps,f1,f2 = powerspectrum(img,ppd,freqrange=[-6,6])
# heatmap(f1,f2,ps)
eori,esf = freqimagestats(ps,f1,f2)
@test sf≈esf
@test ori≈eori

end
