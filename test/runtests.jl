using NeuroAnalysis, Test, BenchmarkTools



@testset begin

spike = sort(rand(1:100000,1000))
on = collect(10:100:95000)
off = collect(90:100:95080)
n1 = subrvr(spike,on,off,israte=false)
n2 = subrvr_ono(spike,on,off,israte=false)
@test n1==n2
# @btime n1 = subrvr($(spike),$(on),$(off),israte=false)
# @btime n2 = subrvr_ono($(spike),$(on),$(off),israte=false)

end
