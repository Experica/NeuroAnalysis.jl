condtable = DataFrame(Ori=repeat([0,90,180,270],2),SpatialFreq=repeat([2,4],inner=4))
ctc = repeat(condtable,8)

@test_throws ErrorException flin(DataFrame(i=1:10))
@test_throws ErrorException flin(DataFrame(n=1:10))
@test_throws ErrorException flin(DataFrame(i=1:10,n=1:10))
@test_throws ErrorException condin(DataFrame(i=1:10))
@test_throws ErrorException condin(DataFrame(n=1:10))
@test_throws ErrorException condin(DataFrame(i=1:10,n=1:10))

cond = condin(ctc)
@test all(cond.n.==8)

fl = flin(ctc)
@test all(fl[:Ori].n.==16)
@test all(fl[:SpatialFreq].n.==32)

[condfactor(r) for r in eachrow(cond)]
@test condfactor(cond) == collect(keys(fl))

finalfactor(cond)
condstring(cond)

condresponse(rand(nrow(ctc)),cond.i)
df = condresponse(rand(nrow(ctc)),cond)
condresponse(Dict(1=>rand(nrow(ctc))),cond)
condresponse(Dict(2=>rand(nrow(ctc))),ctc,[:Ori])

factorresponse(df)

@test isresponsive(randn(nrow(ctc)),randn(nrow(ctc)).+1)
@test isresponsive(randn(nrow(ctc)),randn(nrow(ctc)).+1,cond.i)
ismodulative([DataFrame(Y=randn(nrow(ctc))) ctc])
