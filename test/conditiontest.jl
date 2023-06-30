condtable = DataFrame(Ori=repeat([0,90,180,270],2),SpatialFreq=repeat([2,4],inner=4))
ctc = repeat(condtable,8)

@test_throws ErrorException flin(DataFrame(i=1:5))
@test_throws ErrorException flin(DataFrame(n=1:5))
@test_throws ErrorException flin(DataFrame(i=1:5,n=1:5))
@test_throws ErrorException condin(DataFrame(i=1:5))
@test_throws ErrorException condin(DataFrame(n=1:5))
@test_throws ErrorException condin(DataFrame(i=1:5,n=1:5))

cond = condin(ctc)
@test all(i->i==8, cond.n)

fl = flin(ctc)
@test all(i->i==16, fl.Ori.n)
@test all(i->i==32, fl.SpatialFreq.n)

[condfactor(r) for r in eachrow(cond)]
[finalfactor(r) for r in eachrow(cond)]
finalfactor(cond)
condstring(cond)
@test condfactor(cond) == collect(keys(fl))
findcond(cond,Ori=0,SpatialFreq=2)

rs1 = randn(nrow(ctc))
rs2 = randn(nrow(ctc))
condresponse(rs1,cond.i)
mseuc = condresponse(rs1,cond)
condresponse(Dict(1=>rs1,2=>rs2),cond)

fi,fa = factorspace(cond)
fr = factorresponse(rs1,fi)
fm,fa = factorspace(mseuc;col=:m)
@test fm == fr.m

rs1 = randn(nrow(ctc))
rs2 = randn(nrow(ctc))
@test isresponsive(rs1,rs2.+1)
@test isresponsive(rs1,rs2.+1,cond.i)
# ismodulative([DataFrame(Y=randn(nrow(ctc))) ctc])
ismodulative(rs1,cond.i)
