condtable = DataFrame(Ori=repeat([0,90,180,270],2),SpatialFreq=repeat([2,4],inner=4))
ctc = repeat(condtable,8)

cond = condin(ctc)
@test all(cond.n.==8)

fl = flin(ctc)
@test all(fl[:Ori].n.==16)
@test all(fl[:SpatialFreq].n.==32)

@test condfactor(cond) == collect(keys(fl))

condstr = condstring(cond)
