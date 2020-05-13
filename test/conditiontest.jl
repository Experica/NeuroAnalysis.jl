condtable = DataFrame(Ori=repeat([0,90,180,270],2),SpatialFreq=repeat([2,4],inner=4))
ctc = repeat(condtable,8)

cond = condin(ctc)
@test all(cond.n.==8)
