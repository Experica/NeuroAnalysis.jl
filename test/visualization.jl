using HDF5,JLD
const testpath = splitdir(@__FILE__)[1]
testdata = load(joinpath(testpath,"testdata.jld"))
spike = testdata["spike"];stion = testdata["stion"];stioff = testdata["stioff"]

using NeuroAnalysis.NABase,NeuroAnalysis.NAVisualization
ys,ns,ws,is = subrv(spike,stion,stioff,isminzero=true)
pst = plotspiketrain(ys)
ppsth = plotpsth(ys,0:5:300)
savefig(pst,"spiketrain",path=testpath)
savefig(ppsth,"psth",path=testpath,format="png")
