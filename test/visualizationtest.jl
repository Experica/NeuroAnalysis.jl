x = rand(0:1000,2000)
y = rand(1:200,2000)
plotspiketrain(x,y,timespan=[250,750],timeline=[0,1000])
plotspiketrain(x,y,group=rand(1:5,2000),timespan=[250,750],timeline=[0,1000])

sts = [rand(0:1000,rand(5:10)) for _ in 1:200]
uids = [rand(1:5,length(st)) for st in sts]
plotspiketrain(sts,timespan=[250,750],timeline=[0,1000])
plotspiketrain(sts;uids,timespan=[250,750],timeline=[0,1000])

# single group
rs = rand(0:50,nrow(ctc))
mseuc = condresponse(rs,cond)
plotcondresponse(mseuc)

mseuc = condresponse(rs,condin(ctc[:,[:Ori]]))
plotcondresponse(mseuc)

# multiple groups
mseuc = [DataFrame(m=12 .+rand(10),se=rand(1:0.1:2,10),u="SU1",Ori=range(0,step=36,length=10));
    DataFrame(m=8 .+rand(10),se=rand(1:0.1:2,10),u="SU0",Ori=range(0,step=36,length=10));
    DataFrame(m=4 .+rand(10),se=rand(1:0.1:3,10),u="Pre_SU0",Ori=range(0,step=36,length=10))]

plotcondresponse(mseuc,color=[:black,:gray20,:gray75],linewidth=[3,2,1],grid=true)
plotcondresponse(mseuc,projection=:polar,color=[:black,:gray20,:gray75],linewidth=[3,2,1],grid=true)


plotanalog(randn(100))
plotanalog(randn(50,100),plottype=:line)
plotanalog(rand(100,100),color=:vik,clims=(-3,3),n=rand(100))

## color
huecolors()

range(RGB(0,0,0),RGB(1,1,1))
range(RGB(0,0,1),RGB(1,1,1),RGB(1,0,0))
cgrad(RGB(1,0.0,0),RGB(0,1.0,0))
cgrad(RGB(1,0.0,0),RGB(0,1.0,0),RGB(0,0.0,1))

plotcolormap(ColorMaps["dkl_mcchue_l0"].colors)
plotcolormap(ColorMaps["lidkl_mcchue_l0"].colors)
plotcolormap(ColorMaps["hsl_mshue_l0.4"].colors)

plotcolormap(ColorMaps["dkl_mcclumiso"].colors,shape=:sin)
plotcolormap(ColorMaps["lms_mccliso"].colors,shape=:sin)
plotcolormap(ColorMaps["lms_mccmiso"].colors,shape=:sin)
plotcolormap(ColorMaps["lms_mccsiso"].colors,shape=:sin)

plotcolormap(ColorMaps["dkl_mcclumiso"].colors,shape=:diag,markersize=10)
plotcolormap!(ColorMaps["dkl_mcclmiso"].colors,shape=:hline,markersize=10)
plotcolormap!(ColorMaps["dkl_mccslmiso"].colors,shape=:vline,markersize=10)
plotcolormap!(ColorMaps["lidkl_mcchue_l0"].colors,markersize=12)


# foreach(i->save(joinpath(resultdir,"cm_dkl_mcchue_l$lum$i"),cm_dkl),[".yaml",".mat"])
#
# foreach(i->save(joinpath(resultdir,"cm_lidkl_mcchue_l$lum$i"),cm_lidkl),[".yaml",".mat"])
#
# foreach(i->save(joinpath(resultdir,"cm_hsl_mshue_l$l$i"),cm_hsl),[".yaml",".mat"])
