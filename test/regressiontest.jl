x = 0:0.001:2π
y = sin.(x)

@test_throws Exception searchclosest(0,y;step=0)
@test searchclosest(0,y;start=1) == 1

y = vmf.(x,μ=π)
maxy,maxi = findmax(y)
plot(x,y)
hline!([maxy/2])
vline!([x[maxi]])
hw = halfwidth(y;circ=false,x)

y = vmf.(x,μ=0.25π)
hw = halfwidth(y;circ=false,x)
hw = halfwidth(y;circ=true,x)

y = vmf.(x,μ=1.75π)
hw = halfwidth(y;circ=false,x)
hw = halfwidth(y;circ=true,x)

## Circular Tuning Regression
# save(joinpath(@__DIR__,"circdata.jld2"),"y",y,"x",x,"y1",y1,"x1",x1)
y,x,y1,x1=load(joinpath(@__DIR__,"circdata.jld2"),"y","x","y1","x1")

# Tuning Curve, von Mises and Generalized von Mises
cond = condin(DataFrame(Ori=rad2deg.(x)))
mseuc = condresponse(y,cond,u=0)
plotcondresponse(mseuc,projection=:polar)
cond1 = condin(DataFrame(Ori=rad2deg.(x1)))
mseuc1 = condresponse(y1,cond1,u=1)
plotcondresponse(mseuc1,projection=:polar)

plot([x->vmf(x,κ=200,n=2),x->gvmf(x,β=0.3,μ₁=0.25π,κ₁=1,μ₂=0.75π,κ₂=1)],0:0.01:2π,label=["vmf" "gvmf"],projection=:polar)

# vM and GvM Fitting
mfit = fitmodel(:vm,x,y)
plot(x->mfit.fun(x,mfit.param),0,2π,proj=:polar)
mfit = fitmodel(:vmn2,x1,y1)
plot(x->mfit.fun(x,mfit.param),0,2π,proj=:polar)

mfit = fitmodel(:gvm,x1,y1)
plot(x->mfit.fun(x,mfit.param),0,2π,proj=:polar)
mfit = fitmodel(:gvm,x,y)
plot(x->mfit.fun(x,mfit.param),0,2π,proj=:polar)

# Tuning Properties
rf=factorresponsefeature(cond.Ori,map(i->y[i],cond.i),factor=:Ori)
plot(x->rf.fit.mfit.fun(x,rf.fit.mfit.param),0,2π,lw=2)
scatter!(x,y)

rf1=factorresponsefeature(cond1.Ori,map(i->y1[i],cond1.i),factor=:Angle)
plot(x->rf1.fit.mfit.fun(x,rf1.fit.mfit.param),0,2π,lw=2)
scatter!(x1,y1)


## Spatial Frequency Tuning Regression
# save(joinpath(@__DIR__,"sfdata.jld2"),"y",y,"x",x,"y1",y1,"x1",x1)
y,x,y1,x1=load(joinpath(@__DIR__,"sfdata.jld2"),"y","x","y1","x1")

cond = condin(DataFrame(SpatialFreq=x))
mseuc = condresponse(y,cond,u=0)
plotcondresponse(mseuc)
cond1 = condin(DataFrame(SpatialFreq=x1))
mseuc1 = condresponse(y1,cond1,u=1)
plotcondresponse(mseuc1)

# sfgaussian and sfdog Fitting
mfit = fitmodel(:sfdog,x,y,MinDeltaFitnessTolerance=1e-2)
plot(x->mfit.fun(x,mfit.param),0,8)
mfit = fitmodel(:sfdog,x1,y1,MinDeltaFitnessTolerance=1e-2)
plot(x->mfit.fun(x,mfit.param),0,8)

mfit = fitmodel(:sfgaussian,x,y,MinDeltaFitnessTolerance=1e-2)
plot(x->mfit.fun(x,mfit.param),0,8)
mfit = fitmodel(:sfgaussian,x1,y1,MinDeltaFitnessTolerance=1e-2)
plot(x->mfit.fun(x,mfit.param),0,8)

# Tuning Properties
rf=factorresponsefeature(cond.SpatialFreq,map(i->y[i],cond.i),factor=:SpatialFreq)
plot(x->rf.fit.mfit.fun(x,rf.fit.mfit.param),extrema(x)...,lw=2)
scatter!(x,y)

rf1=factorresponsefeature(cond1.SpatialFreq,map(i->y1[i],cond1.i),factor=:SpatialFreq)
plot(x->rf1.fit.mfit.fun(x,rf1.fit.mfit.param),extrema(x1)...,lw=2)
scatter!(x1,y1)