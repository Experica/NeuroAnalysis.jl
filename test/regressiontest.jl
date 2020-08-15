## Tuning Feature
x = -3:0.002:3
y = gaussianf.(x,σ=1)
hw = halfwidth(y,circ=false,x=x)

x = -5:0.002:1
y = gaussianf.(x,σ=1)
hw = halfwidth(y,circ=false,x=x)
hw = halfwidth(y,circ=true,x=x)

x = -1:0.002:5
y = gaussianf.(x,σ=1)
hw = halfwidth(y,circ=false,x=x)
hw = halfwidth(y,circ=true,x=x)

## Circular Tuning Regression

# save(joinpath(@__DIR__,"circdata.jld2"),"y",y,"x",x,"y1",y1,"x1",x1)
y,x,y1,x1=load(joinpath(@__DIR__,"circdata.jld2"),"y","x","y1","x1")

# Tuning Curve, von Mises and Generalized von Mises
plotcondresponse(y,DataFrame(Ori=rad2deg.(x)),u=0,projection=:polar)
plotcondresponse(y1,DataFrame(Ori=rad2deg.(x1)),u=1,projection=:polar)

plot([vmf,gvmf],0,2π,label=["vmf" "gvmf"],projection=:polar)

# GvM and vM Fitting
mfit = fitmodel(:gvm,x,y)
plot(x->mfit.fun(x,mfit.param),0,2π,proj=:polar)

mfit = fitmodel(:vmn2,x1,y1)
plot(x->mfit.fun(x,mfit.param),0,2π,proj=:polar)

## Tuning Properties
f=factorresponsefeature(rad2deg.(x),y,factor=:Ori)

plot(x->f.fit.gvm.fun(x,f.fit.gvm.param),0,2π,lw=2)
scatter!(mod.(x.+0.5π,2π),y)

plot(x->f.fit.vmn2.fun(x,f.fit.vmn2.param),0,2π,lw=2)
scatter!(x,y)

f1=factorresponsefeature(rad2deg.(x1),y1,factor=:Ori)

plot(x->f1.fit.gvm.fun(x,f1.fit.gvm.param),0,2π,lw=2)
scatter!(mod.(x1.+0.5π,2π),y1)

plot(x->f1.fit.vmn2.fun(x,f1.fit.vmn2.param),0,2π,lw=2)
scatter!(x1,y1)
