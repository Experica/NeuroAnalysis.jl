st1 = vcat([10i .+ (0:9) for i in 0:9]...)
# a epoch of spike train
y,n,w,i = epochspiketrain(st1,10,20)
@test y==10:19
@test n==10
@test w==(10,20)
@test i==11:20
# two epochs of spike train
ys,ns,ws,is = epochspiketrain(st1,[10,20],[20,30])
@test ys[1]==y
@test ns[1]==n
@test ws[1]==w
@test is[1]==i
@test ys[2]==20:29
@test ns[2]==10
@test ws[2]==(20,30)
@test is[2]==21:30
ys1,ns1,ws1,is1 = epochspiketrain(st1,[10 20;20 30])
@test ys1==ys
@test ns1==ns
@test ws1==ws
@test is1==is
ys2,ns2,ws2,is2 = epochspiketrain(st1,[10,20,30])
@test ys2==ys
@test ns2==ns
@test ws2==ws
@test is2==is
# two epochs of two spike trains
st1st1 = [st1,st1]
yss,nss,wss,iss = epochspiketrains(st1st1,[10,20,30])
@test yss[1]==yss[2]
@test nss[1]==nss[2]
@test iss[1]==iss[2]
@test yss[1]==ys
@test nss[1]==ns
@test wss == ws
@test iss[1]==is
# two epoch responses of spike train
r = epochspiketrainresponse(st1,[10,20],[20,30],israte=false)
@test r[1]==r[2]==10
r1 = epochspiketrainresponse(st1,[10 20;20 30],israte=false)
@test r1[1]==r1[2]==10
# epochspiketrainresponse and epochspiketrainresponse_ono
spike = sort(rand(1:100000,1000))
on = collect(0:150:95000)
off = collect(100:150:95100)
n1 = epochspiketrainresponse(spike,on,off,israte=false)
n2 = epochspiketrainresponse_ono(spike,on,off,israte=false)
@test n1==n2
@btime epochspiketrainresponse($(spike),$(on),$(off),israte=false)
@btime epochspiketrainresponse_ono($(spike),$(on),$(off),israte=false)
# mean psth of two spike trains
m,se,x = psthspiketrains(st1st1,[10,20,30],israte=false)
@test m==[10,10]
@test se==[0,0]
@test x==[15,25]
m1,se1,x1 = psthspiketrains(st1st1,10:10:30,israte=false)
@test m1==m
@test se1==se
@test x1==x


## Poisson Model of Spike Generation, [David Heeger(2000), Poisson Model of Spike Generation](http://www.cns.nyu.edu/~david/handouts/poisson.pdf)
rate = 25
duration = 1000
refractoryperiod = 2
ntrial = 50

## Homogeneous Poisson Process, where instantaneous ﬁring rate is constant
# Inter-Spike-Interval Exponential Distribution Method
sts = [poissonspiketrain(rate,duration,rp=refractoryperiod) for _ in 1:ntrial]
plotspiketrain(sts,timeline=[0,duration])

# Instantaneous Firing Rate Method
sts = [poissonspiketrain(t->rate,duration,rp=refractoryperiod) for _ in 1:ntrial]
plotspiketrain(sts,timeline=[0,duration])

## Inhomogeneous Poisson Process, where instantaneous ﬁring rate is changing
sts = [poissonspiketrain(t->rate*(sin(0.05t)+1),duration,rp=refractoryperiod) for _ in 1:ntrial]
plotspiketrain(sts,timeline=[0,duration])


## spike train jitter resampling
st = sort(rand(25))*duration
jst = spikejitter(st;n=ntrial,l=50.5,win=:center)
jst = spikejitter(st;n=ntrial,l=100,win=:fix)
plotspiketrain([st jst]',timeline=0:100:duration)

st = floor.(Int,st)
jst = spikejitter(st;n=ntrial,l=50.5,win=:center)
jst = spikejitter(st;n=ntrial,l=100,win=:fix)
plotspiketrain([st jst]',timeline=0:100:duration)

bst = zeros(duration,ntrial)
foreach(t->bst[clamp!(ceil.(Int,jst[:,t]),1,duration),t] .= 1,1:size(bst,2))
jbst = shufflejitter(bst;l=100)
plotspiketrain(map(t->findall(jbst[:,t].==1),1:size(jbst,2)),timeline=0:100:duration)

bst=psthspiketrains(sts,0:duration,israte=false,ismean=false).mat
heatmap(bst,xlabel="Time (ms)",ylabel="Trials",tickdir=:out)
jbst = shufflejitter(bst',100)
heatmap(jbst',xlabel="Time (ms)",ylabel="Trials",tickdir=:out)
