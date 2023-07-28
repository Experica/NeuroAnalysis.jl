include("CSD.jl")

"""
Epochs of Channel Sample.
"""
function epochsample(x,fs,epochs,chs;whiten=nothing,fun=nothing)
    nepoch = size(epochs,1)
    epochis = round.(Int,epochs.*SecondPerUnit.*fs)
    minepochlength = minimum(diff(epochis,dims=2))
    ys=Array{Float64}(undef,length(chs),minepochlength,nepoch)
    for i in 1:nepoch
        y = x[chs,range(max(1,epochis[i,1]),length=minepochlength)]
        isnothing(whiten) || (y=whiten*y)
        isnothing(fun) || (y=fun(y))
        ys[:,:,i] = y
    end
    return nepoch==1 ? dropdims(ys,dims=3) : ys
end

"fill data in the shape of mask, where masked channels are replaced with local average"
function fill2mask(ys, mask;chmap=1:size(ys,1), replacemask = true,randreplace=false)
    nrow, ncol = size(mask)
    yss = zeros(nrow, ncol, size(ys)[2:end]...)
    for i in eachindex(chmap)
        r,c = chshapenp(chmap[i],nrow,ncol)
        @views yss[r,c,:,:] = ys[i,:,:]
    end
    if replacemask # Local unmasked channels averaging replacement for masked channels
        rmask = copy(mask)
        for (r, c) in Tuple.(findall(mask))
            rd = r - 1;ru = r + 1
            while 1 <= rd && rmask[rd, c]
                rd -= 1
            end
            while ru <= nrow && rmask[ru, c]
                ru += 1
            end
            if rd < 1
                if ru > nrow
                    @error "No Unmasked Channels in Data to Replace Masked Channel[$r, $c]."
                else
                    yss[r, c, :, :] = yss[ru, c, :, :]
                end
            else
                if ru > nrow
                    yss[r, c, :, :] = yss[rd, c, :, :]
                else
                    @views yss[r, c, :, :] = (yss[rd, c, :, :] .+ yss[ru, c, :, :]) / 2
                end
            end
            if randreplace
                @views yss[r,c,:,:] = yss[r,c,:,:] .+ 3*randn(size(yss)[3:end]...)*std(yss[r,c,:,:])
            end
            rmask[r, c] = false
        end
    end
    return yss
end

"Spatial and Temporal filter of data stream"
function stfilter(rm;spatialtype=:none,ir=1,or=2,temporaltype=:none,ti=1:size(rm,2),hbordervalue=nothing)
    nd = ndims(rm)
    if nd==3
        nr,nc,n = size(rm)
        y = Array{Float64}(undef,nr,nc,n)
        for i in 1:n
            y[:,:,i] = stfilter(rm[:,:,i],spatialtype=spatialtype,ir=ir,or=or,temporaltype=temporaltype,ti=ti,hbordervalue=hbordervalue)
        end
        return y
    else
        if spatialtype==:all
            rm .-= mean(rm,dims=1)
        elseif spatialtype==:annulus
            foreach(i->rm[i,:] .-= dropdims(mean(rm[filter(j-> 1<=j<=size(rm,1),union(i-or:i-ir,i+ir:i+or)),:],dims=1),dims=1),1:size(rm,1))
        end
        if temporaltype == :all
            rm .-= mean(rm,dims=2)
        elseif temporaltype == :sub
            rm .-= mean(rm[:,ti],dims=2)
        elseif temporaltype == :rc
            rm = rm ./ mean(rm[:,ti],dims=2) .- 1
        elseif temporaltype == :rcb
            rm = relativechange(rm,mean(rm[:,ti],dims=2);balance=true)
        elseif temporaltype == :log2r
            rm = log2.(rm ./ mean(rm[:,ti],dims=2))
        elseif temporaltype == :z
            rm = (rm .- mean(rm[:,ti],dims=2)) ./ std(rm[:,ti],dims=2)
        end
        if !isnothing(hbordervalue)
            rm[[1,end],:].=hbordervalue
        end
        return rm
    end
end

function relativechange(response,background;balance=false)
    rc = response./background
    if balance
        li = rc .< 1
        hi = .!li
        @views rc[li] = 1 .- 1 ./ rc[li]
        @views rc[hi] = rc[hi] .- 1
    else
        rc .-= 1
    end
    rc
end

"Remove line noise and its harmonics by notch filter"
rmline(y,fs;freq=60,nh=3,bw=3) = rmline!(copy!(similar(y,Float64),y),fs;freq,nh,bw)
function rmline!(y,fs;freq=60,nh=3,bw=3)
    for i=1:nh
        f = iirnotch(freq*i,bw;fs)
        @views for j = 1:size(y,1)
            y[j,:]=filtfilt(f,y[j,:])
        end
    end
    return y
end

"High pass and low pass filtering"
function hlpass(y,fs;high=0,low=Inf)
    fy = similar(y,Float64)
    copy!(fy,y)
    hlpass!(fy,fs;high,low)
end
function hlpass!(y,fs;high=0,low=Inf)
    if high>0
        f = digitalfilter(Highpass(high;fs=fs), Butterworth(4))
        for j=1:size(y,1)
            @views y[j,:]=filtfilt(f,y[j,:])
        end
    end
    if low<Inf
        f = digitalfilter(Lowpass(low;fs=fs), Butterworth(4))
        for j=1:size(y,1)
            @views y[j,:]=filtfilt(f,y[j,:])
        end
    end
    return y
end

"time to sample index"
time2sampleindex(t,fs;secondperunit=1,minsampleindex::Integer=0,maxsampleindex::Integer=typemax(Int64)) = clamp(round(Int,t*secondperunit*fs),minsampleindex,maxsampleindex)
"sample index to time"
sampleindex2time(s,fs;secondperunit=1,mintime::Real=0,maxtime::Real=typemax(Float64)) = clamp(s/fs/secondperunit,mintime,maxtime)

"Convert epochs in time to sample index"
function epoch2sampleindex(epochs,fs;minsampleindex::Integer=1,maxsampleindex::Integer=typemax(Int64))
    nepoch = size(epochs,1)
    epochis = time2sampleindex.(epochs,fs;secondperunit=SecondPerUnit,minsampleindex,maxsampleindex)
    minepochlength = minimum(diff(epochis,dims=2))
    si = [range(epochis[i,1],length=minepochlength) for i in 1:nepoch]
    return nepoch==1 ? si[1] : si
end

"""
Get digital rising and falling edges in a analog stream
return di: edges index
       dv: edges active/inactive state
"""
function parsedigitalinanalog(x,ht,lt=ht,n=length(x);activehigh=true)
    di=Int[];dv=Bool[]
    if activehigh
        cs=x[1]>ht
        for i in 2:n
            t = x[i]
            if t > ht && !cs
                cs=true;push!(di,i);push!(dv,cs)
            elseif t < lt && cs
                cs=false;push!(di,i);push!(dv,cs)
            end
        end
    else
        cs=x[1]<lt
        for i in 2:n
            t = x[i]
            if t < lt && !cs
                cs=true;push!(di,i);push!(dv,cs)
            elseif t > ht && cs
                cs=false;push!(di,i);push!(dv,cs)
            end
        end
    end
    return di,dv
end

"Multi-Taper power spectrum estimation"
function powerspectrum(x,fs;freqrange=[0,100],nw=4)
    nd=ndims(x)
    if nd==3
        n=size(x,1);ne=size(x,3)
        @views ps = mt_pgram(x[1,:,1];fs,nw)
        fi = first(freqrange).<=freq(ps).<=last(freqrange)
        freqs = freq(ps)[fi]
        p = Array{Float64}(undef,n,length(freqs),ne)
        Threads.@threads for i in 1:n
            for j in 1:ne
                @views p[i,:,j] = power(mt_pgram(x[i,:,j];fs,nw))[fi]
            end
        end
    elseif nd == 2
        @views ps = mt_pgram(x[1,:];fs,nw)
        fi = first(freqrange).<=freq(ps).<=last(freqrange)
        freqs = freq(ps)[fi];n = size(x,1)
        p = Array{Float64}(undef,n,length(freqs))
        p[1,:]=power(ps)[fi]
        Threads.@threads for i in 2:n
            @views p[i,:]=power(mt_pgram(x[i,:];fs,nw))[fi]
        end
    else
        ps = mt_pgram(x;fs,nw)
        fi = first(freqrange).<=freq(ps).<=last(freqrange)
        freqs = freq(ps)[fi];p = power(ps)[fi]
    end
    return p,freqs
end

"Multi-Taper channel pair-wise coherence spectrum estimation"
function coherencespectrum(x,fs;freqrange=[0,100],nw=4,ismean=false)
    nd=ndims(x)
    if nd==3
        n=size(x,1);ne=size(x,3)
        @views c1,freqs = coherencespectrum(x[:,:,1],fs;freqrange,nw,ismean)
        c = Array{Float64}(undef,n,n,length(freqs),ne)
        c[:,:,:,1]=c1
        Threads.@threads for i in 2:ne
            @views ci,freqs = coherencespectrum(x[:,:,i],fs;freqrange,nw,ismean)
            c[:,:,:,i]=ci
        end
    elseif nd == 2
        cs = mt_coherence(x;freq_range=freqrange,fs,nw)
        c = coherence(cs);freqs=freq(cs)
        if ismean
            c = mean(c,dims=3)
            freqs=mean(freqs)
        end
    else
        @error "No coherence for one signal."
    end
    length(freqs)==1 && (c=dropdims(c,dims=3))
    return c,freqs
end

"Weighted average of banded masked matrix"
function bandmean(x;r=5,s=1)
    mask = BandedMatrix(Ones(size(x)...),(r,r))

    g = OffsetArray(gaussianf.(-r:r,σ=s),-r:r)
    g[0] = 0
    g ./= sum(g)
    foreach(i->mask[band(i)].=g[i],-r:r)

    bm = dropdims(sum(x.*Matrix(mask),dims=1),dims=1)
end

"local channel pairs and gaussian weights within a circle centered for each channel"
function localpairweight(chpos;lr=55,sigma=25)
    chlocalpair=[];chlocalweight=[];nch=size(chpos,1)
    @views @inbounds for ch in 1:nch
        dist2ch = [norm(chpos[i,:].-chpos[ch,:]) for i in 1:nch]
        lch = findall(0 .< dist2ch .< lr)
        lw = gaussianf.(dist2ch[lch],σ=sigma)
        push!(chlocalpair,[Set((ch,j)) for j in lch])
        push!(chlocalweight,lw./sum(lw))
    end
    allpair = reduce(union,chlocalpair)
    return chlocalpair,chlocalweight,allpair
end

"Weighted average of local channel pair coherences for each channel"
function localcoherence(x,fs,chlocalpair,chlocalweight,allpair = reduce(union,chlocalpair);freqrange=[0,100],nw=4,chgi=[])
    nch=size(x,1)
    @views allpairco = [mt_coherence(x[[p...],:];freq_range=freqrange,fs,nw) for p in allpair]
    freqs = freq(first(allpairco));nfreq = length(freqs)
    allpairc = Array{Float64}(undef,nfreq,length(allpairco))
    @views foreach(i->allpairc[:,i]=coherence(allpairco[i])[1,2,:], eachindex(allpairco))
    chlc = Array{Float64}(undef,nch,nfreq)
    @views foreach(i->chlc[i,:]=allpairc[:,indexin(chlocalpair[i],allpair)]*chlocalweight[i],eachindex(chlocalpair))

    isempty(chgi) && return chlc,freqs

    chglc = Array{Float64}(undef,length(chgi),nfreq)
    @views foreach(i->chglc[i,:] = mapreduce(g->chlc[g,:],.+,chgi[i])/length(chgi[i]),eachindex(chgi))

    return chglc,freqs
end
function localcoherence(x,chpos,fs;freqrange=[0,100],nw=4,lr=55,sigma=25,chgroupdim=2)
    nd = ndims(x)
    chlocalpair,chlocalweight,allpair=localpairweight(chpos;lr,sigma)
    chgi = chgroupdim==0 ? [] : [findall(chpos[:,chgroupdim].==up) for up in unique(chpos[:,chgroupdim])]
    if nd==3
        ne=size(x,3)
        @views lc1,freqs = localcoherence(x[:,:,1],fs,chlocalpair,chlocalweight,allpair;freqrange,nw,chgi)
        chlc = Array{Float64}(undef,size(lc1,1),length(freqs),ne)
        chlc[:,:,1]=lc1
        Threads.@threads for i in 2:ne
            @views lci,freqs = localcoherence(x[:,:,i],fs,chlocalpair,chlocalweight,allpair;freqrange,nw,chgi)
            chlc[:,:,i]=lci
        end
    elseif nd==2
        chlc,freqs = localcoherence(x,fs,chlocalpair,chlocalweight,allpair;freqrange,nw,chgi)
    else
        @error "No coherence for one signal."
    end
    return chlc,freqs
end

@doc raw"""
Discrete Fourier Transform at frequencies

```math
DFT[k] = \sum_{n=0}^{N-1} x[n] e^{\frac{-i2\pi kn}{N}}, k=0:N-1, N=length(x), fₛ=SamplingFreq(x), fₖ=fₛ/N=FreqResolution(DTFT)
```

1. x: signal
2. fs: simpling frequency of signal
3. f...: at which frequencies DFT are directly evaluated
"""
function dft(x,fs,f...)
    N = length(x)
    fₖ = fs/N
    ks = round.(Int,f./fₖ)
    Fs = [zero(ComplexF64) for _ in f]
    Ω = [exp(-im*2π*n/N) for n in 0:(N-1)]
    @inbounds for n in 0:(N-1)
        @inbounds for i in eachindex(f)
            Fs[i] += x[n+1] * Ω[((n*ks[i])%N)+1]
        end
    end
    return Fs
end
