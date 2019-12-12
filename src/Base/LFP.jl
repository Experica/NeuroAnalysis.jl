include("CSD.jl")

export subrm,reshape2mask,rmline!,hlpass,time2sample,sample2time,epoch2samplerange,parsedigitalinanalog,powerspectrum

function subrm(rm,fs,epochs,chs;fun=nothing)
    nepoch = size(epochs,1)
    epochis = floor.(Int,epochs.*SecondPerUnit.*fs)
    minepochlength = minimum(diff(epochis,dims=2))
    ys=Array{Float64}(undef,length(chs),minepochlength,nepoch)
    for i in 1:nepoch
        y = rm[chs,range(max(1,epochis[i,1]),length=minepochlength)]
        if !isnothing(fun)
            y=fun(y)
        end
        ys[:,:,i] = y
    end
    return nepoch==1 ? dropdims(ys,dims=3) : ys
end

function reshape2mask(ys,chmask;replacemask=true)
    nrow,ncol=size(chmask)
    yss=Array{Float64}(undef,nrow,ncol,size(ys)[2:end]...)
    for c in 1:ncol
        yss[:,c,:,:] = ys[c:ncol:end,:,:] # channel counting in chmask is from cols(left -> right), then rows(up -> down)
    end
    if replacemask
        for (r,c) in Tuple.(findall(chmask))
            dr = r-1;ur = r+1;
            while true
                (!chmask[dr,c] || dr<=1) && break
                dr-=1
            end
            while true
                (!chmask[ur,c] || ur>=nrow) && break
                ur+=1
            end
            yss[r,c,:,:] = (yss[dr,c,:,:] .+ yss[ur,c,:,:]) / 2 # Local non-mask channels averaging replacement for masking channels
        end
    end
    return yss
end

"Remove line noise and its harmonics by notch filter"
function rmline!(y,fs;freq=60,nh=3,bw=3)
    for i=1:nh
        f = iirnotch(freq*i,bw;fs=fs)
        for j=1:size(y,1)
            y[j,:]=filtfilt(f,y[j,:])
        end
    end
    return y
end

"High pass and low pass filters"
function hlpass(y,fs;high=0,low=Inf)
    fy=copy(y)
    if high>0
        f = digitalfilter(Highpass(high;fs=fs), Butterworth(4))
        for j=1:size(fy,1)
            fy[j,:]=filtfilt(f,fy[j,:])
        end
    end
    if low<Inf
        f = digitalfilter(Lowpass(low;fs=fs), Butterworth(4))
        for j=1:size(fy,1)
            fy[j,:]=filtfilt(f,fy[j,:])
        end
    end
    return fy
end

time2sample(t,fs;secondperunit=1,t0sample::Integer=1) = round(Int,t*secondperunit*fs)+t0sample
sample2time(s,fs;secondperunit=1,t0sample::Integer=1) = (s-t0sample)/fs/secondperunit

function epoch2samplerange(epochs,fs)
    nepoch = size(epochs,1)
    epochis = floor.(Int,epochs.*SecondPerUnit.*fs)
    minepochlength = minimum(diff(epochis,dims=2))
    sr = [range(max(1,epochis[i,1]),length=minepochlength) for i in 1:nepoch]
    return nepoch==1 ? sr[1] : sr
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

function powerspectrum(x,fs;freqrange=0:100,nw=4)
    nd=ndims(x)
    if nd==3
        n=size(x,1);ne=size(x,3)
        ps = mt_pgram(x[1,:,1],fs=fs,nw=nw)
        fi = in.(freq(ps),[freqrange])
        freqs = freq(ps)[fi]
        p = Array{Float64}(undef,n,length(freqs),ne)
        for i in 1:n,j in 1:ne
            p[i,:,j] = power(mt_pgram(x[i,:,j],fs=fs,nw=nw))[fi]
        end
    else
        ps = mt_pgram(x[1,:],fs=fs,nw=nw)
        fi = in.(freq(ps),[freqrange])
        freqs = freq(ps)[fi];n = size(x,1)
        p = Array{Float64}(undef,n,length(freqs))
        p[1,:]=power(ps)[fi]
        for i in 2:n
            p[i,:]=power(mt_pgram(x[i,:],fs=fs,nw=nw))[fi]
        end
    end
    return p,freqs
end
