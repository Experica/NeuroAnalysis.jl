include("CSD.jl")

"""
Epochs of Channel Sample.
"""
function epochsample(x,fs,epochs,chs;fun=nothing)
    nepoch = size(epochs,1)
    epochis = round.(Int,epochs.*SecondPerUnit.*fs)
    minepochlength = minimum(diff(epochis,dims=2))
    ys=Array{Float64}(undef,length(chs),minepochlength,nepoch)
    for i in 1:nepoch
        y = x[chs,range(max(1,epochis[i,1]),length=minepochlength)]
        if !isnothing(fun)
            y=fun(y)
        end
        ys[:,:,i] = y
    end
    return nepoch==1 ? dropdims(ys,dims=3) : ys
end

"reshape data stream in the same shape of chmask, where masked channels are replaced with local average"
function reshape2mask(ys, chmask; replacemask = true)
    nrow, ncol = size(chmask)
    yss = Array{Float64}(undef, nrow, ncol, size(ys)[2:end]...)
    for c = 1:ncol
        yss[:, c, :, :] = ys[c:ncol:end, :, :] # channel counting in chmask is from cols(left -> right), then rows(up -> down)
    end
    if replacemask # Local unmasked channels averaging replacement for masked channels
        rchmask = copy(chmask)
        for (r, c) in Tuple.(findall(chmask))
            rd = r - 1;ru = r + 1
            while 1 <= rd && rchmask[rd, c]
                rd -= 1
            end
            while ru <= nrow && rchmask[ru, c]
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
                    yss[r, c, :, :] = (yss[rd, c, :, :] .+ yss[ru, c, :, :]) / 2
                end
            end
            rchmask[r, c] = false
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
        end
        if !isnothing(hbordervalue)
            rm[[1,end],:].=hbordervalue
        end
        return rm
    end
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

"High pass and low pass filtering"
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

"time to sample index"
time2sampleindex(t,fs;secondperunit=1,t0sampleindex::Integer=1) = round(Int,t*secondperunit*fs)+t0sampleindex
"sample index to time"
sampleindex2time(s,fs;secondperunit=1,t0sampleindex::Integer=1) = (s-t0sampleindex)/fs/secondperunit

"Convert epochs in time to sample index"
function epoch2sampleindex(epochs,fs;t0sampleindex::Integer=1)
    nepoch = size(epochs,1)
    epochis = time2sampleindex.(epochs,fs;secondperunit=SecondPerUnit,t0sampleindex=t0sampleindex)
    minepochlength = minimum(diff(epochis,dims=2))
    si = [range(max(1,epochis[i,1]),length=minepochlength) for i in 1:nepoch]
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
        ps = mt_pgram(x[1,:,1],fs=fs,nw=nw)
        fi = map(f->freqrange[1]<=f<=freqrange[2],freq(ps))
        freqs = freq(ps)[fi]
        p = Array{Float64}(undef,n,length(freqs),ne)
        for i in 1:n,j in 1:ne
            p[i,:,j] = power(mt_pgram(x[i,:,j],fs=fs,nw=nw))[fi]
        end
    elseif nd == 2
        ps = mt_pgram(x[1,:],fs=fs,nw=nw)
        fi = map(f->freqrange[1]<=f<=freqrange[2],freq(ps))
        freqs = freq(ps)[fi];n = size(x,1)
        p = Array{Float64}(undef,n,length(freqs))
        p[1,:]=power(ps)[fi]
        for i in 2:n
            p[i,:]=power(mt_pgram(x[i,:],fs=fs,nw=nw))[fi]
        end
    else
        ps = mt_pgram(x,fs=fs,nw=nw)
        fi = map(f->freqrange[1]<=f<=freqrange[2],freq(ps))
        freqs = freq(ps)[fi];p = power(ps)[fi]
    end
    return p,freqs
end
