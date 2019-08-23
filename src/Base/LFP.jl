include("CSD.jl")

export rmline!,hlpass,epoch2samplerange,parsedigitalinanalog,powerspectrum

"Remove line noise and its harmonics by notch filter"
function rmline!(y,fs;freq=60,nh=3,bw=3)
    for i=1:nh
        f = iirnotch(freq*i,bw;fs=fs)
        for j=1:size(y,1)
            y[j,:]=filtfilt(f,y[j,:])
        end
    end
end

"High pass and low pass filters"
function hlpass(y,fs;low=Inf,high=0)
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

function epoch2samplerange(epochs,fs)
    nepoch = size(epochs,1)
    epochis = floor.(Int,epochs.*fs)
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
