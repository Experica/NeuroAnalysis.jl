include("CSD.jl")

export rmline!,hlpass,epoch2samplerange,parsedigitalinanalog

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
function hlpass(y;low=Inf,high=0,fs=0)
    fy=copy(y)
    if high!=Inf && high>0 && fs>0
        f = digitalfilter(Highpass(high;fs=fs), Butterworth(4))
        for j=1:size(fy,1)
            fy[j,:]=filtfilt(f,fy[j,:])
        end
    end
    if low!=Inf && low >0 && fs>0
        f = digitalfilter(Lowpass(low;fs=fs), Butterworth(4))
        for j=1:size(fy,1)
            fy[j,:]=filtfilt(f,fy[j,:])
        end
    end
    fy
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
