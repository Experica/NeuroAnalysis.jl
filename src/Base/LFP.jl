include("CSD.jl")

export rmline!,hlpass

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
        f = digitalfilter(Highpass(high;fs=fs), Butterworth(2))
        for j=1:size(fy,1)
            fy[j,:]=filtfilt(f,fy[j,:])
        end
    end
    if low!=Inf && low >0 && fs>0
        f = digitalfilter(Lowpass(low;fs=fs), Butterworth(2))
        for j=1:size(fy,1)
            fy[j,:]=filtfilt(f,fy[j,:])
        end
    end
    fy
end
