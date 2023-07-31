"""
Epoch of a Spike Train `x`, where `min <= x[i] < max`.

kwargs `isminzero`, `ismaxzero` and `shift` set epoch zero value to `min+shift` or `max+shift`.

return:
- y: epoch of the Spike Train
- n: number of elements in the epoch, or divided by duration(max-min) of the epoch, based on kwargs `israte`
- w: epoch window(min,max)
- i: indices of epoch elements in the original Spike Train, such that y = x[i]

See also: [`epochspiketrains`](@ref), [`epochspiketrainresponse`](@ref)
"""
function epochspiketrain(x,min::Real,max::Real;isminzero::Bool=false,ismaxzero::Bool=false,shift::Real=0,israte::Bool=false)
    if ismaxzero && isminzero
        error("Zero setting conflicts, only one of them is allowed to be true.")
    end
    i = findall(min .<= x .< max)
    w = (min, max)
    n = Float64(length(i))
    y = x[i]
    if isminzero
        y .-= min+shift
    end
    if ismaxzero
        y .-= max+shift
    end
    if israte
        n /= ((max-min)*SecondPerUnit)
    end
    return (;y,n,w,i)
end
function epochspiketrain(x,mins,maxs;isminzero::Bool=false,ismaxzero::Bool=false,shift::Real=0,israte::Bool=false)
    yn = length(mins)
    if yn != length(maxs)
        error("Length of mins and maxs do not match.")
    end
    ys = Vector{Vector{Float64}}(undef,yn)
    ns = Vector{Float64}(undef,yn)
    ws = Vector{Tuple{Float64,Float64}}(undef,yn)
    is = Vector{Vector{Int}}(undef,yn)
    for i in 1:yn
        ys[i],ns[i],ws[i],is[i] = epochspiketrain(x,mins[i],maxs[i];isminzero,ismaxzero,shift,israte)
    end
    return (y=ys,n=ns,w=ws,i=is)
end
epochspiketrain(x,minmaxs::AbstractMatrix;isminzero::Bool=false,ismaxzero::Bool=false,shift::Real=0,israte::Bool=false) = @views epochspiketrain(x,minmaxs[:,1],minmaxs[:,2];isminzero,ismaxzero,shift,israte)
"Epochs of a Spike Train in between binedges"
function epochspiketrain(x,binedges;isminzero::Bool=false,ismaxzero::Bool=false,shift::Real=0,israte::Bool=false)
    nbinedges = length(binedges); nbinedges<2 && error("Have $nbinedges binedges, need at least two binedges.")
    epochspiketrain(x,binedges[1:end-1],binedges[2:end];isminzero,ismaxzero,shift,israte)
end

"Epochs for each Spike Trains"
function epochspiketrains(xs,binedges;isminzero::Bool=false,ismaxzero::Bool=false,shift::Real=0,israte::Bool=false)
    tn = length(xs)
    yss = Vector{Vector{Vector{Float64}}}(undef,tn)
    nss = Vector{Vector{Float64}}(undef,tn)
    ws = nothing
    iss = Vector{Vector{Vector{Int}}}(undef,tn)
    for i in 1:tn
        yss[i],nss[i],ws,iss[i] = epochspiketrain(xs[i],binedges;isminzero,ismaxzero,shift,israte)
    end
    return (y=yss,n=nss,w=ws,i=iss)
end

"""
Response of each epoch of a Spike Train, could be mean firing rate or number of spikes based on kwarg `israte`.

See also: [`epochspiketrainresponse_ono`](@ref)
"""
epochspiketrainresponse(x,minmaxs::AbstractMatrix;israte::Bool=true) = epochspiketrain(x,minmaxs;israte)[2]
epochspiketrainresponse(x,mins,maxs;israte::Bool=true) = epochspiketrain(x,mins,maxs;israte)[2]

"""
Response of each epoch of a Spike Train, could be mean firing rate or number of spikes based on kwarg `israte`.

!!! note
    This is a faster(~500x) version compared to [`epochspiketrainresponse`](@ref), but only works when `x`, `mins` and `maxs` are ascendingly ordered, and each `maxs-mins` range are none-overlapping.
"""
epochspiketrainresponse_ono(x,minmaxs::AbstractMatrix;israte::Bool=true,isnan2zero::Bool=true) = @views epochspiketrainresponse_ono(x,minmaxs[:,1],minmaxs[:,2];israte,isnan2zero)
function epochspiketrainresponse_ono(x,mins,maxs;israte::Bool=true,isnan2zero::Bool=true)
    n = length(mins)
    if n != length(maxs)
        error("Length of mins and maxs do not match.")
    end
    ns = zeros(n)
    ni=1
    for v in x
        @label start
        if mins[ni]<=v
            if v<maxs[ni]
                ns[ni] += 1
            else
                ni+=1
                ni>n && break
                @goto start
            end
        end
    end
    if israte
        ns = ns ./ ((maxs.-mins).*SecondPerUnit)
        isnan2zero && replace!(ns,NaN=>0)
    end
    return ns
end

"Generate a Homogeneous Poisson Spike Train"
function poissonspiketrain(r,dur;t0=0.0,rp=2.0)
    isid = Exponential(1000/r)
    st = Float64[-rp]
    while true
        i = rand(isid)
        if i > rp
            s = i + st[end]
            if s <= dur
                push!(st,s)
            else
                break
            end
        end
    end
    deleteat!(st,1)
    return st.+t0
end
"Generate a Inhomogeneous Poisson Spike Train"
function poissonspiketrain(ir::Function,dur;t0=0.0,rp=2.0)
    st=Float64[-rp]
    for t in 0.05:0.1:dur
        rand()<=ir(t)/10000 && (t-st[end])>rp && push!(st,t)
    end
    deleteat!(st,1)
    return st.+t0
end

"Flat Spike Trains to Vector of SpikeTimes and Trials"
function flatspiketrains(sts::AbstractVector;trialorder=[])
    nt = length(sts)
    ntr = length(trialorder)
    if ntr==nt
        ti = sortperm(trialorder)
        order = sort(trialorder)
    else
        ntr>0 && @warn "Number of trials do not match, trial order ignored."
        ti = 1:nt
        order=[]
    end

    spike=Float64[];trial=Int[]
    for i in 1:nt
        st = sts[i];n=length(st)
        n==0 && continue
        append!(spike,st);append!(trial,fill(ti[i],n))
    end
    isempty(order) || (order = order[trial])
    return spike,trial,order
end
"Flat Spike Trains to Vector of SpikeTimes and Trials"
function flatspiketrains(sts::AbstractMatrix;trialorder=[])
    nt,n = size(sts)
    ntr = length(trialorder)
    if ntr==nt
        ti = sortperm(trialorder)
        order = sort(trialorder)
    else
        ntr>0 && @warn "Number of trials do not match, trial order ignored."
        ti = 1:nt
        order=[]
    end

    spike = sts[:]
    trial = repeat(ti,outer=n)
    isempty(order) || (order = order[trial])
    return spike,trial,order
end

"""
`Mean` and `SEM` of an Array along `dims`, 
optionally apply `sfun` on each slice along `dims` before and `mfun` on mean after.
"""
function meanse(x;dims=1,sfun=nothing,mfun=nothing)
    n = size(x,dims)
    if !isnothing(sfun)
        x = stack(sfun,eachslice(x;dims);dims)
    end
    m,se=dropdims.(mean_and_std(x,dims);dims)
    se /= sqrt(n)
    if !isnothing(mfun)
        m = mfun(m)
    end
    return (;m,se)
end

"PSTH of Spike Trains"
function psthspiketrains(xs,binedges;israte::Bool=true,ismean::Bool=true,sfun=nothing,mfun=nothing)
    nss = epochspiketrains(xs,binedges;israte)[2]
    halfbinwidth = (binedges[2]-binedges[1])/2
    x = binedges[1:end-1].+halfbinwidth
    mat = stack(nss,dims=1)
    if ismean
        m,se = meanse(mat;dims=1,sfun,mfun)
        return (;m,se,x)
    else
        return (;mat,x)
    end
end
