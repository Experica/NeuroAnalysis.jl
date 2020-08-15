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
    n = length(i)
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
        ys[i],ns[i],ws[i],is[i] = epochspiketrain(x,mins[i],maxs[i],isminzero=isminzero,ismaxzero=ismaxzero,shift=shift,israte=israte)
    end
    return (y=ys,n=ns,w=ws,i=is)
end
epochspiketrain(x,minmaxs::AbstractMatrix;isminzero::Bool=false,ismaxzero::Bool=false,shift::Real=0,israte::Bool=false) = epochspiketrain(x,minmaxs[:,1],minmaxs[:,2],isminzero=isminzero,ismaxzero=ismaxzero,shift=shift,israte=israte)
"Epochs of a Spike Train in between binedges"
function epochspiketrain(x,binedges;isminzero::Bool=false,ismaxzero::Bool=false,shift::Real=0,israte::Bool=false)
    nbinedges = length(binedges); nbinedges<2 && error("Have $nbinedges binedges, need at least two binedges.")
    epochspiketrain(x,binedges[1:end-1],binedges[2:end],isminzero=isminzero,ismaxzero=ismaxzero,shift=shift,israte=israte)
end

"Epochs for each Spike Trains"
function epochspiketrains(xs,binedges;isminzero::Bool=false,ismaxzero::Bool=false,shift::Real=0,israte::Bool=false)
    tn = length(xs)
    yss = Vector{Vector{Vector{Float64}}}(undef,tn)
    nss = Vector{Vector{Float64}}(undef,tn)
    ws = nothing
    iss = Vector{Vector{Vector{Int}}}(undef,tn)
    for i in 1:tn
        yss[i],nss[i],ws,iss[i] = epochspiketrain(xs[i],binedges,isminzero=isminzero,ismaxzero=ismaxzero,shift=shift,israte=israte)
    end
    return (y=yss,n=nss,w=ws,i=iss)
end

"""
Response of each epoch of a Spike Train, could be mean firing rate or number of spikes based on kwarg `israte`.

See also: [`epochspiketrainresponse_ono`](@ref)
"""
epochspiketrainresponse(x,minmaxs::AbstractMatrix;israte::Bool=true) = epochspiketrain(x,minmaxs,israte=israte)[2]
epochspiketrainresponse(x,mins,maxs;israte::Bool=true) = epochspiketrain(x,mins,maxs,israte=israte)[2]

"""
Response of each epoch of a Spike Train, could be mean firing rate or number of spikes based on kwarg `israte`.

!!! note
    This is a faster(~500x) version compared to [`epochspiketrainresponse`](@ref), but only works when `x`, `mins` and `maxs` are ascendingly ordered, and each `maxs-mins` range are none-overlapping.
"""
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

function subrvr_onotest(rv::RealVector,mins::RealVector,maxs::RealVector;israte::Bool=true,isnan2zero::Bool=true,isbin::Bool=true)
    n = length(mins)
    if n != length(maxs)
        error("Length of mins and maxs do not match.")
    end
    ns = zeros(n)  # store FR in condition on/off bins
    ni=1  # bin counter
    sp=0
    for v in rv
        @label start
        if (mins[ni]<=v) && (v<maxs[ni])

        else
            ni+=1
            ni>n && break
            @goto start
        end
        ns[ni] += 1
    end
    if israte
        ns = ns ./ ((maxs.-mins).*SecondPerUnit)
        isnan2zero && replace!(ns,NaN=>0)
    end
    return ns
end

function subrmr_ono(rv::RealVector,mins::RealVector,maxs::RealVector;israte::Bool=true,isnan2zero::Bool=true)
end

"Generate a Poisson Spike Train"
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
function poissonspiketrain(ir::Function,dur;t0=0.0,rp=2.0)
    st=Float64[-rp]
    for t in 0.05:0.1:dur
        rand()<=ir(t)/10000 && (t-st[end])>rp && push!(st,t)
    end
    deleteat!(st,1)
    return st.+t0
end

"Flat Spike Trains to SpikeTimes and Trials, optionally sort trial based on `sv`"
function flatspiketrains(xs,sv=[])
    tn = length(xs)
    if isempty(sv)
        issort=false
    elseif length(sv)==tn
        issort=true
    else
        @warn """Length of "xs" and "sv" do not match, sorting ignored."""
        issort=false
    end
    if issort
        sxs=xs[sortperm(sv)]
        ssv=sort(sv)
    else
        sxs=xs
        ssv=sv
    end
    x=Float64[];y=Float64[];s=[]
    for i in 1:tn
        v = sxs[i];n=length(v)
        n==0 && continue
        append!(x,v);append!(y,fill(i,n))
        issort && append!(s,fill(ssv[i],n))
    end
    return x,y,s
end

"Vertical stack same length vectors to matrix"
function vstack(xs)
    tn = length(xs)
    n = length(xs[1])
    mat = Matrix{Float64}(undef,tn,n)
    for i in 1:tn
        mat[i,:] = xs[i]
    end
    return mat
end

"Vertical Mean and SEM of a matrix"
function vmeanse(mat::AbstractMatrix;normfun=nothing)
    n = size(mat,1)
    if !isnothing(normfun)
        for i=1:n
            mat[i,:]=normfun(mat[i,:])
        end
    end
    m = dropdims(mean(mat,dims=1),dims=1)
    se = dropdims(std(mat,dims=1),dims=1)/sqrt(n)
    return (;m,se)
end

"PSTH of Spike Trains"
function psthspiketrains(xs,binedges;israte::Bool=true,ismean::Bool=true,normfun=nothing)
    nss = epochspiketrains(xs,binedges,israte=israte)[2]
    halfbinwidth = (binedges[2]-binedges[1])/2
    x = binedges[1:end-1].+halfbinwidth
    mat = vstack(nss)
    if ismean
        m,se = vmeanse(mat,normfun=normfun)
        return (;m,se,x)
    else
        return (;mat,x)
    end
end
