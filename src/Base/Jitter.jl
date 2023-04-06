"""
Find and extract the maximal r-patterns from the sorted spike train
_____       ___     _       _
| | |       | |     |       |
 0 0    1    0   1      1 
"""
function findpattern(st;r=50)
    n=length(st)
    pbreak = findall(diff(st) .> r)
    prange = map((i,j)->i:j,[1;pbreak.+1],[pbreak;n])
end

"""
A single step in the (reverse) dynamic programming algorithm for sampler
"""
function pjdp(aj,bj,rk,αk,nk,βk)

    # change bj to the largest possible value
    bj = min(bj,αk+nk-rk-2)

    # check to see if impossible
    if isnan(αk) || bj < aj
        return (α = NaN, n = NaN, β = NaN)
    end

    ## SPECIAL CASES

    # check for fixed output spike
    if aj == bj
        return (α = aj, n = 1, β = 1)
    end

    # now the interval is determined
    α = aj
    n = bj-aj+1

    # check to see if no interaction (case E and fixed input spike)
    if bj+rk < αk
        # return (;α,n,β=(n:-1:1).'/n)
    end

    ## MAIN ALGORITHM

    # compute betaj

    # initialize
    β = zeros(n)

    betasum = 0;
    j2k = aj+rk-αk+1

    # just in case for some reason betak is not normalized
    betak1 = βk(1)  # this should be 1

    # loop in reverse
    for j = n:-1:1
        # get the index into betak
        k = j + j2k; # should always be <= nk
        if k <= 1
            betasum += betak1
        else
            betasum += betak(k)
        end
        betaj[j] = betasum
    end

    # normalize
    betaj ./= betasum
    return (;α, n, β)
end

"""
Draws n samples from the description in Sampler 
"""
function pjsample(n,S)
    ## precomputations

    # number of patterns
    m = length(S.p)

    # extract some useful info
    patn = zeros(m,1);
    for k = 1:m
        patn[k] = SS(k).patn;
    end
    patn = S.pn

    # where the patterns start and end
    ndx = [0 ; cumsum(patn)]

    # total number of spikes
    nt = ndx[end]

    nt == 0 && return []

    # initialize
    t = zeros(nt,n)

    ## sampling

    # loop over spikes
    for k = 1:m
        t0 = pjrng(-Inf,S(k))
        # add the pattern
        pat = SS(k).pat
        for j = 1:patn(k)
            t[ndx(k)+j,:] = t0+pat(j)
        end
    end
    t
end

"""
Samples the next spike given the previous spike time using the pattern jitter sampling
"""
function pjrng(tj,rk,αk,nk,βk)
    
    # check for a priori impossible
    if isnan(alphak)
        tk = nan(size(tj));
        return
    end

    # check for fixed spike time
    if nk == 1
        tk = tj + rk
        tk[tk >= alphak] = NaN
        tk[!isnan(tk)] = alphak
        return
    end

    # get the earliest next spike time
    # map to indices
    tj = tj+(rk+2-alphak);

    # deal with impossible given this spike time
    tj[tj > nk] = NaN

    # store the uniform random numbers in the output
    tk = rand(size(tj))

    # loop over spikes
    alphak1 = alphak-1
    for h = 1:numel(tj)

        tjh = tj(h);
        if isnan(tjh)== tk(h) == NaN
            continue
        end
        
        # rescale to account for previous spike
        if tjh > 1
            tkh = tk(h)*betak(tjh);
        else
            tkh = tk(h);
        end
        
        # figure out which segment
        for i = nk:-1:1
            if tkh <= betak(i)
                break
            end
        end
        if i < tjh
            [h tjh i]
            error("algorithm error 1")
        end
        
        # store the sample
        tk[h] = i + alphak1
    end
end

"""
"""
function patternjitter(st;r=50,l=25,n=5,win=:fix)
    prange = findpattern(st;r)
    np = length(prange)
    pdur = [st[last(i)]-st[first(i)] for i in prange]
    pn = map(length,prange)
    pstart = [st[first(i)] for i in prange]
    p = map((i,s)->st[i].-s,prange,pstart)

    # get the appropriate r's
    pr = [r; pdur[1:end-1].+r]

    # compute the appropriate intervals
    if win==:center
        # if L ~= round(L) || L < 1, error('bad L'), end
        #         if any(t ~= round(t)) || any(t0 ~= round(t0)), error('non integer t'), end
        #         if r ~= round(r), error('bad r'), end
        u = round(pstart - floor(l/2))
    elseif win==:fix
        u = round(floor(pstart/l)*l)
    else
        error("Unknown win=$win")
    end

    # sampler
    S = (r=pr,α=[],n=[],β=[],p=p,pn=pn)
    # initialize
    S.α[np], S.n[np], S.β[np] = pjdp(u(n),u(n)+L-1,-inf,inf,1,1)
    # dynamic programming
    for k = n-1:-1:1
        S.α[k], S.n[k], S.β[k] = pjdp(u(k),u(k)+L-1,S.r[k+1],S.α[k+1],S.n[k+1],S.β[k+1])
    end

    # sampling
    sst = pjsample(n,S)
end

"""
Independently and uniformly jitters the spike times in `st` over jitter
windows of length `l`(>0). This is repeated `n` times with the `kth` jittered 
spike train stored in `jst[:,k]`.

- win: 
    - :center - the jitter window for st[k] is st[k] + [-l/2, l/2]
    - :fix - the jitter window for st[k] is floor(st[k]/l)*l + [0, l]
- sort: if sorting jittered spike trains
"""
function spikejitter(st::AbstractVector{T};n::Signed=10,l::Tl=25,win::Symbol=:fix,sort::Bool=false) where {T<:Real,Tl<:Real}
    if T <: Integer && !(Tl <: Integer)
        @warn "Round l=$l to $T"
        l = round(T,l)
    end
    l <= 0 && error("l <= 0")

    if win==:center
        hl = l/2
        T <: Integer && (hl=floor(T,hl))
        u = st .- hl
    elseif win==:fix
        u = floor.(st./l).*l
    else
        error("Unknown win=$win")
    end

    jst = u .+ l*rand(length(u),n)
    sort && sort!(jst,dims=1)
    n==1 && (jst = dropdims(jst,dims=2))
    T <: Integer ? floor.(T,jst) : jst
end

"""
shuffle spikes in bined(1ms) spike trains `bst`(nbin x ntrial) between trials 
in fixed jitter windows of length `l`(>0ms)

(Smith, Matthew A., and Adam Kohn (2008). Spatial and Temporal Scales of Neuronal Correlation in Primary Visual Cortex. J. Neurosci. 28.48 : 12591-12603.)
"""
function shufflejitter(bst;l::Integer=25)
    l <= 0 && error("l <= 0")

    nbin,ntrial = size(bst)
    nl,rl = fldmod(nbin,l)
    lbi=[]
    nl>0 && append!(lbi,[(1:l) .+ l*i for i in 0:nl-1])
    rl>0 && push!(lbi,(nl*l + 1):nbin)

    # psth in jitter windows as resample distribution
    psth = dropdims(mean(bst,dims=2),dims=2)
    @views lbp = map(i->psth[i],lbi)
    lbp = map(i->weights(i./sum(i)),lbp) 
    # n spikes in jitter windows to resample
    @views lns = map(i->dropdims(sum(bst[i,:],dims=1),dims=1),lbi)
    
    jbst = zeros(nbin,ntrial)
    jbi = map((i,p,ns)->map(n->sample(i,p,Int(n),replace=false),ns),lbi,lbp,lns)
    for i in jbi
        for t in 1:ntrial
            jbst[i[t],t] .= 1
        end
    end
    jbst
end
