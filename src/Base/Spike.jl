export subrv,subrvr,isi,poissonspike,flatrvv,histrv,histmatrix,psth

"Sub set of `RealVector` V, where `min <= V[i] < max`. KwArgs `isminzero` and `ismaxzero` set subrv zero to min and max."
function subrv(rv::RealVector,min::Real,max::Real;isminzero::Bool=false,ismaxzero::Bool=false,israte::Bool=false)
    if ismaxzero && isminzero
        error("Zero setting conflicts, only one of them is allowed to be true.")
    end
    i = findall(min .<= rv .< max)
    n = length(i)
    y = rv[i]
    if isminzero
        y .-= min
    end
    if ismaxzero
        y .-= max
    end
    if israte
        n /= ((max-min)*SecondPerUnit)
    end
    w = (min, max)
    return y,n,w,i
end
function subrv(rv::RealVector,mins::RealVector,maxs::RealVector;isminzero::Bool=false,ismaxzero::Bool=false,israte::Bool=false)
    yn = length(mins)
    if yn != length(maxs)
        error("Length of mins and maxs do not match.")
    end
    ys = Vector{RealVector}(undef,yn)
    ns = Vector{Real}(undef,yn)
    ws = Vector{Tuple{Real,Real}}(undef,yn)
    is = Vector{Vector{Int}}(undef,yn)
    for i in 1:yn
        ys[i],ns[i],ws[i],is[i] = subrv(rv,mins[i],maxs[i],isminzero=isminzero,ismaxzero=ismaxzero,israte=israte)
    end
    return ys,ns,ws,is
end
"Sub sets of `RealVector` in between binedges"
function subrv(rv::RealVector,binedges::RealVector;isminzero::Bool=false,ismaxzero::Bool=false,israte::Bool=false)
    nbinedges = length(binedges);nbinedges<2 && error("Have $nbinedges binedges, need at least two binedges.")
    subrv(rv,binedges[1:end-1],binedges[2:end],isminzero=isminzero,ismaxzero=ismaxzero,israte=israte)
end
function subrv(rvv::RVVector,binedges::RealVector;isminzero::Bool=false,ismaxzero::Bool=false,israte::Bool=false)
    yn = length(rvv)
    ys = Vector{RVVector}(undef,yn)
    ns = Vector{Vector{Real}}(undef,yn)
    ws = nothing
    is = Vector{Vector{Vector{Int}}}(undef,yn)
    for i in 1:yn
        ys[i],ns[i],ws,is[i] = subrv(rvv[i],binedges,isminzero=isminzero,ismaxzero=ismaxzero,israte=israte)
    end
    return ys,ns,ws,is
end

function subrvr(rv::RealVector,mins::RealVector,maxs::RealVector)
    _,ns,_,_ = subrv(rv,mins,maxs,israte=true)
    return ns
end

function isi(rv::RealVector)
    diff(sort(rv))
end

"Generate Poisson Spike Train"
function poissonspike(r::Real,dur::Real;t0::Real=0.0,rp::Real=2.0)
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
function poissonspike(ir::Function,dur::Real;t0::Real=0.0,rp::Real=2.0)
    st=Float64[-rp]
    for t in 0.05:0.1:dur
        rand()<=ir(t)/10000 && (t-st[end])>rp && push!(st,t)
    end
    deleteat!(st,1)
    return st.+t0
end

function flatrvv(rvv::RVVector,sv=[])
    nrv = length(rvv)
    if isempty(sv)
        issort=false
    elseif length(sv)==nrv
        issort=true
    else
        @warn """Length of "rvv" and "sv" do not match, sorting ignored."""
        issort=false
    end
    if issort
        srvv=rvv[sortperm(sv)]
        ssv=sort(sv)
    else
        srvv=rvv
        ssv=sv
    end
    x=[];y=[];s=[]
    for i in 1:nrv
        rv = srvv[i];n=length(rv)
        n==0 && continue
        append!(x,rv);append!(y,fill(i,n))
        issort && append!(s,fill(ssv[i],n))
    end
    return x,y,s
end

function histrv(rv::RealVector,min::Real=minimum(rv),max::Real=maximum(rv);nbins::Integer=10,binwidth::Real=0.0,isminzero::Bool=false,ismaxzero::Bool=false,israte::Bool=false)
    if binwidth <= 0.0
        binwidth = (max-min)/nbins
    end
    subrv(rv,min:binwidth:max,isminzero=isminzero,ismaxzero=ismaxzero,israte=israte)
end
function histrv(rvv::RVVector,min::Real,max::Real;nbins::Integer=10,binwidth::Real=0.0,isminzero::Bool=false,ismaxzero::Bool=false,israte::Bool=false)
    if binwidth <= 0.0
        binwidth = (max-min)/nbins
    end
    subrv(rvv,min:binwidth:max,isminzero=isminzero,ismaxzero=ismaxzero,israte=israte)
end

function histmatrix(hs::Vector{Vector{Real}},ws::Vector{Tuple{Real,Real}})
    hn = length(hs)
    nbins = length(ws)
    ((hn == 0) || (nbins == 0)) && error("Arguments Empty.")
    nbins!=length(hs[1]) && error("nbins does not match.")
    hm = Matrix{Real}(undef,hn,nbins)
    for i in 1:hn
        hm[i,:] = hs[i]
    end
    binwidth = ws[1][2]-ws[1][1]
    bincenters = [ws[i][1]+binwidth/2.0 for i=1:nbins]
    return hm,bincenters
end
function histmatrix(rvv::RVVector,binedges::RealVector;israte::Bool=false)
    ys,ns,ws,is = subrv(rvv,binedges,israte=israte)
    hm,x = histmatrix(ns,ws)
end

function psth(hm::Matrix{Real},x::RealVector;normfun=nothing)
    n = size(hm,1)
    if normfun!=nothing
        for i=1:n
            hm[i,:]=normfun(hm[i,:])
        end
    end
    m = dropdims(mean(hm,dims=1),dims=1)
    se = dropdims(std(hm,dims=1),dims=1)/sqrt(n)
    return m,se,x
end
function psth(hs::Vector{Vector{Real}},ws::Vector{Tuple{Real,Real}};normfun=nothing)
    hm,x = histmatrix(hs,ws)
    psth(hm,x,normfun=normfun)
end
function psth(rvv::RVVector,binedges::RealVector;israte::Bool=true,normfun=nothing)
    hm,x = histmatrix(rvv,binedges,israte=israte)
    psth(hm,x,normfun=normfun)
end
