export subrv,subrvr,isi,poissonspike,flatrvv,histrv,histmatrix,psth

function subrv(rv::RealVector,min::Real,max::Real;isminzero::Bool=false,ismaxzero::Bool=false)
    if ismaxzero && isminzero
        error("Zero Setting Conflicts.")
    end
    i = find(min .<= rv .< max)
    n = length(i)
    y = rv[i]
    if isminzero
        y -= min
    end
    if ismaxzero
        y -= max
    end
    w = [min, max]
    return y,n,w,i
end
function subrv(rv::RealVector,mins::RealVector,maxs::RealVector;isminzero::Bool=false,ismaxzero::Bool=false,israte::Bool=false)
    yn = length(mins)
    if yn != length(maxs)
        error("Length of mins and maxs do not match.")
    end
    ys = Array{RealVector}(yn)
    ns = Array{Int}(yn)
    ws = Array{RealVector}(yn)
    is = Array{Vector{Int}}(yn)
    for i in 1:yn
        ys[i],ns[i],ws[i],is[i] = subrv(rv,mins[i],maxs[i],isminzero=isminzero,ismaxzero=ismaxzero)
    end
    if israte
        ns = ns./((maxs-mins)*SecondPerUnit)
    end
    return ys,ns,ws,is
end
function subrvr(rv::RealVector,mins::RealVector,maxs::RealVector)
    ys,ns,ws,is = subrv(rv,mins,maxs,israte=true)
    return ns
end

function isi(rv::RealVector)
    diff(sort(rv))
end

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
    return length(st)>1?st[2:end]+t0:Float64[]
end

function poissonspike(r::Function,dur::Real;t0::Real=0.0,rp::Real=2.0)
    st=Float64[-rp]
    for t in 0.05:0.1:dur
        rand()<=r(t)/10000 && (t-st[end])>rp && push!(st,t)
    end
    return length(st)>1?st[2:end]+t0:Float64[]
end

function flatrvv(rvv::RVVector,sv=[])
    nrv = length(rvv)
    if isempty(sv)
        issort=false
    elseif length(sv)==nrv
        issort=true
    else
        warn("Length of rvv and sv do not match, sorting ignored.")
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
        append!(x,rv);append!(y,ones(n)*i)
        issort && append!(s,fill(ssv[i],n))
    end
    return x,y,s
end

function histrv(rv::RealVector,binedges::RealVector;isminzero::Bool=false,ismaxzero::Bool=false)
    nbinedges = length(binedges)
    nbinedges<2 && error("Have $nbinedges binedges, Need at least two binedges.")
    subrv(rv,binedges[1:end-1],binedges[2:end],isminzero=isminzero,ismaxzero=ismaxzero)
end
function histrv(rv::RealVector,min::Real,max::Real;nbins::Integer=10,binwidth::Real=0.0,isminzero::Bool=false,ismaxzero::Bool=false)
    if binwidth <= 0.0
        binwidth = (max-min)/nbins
    end
    histrv(rv,min:binwidth:max,isminzero=isminzero,ismaxzero=ismaxzero)
end
histrv(rv::RealVector;nbins::Integer=10,binwidth::Real=0.0,isminzero::Bool=false,ismaxzero::Bool=false) = histrv(rv,minimum(rv),maximum(rv),nbins=nbins,binwidth=binwidth,isminzero=isminzero,ismaxzero=ismaxzero)
function histrv(rvs::RVVector,binedges::RealVector;isminzero::Bool=false,ismaxzero::Bool=false)
    yn = length(rvs)
    yn == 0 && error("Empty RVVector")
    ys = Array{RVVector}(yn)
    ns = Array{Vector{Int}}(yn)
    ws = []
    is = Array{Vector{Vector{Int}}}(yn)
    for i in 1:yn
        ys[i],ns[i],ws,is[i] = histrv(rvs[i],binedges,isminzero=isminzero,ismaxzero=ismaxzero)
    end
    return ys,ns,ws,is
end
function histrv(rvs::RVVector,min::Real,max::Real;nbins::Integer=10,binwidth::Real=0.0,isminzero::Bool=false,ismaxzero::Bool=false)
    if binwidth <= 0.0
        binwidth = (max-min)/nbins
    end
    histrv(rvs,min:binwidth:max,isminzero=isminzero,ismaxzero=ismaxzero)
end

function histmatrix(hv::Vector{Vector{Int}},ws::RVVector)
    hn = length(hv)
    nbins = length(ws)
    ((hn == 0) || (nbins == 0)) && error("Arguments Empty.")
    hn>0 && nbins>0 && nbins!=length(hv[1]) && error("nbins does not match.")
    binwidth = ws[1][2]-ws[1][1]
    hm = Array{Int}(hn,nbins)
    for i in 1:hn
        hm[i,:] = hv[i]
    end
    bincenters = Float64[ws[i][1]+binwidth/2.0 for i=1:nbins]
    return hm,bincenters
end
function histmatrix(rvs::RVVector,binedges::RealVector)
    ys,ns,ws,is = histrv(rvs,binedges)
    hm,x = histmatrix(ns,ws)
end

function psth(hm::Matrix{Int},x::RealVector;normfun=nothing)
    binwidth = x[2]-x[1]
    hmr = hm / (binwidth*SecondPerUnit)
    n = size(hmr,1)
    if normfun!=nothing
        for i=1:n
            hmr[i,:]=normfun(hmr[i,:])
        end
    end
    m = mean(hmr,1)[:]
    se = std(hmr,1)[:]/sqrt(n)
    return m,se,x
end
function psth(hv::Vector{Vector{Int}},ws::RVVector;normfun=nothing)
    hm,x = histmatrix(hv,ws)
    psth(hm,x,normfun=normfun)
end
function psth(rvs::RVVector,binedges::RealVector;normfun=nothing)
    hm,x = histmatrix(rvs,binedges)
    psth(hm,x,normfun=normfun)
end
