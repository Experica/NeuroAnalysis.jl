function alphablend(src,dst,srcfactor,dstfactor=1-srcfactor)
    srcfactor.*src+dstfactor.*dst
end
function alphablend(src::Colorant,dst::Colorant)
    srcfactor = alpha(src)
    srcfactor.*src+(1-srcfactor).*dst
end
function alphamask(src;radius=0.5,sigma=0.15,masktype="Disk")
    if masktype=="Disk"
        alphamask_disk(src,radius)
    elseif masktype=="Gaussian"
        alphamask_gaussian(src,sigma)
    elseif masktype=="DiskFade"
        alphamask_diskfade(src,radius,sigma)
    else
        return copy(src),Int[]
    end
end
function alphamask_disk(src,radius)
    dims = size(src);dim1=dims[1];dim2=dims[2];mindim=min(dim1,dim2)
    hh = dim1/2;hw = dim2/2;dst = copy(src);unmaskidx=Int[];li = LinearIndices(dims)
    for i=1:dim1,j=1:dim2
        d = sqrt((i-hh)^2+(j-hw)^2)-radius*mindim
        if d>0
            dst[i,j]=coloralpha(color(dst[i,j]),0)
        else
            push!(unmaskidx,li[i,j])
        end
    end
    return dst,unmaskidx
end
function alphamask_gaussian(src,sigma)
    dims = size(src);dim1=dims[1];dim2=dims[2];mindim=min(dim1,dim2)
    hh = dim1/2;hw = dim2/2;dst = copy(src);unmaskidx=Int[];li = LinearIndices(dims)
    for i=1:dim1,j=1:dim2
        d = ((i-hh)^2+(j-hw)^2)/(0.5*mindim)^2
        dst[i,j]=coloralpha(color(dst[i,j]),alpha(dst[i,j])*exp(-d/(2*sigma^2)))
        push!(unmaskidx,li[i,j])
    end
    return dst,unmaskidx
end
function alphamask_diskfade(src,radius,sigma)
    dims = size(src);dim1=dims[1];dim2=dims[2];mindim=min(dim1,dim2)
    hh = dim1/2;hw = dim2/2;dst = copy(src);unmaskidx=Int[];li = LinearIndices(dims)
    for i=1:dim1,j=1:dim2
        d = sqrt((i-hh)^2+(j-hw)^2)/mindim-radius
        if d>0
            dst[i,j] = coloralpha(color(dst[i,j]),alpha(dst[i,j])*erfc(sigma*d))
        else
            push!(unmaskidx,li[i,j])
        end
    end
    return dst,unmaskidx
end

function clampscale(x,min::Real,max::Real)
    scaleminmax(min,max).(x)
end
clampscale(x) = clampscale(x,extrema(x)...)
function clampscale(x,sdfactor)
    m=mean(x);sd=std(x)
    clampscale(x,m-sdfactor*sd,m+sdfactor*sd)
end
function oiframeresponse(frames;frameindex=nothing,baseframeindex=nothing)
    if frameindex==nothing
        r = dropdims(sum(frames,dims=3),dims=3)
    else
        r = dropdims(sum(frames[:,:,frameindex],dims=3),dims=3)
    end
    if baseframeindex!=nothing
        r./=dropdims(sum(frames[:,:,baseframeindex],dims=3),dims=3)
        r.-=1
    end
    return r
end
function oiresponse(response,stimuli;ustimuli=sort(unique(stimuli)),blankstimuli=0,
    stimuligroup=Any[1:length(findall(ustimuli.!=blankstimuli))],filter=nothing,sdfactor=nothing)
    if filter==nothing
        if sdfactor==nothing
            rs = map(i->cat(response[stimuli.==i]...,dims=3),ustimuli)
        else
            rs = map(i->cat(clampscale.(response[stimuli.==i],sdfactor)...,dims=3),ustimuli)
        end
    else
        if sdfactor==nothing
            rs = map(i->cat(imfilter.(response[stimuli.==i],[filter])...,dims=3),ustimuli)
        else
            rs = map(i->cat(clampscale.(imfilter.(response[stimuli.==i],[filter]),sdfactor)...,dims=3),ustimuli)
        end
    end
    responsemean = map(i->dropdims(mean(i,dims=3),dims=3),rs)
    responsesd = map(i->dropdims(std(i,dims=3),dims=3),rs)
    responsen = map(i->size(i,3),rs)

    blank = responsemean[findfirst(ustimuli.==blankstimuli)]
    rindex = ustimuli.!=blankstimuli
    rmap=responsemean[rindex]
    cocktail=Any[];cocktailmap=Any[]
    for ig in stimuligroup
        c = dropdims(mean(cat(rmap[ig]...,dims=3),dims=3),dims=3)
        cm = map(i->i./c,rmap[ig])
        push!(cocktail,c)
        append!(cocktailmap,cm)
    end
    return blank,cocktail,DataFrame(stimuli=ustimuli[rindex],map=rmap,blankmap=map(i->i./blank,rmap),cocktailmap=cocktailmap),
    DataFrame(stimuli=ustimuli,map=responsemean,mapsd=responsesd,mapn=responsen)
end
function oicomplexmap(maps,angles;isangledegree=true,isangleinpi=true,presdfactor=3,filter=Kernel.DoG((3,3),(30,30)),sufsdfactor=3)
    if isangledegree
        angledegree = sort(angles)
        angles = deg2rad.(angles)
    else
        angledegree = sort(rad2deg.(angles))
    end
    if isangleinpi
        angles *=2
    end
    if presdfactor!=nothing
        maps=map(i->clampscale(i,presdfactor),maps)
    end
    if filter != nothing
        maps=map(i->imfilter(i,filter),maps)
    end
    if sufsdfactor!=nothing
        maps=map(i->clampscale(i,sufsdfactor),maps)
    end

    cmap=dropdims(sum(cat(map((m,a)->Complex(cos(a),sin(a)).*-m,maps,angles)...,dims=3),dims=3),dims=3)
    amap,mmap = angleabs(cmap)
    return Dict("complex"=>cmap,"angle"=>amap,"abs"=>mmap,"rad"=>sort(angles),"deg"=>angledegree)
end
function angleabs(cmap)
    amap = angle.(cmap);amap[amap.<0]=amap[amap.<0] .+ 2pi
    mmap = clampscale(abs.(cmap))
    return amap,mmap
end
anglemode(a,theta) = theta[findclosestangle(a,theta)]
findclosestangle(a,theta) = argmin(abs.(angle.(Complex(cos(a),sin(a))./Complex.(cos.(theta),sin.(theta)))))

"""
Generate Grating Image.

- θ: Orientation (radius)
- sf: SpatialFreq (cycle/deg)
- phase: Sin wave phase in [0,1]
- stisize: Tuple of image size in visual degree
- ppd: Pixel per degree

return image in [0,1]
"""
function grating(;θ=0,sf=1,phase=0,tf=1,t=0,stisize=(10,10),ppd=50)
    pc = round.(Int,stisize.*ppd./2)
    psize = pc.*2
    g = fill(0.5,psize)
    isnan(θ) && return g
    sinv,cosv = sincos(θ)
    for i in 1:psize[1], j in 1:psize[2]
        u = (j-pc[2])/pc[2]/2
        v = (pc[1]-i)/pc[1]/2
        d = cosv * v * stisize[1] - sinv * u * stisize[2]
        g[i,j] = (sin(2π * (sf * d - tf * t + phase)) + 1) / 2
    end
    return g
end

function powerspectrum(x::AbstractMatrix,fs;freqrange=[-15,15])
    ps = periodogram(x,fs=fs)
    p = power(ps)
    freq1,freq2 = freq(ps)
    fi1 = map(f->freqrange[1]<=f<=freqrange[2],freq1)
    fi2 = map(f->freqrange[1]<=f<=freqrange[2],freq2)
    p = p[fi1,fi2];freq1 = freq1[fi1];freq2 = freq2[fi2]
    si1=sortperm(freq1);si2=sortperm(freq2)
    return p[si1,si2],freq1[si1],freq2[si2]
end

function freqimagestats(x,f1,f2)
    f0i = (findfirst(f1.==0),findfirst(f2.==0))
    p = copy(x);p[f0i...]=0
    mi = argmax(p)
    mf = [f2[mi[2]],f1[mi[1]]]
    sf = norm(mf)
    ori = mod(atan(mf...),π)
    return ori,sf
end
