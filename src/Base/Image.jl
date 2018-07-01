export alphablend,alphamask,alphamask_disk,alphamask_gaussian,alphamask_diskfade,
clampscale,oiframeresponse,oiresponse,oicomplexmap,anglemode,angleabs

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
    dim = size(src);dim1=dim[1];dim2=dim[2];mindim=min(dim1,dim2)
    mh = dim1/2;mw = dim2/2;dst = copy(src);didx=Int[]
    for i=1:dim1,j=1:dim2
        d = sqrt((i-mh)^2+(j-mw)^2)-radius*mindim
        if d>0
            dst[i,j]=coloralpha(color(dst[i,j]),0)
        else
            push!(didx,sub2ind(dim,i,j))
        end
    end
    return dst,didx
end
function alphamask_gaussian(src,sigma)
    dim = size(src);dim1=dim[1];dim2=dim[2];mindim=min(dim1,dim2)
    mh = dim1/2;mw = dim2/2;dst = copy(src);didx=Int[]
    for i=1:dim1,j=1:dim2
        d = ((i-mh)^2+(j-mw)^2)/(0.5*mindim)^2
        dst[i,j]=coloralpha(color(dst[i,j]),alpha(dst[i,j])*exp(-d/(2*sigma^2)))
        push!(didx,sub2ind(dim,i,j))
    end
    return dst,didx
end
function alphamask_diskfade(src,radius,sigma)
    dim = size(src);dim1=dim[1];dim2=dim[2];mindim=min(dim1,dim2)
    mh = dim1/2;mw = dim2/2;dst = copy(src);didx=Int[]
    for i=1:dim1,j=1:dim2
        d = sqrt((i-mh)^2+(j-mw)^2)/mindim-radius
        if d>0
            dst[i,j]=coloralpha(color(dst[i,j]),alpha(dst[i,j])*erfc(sigma*d))
        else
            push!(didx,sub2ind(dim,i,j))
        end
    end
    return dst,didx
end

function clampscale(x,min::Real,max::Real)
    scaleminmax(min,max).(x)
end
clampscale(x) = clampscale(x,extrema(x)...)
function clampscale(x,sdfactor)
    sdfactor <= 0 && return clampscale(x)
    m=mean(x);sd=std(x)
    clampscale(x,m-sdfactor*sd,m+sdfactor*sd)
end
function oiframeresponse(frames;frameindex=nothing,baseframeindex=nothing)
    if frameindex==nothing
        r = squeeze(sum(frames,3),3)
    else
        r = squeeze(sum(frames[:,:,frameindex],3),3)
    end
    if baseframeindex!=nothing
        r./=squeeze(sum(frames[:,:,baseframeindex],3),3)
        r-=1
    end
    return r
end
function oiresponse(response,stimuli;ustimuli=sort(unique(stimuli)),blankstimuli=0,
    stimuligroup=Any[find(ustimuli[ustimuli.!=blankstimuli])],filter=nothing,sdfactor=nothing)
    if filter==nothing
        if sdfactor==nothing
            rs = map(i->cat(3,response[stimuli.==i]...),ustimuli)
        else
            rs = map(i->cat(3,clampscale.(response[stimuli.==i],sdfactor)...),ustimuli)
        end
    else
        if sdfactor==nothing
            rs = map(i->cat(3,imfilter.(response[stimuli.==i],[filter])...),ustimuli)
        else
            rs = map(i->cat(3,clampscale.(imfilter.(response[stimuli.==i],[filter]),sdfactor)...),ustimuli)
        end
    end
    responsemean = map(i->squeeze(mean(i,3),3),rs)
    responsesd = map(i->squeeze(std(i,3),3),rs)
    responsen = map(i->size(i,3),rs)

    blank = responsemean[ustimuli.==blankstimuli][1]
    rindex = ustimuli.!=blankstimuli
    rmap=responsemean[rindex]
    cocktail=Any[];cocktailmap=Any[]
    for ig in stimuligroup
        c = squeeze(mean(cat(3,rmap[ig]...),3),3)
        cm = map(i->i./c,rmap[ig])
        cocktail=cat(1,cocktail,Any[c])
        cocktailmap=cat(1,cocktailmap,cm)
    end
    return blank,cocktail,DataFrame(stimuli=ustimuli[rindex],map=rmap,blankmap=map(i->i./blank,rmap),cocktailmap=cocktailmap),
    DataFrame(stimuli=ustimuli,map=responsemean,mapsd=responsesd,mapn=responsen)
end
function oicomplexmap(maps,angles;isangledegree=true,isangleinpi=true,presdfactor=3,filter=Kernel.DoG((3,30)),sufsdfactor=3)
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
        map!(i->clampscale(i,presdfactor),maps,maps)
    end
    if filter != nothing
        map!(i->imfilter(i,filter),maps,maps)
    end
    if sufsdfactor!=nothing
        map!(i->clampscale(i,sufsdfactor),maps,maps)
    end

    cmap=squeeze(sum(cat(3,map((m,a)->Complex(cos(a),sin(a)).*-m,maps,angles)...),3),3)
    amap,mmap = angleabs(cmap)
    return Dict("complex"=>cmap,"angle"=>amap,"abs"=>mmap,"rad"=>sort(angles),"deg"=>angledegree)
end
function angleabs(cmap)
    amap = angle.(cmap);amap[amap.<0]=amap[amap.<0]+2pi
    mmap = clampscale(abs.(cmap))
    return amap,mmap
end
function anglemode(a,theta)
    theta[findmin(abs.(angle.(Complex(cos(a),sin(a))./Complex.(cos.(theta),sin.(theta)))))[2]]
end
