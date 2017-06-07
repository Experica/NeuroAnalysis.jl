export alphablend,alphamask,alphamask_disk,alphamask_gaussian,alphamask_diskfade,
clampscale,oiframeresponse,oiresponse,oicomplexmap,anglemode

using Colors,Images,ImageFiltering,DataFrames

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
    end
end
function alphamask_disk(src,radius)
    dim = size(src);dim1=dim[1];dim2=dim[2];mindim=min(dim1,dim2)
    mh = dim1/2;mw = dim2/2;dst = copy(src);didx=[]
    for i=1:dim1,j=1:dim2
        d = sqrt((i-mh)^2+(j-mw)^2)-radius*mindim
        if d>=0
            dst[i,j]=coloralpha(color(dst[i,j]),0)
        else
          didx=[didx;sub2ind(dim,i,j)]
        end
    end
    return dst,didx
end
function alphamask_gaussian(src,sigma)
    dim = size(src);dim1=dim[1];dim2=dim[2];mindim=min(dim1,dim2)
    mh = dim1/2;mw = dim2/2;dst = copy(src);didx=[]
    for i=1:dim1,j=1:dim2
        d = ((i-mh)^2+(j-mw)^2)/(0.5*mindim)^2
        dst[i,j]=coloralpha(color(dst[i,j]),alpha(dst[i,j])*exp(-d/(2*sigma^2)))
        didx=[didx;sub2ind(dim,i,j)]
    end
    return dst,didx
end
function alphamask_diskfade(src,radius,sigma)
    dim = size(src);dim1=dim[1];dim2=dim[2];mindim=min(dim1,dim2)
    mh = dim1/2;mw = dim2/2;dst = copy(src);didx=[]
    for i=1:dim1,j=1:dim2
        d = sqrt((i-mh)^2+(j-mw)^2)/mindim-radius
        if d>=0
            dst[i,j]=coloralpha(color(dst[i,j]),alpha(dst[i,j])*erfc(sigma*d))
        else
          didx=[didx;sub2ind(dim,i,j)]
        end
    end
    return dst,didx
end

function clampscale(x,min::Real,max::Real)
    scaleminmax(min,max).(x)
end
clampscale(x) = clampscale(x,extrema(x)...)
function clampscale(x,sdfactor)
    m=mean(x);sd=std(x)
    clampscale(x,m-sdfactor*sd,m+sdfactor*sd)
end
function oiframeresponse(frames;filter=Kernel.gaussian(2),frameindex=nothing)
    if frameindex==nothing
    imfilter(squeeze(sum(frames,3),3),filter)
    else
        b = squeeze(sum(frames[:,:,frameindex[1]],3),3)
        r = squeeze(sum(frames[:,:,frameindex[2]],3),3)
        imfilter((r-b)./b,filter)
    end
end
function oiresponse(response,stimuli;ustimuli=sort(unique(stimuli)),blankstimuli=0,ustimuliindexgroup=Any[find(ustimuli[ustimuli.!=blankstimuli])],sdfactor=3)
    if sdfactor==nothing
        meanresponse = map(i->squeeze(mean(cat(3,response[stimuli.==i]...),3),3),ustimuli)
    else
        meanresponse = map(i->squeeze(mean(cat(3,clampscale.(response[stimuli.==i],sdfactor)...),3),3),ustimuli)
    end
    blank = meanresponse[ustimuli.==blankstimuli][1]
    oimap = DataFrame(stimuli=ustimuli[ustimuli.!=blankstimuli])
    oimap[:map] = meanresponse[ustimuli.!=blankstimuli]
    oimap[:blankmap]=map(i->i./blank,oimap[:map])
    cocktail=Any[];cocktailmap=Any[]
    for ig in ustimuliindexgroup
        c = squeeze(mean(cat(3,oimap[:map][ig]...),3),3)
        cm = map(i->i./c,oimap[:map][ig])
        cocktail=cat(1,cocktail,Any[c])
        cocktailmap=cat(1,cocktailmap,cm)
    end
    oimap[:cocktailmap]=cocktailmap
    return blank,cocktail,oimap
end
function oicomplexmap(maps,stimuli,cond;filter=Kernel.gaussian(4),sdfactor=3)
    if haskey(cond,"Ori")
        fl = cond["Ori"]
        theta=deg2rad.(fl[stimuli])*2
    end
    cmap=squeeze(sum(cat(3,map((m,a)->Complex(cos(a),sin(a)).*(1-clampscale(imfilter(m,filter),sdfactor)),maps,theta)...),3),3)
    amap = angle.(cmap);amap[amap.<0]=amap[amap.<0]+2pi
    mmap = clampscale(abs.(cmap))
    return amap,mmap,theta[sortperm(stimuli)],fl[sort(stimuli)]
end
function anglemode(a,theta,ctheta)
    theta[findmin(abs(angle.(Complex(cos(a),sin(a))./ctheta)))[2]]
end
