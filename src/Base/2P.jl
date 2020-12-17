

function sbxsubrm(rm,epochs,rois;fun=nothing)
    nepoch = size(epochs,1)   # number of trial
    minepochlength = floor(Int,minimum(diff(epochs,dims=2)))
    ys=Array{Float64}(undef,length(rois),minepochlength,nepoch)
    for i in 1:nepoch
        y = rm[rois,range(max(1,epochs[i,1]),length=minepochlength)]
        if !isnothing(fun)
            y=fun(y)
        end
        ys[:,:,i] = y
    end
    return nepoch==1 ? dropdims(ys,dims=3) : ys
end

function dFoF(bi)
        return x->begin
        bl = mean(x[:,bi],dims=2)
        x = x./bl .-1
        x
        end
    end

sbxcondin(ctc::Dict)=sbxcondin(DataFrame(ctc))
function sbxcondin(ctc::DataFrame)
    t = [ctc DataFrame(i=1:nrow(ctc))]
    t = by(t, names(ctc),g->DataFrame(n=nrow(g), i=[g[:,:i]]))
    sort!(t.n); sort!(t[1:end-1]); return t
end


"""
# Normalize values of an array to be between -1 and 1
# original sign of the array values is maintained.
"""
function normalized(data;equal=true)

    if isequal(equal,true)
        if abs(min(data...)) > max(data...)
            max_range_value = abs(min(data...))
            min_range_value = min(data...);
        else
            max_range_value = max(data...)
            min_range_value = -max(data...)
        end
        norm_value = 2 .* data / (max_range_value - min_range_value)
    elseif isequal(equal,false)
        max_range_value = max(data...)
        min_range_value = abs(min(data...))
        b=zeros(size(data))
        bw = data.> 0
        idx = findall(x->x==1,bw)
        b[idx] = data[idx] / max_range_value
        bw = data.< 0
        idx = findall(x->x==1,bw)
        b[idx] = data[idx] / min_range_value
        norm_value = b
    end
    return norm_value
end


"""
# Find the extreme value
"""
function sbxexd(stas,stamag,goodTimeidx;sig=true)
    bstdi = argmax(stamag)
    exi = [argmin(stas[:,:,bstdi]),argmax(stas[:,:,bstdi])]
    ex = stas[exi,bstdi]
    absexi = argmax(abs.(ex))
    (ex=sig ? ex[absexi] : 0, mag=stamag[bstdi], pd=goodTimeidx[bstdi])
end

sfdog(x,p...) = dogf(x;aₑ=p[1],μₑ=p[2],σₑ=p[3],aᵢ=p[4],μᵢ=p[5],σᵢ=p[6])


"""
# Peichao:
# Join the results from differernt Hartleys after hartleyFourier analysis
"""
function sbxjoinhartleyFourier(harts)
    cn = length(harts)
    uids = length(harts[1]["result"].cellId)
    kernsize = size(harts[1]["result"].kernnor[1])
    dataset = Dict()
    signifs=Dict();taumaxs=Dict();kstdmaxs=Dict();kdeltas=Dict();slambdas=Dict();
    orimaxs=Dict();sfmaxs=Dict();orimeans=Dict();sfmeans=Dict();
    kerns=Dict();kernraws=Dict();kernests=Dict();oris=Dict();oriidxs=Dict();
    orifits=Dict();sfidxs=Dict(); sffits=Dict();sfs=Dict();

    for u in 1:uids
        kern = Array{Float64}(undef,kernsize...,cn)  # Normlized kernels of all hartley
        kernraw = Array{Float64}(undef,kernsize...,cn)  # Raw kernels of all hartley
        kernest = Array{Float64}(undef,kernsize...,cn)  # Estimated kernels of all hartley
        signif = zeros(1,cn);taumax = zeros(1,cn);kstdmax = zeros(1,cn);
        kdelta = zeros(1,cn);slambda = zeros(1,cn);orimax = zeros(1,cn);
        sfmax = zeros(1,cn);orimean = zeros(1,cn);sfmean = zeros(1,cn);
        ori=[]; oriidx=[]; sf=[]; sfidx=[]; orifit=[]; sffit=[];

        for c in 1:cn
            signif[1,c] = harts[c]["result"].signif[u]
            taumax[1,c] = harts[c]["result"].taumax[u]
            kstdmax[1,c] = harts[c]["result"].kstdmax[u]
            kdelta[1,c] = harts[c]["result"].kdelta[u]
            slambda[1,c] = harts[c]["result"].slambda[u]
            orimax[1,c] = harts[c]["result"].orimax[u]
            sfmax[1,c] = harts[c]["result"].sfmax[u]
            orimean[1,c] = harts[c]["result"].orimean[u]
            sfmean[1,c] = harts[c]["result"].sfmean[u]

            kern[:,:,c] = harts[c]["result"].kernnor[u]
            kernraw[:,:,c] = harts[c]["result"].kernraw[u]
            kernest[:,:,c] = harts[c]["result"].kernest[u]

            θ = deg2rad.(harts[c]["result"].oriidx[u])
            fr = harts[c]["result"].oricurve[u]

            push!(oriidx,θ);
            push!(ori,fr);
            # fit von Mises for orientation
            fit=()
            try
                # vmfit = curve_fit((x,p)->vmf.(x,p...,n=2),θ,fr,[1.0,0,1])
                # if vmfit.converged
                #     x = 0:0.004:2π
                #     y = vmf.(x,vmfit.param...,n=2)
                #     fit = (fit...,circtuningstats(x,y,od=0.5π,s=:o)...,vm=vmfit)
                # end
                mfit = fitmodel(:vmn2,θ,fr)
                fit = (fit...,circtuningfeature(mfit,od=0.5π,fn=:o)...,vmn2=mfit)
            catch
            end
            push!(orifit,fit)

            fl = harts[c]["result"].sfidx[u]
            fr = harts[c]["result"].sfcurve[u]

            push!(sfidx,fl);
            push!(sf,fr);
            # fit difference of gaussians
            fit=()
            try
                # dogfit = curve_fit((x,p)->sfdog.(x,p...),fl,fr,[1.0,0,1,1,0,1])
                # if dogfit.converged
                #     fit = (dog=dogfit)
                # end
                mfit = fitmodel(:dog,fl,fr)
                fit = (sftuningfeature(mfit)...,dog=mfit)
            catch
            end
            push!(sffit,fit)

        end
        signifs[u] = signif
        taumaxs[u] = taumax
        kstdmaxs[u] = kstdmax
        kdeltas[u] = kdelta
        slambdas[u] = slambda
        orimaxs[u] = orimax
        sfmaxs[u] = sfmax
        orimeans[u] = orimean
        sfmeans[u] = sfmean

        kerns[u] = kern
        kernraws[u] = kernraw
        kernests[u] = kernest

        oriidxs[u] = oriidx
        oris[u] = ori
        orifits[u] = orifit
        sfidxs[u] = sfidx
        sfs[u] = sf
        sffits[u] = sffit
    end
    dataset["signif"] = signifs
    dataset["taumax"] = taumaxs
    dataset["kstdmax"] = kstdmaxs
    dataset["kdelta"] = kdeltas
    dataset["slambda"] = slambdas
    dataset["orimax"] = orimaxs
    dataset["sfmax"] = sfmaxs
    dataset["orimean"] = orimeans
    dataset["sfmean"] = sfmeans

    dataset["kern"] = kerns
    dataset["kernraw"] = kernraws
    dataset["kernest"] = kernests

    dataset["oriraw"] = oris
    dataset["oriidx"] = oriidxs
    dataset["orifit"] = orifits
    dataset["sfraw"] = sfs
    dataset["sfidx"] = sfidxs
    dataset["sffit"] = sffits

    return dataset
end


function sbxjoinsta(stas,lbTime,ubTime,blkTime)
    cn = length(stas)
    imagesize = stas[1]["imagesize"]
    maskradius = stas[1]["maskradius"]
    delays = stas[1]["delays"]
    stisize = stas[1]["stisize"]
    ppd = imagesize[1]/stisize[1]
    xi = stas[1]["xi"]     # pixels in stimulated region
    nxi = setdiff(1:prod(imagesize),xi)  # pxiels in background
    cii=CartesianIndices(imagesize)

    goodTimeidx=findall(x -> lbTime .<= x .<= ubTime, delays)
    bdi=findall(x -> (x .<=blkTime), delays)   # select as blanks

    dataset = Dict("stisize"=>stisize,"imagesize"=>imagesize,"maskradius"=>maskradius,"delays"=>delays,"ppd"=>ppd,"bdi"=>bdi,"xi"=>xi)
    # dataset["log"] = [i["log"] for i in stas]
    dataset["color"]=[i["color"] for i in stas]
    # dataset["minmaxcolor"] = [(i["mincolor"],i["maxcolor"]) for i in stas]
    # dataset["minmaxcg"] = map(i->minmaxcolorgradient(i...,n=256),dataset["minmaxcolor"])
    rawsta = Dict();usta = Dict();uzsta = Dict();ustamag = Dict();ustavec = Dict(); ucex=Dict();uresponsive=Dict()
    # uids = mapreduce(c->size(stas[c]["usta"],1),intersect,1:cn)
    uids = size(stas[1]["usta"],1)
    # uids = mapreduce(c->keys(stas[c]["usta"]),intersect,1:cn)
    for u in 1:uids
        rsta = Array{Float64}(undef,imagesize...,length(delays),cn)  # raw sta
        csta = Array{Float64}(undef,imagesize...,length(delays),cn)  # unit sta
        zsta = Array{Float64}(undef,imagesize...,length(delays),cn)  # z-scored sta
        # bsta = Array{Float64}(undef,imagesize...,length(bdi))  # blank sta
        cstamag = Array{Float64}(undef,length(delays),cn)
        cstavec = Array{Float64}(undef,length(xi),length(delays),cn)

        cex = []
        for c in 1:cn
            bsta = stas[c]["usta"][u,bdi,:]
            bm = mean(bsta);bsd = std(bsta)
            zzsta = zscore(stas[c]["usta"][u,:,:],bm,bsd)
            for d in eachindex(delays)
                rawstavec = stas[c]["usta"][u,d,:]  # raw sta of grating stim
                mag = norm(rawstavec)  # magnitude of the raw sta vector
                ucsta = rawstavec/mag   # unit vector sta
                cstamag[d,c] = mag
                csta[cii[nxi],d,c].= 0  # NonS-stimulated region, background(0 0 0)
                csta[cii[xi],d,c] = ucsta   # unit vector STA; Stimulated region,RF  # now STA is 2D matrix, not a vector
                cstavec[:,d,c] = ucsta
                rsta[cii[nxi],d,c].= 0
                rsta[cii[xi],d,c] = rawstavec  # raw vector sta
                zsta[cii[nxi],d,c].= 0
                zsta[cii[xi],d,c] = zzsta[d,:]  # z-scored sta
            end
            push!(cex,sbxexd(csta[:,:,goodTimeidx,c], cstamag[goodTimeidx,c],goodTimeidx))   #PL: To check the extrema of STA, I set the time window constrain here
        end
        rawsta[u] = rsta
        usta[u] = csta
        uzsta[u] = zsta
        ustamag[u] = cstamag
        ustavec[u] = cstavec
        ucex[u] = cex
        uresponsive[u] = false
    end
    dataset["rawsta"] = rawsta
    dataset["usta"] = usta
    dataset["zsta"] = uzsta
    dataset["ustamag"] = ustamag
    dataset["ustavec"] = ustavec
    dataset["ucex"] = ucex
    dataset["uresponsive"] = uresponsive
    return dataset
end


"peak ROI region and its delay"
function sbxpeakroi(clc)
    pi = [Tuple(argmax(clc))...]  # peak index
    plc = clc[:,:,pi[3:end]...]   # peak local contrast at best delay
    segs = seeded_region_growing(plc,[(CartesianIndex(1,1),1),(CartesianIndex(pi[1:2]...),2)])  # find ROI in STA image at best delay
    idx = findall(labels_map(segs).==2)   # find pixels in ROI
    idxlims = dropdims(extrema(mapreduce(i->[Tuple(i)...],hcat,idx),dims=2),dims=2)
    hw = map(i->i[2]-i[1],idxlims)
    center = round.(Int,mean.(idxlims))
    radius = round(Int,maximum(hw)/2)
    return (i=idx,center=center,radius=radius,pd=pi[3])
end
"check local mean or contrast in a region at a delay significently higher than baseline"
# function sbxisresponsive(sta,idx,d,bdi,badTimeidx;msdfactor=3.0,csdfactor=3.0)
#     d in badTimeidx && return false
#     blm = [mean(sta[idx,j]) for j in bdi]
#     bmm = mean(blm);bmsd=std(blm)
#     blc = [std(sta[idx,j]) for j in bdi]
#     bcm = mean(blc);bcsd=std(blc)
#     (std(sta[idx,d]) > bcm+csdfactor*bcsd) || (mean(sta[idx,d]) > bmm+msdfactor*bmsd)
# end

function sbxisresponsive(sta,idx,bdi,badTimeidx;mfactor=3.0,cfactor=3.0)
    lc = [std(sta[idx,j]) for j in 1:size(sta,3)]   # local variance of ROI for each delay
    lm = [mean(sta[idx,j]) for j in 1:size(sta,3)]  # local mean of ROI for each delay
    lcmaxd = argmax(lc); lcmax = lc[lcmaxd]  # maxi local variance
    lmmaxd = argmax(lm); lmmax = lm[lmmaxd]  # maxi local mean
    lmmind = argmin(lm); lmmin = lm[lmmind]  # min local mean
    bmm = mean(lm[bdi]);bmsd=std(lm[bdi])  # mean and std of local mean of blank condition
    bcm = mean(lc[bdi]);bcsd=std(lc[bdi])  # mean and std of local variance of blank condition

    (!(lcmaxd in badTimeidx) && lcmax > bcm+cfactor*bcsd) ||
    (!(lmmaxd in badTimeidx) && lmmax > bmm+mfactor*bmsd) ||
    (!(lmmind in badTimeidx) && lmmin < bmm-mfactor*bmsd)
end


function sbxresponsivesta!(dataset,lbTime,ubTime,respThres;ws=0.5,msdfactor=3.5,csdfactor=3.5)
    imagesize = dataset["imagesize"]
    maskradius = dataset["maskradius"]
    ppd=dataset["ppd"]
    xi = dataset["xi"]     # pixels in stimulated region
    bdi = dataset["bdi"]
    usta = dataset["usta"]
    uzsta = dataset["zsta"]
    ustamag = dataset["ustamag"]
    ucex = dataset["ucex"]
    delays = dataset["delays"]
    uresponsive=dataset["uresponsive"]

    center = [imagesize[1]÷2,imagesize[2]÷2]
    radius = Int64(ceil(max(imagesize[1],imagesize[2])*maskradius))
    goodTimeidx=findall(x -> lbTime .<= x .<= ubTime, delays)
    badTimeidx=findall(x -> (x .< lbTime || x .> ubTime), delays)

    ulsta=Dict();ulroi=Dict();ulcex=Dict();uconeresponsive=Dict();
    p = ProgressMeter.Progress(length(usta),desc="Test STAs ... ")
    for u in keys(usta)
        # u=671
        clc = localcontrast(usta[u],round(Int,ws*ppd))
        # cproi = map(c->sbxpeakroi(clc[:,:,:,c]),1:size(clc,4))
        cproi = map(c->sbxpeakroi(clc[:,:,:,c]),1:size(clc,4))
        # idx,c,cpd = sbxpeakroi(clc)

        # PL: Check responsiveness of each cone type and achromatic
        # PL: based on std and mean of STA image (compare with blank, which are STAs from the stim onset or before that), and delay range
        rs= [sbxisresponsive(uzsta[u][:,:,:,c],cproi[c].i,bdi,badTimeidx,mfactor=msdfactor,cfactor=csdfactor) for c in 1:size(clc,4)]
        # PL: based on threshold and delay range, to filter out low-response and cell 'response' too early or too late
        # rs2 = [((ustamag[u][cproi[i].pd,i]-mean(ustamag[u][bdi,i]))/mean(ustamag[u][bdi,i])>respThres) && (ucex[u][i].pd in goodTimeidx) for i in 1:size(clc,4)]
        rs2 = [(ustamag[u][ucex[u][i].pd,i]-mean(ustamag[u][bdi,i]))/mean(ustamag[u][bdi,i])>respThres  for i in 1:size(clc,4)]
        # rs2 = [(abs(ucex[u][i].ex)>respThres) && (ucex[u][i].pd in goodTimeidx) for i in 1:size(clc,4)]

        uconeresponsive[u] = rs .& rs2  # save cell's responesiveness to each cone & achromatic stimuli

        # PL: if cell response to one of cone stim and pass the threshold, it is responsive.
        if any(rs .& rs2)
            uresponsive[u] = true
            # vr = map(i->intersect(i.+(-radius:radius),1:imagesize[1]),center)
            # radius = (minimum ∘ mapreduce)((r,c)->abs.([r[1],r[end]].-c),append!,vr,center)
            # ulsta[u] = usta[u][xi]
            ulsta[u] = usta[u][map(i->i.+(-radius:radius),center)...,:,:]  # PL: I cut the STA image according to stim radius, not radius detected by localcontrast
            # ulroi[u] = (center=center,stisize=(2*radius+1)/ppd,d=cproi[i].pd,c=c)
            ulcex[u]=map(i->sbxexd(ulsta[u][:,:,goodTimeidx,i],ustamag[u][goodTimeidx,i],goodTimeidx,sig=uconeresponsive[u][i]),1:size(ulsta[u],4))
        else
            uresponsive[u] = false
        end
        next!(p)
    end
    dataset["ulsta"]=ulsta
    # dataset["ulroi"]=ulroi
    dataset["ulcex"]=ulcex
    dataset["uconeresponsive"] = uconeresponsive
    return dataset
end

rfgabor(x,y,p...) = gaborf(x,y,a=p[1],μ₁=p[2],σ₁=p[3],μ₂=p[4],σ₂=p[3]*p[5],θ=p[6],f=p[7],phase=p[8])
# rfdog(x,y,p...) = dogf(x,y,aₑ=p[1],μₑ₁=p[2],σₑ₁=p[3],μₑ₂=p[4],σₑ₂=p[3]*p[5],θₑ=p[6],aᵢ=p[7],μᵢ₁=p[2]+p[8],σᵢ₁=p[3]*p[9],μᵢ₂=p[4]+p[10],σᵢ₂=p[3]*p[9]*p[5],θᵢ=p[6])
rfdog(x,y,p...) = dogf(x,y,aₑ=p[1],μₑ₁=p[2],σₑ₁=p[3],μₑ₂=p[4],σₑ₂=p[3],θₑ=0,aᵢ=p[5],μᵢ₁=p[2],σᵢ₁=p[3]*p[6],μᵢ₂=p[4],σᵢ₂=p[3]*p[6],θᵢ=0)

function sbxmodelfit(data,ppd;model=:gabor)
    alb,aub = abs.(extrema(data))
    ab = max(alb,aub)
    rpx = (size(data)[1]-1)/2
    r = rpx/ppd

    x = (mapreduce(i->[i[2] -i[1]],vcat,CartesianIndices(data)) .+ [-(rpx+1) rpx+1])/ppd
    y = vec(data)
    rlt = mfun = missing
    try
        if model == :dog
            if aub >= alb
                ai = 3.5alb
                ae = aub + ai
                es = 0.2r;esl=0.15r;esu=0.3r
                ier=2;ierl = 1.1;ieru = 3
            else
                ae = 3.5aub
                ai = alb + ae
                es = 0.4r;esl=0.16r;esu=0.6r
                ier=0.5;ierl = 0.3;ieru = 0.9
            end
            # lb=[0,          -0.4sr,    0.1sr,   -0.4sr,    0.5,    0,     0,       -0.1sr,     0.1,    -0.1sr]
            # ub=[10,         0.4sr,    0.5sr,    0.4sr,    2,      π,     Inf,      0.1sr,     10,       0.1sr]
            # p0=[0,       0,        0.3sr,    0,        1,      π/4,   aei[2],    0,         0.25,       0]
            ub=[1.5ae,    0.6r,    esu,    0.6r,     1.5ai,    ieru]
            lb=[0.2ae,   -0.6r,    esl,   -0.6r,     0.5ai,    ierl]
            p0=[ae,       0,        es,     0,          ai,    ier]
            mfun = (x,p) -> rfdog.(x[:,1],x[:,2],p...)
        elseif model == :gabor
            ori,sf = freqimagestats(powerspectrum2(data,ppd)...)
            fub = min(1.5sf,8);flb=max(0.5sf,0.2)
            ub=[1.5ab,   0.6r,   1.0r,   0.6r,   2.3r,    π,    fub,   1]
            lb=[0.2ab,  -0.6r,   0.1r,  -0.6r,   0.3r,    0,    flb,   0]
            p0=[ab,      0,      0.2r,    0,     0.2r,  ori,     sf,   0]
            mfun = (x,p) -> rfgabor.(x[:,1],x[:,2],p...)
        end
        if !ismissing(mfun)
            mfit = curve_fit(mfun,x,y,p0,lower=lb,upper=ub,
            maxIter=3000,x_tol=1e-11,g_tol=1e-15,min_step_quality=1e-4,good_step_quality=0.25,lambda_increase=5,lambda_decrease=0.2)
            rlt = (model=model,radius=r,param=mfit.param,converged=mfit.converged,resid=mfit.resid,r=cor(y,mfun(x,mfit.param)))
        end
    catch
    end
    return rlt
end
"fit responsive sta to each type of models"
function sbxrffit!(dataset;model=[:gabor,:dog])
    ulsta = dataset["ulsta"]
    ulroi = dataset["ulroi"]
    ppd = dataset["ppd"]
    uconeresponsive = dataset["uconeresponsive"]

    if !haskey(dataset,"urf")
        dataset["urf"]=Dict()
    end
    urf = dataset["urf"]
    p = ProgressMeter.Progress(length(ulsta),desc="Fit RFs ... ")
    for u in keys(ulsta)
        if !haskey(urf,u)
            urf[u] = Dict()
        end
        exs=map(i->i.ex,dataset["ucex"][u])
        ds=map(i->i.d,dataset["ucex"][u])
        for m in model
            urf[u][m] = map(i->isequal(uconeresponsive[u][i],false) ? nothing : sbxmodelfit(ulsta[u][:,:,ds[i],i],ppd,model=m),1:length(exs))
        end
        next!(p)
    end
    return dataset
end

function sbxgetbestconesta(dataset)
    uconeresponsive = dataset["uconeresponsive"]
    delays = dataset["delays"]

    ublsta=Dict(); ubusta=Dict(); ubrawsta=Dict(); ubstamag=Dict(); ubcone=Dict(); ubmag=Dict();ustasign=Dict(); uachro=[];cellId=[];achroResp=DataFrame(); # save Cone Weights
    for u in keys(dataset["ulsta"])
        # check on/off dominance of a cell for each cone (inculde of achromatic), not restrict the delay is the same.
        # if any(uconeresponsive[u][1:4])
        exs=map(i->i.ex,dataset["ulcex"][u][1:4])
            # a = uconeresponsive[u]
            # b = convert(Array{Float64},a)
            # c = sign.(exs)
            # ustasign[u] = b.*c
        ustasign[u] = sign.(exs)
        # end
        # for calculating cone weight, so the time delay is restricted to the same delay.
        if any(uconeresponsive[u][1:3])
            exs=map(i->i.ex,dataset["ulcex"][u][1:3])
            mags=map(i->i.mag,dataset["ulcex"][u][1:3])
            ds=map(i->i.pd,dataset["ulcex"][u][1:3])

            # bestconeidx = findmax(abs.(exs))[2]  # define by max exs(extrem value)
            bestconeidx = findmax(mags)[2]  # define dominant cone by magnitude.
            bestdly = ds[bestconeidx]  # the best deley is just for best cone, not best deley for each cone.

            ulsta = dataset["ulsta"][u]
            usta = dataset["usta"][u]
            urawsta = dataset["rawsta"][u]
            ustamag = dataset["ustamag"][u]
            ublsta[u] = ulsta[:,:,bestdly,:]        # for ploting
            ubusta[u] = usta[:,:,bestdly,:]        # for cone weight calculation
            ubrawsta[u] = urawsta[:,:,bestdly,:]   # for cone weight calculation
            ubstamag[u] = ustamag[bestdly,:]
            ubcone[u] = (resp=sign.(exs), bstidx=bestconeidx, bstdly=delays[bestdly],bstex=exs[bestconeidx],bstmag=mags[bestconeidx])
        end
        if uconeresponsive[u][4]
            exs=map(i->i.ex,dataset["ulcex"][u][4:4])
            push!(uachro, exs[1]);push!(cellId,u)
        end
    end
    cellNum = length(cellId)
    achroResp.py = 0:cellNum-1
    achroResp.cellId=cellId
    achroResp.uachro=uachro
    dataset["ublsta"] = ublsta
    dataset["ubusta"] = ubusta
    dataset["ubrawsta"] = ubrawsta
    dataset["ubstamag"] = ubstamag
    dataset["ubcone"] = ubcone
    dataset["ustasign"] = ustasign
    dataset["achroResp"] = achroResp
    return dataset
end


function coneweight(lcone,mcone,scone)
    tcw = abs(lcone)+abs(mcone)+abs(scone)
    lcw = lcone/tcw
    mcw = mcone/tcw
    scw = scone/tcw
    return lcw,mcw,scw
end

function findextreme(a)
    b=findmax(a)
    c=findmin(a)
    if abs(b[1]) >= abs(c[1])
        d=a[b[2]]
    elseif abs(b[1]) < abs(c[1])
        d=a[c[2]]
    end
    return d
end

function sbxgetconeweight(dataset)
    ubrawsta = dataset["ubrawsta"]
    ubstamag = dataset["ubstamag"]
    ubusta = dataset["ubusta"]
    ubcone = dataset["ubcone"]


    cw = DataFrame()
    cellId=[];Lconemn=[];Mconemn=[];Sconemn=[];Lconemg=[];Mconemg=[];Sconemg=[];Lconemgall=[];Mconemgall=[];Sconemgall=[];

    for u in keys(ubcone)
        ## make STA mask
        if ubcone[u].bstex>=0
            img = ubusta[u][:,:,ubcone[u].bstidx]
        elseif ubcone[u].bstex<0
            img = -ubusta[u][:,:,ubcone[u].bstidx]
        end
        # smoothing
        imgf = imfilter(img, Kernel.gaussian(3))
        mask = Gray.(imgf).> max(imgf...)*0.95
        # calculate relative cone weights
        # coneresp=ubcone[u].resp[1:3]

        lconemn = sum(ubusta[u][:,:,1] .* mask)/sum(mask) #* abs(coneresp[1])
        mconemn = sum(ubusta[u][:,:,2] .* mask)/sum(mask) #* abs(coneresp[2])
        sconemn = sum(ubusta[u][:,:,3] .* mask)/sum(mask) #* abs(coneresp[3])
        lcwmn,mcwmn,scwmn = coneweight(lconemn,mconemn,sconemn)

        idx=findall(mask.==1)
        lconemg = norm(ubrawsta[u][:,:,1][idx]) * sign(findextreme(ubrawsta[u][:,:,1][idx]))
        mconemg = norm(ubrawsta[u][:,:,2][idx]) * sign(findextreme(ubrawsta[u][:,:,2][idx]))
        sconemg = norm(ubrawsta[u][:,:,3][idx]) * sign(findextreme(ubrawsta[u][:,:,3][idx]))
        lcwmg,mcwmg,scwmg = coneweight(lconemg,mconemg,sconemg)

        lconemgall = ubstamag[u][1] * sign(findextreme(ubrawsta[u][:,:,1][idx]))
        mconemgall = ubstamag[u][2] * sign(findextreme(ubrawsta[u][:,:,2][idx]))
        sconemgall = ubstamag[u][3] * sign(findextreme(ubrawsta[u][:,:,3][idx]))
        lcwmgall,mcwmgall,scwmgall = coneweight(lconemgall,mconemgall,sconemgall)

        push!(cellId, u); push!(Lconemn,lcwmn);push!(Mconemn,mcwmn);push!(Sconemn,scwmn);
        push!(Lconemg,lcwmg);push!(Mconemg,mcwmg);push!(Sconemg,scwmg);
        push!(Lconemgall,lcwmgall);push!(Mconemgall,mcwmgall);push!(Sconemgall,scwmgall);
    end
    cellNum = length(cellId)
    cw.py = 0:cellNum-1
    cw.cellId=cellId
    cw.Lconemn=Lconemn
    cw.Mconemn=Mconemn
    cw.Sconemn=Sconemn
    cw.Lconemg=Lconemg
    cw.Mconemg=Mconemg
    cw.Sconemg=Sconemg
    cw.Lconemgall=Lconemgall
    cw.Mconemgall=Mconemgall
    cw.Sconemgall=Sconemgall
    dataset["coneweight"]=cw

    return dataset
end
