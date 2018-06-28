export batchtests,processtest,processori,processlaser,processimage,processlaserimage

using ProgressMeter,ePPR

function batchtests(tests::DataFrame,dataroot,resultroot,datatype...;condroot=Dict{String,Any}("rootdir"=>""),delay=20,binwidth=10,isplot=true)
    udf=[];cdf=[]
    p = ProgressMeter.Progress(nrow(tests),1,"Batch Tests ",50)
    for t in eachrow(tests)
        try
            u,c=processtest(prepare(abspath(dataroot,t[:files]),datatype...),resultroot,
            uuid=t[:UUID],condroot=condroot,delay=delay,binwidth=binwidth,isplot=isplot)
            push!(udf,u)
            push!(cdf,c)
        catch exc
            display.(catch_stacktrace())
        end
        next!(p)
    end
    return vcat(udf...),vcat(cdf...)
end

function processtest(dataset::Dict,resultroot;uuid="",condroot=Dict{String,Any}("rootdir"=>""),delay=20,binwidth=10,isplot=true)
    if haskey(dataset,"ex")
        testid = dataset["ex"]["ID"]
        if testid=="OriGrating"
            processori(dataset,resultroot,uuid=uuid,delay=delay,binwidth=binwidth,isplot=isplot)
        elseif testid=="Laser"
                processlaser(dataset,resultroot,uuid=uuid,delay=delay,binwidth=binwidth,isplot=isplot)
        elseif testid=="Image"
            processimage(dataset,condroot,resultroot,uuid=uuid,delay=delay,binwidth=binwidth,isplot=isplot)
        elseif testid=="LaserImage"
            processlaserimage(dataset,condroot,resultroot,uuid=uuid,delay=delay,binwidth=binwidth,isplot=isplot)
        end
    end
end

function processlaserimage(dataset::Dict,condroot::Dict{String,Any},resultroot;uuid="",delay=20,binwidth=10,minpredur=10,mincondtest=12000,
    nscale=2,downsample=2,sigma=1.5,pixelscale=255,isplot=true,forceprocess=true)
    ex = dataset["ex"];envparam = ex["EnvParam"];preicidur = ex["PreICI"];conddur = ex["CondDur"];suficidur = ex["SufICI"]
    bgcolor=RGBA(getparam(envparam,"BGColor")...)
    imagesetname = replace(getparam(envparam,"ImageSet","ImageQuad"),"Ã—","_")
    imagemasktype = getparam(envparam,"MaskType","ImageQuad")
    imagemaskradius = getparam(envparam,"MaskRadius","ImageQuad")
    imagemasksigma = getparam(envparam,"Sigma","ImageQuad")
    ct,ctc=ctctc(ex)

    lfactor = filter(f->f!=:Image,names(ctc))
    lcond = condin(ctc[:,lfactor])
    all(lcond[:n] .< mincondtest) && return [],[]

    if !haskey(condroot,imagesetname) && haskey(condroot,"rootdir") && isdir(condroot["rootdir"])
        pyramid = Dict{Symbol,Any}(:pyramid => map(i->gaussian_pyramid(i, nscale-1, downsample, sigma),
        loadimageset(joinpath(condroot["rootdir"],imagesetname),alpha=true)))
        pyramid[:size] = map(i->size(i),pyramid[:pyramid][1])
        condroot[imagesetname] = pyramid
    end
    imageset = condroot[imagesetname]
    bgimagecolor = oftype(imageset[:pyramid][1][1][1],bgcolor)
    unmaskindex = map(i->alphamask(i,radius=imagemaskradius,sigma=imagemasksigma,masktype=imagemasktype)[2],imageset[:pyramid][1])
    imagestimuli = map(s->map(i->alphablend.(alphamask(i[s],radius=imagemaskradius,sigma=imagemasksigma,masktype=imagemasktype)[1],[bgimagecolor]),imageset[:pyramid]),1:nscale)
    
    s = 2
    ximagesize = imageset[:size][s]
    xi = unmaskindex[s]
    imagestimulimatrix = Array{Float64}(length(imagestimuli[s]),prod(ximagesize))
    for i in 1:size(imagestimulimatrix,1)
        imagestimulimatrix[i,:] = vec(gray.(imagestimuli[s][i]))
    end

    predur = max(preicidur,minpredur)
    resultroot=abspath(resultroot,ex["ID"])
    !isdir(resultroot) && mkpath(resultroot)
    
    udf = []
    if haskey(dataset,"spike")
        spike = dataset["spike"];spikeeid = spike["electrodeid"];spikeuuid = spike["uuid"];spikeuid = spike["unitid"];spiketime = spike["time"]
        for e in spikeeid
            ei = findfirst(spikeeid.==e);est = spiketime[ei];esu = spikeuid[ei];euuid = spikeuuid[ei]
            for u in euuid
                preurs = subrvr(est[esu.==u],ct[:CondOn]+delay-predur,ct[:CondOn]+delay)
                y = subrvr(est[esu.==u],ct[:CondOn]+delay,ct[:CondOff]+delay)
                !isresponsive(preurs,y,lcond[:i]) && continue

                udir = joinpath(resultroot,"$(uuid)_E$(e)_U$(u)")
                !forceprocess && isdir(udir) && continue
                plotdir = isplot?udir:nothing
                mrs = []
                for l in eachrow(lcond)
                    debug = isplot?ePPRDebugOptions(level=DebugVisual,logdir=joinpath(plotdir,condstring(l,lfactor))):ePPRDebugOptions()
                    hp = ePPRHyperParams(ximagesize...,xindex=xi,ndelay=1,blankcolor=gray(bgimagecolor)*pixelscale)
                    hp.nft = [6]
                    hp.lambda = 30000
                    model,models = epprcv(imagestimulimatrix[Int.(ctc[:Image][l[:i]]),:]*pixelscale,y[l[:i]],hp,debug)

                    if isplot && model!=nothing
                        debug(plotmodel(model,hp),log="Model_Final")
                    end
                    push!(mrs,(clean(model),clean.(models),hp,l))
                end

                push!(udf,DataFrame(UUID=uuid,e=e,u=u,modelresults=mrs))
            end
        end
    end

    return vcat(udf...),[]
end

function processimage(dataset::Dict,condroot::Dict{String,Any},resultroot;uuid="",delay=20,binwidth=10,minpredur=10,mincondtest=12000,
    nscale=2,downsample=2,sigma=1.5,pixelscale=255,isplot=true,forceprocess=true)
    ex = dataset["ex"];envparam = ex["EnvParam"];preicidur = ex["PreICI"];conddur = ex["CondDur"];suficidur = ex["SufICI"]
    bgcolor=RGBA(getparam(envparam,"BGColor")...)
    imagesetname = replace(getparam(envparam,"ImageSet","ImageQuad"),"Ã—","_")
    imagemasktype = getparam(envparam,"MaskType","ImageQuad")
    imagemaskradius = getparam(envparam,"MaskRadius","ImageQuad")
    imagemasksigma = getparam(envparam,"Sigma","ImageQuad")
    ct,ctc=ctctc(ex)
    nrow(ct) < mincondtest && return [],[]

    if !haskey(condroot,imagesetname) && haskey(condroot,"rootdir") && isdir(condroot["rootdir"])
        pyramid = Dict{Symbol,Any}(:pyramid => map(i->gaussian_pyramid(i, nscale-1, downsample, sigma),
        loadimageset(joinpath(condroot["rootdir"],imagesetname),alpha=true)))
        pyramid[:size] = map(i->size(i),pyramid[:pyramid][1])
        condroot[imagesetname] = pyramid
    end
    imageset = condroot[imagesetname]
    bgimagecolor = oftype(imageset[:pyramid][1][1][1],bgcolor)
    unmaskindex = map(i->alphamask(i,radius=imagemaskradius,sigma=imagemasksigma,masktype=imagemasktype)[2],imageset[:pyramid][1])
    imagestimuli = map(s->map(i->alphablend.(alphamask(i[s],radius=imagemaskradius,sigma=imagemasksigma,masktype=imagemasktype)[1],[bgimagecolor]),imageset[:pyramid]),1:nscale)
    
    s = 2
    ximagesize = imageset[:size][s]
    xi = unmaskindex[s]
    imagestimulimatrix = Array{Float64}(length(imagestimuli[s]),prod(ximagesize))
    for i in 1:size(imagestimulimatrix,1)
        imagestimulimatrix[i,:] = vec(gray.(imagestimuli[s][i]))
    end
    x = imagestimulimatrix[Int.(ctc[:Image]),:]*pixelscale;

    predur = max(preicidur,minpredur)
    resultroot=abspath(resultroot,ex["ID"])
    !isdir(resultroot) && mkpath(resultroot)
    
    udf = []
    if haskey(dataset,"spike")
        spike = dataset["spike"];spikeeid = spike["electrodeid"];spikeuuid = spike["uuid"];spikeuid = spike["unitid"];spiketime = spike["time"]
        for e in spikeeid
            ei = findfirst(spikeeid.==e);est = spiketime[ei];esu = spikeuid[ei];euuid = spikeuuid[ei]
            for u in euuid
                preurs = subrvr(est[esu.==u],ct[:CondOn]+delay-predur,ct[:CondOn]+delay)
                y = subrvr(est[esu.==u],ct[:CondOn]+delay,ct[:CondOff]+delay)
                !isresponsive(preurs,y) && continue

                udir = joinpath(resultroot,"$(uuid)_E$(e)_U$(u)")
                !forceprocess && isdir(udir) && continue
                plotdir = isplot?udir:nothing
                debug = isplot?ePPRDebugOptions(level=DebugVisual,logdir=plotdir):ePPRDebugOptions()
                hp = ePPRHyperParams(ximagesize...,xindex=xi,ndelay=1,blankcolor=gray(bgimagecolor)*pixelscale)
                hp.nft = [3,3,3]
                hp.lambda = 30000
                model,models = epprcv(x,y,hp,debug)

                if isplot && model!=nothing
                    debug(plotmodel(model,hp),log="Model_Final")
                end

                push!(udf,DataFrame(UUID=uuid,e=e,u=u,model=clean(model),models=clean.(models),hp=hp))
            end
        end
    end

    return vcat(udf...),[]
end

function processori(dataset::Dict,resultroot;uuid="",delay=20,binwidth=10,minpredur=100,isplot=true)
    ex = dataset["ex"];envparam = ex["EnvParam"];preicidur = ex["PreICI"];conddur = ex["CondDur"];suficidur = ex["SufICI"]
    ct,ctc=ctctc(ex)
    fl,fln,fli=flni(ctc)
    cond=condni(ctc)

    ff = finalfactor(cond)[1]
    oris = ctc[ff]
    predur = max(preicidur,minpredur)
    resultroot=abspath(resultroot,ex["ID"])
    !isdir(resultroot) && mkpath(resultroot)
    
    udf = []
    if haskey(dataset,"spike")
        spike = dataset["spike"];spikeeid = spike["electrodeid"];spikeuuid = spike["uuid"];spikeuid = spike["unitid"];spiketime = spike["time"]
        for e in spikeeid
            ei = findfirst(spikeeid.==e);est = spiketime[ei];esu = spikeuid[ei];euuid = spikeuuid[ei]
            for u in euuid
                preurs = subrvr(est[esu.==u],ct[:CondOn]+delay-predur,ct[:CondOn]+delay)
                urs = subrvr(est[esu.==u],ct[:CondOn]+delay,ct[:CondOff]+delay)
                !isresponsive(preurs,urs,cond[:i]) && continue

                if isplot
                    plotname = "$(uuid)_E$(e)_U$(u)"
                    plotcondresponse(urs,cond,u,title=plotname,legend=:none)
                    png(joinpath(resultroot,plotname))
                end

                stats = statsori(urs,oris)

                push!(udf,DataFrame(UUID=uuid,e=e,u=u,dcv=stats[:dcv],pdir=stats[:pdir],ocv=stats[:ocv],pori=stats[:pori]))
            end
        end
    end

    cond[:UUID]=uuid
    return vcat(udf...),cond
end

function processlaser(dataset::Dict,resultroot;uuid="",delay=20,binwidth=10,minpredur=100,isplot=true)
    ex = dataset["ex"];envparam = ex["EnvParam"];preicidur = ex["PreICI"];conddur = ex["CondDur"];suficidur = ex["SufICI"]
    ct,ctc=ctctc(ex)
    cond=condin(ctc)

    predur = max(preicidur,minpredur)
    resultroot=abspath(resultroot,ex["ID"])
    !isdir(resultroot) && mkpath(resultroot)
    
    udf = []
    if haskey(dataset,"spike")
        spike = dataset["spike"];spikeeid = spike["electrodeid"];spikeuuid = spike["uuid"];spikeuid = spike["unitid"];spiketime = spike["time"]
        for e in spikeeid
            ei = findfirst(spikeeid.==e);est = spiketime[ei];esu = spikeuid[ei];euuid = spikeuuid[ei]
            for u in euuid
                preurs = subrvr(est[esu.==u],ct[:CondOn]+delay-predur,ct[:CondOn]+delay)
                urs = subrvr(est[esu.==u],ct[:CondOn]+delay,ct[:CondOff]+delay)
                !isresponsive(preurs,urs,cond[:i]) && continue
                ksadp = pvalue(KSampleADTest(getindex.([urs],cond[:i])...))

                if isplot
                    for f in finalfactor(ctc)
                        plotname = "$(uuid)_E$(e)_U$(u)_$f"
                        plotcondresponse(urs,ctc,u,factor=f,title=plotname,legend=:none)
                        png(joinpath(resultroot,plotname))
                    end
                end

                push!(udf,DataFrame(UUID=uuid,e=e,u=u,ksadp=ksadp))
            end
        end
    end

    cond[:UUID]=uuid
    return vcat(udf...),cond
end