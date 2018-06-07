export batchtests,processtest,processori

using ProgressMeter,ePPR

function batchtests(tests::DataFrame,dataroot,resultroot,datatype...;condroot=Dict{String,Any}("rootdir"=>""),delay=20,binwidth=10,isplot=true)
    udf=[];cdf=[]
    p = ProgressMeter.Progress(nrow(tests),1,"Batch Tests ",50)
    for t in eachrow(tests)
        u,c=processtest(prepare(abspath(dataroot,t[:files]),datatype...),resultroot,
        uuid=t[:UUID],condroot=condroot,delay=delay,binwidth=binwidth,isplot=isplot)
        push!(udf,u)
        push!(cdf,c)
        next!(p)
    end
    return vcat(udf...),vcat(cdf...)
end

function processtest(dataset::Dict,resultroot;uuid="",condroot=Dict{String,Any}("rootdir"=>""),delay=20,binwidth=10,isplot=true)
    if haskey(dataset,"ex")
        testid = dataset["ex"]["ID"]
        if testid=="OriGrating"
            processori(dataset,resultroot,uuid=uuid,delay=delay,binwidth=binwidth,isplot=isplot)
        elseif testid=="Image"
            processimage(dataset,condroot,resultroot,uuid=uuid,delay=delay,binwidth=binwidth,isplot=isplot)
        end
    end
end

function processimage(dataset::Dict,condroot::Dict{String,Any},resultroot;uuid="",delay=20,binwidth=10,minpredur=10,
    nscale=2,downsample=2,sigma=1.5,pixelscale=255,isplot=true)
    ex = dataset["ex"];envparam = ex["EnvParam"];preicidur = ex["PreICI"];conddur = ex["CondDur"];suficidur = ex["SufICI"]
    bgcolor=RGBA(getparam(envparam,"BGColor")...)
    imagesetname = replace(getparam(envparam,"ImageSet","ImageQuad"),"Ã—","_")
    imagemasktype = getparam(envparam,"MaskType","ImageQuad")
    imagemaskradius = getparam(envparam,"MaskRadius","ImageQuad")
    imagemasksigma = getparam(envparam,"Sigma","ImageQuad")
    ct,ctc=ctctc(ex)
    fl,fln,fli=flni(ctc)
    cond=condni(ctc)

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
    imagestimulimatrix=cat(2,map(i->vec(gray.(i)),imagestimuli[s])...)'
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

                plotdir = isplot?joinpath(resultroot,"$(uuid)_E$(e)_U$(u)"):nothing
                debug = isplot?ePPRDebugOptions(level=DebugVisual,logdir=plotdir):ePPRDebugOptions()
                hp = ePPRHyperParams(ximagesize...,xindex=xi,ndelay=3,blankcolor=gray(bgimagecolor)*pixelscale)
                hp.nft = [6]
                hp.lambda = 30000
                model,models = epprcv(x,y,hp,debug)

                if isplot && model!=nothing
                    debug(plotmodel(model,hp),log="Model_Final")
                end

                push!(udf,DataFrame(UUID=uuid,e=e,u=u,model=model,models=models,hp=hp))
            end
        end
    end

    cond[:UUID]=uuid
    return vcat(udf...),cond
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