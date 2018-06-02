export batchtests,processtest,processori

function batchtests(tests::DataFrame,dataroot,resultroot,datatype...;delay=20,binwidth=10,isplot=true)
    udf=[];cdf=[]
    for t in eachrow(tests)
        u,c=processtest(prepare(abspath(dataroot,t[:files]),datatype...),resultroot,
        uuid=t[:UUID],delay=delay,binwidth=binwidth,isplot=isplot)
        push!(udf,u)
        push!(cdf,c)
    end
    return vcat(udf...),vcat(cdf...)
end

function processtest(dataset::Dict,resultroot;uuid="",delay=20,binwidth=10,isplot=true)
    if haskey(dataset,"ex")
        testid = dataset["ex"]["ID"]
        if testid=="OriGrating"
            processori(dataset,resultroot,uuid=uuid,delay=delay,binwidth=binwidth,isplot=isplot)
        end
    end
end

function processori(dataset::Dict,resultroot;uuid="",delay=20,binwidth=10,minpredur=100,isplot=true)
    ex = dataset["ex"];envparam = ex["EnvParam"];preicidur = ex["PreICI"];conddur = ex["CondDur"];suficidur = ex["SufICI"]
    ct,ctc=ctctc(ex)
    fl,fln,fli=flni(ctc)
    cond=condni(ctc)

    ff = finalfactor(cond)[1]
    oris = ctc[ff]
    predur = max(preicidur,minpredur)
    resultroot=abspath(resultroot,"LZ_Ori")
    !isdir(resultroot) && mkdir(resultroot)
    
    udf = []
    if haskey(dataset,"spike")
        spike = dataset["spike"];spikeeid = spike["electrodeid"];spikeuuid = spike["uuid"];spikeuid = spike["unitid"];spiketime = spike["time"]
        for e in spikeeid
            ei = findfirst(spikeeid.==e);est = spiketime[ei];esu = spikeuid[ei];euuid = spikeuuid[ei]
            for u in euuid
                preurs = subrvr(est[esu.==u],ct[:CondOn]-predur,ct[:CondOn])
                urs = subrvr(est[esu.==u],ct[:CondOn]+delay,ct[:CondOff]+delay)
                !isresponsive(preurs,urs,cond[:i]) && continue

                if isplot
                    plotname = "$(uuid)_E$(e)_U$(u)"
                    plotcondresponse(urs,cond,u,title=plotname)
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

isresponsive(baseline,response;alpha=0.05) = pvalue(SignedRankTest(baseline,response)) < alpha
isresponsive(baseline,response,is;alpha=0.05) = any(map(i->isresponsive(baseline[i],response[i],alpha=alpha),is))