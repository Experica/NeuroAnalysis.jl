export batchtests,processtest,processori

function batchtests(tests::DataFrame,dataroot,datatype...;delay=20,binwidth=10)
    udf=DataFrame();cdf=DataFrame()
    for t in eachrow(tests)
        u,c=processtest(prepare(abspath(dataroot,t[:files]),datatype...),uuid=t[:UUID],delay=delay,binwidth=binwidth)
        udf=vcat(udf,u)
        cdf=vcat(cdf,c)
    end
    return udf,cdf
end

function processtest(dataset::Dict;uuid="",delay=20,binwidth=10)
    if haskey(dataset,"ex")
        testid = dataset["ex"]["ID"]
        if testid=="OriGrating"
            processori(dataset,uuid=uuid,delay=delay,binwidth=binwidth)
        end
    end
end

function processori(dataset::Dict;uuid="",delay=20,binwidth=10)
    ex = dataset["ex"];envparam = ex["EnvParam"];preicidur = ex["PreICI"];conddur = ex["CondDur"];suficidur = ex["SufICI"]
    ct,ctc=ctctc(ex)
    fl,fln,fli=flni(ctc)
    cond=condni(ctc)
    
    udf = DataFrame()
    if haskey(dataset,"spike")
        spike = dataset["spike"];spikeeid = spike["electrodeid"];spikeuuid = spike["uuid"];spikeuid = spike["unitid"];spiketime = spike["time"]
        for e in spikeeid
            ei = findfirst(spikeeid.==e);est = spiketime[ei];esu = spikeuid[ei];euuid = spikeuuid[ei]
            for u in euuid
                urs = subrvr(est[esu.==u],ct[:CondOn]+delay,ct[:CondOff]+delay)
                umse = condresponse(urs,cond[:i])
                stats = statsori(umse[:m],cond[finalfactor(cond)][1])

                udf = vcat(udf,DataFrame(UUID=uuid,e=e,u=u,r=umse,cv=stats[:cv]))
            end
        end
    end

    cond[:UUID]=uuid
    return udf,cond
end