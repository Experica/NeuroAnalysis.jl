function load(source::Val{Experica},filepath;ext=".yaml",prepare=true)
    filedir,filename = splitdir(filepath)
    filename = splitext(filename)[1]
    fullfilename = filename*ext
    fullfilepath =joinpath(filedir,fullfilename)
    if !isfile(fullfilepath)
        @error "No Experica File:    $fullfilepath"
        return nothing
    end

    @info "Reading Experica File:    $fullfilename    ...."
    dataset = Dict{String,Any}()
    dataset["ex"] = YAML.load_file(fullfilepath)
    dataset["filename"] = filename
    dataset["filetime"] = ctime(fullfilepath)
    dataset["sourceformat"] = Experica
    @info "Reading Experica File:    $fullfilename    Done"

    if prepare
        v = dataset["ex"]["Version"]
        @info "Preparing Experica Data(v$v):    $fullfilename    ...."
        prepare!(Val(Experica),Val(v),dataset["ex"])
        @info "Preparing Experica Data(v$v):    $fullfilename    Done"
    end
    return dataset
end

tryparseparam(v)=v
function tryparseparam(v::AbstractString)
    if contains(v,' ')
        try
            return parse.(Float64,split(v,' '))
        catch
        end
    end
    v
end
function tryshortenparamname!(d::Dict)
    ks = filter(k->contains(k,'@'),collect(keys(d)))
    sks = map(k->split(k,'@')[1],ks)
    cm = countmap(sks)
    vsks = [sk for (sk,c) in cm if c==1]
    vks = ks[indexin(vsks,sks)]
    foreach((sk,k)->d[sk]=d[k],vsks,vks)
    foreach(k->delete!(d,k),vks)
    d
end
function trygetparam(d::Dict,name,fromobject="")
    ks = filter(k->startswith(k,name) && endswith(k,fromobject),collect(keys(d)))
    if length(ks)==0
        @warn "No Matching Param Found: $name*$fromobject"
        return nothing
    elseif length(ks)>1
        @warn "Multiple Matching Param Found, Need To Be More Specific: $ks"
        return nothing
    else
        return d[ks[1]]
    end
end
function trygetscreensize(ps::Dict)
    d = trygetparam(ps,"ScreenToEye")
    sh = trygetparam(ps,"ScreenHeight")
    sr = trygetparam(ps,"ScreenAspect")
    if isnothing(d) || isnothing(sh) || isnothing(sr)
        return nothing
    else
        h = atand(sh/2,d)*2
        w = h * sr
        return (w,h)
    end
end

function prepare!(source::Val{Experica},version::Val{2},ex)
    envparam = ex["EnvParam"]
    tryshortenparamname!(envparam)
    foreach(k->envparam[k]=tryparseparam(envparam[k]), keys(envparam))
    exparam = ex["Param"]
    foreach(k->exparam[k]=tryparseparam(exparam[k]), keys(exparam))
    cond = ex["Cond"]
    foreach(k->cond[k]=tryparseparam.(cond[k]), keys(cond))
end

function prepare!(source::Val{Experica},version::Val{3},ex)
    envparam = ex["EnvParam"]
    tryshortenparamname!(envparam)
    foreach(k->envparam[k]=tryparseparam(envparam[k]), keys(envparam))
    extparam = ex["ExtendParam"]
    foreach(k->extparam[k]=tryparseparam(extparam[k]), keys(extparam))
    cond = ex["Cond"]
    foreach(k->cond[k]=tryparseparam.(cond[k]), keys(cond))
    cond = cond |> DataFrame
    
    ct = ex["CondTest"]
    nct = mapreduce(length, max, values(ct))
    foreach(k->begin
        v = eltype(ct[k])==Any ? ct[k] : Any[ct[k]...]
        l = length(v)

        for i in eachindex(v)
            if isnothing(v[i]) || isempty(v[i])
                v[i] = missing
            end
        end
        if endswith(k,"Index")
            @warn "Change CondTest[$k]: 0-based to 1-based"
            v.+=1            
        end

        if l < nct
            @warn "Pad CondTest[$k] Length: $l to $nct"
            append!( v, fill(missing, nct - l) )
        end

        ct[k] = v
    end, keys(ct))
    ct = ct |> DataFrame

    ex["Cond"] = cond
    ex["CondTest"] = ct

    if "CondIndex" in names(ct)
        ctc = transform(ct,:CondIndex => ByRow(i-> NamedTuple([f=>ismissing(i) ? missing : cond[i,f] for f in propertynames(cond)])) => AsTable)
        ex["CondTestCond"] = ctc
    end    

end