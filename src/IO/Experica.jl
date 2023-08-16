function load(source::Val{Experica},filepath;ext=".yaml",prepare=true)
    filedir,filename = splitdir(filepath)
    filename = splitext(filename)[1]
    fullfilename = filename*ext
    fullfilepath =joinpath(filedir,fullfilename)
    if !isfile(fullfilepath)
        @warn "No Experica File:    $fullfilename"
        return nothing
    end

    @info "Reading Experica File:    $fullfilename    ...."
    dataset = Dict{String,Any}()
    dataset["ex"] = YAML.load_file(fullfilepath)
    dataset["filename"] = filename
    dataset["filetime"] = ctime(fullfilepath)
    dataset["datasource"] = Experica
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
    if allunique(sks)
        foreach((sk,k)->d[sk]=d[k],sks,ks)
        foreach(k->delete!(d,k),ks)
    end
    d
end

function prepare!(source::Val{Experica},version::Val{2},ex)
    envparam = ex["EnvParam"]
    tryshortenparamname!(envparam)
    foreach(k->envparam[k]=tryparseparam(envparam[k]), keys(envparam))
    exparam = ex["Param"]
    foreach(k->exparam[k]=tryparseparam(exparam[k]), keys(exparam))
    cond = ex["Cond"]
    foreach(k->cond[k]=tryparseparam.(cond[k]), keys(cond))
    ct = ex["CondTest"]
end