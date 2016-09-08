module NAIO
using MAT,DataFrames

export readmat,prepare,prepare!,prepare_ripple!,prepare_vlab!,statetime,matchfile

"Read `Matlab` exported data"
function readmat(f::String;v="dataset")
  d = matread(f)[v]
end

prepare(f::String)=prepare!(readmat(f))
function prepare!(d::Dict)
    if(haskey(d,"sourceformat"))
        sf = d["sourceformat"]
        if(sf=="Ripple")
            d=prepare_ripple!(d)
        end
    end
    return d
end
function prepare_ripple!(d::Dict)
    if(haskey(d,"spike"))
        d["spike"]["electrodeid"]=collect(Int,d["spike"]["electrodeid"])
        d["spike"]["time"]=map(i->collect(Float64,i),collect(d["spike"]["time"]))
        d["spike"]["unitid"]=map(i->collect(Int,i),collect(d["spike"]["unitid"]))
    end
    if(haskey(d,"digital"))
        dc = i-> begin
            s = split(i)
            c = s[1]=="SMA"?parse(Int,s[2]):s
        end
        d["digital"]["channel"]=map(i->dc(i[1,1]),collect(d["digital"]["channel"]))
        d["digital"]["time"]=map(i->collect(Float64,i),collect(d["digital"]["time"]))
        d["digital"]["data"]=map(i->collect(Int,i),collect(d["digital"]["data"]))
    end
    if(haskey(d,"analog1k"))
        d["analog1k"]["electrodeid"]=collect(Int,d["analog1k"]["electrodeid"])
        d["analog1k"]["time"]=collect(Float64,d["analog1k"]["time"])
    end
    if(haskey(d,"ex"))
        d["ex"]=prepare_vlab!(d["ex"])
    end
    return d
end
function prepare_vlab!(d::Dict)
    if(haskey(d["CondTest"],"CONDSTATE"))
        d["CondTest"]["CONDSTATE"]= map(i->collect(i),collect(d["CondTest"]["CONDSTATE"]))
    end
    return d
end
function statetime(ct::Dict;statetype::AbstractString="CONDSTATE",state::AbstractString="COND")
    if(haskey(ct,statetype))
        filter!(l->!isempty(l),map(i->begin
        t = filter!(k->!isempty(k),map(j->haskey(j,state)?j[state]:[],i))
            length(t)==1?t[1]:t
            end,ct[statetype]))
    else
        []
    end
end

"convert `Matlab` struct of array to `DataFrame`"
matdictarray2df(d) = DataFrame(Any[squeeze(v,2) for v in values(d)],[Symbol(k) for k in keys(d)])

"Regular expression to match `VLab` data file names"
function vlabregex(testtype;subject="[0-9]",maxch=1,cell="[a-z]",maxrepeat=3,format="mat")
  mr = lpad(maxrepeat,2,0)
  Regex("^$subject+[A-Za-z]+[1-$maxch]$cell$testtype[0-$(mr[1])][0-$(mr[2])]\\.$format")
end

"Get matched file names in path"
function getdatafile(testtype;subject="[0-9]",path="./data")
  matchfile(vlabregex(testtype,subject=subject),path=path)
end

function matchfile(pattern::Regex;path="")
    fs = readdir(path)
    mi = map(f->ismatch(pattern,f),fs)
    mfs = fs[mi]
end
