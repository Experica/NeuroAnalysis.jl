module NAIO

include("TDTIO.jl")
include("NeuroShareIO.jl")

export matchfile


function matchfile(pattern::Regex;path="")
    fs = readdir(path)
    mi = map(f->ismatch(pattern,f),fs)
    mfs = fs[mi]
end

end # module
