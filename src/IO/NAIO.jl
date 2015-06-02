module NAIO

include("TDTIO.jl")
include("NeuroShareIO.jl")

export matchfile


function matchfile(pattern::Regex;path="")
    fs = readdir(path)
    mi = [ismatch(pattern,fs[i]) for i=1:length(fs)]
    mfs = fs[mi]
end

end # module
