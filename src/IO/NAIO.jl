module NAIO

include("TDTIO.jl")
include("NeuroShareIO.jl")

export matchfile


function matchfile(pattern::Regex;path="")
    fs = readdir(path)
<<<<<<< HEAD
    mi = map(f->ismatch(pattern,f),fs)
=======
    mi = [ismatch(pattern,fs[i]) for i=1:length(fs)]
>>>>>>> 107f77da1a310ca40b52df72ae4d71ec825ec5d1
    mfs = fs[mi]
end

end # module
