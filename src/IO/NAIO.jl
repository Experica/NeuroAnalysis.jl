module NAIO

include("TDTIO.jl")
include("NeuroShareIO.jl")

using DataFrames
export squeezedictdf,matchfiles

function squeezedictdf(d)
  for k in keys(d)
    d[k] = squeeze(d[k],2)
  end
  convert(DataFrame,d)
end

function matchfiles(pattern::Regex;path="")
    fs = readdir(path)
    midx = [ismatch(pattern,fs[i]) for i=1:length(fs)]
    mfs = fs[midx]
end

end # module
