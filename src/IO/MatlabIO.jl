module MatlabIO

using MAT

export readdataset

function readdataset(filename)
  vars = matread(filename)
  vars = vars["dataset"]
end

end # module
