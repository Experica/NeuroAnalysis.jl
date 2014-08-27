module MatlabIO

import MAT: matread

export readdataset,matread

function readdataset(filename)
  vars = matread(filename)
  vars = vars["dataset"]
end

end # module
