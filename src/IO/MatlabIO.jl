module MatlabIO

using MAT

function read(filename)
   vars = matread(filename);
end

end # module
