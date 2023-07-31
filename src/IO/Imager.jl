"Read list of files encoding raw mono 8bit image"
function readrawim_Mono8(file::Array,width,height;n=length(file),T=Float64)
    n == 0 && error("No files to read.")
    imgs = Array{T}(undef,height,width,n)
    @inbounds for i in eachindex(file)
        imgs[:,:,i] = readrawim_Mono8(file[i],width,height)
    end
    imgs
end
function readrawim_Mono12Packed(file::Array,width,height;n=length(file),T=Float64)
    n == 0 && error("No files to read.")
    imgs = Array{T}(undef,height,width,n)
    raw = Vector{UInt8}(undef,Int(width*height*1.5))
    npack = Int(length(raw)/3)
    img = Vector{UInt8}(undef,4npack)
    @inbounds for i in eachindex(file)
        imgs[:,:,i] = readrawim_Mono12Packed(file[i],width,height;raw,npack,img)
    end
    imgs
end

readrawim_Mono8(file::String,width,height) = permutedims(Mmap.mmap(file,Matrix{UInt8},(width,height)))
function readrawim_Mono12Packed(file::String,width,height;raw = Vector{UInt8}(undef,Int(width*height*1.5)),
                                npack = Int(length(raw)/3), img = Vector{UInt8}(undef,4npack))
    # Unpack Mono12Packed 3 bytes to 2 UInt16
    # The two 12bits pixels | A | B | are packed in | A11 A10 A9 A8 A7 A6 A5 A4 | B3 B2 B1 B0 A3 A2 A1 A0 | B11 B10 B9 B8 B7 B6 B5 B4 | 3 bytes
    read!(file,raw)
    @turbo unroll=1 thread=true for i in 0:npack-1
        img[1+4i] = raw[2+3i] << 4 # second byte shift 4 bits left | A3 A2 A1 A0 0 0 0 0 |
        img[2+4i] = raw[1+3i]
        img[3+4i] = raw[2+3i]
        img[4+4i] = raw[3+3i]
    end
    img = reinterpret(UInt16,img)
    @turbo unroll=1 thread=true for i in 1:2npack
        img[i] >>>= 4 # shift UInt16s 4 bits right
    end
    permutedims(reshape(img,width,height))
end

function load(source::Val{Imager},filepath;ext=".meta")
    filedir,filename = splitdir(filepath)
    filename = splitext(filename)[1]
    fullfilename = filename*ext
    fullfilepath =joinpath(filedir,fullfilename)
    if !isfile(fullfilepath)
        @warn "No Imager File:    $fullfilename"
        return nothing
    end

    @info "Reading Imager File:    $fullfilename    ...."
    dataset = Dict{String,Any}()
    dataset["meta"] = YAML.load_file(fullfilepath)
    dataset["filename"] = filename
    dataset["filetime"] = ctime(fullfilepath)
    dataset["datasource"] = Imager
    @info "Reading Imager File:    $fullfilename    Done"
    return dataset
end
