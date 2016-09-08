export matchfile

function matchfile(pattern::Regex;path="")
    fs = readdir(path)
    mi = map(f->ismatch(pattern,f),fs)
    mfs = fs[mi]
end
