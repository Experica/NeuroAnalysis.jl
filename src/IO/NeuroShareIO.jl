module Analysis.IO
# Read Neural Data via NeuroShare API
using Base NeuroShare

function nssetlib(libfile)
	ns_SetLibrary(libfile)
	ns_GetLibraryInfo()
end

function nsopenfile(file)
	nssetlib(libfile)
	ns_OpenFile(file)
	ns_GetFileInfo()

end

function nsreadall(file)

end

function nsreadanalog(file)

end

function nsreadevent(file)
end

function nsreadepoch(file)
end

function nsreadsegment(file)

end

function nsreadcell(file)
end

function nsreadcellassemble(file)
end

end # module