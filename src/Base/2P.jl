

export sbxsubrm, dPrime

function sbxsubrm(rm,epochs,rois;fun=nothing)
    nepoch = size(epochs,1)   # number of trial
    minepochlength = floor(Int,minimum(diff(epochs,dims=2)))
    ys=Array{Float64}(undef,length(rois),minepochlength,nepoch)
    for i in 1:nepoch
        y = rm[rois,range(max(1,epochs[i,1]),length=minepochlength)]
        if !isnothing(fun)
            y=fun(y)
        end
        ys[:,:,i] = y
    end
    return nepoch==1 ? dropdims(ys,dims=3) : ys
end

function dPrime(bi)
        return x->begin
        bl = mean(x[:,bi],dims=2)
        x = x./bl .-1
        x
        end
    end

sbxcondin(ctc::Dict)=sbxcondin(DataFrame(ctc))
function sbxcondin(ctc::DataFrame)
    t = [ctc DataFrame(i=1:nrow(ctc))]
    t = by(t, names(ctc),g->DataFrame(n=nrow(g), i=[g[:,:i]]))
    sort!(t.n); sort!(t[1:end-1]); return t
end
