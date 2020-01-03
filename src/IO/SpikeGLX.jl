export gaincorrectim,chmaskim,refchmaskim,badchmaskim

"""
Convert raw imec data to gain-corrected voltages. The saved-channel id subset in data
will be used to get corresponding gain.

The full conversion with gain is:
   dataVolts = dataInt * fI2V / gain.
Note that each channel may have its own gain.

y: the saved epoch of all imec channels, excluding the last Sync channel.
"""
function gaincorrectim(y,meta)
    chids = meta["savedchans"]
    fi2v = meta["fi2v"]
    apgain = meta["roapgain"]
    lfgain = meta["rolfgain"]
    nch=length(apgain)

    cy=Array{Float64}(undef,size(y))
    for i in 1:size(y,1)
        chid=chids[i]
        if chid <= nch
            f = fi2v/apgain[chid]
        elseif chid <= 2nch
            f = fi2v/lfgain[chid-nch]
        end
        cy[i,:]=y[i,:]*f
    end
    return cy
end
"Logical mask for `IMEC` channels in probe shape"
function chmaskim(nch,chs,nrow,ncol)
    mask = falses(nch)
    mask[chs].=true
    mask = reshape(mask,ncol,nrow) # Neuropixel Probe channel counting is from cols(left -> right), then rows(tip -> tail)
    return mask' # rows are reverted and will be reverted back when XY plotting
end
"Logical mask for `IMEC` reference channels in probe shape"
function refchmaskim(dataset;type="lf")
    nch = dataset["ap"]["meta"]["acqApLfSy"][1]
    rorefch = dataset["ap"]["meta"]["rorefch"][1]
    if rorefch==0
        refch = dataset[type]["meta"]["refch"]
    else
        refch = [rorefch]
    end
    chmaskim(nch,refch,dataset["ap"]["meta"]["nrow"],dataset["ap"]["meta"]["ncol"])
end
"Logical mask for `IMEC` bad and reference channels in probe shape"
function badchmaskim(dataset;type="lf",badch::Vector{Int}=Int[])
    nch = dataset["ap"]["meta"]["acqApLfSy"][1]
    rorefch = dataset["ap"]["meta"]["rorefch"][1]
    if rorefch==0
        refch = dataset[type]["meta"]["refch"]
    else
        refch = [rorefch]
    end
    if haskey(dataset[type]["meta"],"badch")
        badch = union(badch,dataset[type]["meta"]["badch"])
    end
    chmaskim(nch,union(refch,badch),dataset["ap"]["meta"]["nrow"],dataset["ap"]["meta"]["ncol"])
end
