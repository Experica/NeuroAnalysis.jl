export gaincorrectim,refchmaskim

"""
Convert raw imec data to gain-corrected voltages. The saved-channel range `1:nSavedChans`
will convert to the original channel id, then corresponding gain will be applied.

The full conversion with gain is:
   dataVolts = dataInt * fI2V / gain.
Note that each channel may have its own gain.

y: the epoch of all imec channels, excluding the last Sync channel.
"""
function gaincorrectim(y,meta)
    ochans = meta["chans"]
    apgain = meta["apgain"]
    lfgain = meta["lfgain"]
    nap=length(apgain)
    fi2v = meta["fi2v"]

    cy=Array{Float64}(undef,size(y))
    for i in 1:size(y,1)
        j=ochans[i]
        if j <= nap
            f = fi2v/apgain[j]
        elseif j <= 2nap
            f = fi2v/lfgain[j-nap]
        else
            f = 1
        end
        cy[i,:]=y[i,:]*f
    end
    return cy
end
"Logical mask for IMEC reference channels in Probe shape"
function refchmaskim(refs,nrow,ncol)
    mask = falses(ncol,nrow) # channel counting from cols(left -> right), then rows(tip -> tail)
    mask[refs].=true # channel linear indexing to probe shape
    return mask'
end
refchmaskim(meta) = refchmaskim(meta["refch"],meta["nrow"],meta["ncol"])
