export gaincorrectim,refchmaskim

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
"Logical mask for `IMEC` reference channels in probe shape"
function refchmaskim(nch,refs,nrow,ncol)
    mask = falses(nch)
    mask[refs].=true
    mask = reshape(mask,ncol,nrow) # Neuropixel Probe channel counting is from cols(left -> right), then rows(tip -> tail)
    return mask' # rows are reverted and will be reverted back when XY plotting
end
refchmaskim(meta) = refchmaskim(meta["acqApLfSy"][1],meta["refch"],meta["nrow"],meta["ncol"])
