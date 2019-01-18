export circmean,circvar,circr,circaxial,circaxialmean

"Mean direction for circular data"
circmean(α::AbstractVector,w=ones(length(α)))=angle(sum(w.*exp.(im*α)))

"""
Mean Resultant Vector Length for circular data
    α: sample of angles in radians
    w: number of incidences in case of binned angle data
    d: spacing of bin centers for binned data, if supplied
       correction factor is used to correct for bias
"""
function circr(α::AbstractVector,w=ones(length(α)),d=nothing)
    r = abs(sum(w.*exp.(im*α)))/sum(w)
    if d!=nothing
        r *= d/2/sin(d/2)
    end
    r
end

"Circular Variance"
circvar(α::AbstractVector,w=ones(length(α)),d=nothing) = 1-circr(α,w,d)

"""
Transforms p-axial data to a common scale
  α: sample of angles in radians
  p: number of modes
"""
circaxial(α,p=2)=(p*α)%(2pi)

"Mean direction for circular data with axial"
circaxialmean(α::AbstractVector,p=2)=angle(sum(exp.(im*α*p)))/p
