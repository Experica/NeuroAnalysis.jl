"""
1D Current Source Density for voltage traces from a equidistant linear array of electrodes.

1. data: ch x sample or ch x sample x epoch(volts)

2. method:
  - CSD: finite difference approximation of second spatial derivatives of potentials, C = - ∇⋅σ∇ϕ,
    σ is conductivity tensor, ϕ is potential.
      (Nicholson & Freeman, 1975, J Neurophysiol, 38(2): 356-68)

  - iCSD: solve C -> ϕ forward model, then get Ĉ by the inverse transformation matrix,
    here C is modeled by δ-source, i.e. infinitely thin current source disc of radius r.
      (Petterson et al., 2006, J Neurosci Methods, 154(1-2):116-33)

  - kCSD: kernal CSD
      ()

- h: channel spacing(micrometer)
- c: conductivity of extracellular medium(siemans/meter)
- r: radius(micrometer) of δ-source

- return: CSD(amps/meter^3)
"""
function csd(data,method=:iCSD;h=20,c=0.3,r=500)
  nd=ndims(data)
  if nd==3
    nr,nc,n = size(data)
    y = Array{Float64}(undef,nr,nc,n)
    for i in 1:n
      @views y[:,:,i] = csd(data[:,:,i],method;h,c,r)
    end
  else
    nr,nc = size(data)
    𝒉 = h*1e-6
    𝒓 = r*1e-6
    if method==:iCSD
      F = [𝒉^2/(2c) * ( sqrt((j-i)^2 + (𝒓/𝒉)^2) - abs(j-i) ) for j=1:nr,i=1:nr]
      y = F\data
    elseif method==:kCSD

    else
      y = zeros(nr,nc)
      for i in 2:nr-1
        @views y[i,:] = -c * (data[i-1,:] + data[i+1,:] - 2*data[i,:]) / 𝒉^2
      end
    end
  end
  return y
end
