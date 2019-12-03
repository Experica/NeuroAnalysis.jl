export csd

"""
1D Current Source Density for a set of voltage traces from a linear array of electrodes.

CSD: finite difference method of second spatial derivatives of potentials, I = - ∇σ⋅∇ϕ,
      σ is conductivities, ϕ is potentials, on the continuous assumption when spacing is around 20-100 um.
     (Nicholson & Freeman, 1975, J Neurophysiol, 38(2): 356-68)

iCSD: solve I -> ϕ forward problem, then get the inverse transformation matrix.
      (Petterson et al., 2006, J Neurosci Methods, 154(1-2):116-33)

      iCSDdelta: δ method where infinitely thin current source discs of radius r.
      iCSDspline: cubic spline method where current source between electrode be modeled by smooth spline

kCSD: kernal CSD
      ()

data: ch x sample matrix(volts)
h: channel spacing(micrometers)
c: conductivity of the extracellular medium in siemans/meter
r: radius(micrometers) of the model current source/sink in reverse method

return: CSD(amps/meters^3)
"""
function csd(data;filter=nothing,method=:iCSDdelta,h=20,c=0.3,r=500)
  nd=ndims(data)
  if nd==3
    nr,nc,n = size(data)
    y = Array{Float64}(undef,nr,nc,n)
    for i in 1:n
      y[:,:,i] = csd(data[:,:,i],filter=filter,method=method,h=h,c=c,r=r)
    end
  else
    fdata = isnothing(filter) ? data : imfilter(data,filter)
    nr,nc = size(fdata)
    ħ = h*1e-6
    ɽ = r*1e-6
    if method==:iCSDdelta
      f = [ħ^2/(2c) * ( sqrt((j-i)^2 + (ɽ/ħ)^2) - abs(j-i) ) for j=1:nr,i=1:nr]
      y = f\fdata
    elseif method==:iCSDspline # have wired problems, doesn't replicate the original code(CSDPlotter) results
      ep = collect(1:nr)*h
      y,zs = csd_cubicspline(ep,data,f_cubicspline(ep,r,c)...)
    elseif method==:kCSD

    else # traditional CSD
      y=zeros(Float64,nr,nc)
      for i in 2:nr-1
        y[i,:] = -c*(fdata[i-1,:] + fdata[i+1,:] - 2*fdata[i,:]) / ħ^2
      end
    end
  end
  return y
end

"""
cubic spline inverse CSD

contact_positions: contact positions
data: measured potentials
Fcs: cubic spline transformation matrix corresponding to the given contact positions
E0,E1,E2,E3: matrixes containing the "recursive rules"
"""
function csd_cubicspline(contact_positions,data,Fcs,E0,E1,E2,E3,num_out_zs=size(data,1))
  N,num_of_timesteps = size(data)
  cs_data = zeros(N+2,num_of_timesteps)
  cs_data[2:N+1,:] = data # ϕ_1 = ϕ_N+2 = 0
  h = mean(diff(contact_positions))
  el_pos_with_ends = vcat(0,contact_positions,contact_positions[end]+h)

  CSD_coeff = Fcs\cs_data
  # The cubic spline polynomial coeffescients
  A0 = E0*CSD_coeff
  A1 = E1*CSD_coeff
  A2 = E2*CSD_coeff
  A3 = E3*CSD_coeff

  out_zs = el_pos_with_ends[1]:(el_pos_with_ends[N+2]-el_pos_with_ends[1])/(num_out_zs-1):el_pos_with_ends[N+2]
  CSD = zeros(length(out_zs),num_of_timesteps)

  i = 1
  for j=1:length(out_zs)
    if out_zs[j]>el_pos_with_ends[i+1]
      i+=1
    end
    CSD[j,:] = A0[i,:] + A1[i,:]*(out_zs[j]-el_pos_with_ends[i]) +
    A2[i,:]*(out_zs[j]-el_pos_with_ends[i])^2 + A3[i,:]*(out_zs[j]-el_pos_with_ends[i])^3
  end
  return CSD,out_zs
end

"""
Computes the E0, E1, E2 and E3 matrixes used in the cubic spline iCSD method.
These matrixes contains the recursive formulas for finding the F matrix
"""
function e_matrixes(el_pos)
  N = length(el_pos)
  h = mean(diff(el_pos))
  z_js = vcat(0,el_pos,el_pos[end]+h) # electrode positions with two imaginary

  C_vec = 1/diff(z_js) # length: N+1
  # Define transformation matrixes
  C_jm1 = zeros(N+2,N+2)
  C_j0 = zeros(N+2,N+2)
  C_mat3 = zeros(N+1,N+1)

  for i=1:N+1
    for j=1:N+1
      if i == j
        C_jm1[i+1,j+1] = C_vec[i]
        C_j0[i,j] = C_jm1[i+1,j+1]
        C_mat3[i,j] = C_vec[i]
      end
    end
  end
  C_jm1[N+2,N+2] = 0
  C_j0[1,1] = 0

  C_jall = copy(C_j0)
  C_jall[1,1] = 1
  C_jall[N+2,N+2] = 1

  Tjp1 = zeros(N+2,N+2) # converting an element k_j to k_j+1
  Tjm1 = zeros(N+2,N+2) # converting an element k_j to k_j-1
  Tj0  = Matrix{Float64}(I,N+2,N+2)
  Tj0[1,1] = 0
  Tj0[N+2,N+2] = 0

  # C to K
  for i=2:N+2
    for j=1:N+2
      if i==j-1
        Tjp1[i,j] = 1
      end
      if i==j+1
        Tjm1[i,j] = 1
      end
    end
  end

  # C to K transformation matrix
  K = ( C_jm1*Tjm1 + 2*C_jm1*Tj0 + 2*C_jall + C_j0*Tjp1 ) \ (3*( C_jm1^2*Tj0 - C_jm1^2*Tjm1 + C_j0^2*Tjp1 - C_j0^2*Tj0 ))

  # Define matrixes for C to A transformation
  Tja  = zeros(N+1,N+2)      # identity matrix except that it cuts off last element
  Tjp1a  = zeros(N+1,N+2)    # converting k_j to k_j+1 and cutting off last element

  # C to A
  for i=1:N+1
    for j=1:N+2
      if i==j-1
        Tjp1a[i,j] = 1
      end
      if i==j
        Tja[i,j] = 1
      end
    end
  end

  # Define spline coeffiscients
  E0  = Tja
  E1  = Tja*K
  E2  = 3*C_mat3^2*(Tjp1a-Tja) - C_mat3*(Tjp1a+2*Tja)*K
  E3  = 2*C_mat3^3*(Tja-Tjp1a) + C_mat3^2*(Tjp1a+Tja)*K
  return E0,E1,E2,E3
end

"F and E matrix of the cubic spline inverse CSD"
function f_cubicspline(el_pos,r,cond=0.3,cond_top=cond)
  N = length(el_pos)
  h = mean(diff(el_pos))
  z_js = vcat(0,el_pos,el_pos[end]+h) # electrode positions with two imaginary

  E0,E1,E2,E3 = e_matrixes(el_pos)
  # Potential functions
  f0(zeta,zj,diam,cond) = 1/(2cond)*(sqrt((diam/2)^2 + (zj-zeta)^2) - abs(zj-zeta))
  f1(zeta,zj,zi,diam,cond) = (zeta-zi)*f0(zeta,zj,diam,cond)
  f2(zeta,zj,zi,diam,cond) = (zeta-zi)^2*f0(zeta,zj,diam,cond)
  f3(zeta,zj,zi,diam,cond) = (zeta-zi)^3*f0(zeta,zj,diam,cond)
  # Define integration matrixes
  F0  = zeros(N,N+1)
  F1  = zeros(N,N+1)
  F2  = zeros(N,N+1)
  F3  = zeros(N,N+1)

  for j = 1:N
    for i = 1:N+1
      F0[j,i] = hquadrature(x->f0(x,z_js[j+1],r,cond),z_js[i],z_js[i+1])[1]
      F1[j,i] = hquadrature(x->f1(x,z_js[j+1],z_js[i],r,cond),z_js[i],z_js[i+1])[1]
      F2[j,i] = hquadrature(x->f2(x,z_js[j+1],z_js[i],r,cond),z_js[i],z_js[i+1])[1]
      F3[j,i] = hquadrature(x->f3(x,z_js[j+1],z_js[i],r,cond),z_js[i],z_js[i+1])[1]
      if cond != cond_top # image technique if conductivity not constant
        F0[j,i] = F0[j,i] + (cond-cond_top)/(cond+cond_top)*hquadrature(x->f0(x,-z_js[j+1],r,cond),z_js[i],z_js[i+1])[1]
        F1[j,i] = F1[j,i] + (cond-cond_top)/(cond+cond_top)*hquadrature(x->f1(x,-z_js[j+1],z_js[i],r,cond),z_js[i],z_js[i+1])[1]
        F2[j,i] = F2[j,i] + (cond-cond_top)/(cond+cond_top)*hquadrature(x->f2(x,-z_js[j+1],z_js[i],r,cond),z_js[i],z_js[i+1])[1]
        F3[j,i] = F3[j,i] + (cond-cond_top)/(cond+cond_top)*hquadrature(x->f3(x,-z_js[j+1],z_js[i],r,cond),z_js[i],z_js[i+1])[1]
      end
    end
  end

  temp_F = F0*E0+F1*E1+F2*E2+F3*E3  # F matrix(N x N+2)
  # Convert to (N+2 x N+2) matrixes by applying the boundary I_0 = I_N+1 = 0.
  F = zeros(N+2,N+2)
  F[2:N+1,:] = temp_F
  F[1,1] = 1         # implies I_1 = ϕ_1
  F[N+2,N+2] = 1     # implies I_N+2 = ϕ_N+2
  return F,E0,E1,E2,E3
end
