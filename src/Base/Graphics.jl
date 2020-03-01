import Base: convert,getindex,setindex!,+,-,*,/,==,length,norm

# Vector3
struct Vec3{T}
  x::T
  y::T
  z::T
end
Vec3(x,y)=Vec3(x,y,y)
Vec3(x)=Vec3(x,x,x)
Vec3()=Vec3(0.0)
function getindex(V::Vec3,i::Int)
  if i==1
    X=V.x
  elseif i==2
    X=V.y
  elseif i==3
    X=V.z
  else
    error("Out of Bound of Vector3.")
  end
end
function setindex!(V::Vec3,X,i::Int)
  if i==1
    V.x=X
  elseif i==2
    V.y=X
  elseif i==3
    V.z=X
  else
    error("Out of Bound of Vector3.")
  end
end
function convert(::Type{Vec3},a::AbstractVector)
  l=length(a)
  if l==0
    error("Empty Vector.")
  elseif l==1
    v = Vec3(a[1])
  elseif l==2
    v = Vec3(a[1],a[2])
  else
    v = Vec3(a[1],a[2],a[3])
  end
end
function convert(::Type{Vector},v::Vec3)
  a = [v.x,v.y,v.z]
end
function convert{T<:Vec3}(::Type{Matrix},vs::AbstractVector{T})
  l = length(vs)
  if l==0
    error("Empty Vector.")
  else
    a = Array(typeof(vs[1].x),3,l)
    for i in 1:l
      a[1,i]=vs[i].x
      a[2,i]=vs[i].y
      a[3,i]=vs[i].z
    end
  end
  return a
end
function convert{T<:Vec3}(::Type{Vector{T}},a::Matrix)
  s1,s2 = size(a)
  if s1 != 4
    error("Size of Dim 1 of Matrix Doesn't Match Vec4")
  elseif s2 == 0
    return []
  else
    vs = Array(T,s2)
    for i in 1:s2
      vs[i] = Vec3(a[1,i],a[2,i],a[3,i])
    end
    return vs
  end
end
+(a::Vec3,b::Vec3)=Vec3(a.x+b.x,a.y+b.y,a.z+b.z)
-(a::Vec3,b::Vec3)=Vec3(a.x-b.x,a.y-b.y,a.z-b.z)
*(a::Vec3,b::Real)=Vec3(a.x*b,a.y*b,a.z*b)
*(b::Real,a::Vec3)=*(a,b)
/(a::Vec3,b::Real)=Vec3(a.x/b,a.y/b,a.z/b)
/(b::Real,a::Vec3)=/(a,b)
==(a::Vec3,b::Vec3)=(a.x==b.x) && (a.y==b.y) && (a.z==b.z)

# Vector4
struct Vec4{T}
  x::T
  y::T
  z::T
  w::T
end
Vec4(x,y,z)=Vec4(x,y,z,z)
Vec4(x,y)=Vec4(x,y,y,y)
Vec4(x)=Vec4(x,x,x,x)
Vec4()=Vec4(0.0,0.0,0.0,1.0)
function getindex(V::Vec4,i::Int)
  if i==1
    X=V.x
  elseif i==2
    X=V.y
  elseif i==3
    X=V.z
  elseif i==4
    X=V.w
  else
    error("Out of Bound of Vector4.")
  end
end
function setindex!(V::Vec4,X,i::Int)
  if i==1
    V.x=X
  elseif i==2
    V.y=X
  elseif i==3
    V.z=X
  elseif i==4
    V.w=X
  else
    error("Out of Bound of Vector4.")
  end
end
function convert(::Type{Vec4},a::AbstractVector)
  l=length(a)
  if l==0
    error("Empty Vector.")
  elseif l==1
    v = Vec4(a[1])
  elseif l==2
    v = Vec4(a[1],a[2])
  elseif l==3
    v = Vec4(a[1],a[2],a[3])
  else
    v = Vec4(a[1],a[2],a[3],a[4])
  end
end
function convert(::Type{Vector},v::Vec4)
  a = [v.x,v.y,v.z,v.w]
end
function convert{T<:Vec4}(::Type{Matrix},vs::AbstractVector{T})
  l = length(vs)
  if l==0
    error("Empty Vector.")
  else
    a = Array(typeof(vs[1].x),4,l)
    for i in 1:l
      a[1,i]=vs[i].x
      a[2,i]=vs[i].y
      a[3,i]=vs[i].z
      a[4,i]=vs[i].w
    end
  end
  return a
end
function convert{T<:Vec4}(::Type{Vector{T}},a::Matrix)
  s1,s2 = size(a)
  if s1 != 4
    error("Size of Dim 1 of Matrix Doesn't Match Vec4")
  elseif s2 == 0
    return []
  else
    vs = Array(T,s2)
    for i in 1:s2
      vs[i] = Vec4(a[1,i],a[2,i],a[3,i],a[4,i])
    end
    return vs
  end
end
convert(::Type{Vec3},v::Vec4)=Vec3(v.x,v.y,v.z)
convert(::Type{Vec4},v::Vec3)=Vec4(v.x,v.y,v.z,1.0)
+(a::Vec4,b::Vec4)=Vec4(a.x+b.x,a.y+b.y,a.z+b.z,a.w+b.w)
-(a::Vec4,b::Vec4)=Vec4(a.x-b.x,a.y-b.y,a.z-b.z,a.w-b.w)
*(a::Vec4,b::Real)=Vec4(a.x*b,a.y*b,a.z*b,a.w)
*(b::Real,a::Vec4)=*(a,b)
/(a::Vec4,b::Real)=Vec4(a.x/b,a.y/b,a.z/b,a.w)
/(b::Real,a::Vec4)=/(a,b)
==(a::Vec4,b::Vec4)=(a.x==b.x) && (a.y==b.y) && (a.z==b.z) && (a.w==b.w)

function length(V::Union(Vec3,Vec4))
  sqrt(V.x^2+V.y^2+V.z^2)
end
function norm(V::Union(Vec3,Vec4))
  V/length(V)
end

# Transformation
function translation(x,y,z)
  tm =  [1.0 0.0 0.0   x;
         0.0 1.0 0.0   y;
         0.0 0.0 1.0   z;
         0.0 0.0 0.0 1.0]
end
translation(t::Union(Vec3,Vec4))=translation(t.x,t.y,t.z)
translation(;x=0.0,y=0.0,z=0.0)=translation(x,y,z)

function rotationxyz(theta::Real,axis::Int)
  if (3<axis) || (axis <1)
    error("Only X:1, Y:2, Z:3 Axis Allowed.")
  else
    r=Vec3()
    r[axis]=theta
    rotationxyz(r)
  end
end
function rotationxyz(r::Vec3)
  ar = convert(Vector,r)
  ri = find(ar)
  rn = length(ri)
  if rn == 0
    return eye(4)
  elseif rn == 1
    axis = ri[1]
  else
    error("Only One Axis Allowed.")
  end

  c = cos(ar[axis])
  s = sin(ar[axis])

  if axis==1
    rm = [1.0 0.0 0.0 0.0;
          0.0   c  -s 0.0;
          0.0   s   c 0.0;
          0.0 0.0 0.0 1.0]
  elseif axis==2
    rm = [  c 0.0   s 0.0;
          0.0 1.0 0.0 0.0;
          -s 0.0   c 0.0;
          0.0 0.0 0.0 1.0]
  else
    rm = [  c  -s 0.0 0.0;
          s   c 0.0 0.0;
          0.0 0.0 1.0 0.0;
          0.0 0.0 0.0 1.0]
  end
end
rotationxyz(ar::AbstractVector)=rotationxyz(convert(Vec3,ar))

function reflection(r::Vec3)
  ar = convert(Vector,r)
  ri = find(ar)
  rn = length(ri)
  if rn != 2
    error("Only XY, XZ, YZ Plane Allowed.")
  end
  pi = sum(ri)
  if pi==3 # xy
    rm = [1.0  0.0  0.0  0.0;
          0.0  1.0  0.0  0.0;
          0.0  0.0 -1.0  0.0;
          0.0  0.0  0.0  1.0]
  elseif pi==4 # xz
    rm = [1.0  0.0  0.0  0.0;
          0.0 -1.0  0.0  0.0;
          0.0  0.0  1.0  0.0;
          0.0  0.0  0.0  1.0]
  else # yz
    rm = [-1.0  0.0  0.0  0.0;
          0.0  1.0  0.0  0.0;
          0.0  0.0  1.0  0.0;
          0.0  0.0  0.0  1.0]
  end
end
reflection(ar::AbstractVector)=reflection(convert(Vec3,ar))

function transform(v::Vec4,m::Matrix)
  s1,s2 = size(m)
  if (s1 != s2) || (s1 != 4)
    error("Invalid Transform Matrix.")
  end
  convert(Vec4,m*convert(Vector,v))
end
function transform{T<:Vec4}(vs::AbstractVector{T},m::Matrix)
  s1,s2 = size(m)
  if (s1 != s2) || (s1 != 4)
    error("Invalid Transform Matrix.")
  end
  convert(Vector{T},m*convert(Matrix,vs))
end
translate(v::Vec4,t::Vec3)=transform(v,translation(t))
translate{T<:Vec4}(vs::AbstractVector{T},t::Vec3)=transform(vs,translation(t))
rotatexyz(v::Vec4,r::Vec3)=transform(v,rotationxyz(r))
rotatexyz{T<:Vec4}(vs::AbstractVector{T},r::Vec3)=transform(vs,rotationxyz(r))
function transrotatez(v::Vec4,t::Vec3,theta::Real)
  tm = translation(t)
  rm = rotationxyz(Vec3(0.0,0.0,theta))
  convert(Vec4,rm*tm*convert(Vector,v))
end
function transrotatez{T<:Vec4}(vs::AbstractVector{T},t::Vec3,theta::Real)
  tm = translation(t)
  rm = rotationxyz(Vec3(0.0,0.0,theta))
  convert(Vector{T},rm*tm*convert(Matrix,vs))
end
