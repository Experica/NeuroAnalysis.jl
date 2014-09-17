import Base: convert,getindex,setindex!,+,-
export Vec3,Vec4,translation,rotationxyz,
transform,transrotatez,ang2rad,rad2ang

function ang2rad(a)
  r = a * pi/180
end
function rad2ang(r)
  a = r * 180/pi
end

# Vector3
type Vec3{T}
  x::T
  y::T
  z::T
end
Vec3{T}(x::T,y::T)=Vec3{T}(x,y,y)
Vec3{T}(x::T)=Vec3{T}(x,x,x)
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
function convert(::Type{Vec3},a::Vector)
    l=length(a)
    if l==0
        v = Vec3()
    elseif l==1
        v = Vec3(a[1])
    elseif l==2
        v = Vec3(a[1],a[2])
    else l==3
        v = Vec3(a[1],a[2],a[3])
    end
end
function convert(::Type{Vector},v::Vec3)
    a = Array(typeof(v.x),3)
    a[1]=v.x
    a[2]=v.y
    a[3]=v.z
    return a
end
+{T}(a::Vec3{T},b::Vec3{T})=Vec3{T}(a.x+b.x,a.y+b.y,a.z+b.z)
-{T}(a::Vec3{T},b::Vec3{T})=Vec3{T}(a.x-b.x,a.y-b.y,a.z-b.z)

# Vector4
type Vec4{T}
  x::T
  y::T
  z::T
  w::T
end
Vec4{T}(x::T,y::T,z::T)=Vec4{T}(x,y,z,z)
Vec4{T}(x::T,y::T)=Vec4{T}(x,y,y,y)
Vec4{T}(x::T)=Vec4{T}(x,x,x,x)
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
function convert(::Type{Vec4},a::Vector)
    l=length(a)
    if l==0
        v = Vec4()
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
  a = Array(typeof(v.x),4)
    a[1]=v.x
    a[2]=v.y
    a[3]=v.z
    a[4]=v.w
    return a
end
function convert{T}(::Type{Matrix{T}},vs::Vector{Vec4{T}})
    l = length(vs)
    if l==0
        return Array(None,4,l)
    else
        a = Array(T,4,l)
    for i in 1:l
        a[1,i]=vs[i].x
        a[2,i]=vs[i].y
        a[3,i]=vs[i].z
        a[4,i]=vs[i].w
    end
    end
    return a
end
function convert{T}(::Type{Vector{Vec4{T}}},a::Matrix{T})
    s1,s2 = size(a)
    if s1 != 4
        error("Size of Dim1 of Matrix Doesn't Match Vector4")
        elseif s2 == 0
        return []
    else
        vs = Array(Vec4{T},s2)
        for i in 1:s2
            vs[i] = Vec4(a[1,i],a[2,i],a[3,i],a[4,i])
        end
        return vs
    end
end
convert(::Type{Vec3},v::Vec4)=Vec3(v.x,v.y,v.z)
convert(::Type{Vec4},v::Vec3)=Vec4(v.x,v.y,v.z,1.0)
+{T}(a::Vec4{T},b::Vec4{T})=Vec4{T}(a.x+b.x,a.y+b.y,a.z+b.z,a.w+b.w)
-{T}(a::Vec4{T},b::Vec4{T})=Vec4{T}(a.x-b.x,a.y-b.y,a.z-b.z,a.w-b.w)

# Transformation
function translation(x,y,z)
  tm =  [1.0 0.0 0.0   x;
         0.0 1.0 0.0   y;
         0.0 0.0 1.0   z;
         0.0 0.0 0.0 1.0]
end
translation(t::Vec3)=translation(t.x,t.y,t.z)
translation(;x=0.0,y=0.0,z=0.0)=translation(x,y,z)

function rotationxyz(axis::Int,theta::Real)
  if (3<axis) || (axis <1)
    error("Only 1:x, 2:y, 3:z Axis Allowed.")
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
rotationxyz(ar::Vector)=rotationxyz(convert(Vec3,ar))

function transform(v::Vec4,m::Matrix)
  s1,s2 = size(m)
  if (s1 != s2) || (s1 != 4)
    error("Invalid Transform Matrix.")
  end
  convert(Vec4,m*convert(Vector,v))
end
function transrotatez(v::Vec4,t::Vec3,theta::Real)
    tm = translation(t)
    rm = rotationxyz(Vec3(0.0,0.0,theta))
    convert(Vec4,rm*tm*convert(Vector,v))
end
function transrotatez{T}(vs::Vector{Vec4{T}},t::Vec3,theta::Real)
    tm = translation(t)
    rm = rotationxyz(Vec3(0.0,0.0,theta))
    convert(Vector{Vec4{T}},rm*tm*convert(Matrix{T},vs))
end
