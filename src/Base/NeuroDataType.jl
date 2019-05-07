export RealVector,RealMatrix,RVVector,RVVVector,SecondPerUnit

RealVector{T} = AbstractArray{T,1} where T<:Real
RealMatrix{T} = AbstractArray{T,2} where T<:Real
RVVector{T} = AbstractArray{T,1} where T<:RealVector
RVVVector{T} = AbstractArray{T,1} where T<:RVVector

"Time Unit set to millisecond by default."
SecondPerUnit = 0.001
