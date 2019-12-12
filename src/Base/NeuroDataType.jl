export RealVector,RealMatrix,RVVector,RVVVector,SecondPerUnit,settimeunit,timetounit

RealVector{T} = AbstractArray{T,1} where T<:Real
RealMatrix{T} = AbstractArray{T,2} where T<:Real
RVVector{T} = AbstractArray{T,1} where T<:RealVector
RVVVector{T} = AbstractArray{T,1} where T<:RVVector

"Time Unit set to millisecond by default."
SecondPerUnit = 0.001
function settimeunit(secondperunit)
    global SecondPerUnit = secondperunit
end
function timetounit(t,tunit=:ms)
    if tunit==:ms
        ts = t*0.001
    else
        ts = t
    end
    ts/SecondPerUnit
end
