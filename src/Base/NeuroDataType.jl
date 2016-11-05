export RealVector,RVVector,RVVVector,SecondPerUnit,
Spike,SpikeTrain,Cell,CellAssembly,Block,Segment,Subject,Experiment

typealias RealVector{T<:Real} AbstractArray{T,1}
typealias RVVector{T<:RealVector} AbstractArray{T,1}
typealias RVVVector{T<:RVVector} AbstractArray{T,1}

"Time Unit is millisecond by default."
const SecondPerUnit = 0.001

type Spike
  value::RealVector
  time::Real
  delay
  sort
end

type SpikeTrain
  name::AbstractString
  id
  channel::UInt
  fs::Real
  spikes::AbstractVector{Spike}
end

type Channel
  name::AbstractString
  index::UInt64
  coordinate::AbstractVector{Float64}
  signal
  spiketrains::AbstractVector{SpikeTrain}
end

type ChannelGroup
  name::AbstractString
  channels::AbstractVector{Channel}
end

type AnalogSignal
  name::AbstractString
  description
  channel::UInt64
  fs::Float64
  value::AbstractVector{Float64}
  startime::Float64
end

type AnalogSignalArray
  name::AbstractString
  signals::AbstractVector{AnalogSignal}
end

type Event
  name::AbstractString
  value
  time::Float64
end

type EventArray
  name::AbstractString
  events::AbstractVector{Event}
end

type Epoch
  name::AbstractString
  time::Float64
  duration::Float64
  value
end

type EpochArray
  name::AbstractString
  epochs::AbstractVector{Epoch}
end

type Segment
  id
  starttime
  endtime
  duration
  param::Dict
  data
end

type Block
  name::AbstractString
  id
  description
  source
  starttime
  endtime
  duration
  setting::Dict
  segments::AbstractVector{Segment}
end

type Cell
  name::AbstractString
  id
  description
  spiketrain
  tests::AbstractVector{Block}
end

type CellAssembly
  name::AbstractString
  id
  cells::AbstractVector{Cell}
  projections
  tests::AbstractVector{Block}
end

type Subject
  name::AbstractString
  description
  contact
  cellassemblies::AbstractVector{CellAssembly}
end

type Experiment
  name::AbstractString
  description
  designers::AbstractVector{AbstractString}
  experimenters::AbstractVector{AbstractString}
  subjects::AbstractVector{Subject}
end
