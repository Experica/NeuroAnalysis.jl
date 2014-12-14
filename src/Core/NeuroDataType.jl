export TimePoints,TPsVector,TPVVector,Spike,SpikeTrain,Cell,CellAssembly,
Block,Segment,Subject,Experiment

# Time Unit is Millisecond
typealias TimePoints{T<:Real} AbstractArray{T,1}
typealias TPsVector{T<:TimePoints} AbstractArray{T,1}
typealias TPVVector{T<:TPsVector} AbstractArray{T,1}

type Spike
  value::AbstractVector{Real}
  time::Real
  delay::Real
  sort
end

type SpikeTrain
  name::String
  id
  channel::Uint
  fs::Real
  spikes::AbstractVector{Spike}
end

type Channel
  name::String
  index::Uint64
  coordinate::AbstractVector{Float64}
  signal
  spiketrains::AbstractVector{SpikeTrain}
end

type ChannelGroup
  name::String
  channels::AbstractVector{Channel}
end

type AnalogSignal
  name::String
  description
  channel::Uint64
  fs::Float64
  value::AbstractVector{Float64}
  startime::Float64
end

type AnalogSignalArray
  name::String
  signals::AbstractVector{AnalogSignal}
end

type Event
  name::String
  value
  time::Float64
end

type EventArray
  name::String
  events::AbstractVector{Event}
end

type Epoch
  name::String
  time::Float64
  duration::Float64
  value
end

type EpochArray
  name::String
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
  name::String
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
  name::String
  id
  description
  spiketrain
  tests::AbstractVector{Block}
end

type CellAssembly
  name::String
  id
  cells::AbstractVector{Cell}
  projections
  tests::AbstractVector{Block}
end

type Subject
  name::String
  description
  contact
  cellassemblies::AbstractVector{CellAssembly}
end

type Experiment
  name::String
  description
  designers::AbstractVector{String}
  experimenters::AbstractVector{String}
  subjects::AbstractVector{Subject}
end
