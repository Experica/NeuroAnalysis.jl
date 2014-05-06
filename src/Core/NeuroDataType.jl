export TimePoints

##############################
type AbstractVector3{T}
	x::T
	y::T
	z::T
end

type CellError <: Exception

end

typealias TimePoints{T<:Real} AbstractArray{T,1}

##############################
type Spike
	channel::Uint64
	fs::Float64
	value::AbstractVector{Float64}
	time::Float64
	delay::Float64
	sort
end

type SpikeTrain
	name::String
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

type Cell
	name::String
	celltype
	coordinate::AbstractVector{Float64}
	spiketrain::AbstractVector{Float64}
end

type CellAssemble
	name::String
	cells::AbstractVector{Cell}
end

type Segment
  name::String
	description
	starttime
  endtime
  duration
	settings::Dict
  eventarrays::AbstractVector{EventArray}
  epocharrays::AbstractVector{EpochArray}
  cellassembles::AbstractVector{CellAssemble}
  channelgroups::AbstractVector{ChannelGroup}
end

type Block
	name::String
	description
	source
	starttime
  endtime
	duration
	settings::Dict
  segments::AbstractVector{Segment}
end

type RecordSession
	name::String
	description
	region
	date
	experimenters::AbstractVector{String}
	blocks::AbstractVector{Block}
end

type Subject
	name::String
	description
	contact
	gender
	age
	height
	weight
	recordsessions::AbstractVector{RecordSession}
end

type Experiment
	name::String
	description
	designers::AbstractVector{String}
	subjects::AbstractVector{Subject}
end
