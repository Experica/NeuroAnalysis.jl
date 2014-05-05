type AbstratVector3{T}
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
	value::AbstratVector{Float64}
	time::Float64
	delay::Float64
	sort
end

type SpikeTrain
	name::String
	spikes::AbstratVector{Spike}
end

type Channel
	name::String
	index::Uint64
	coordinate::AbstratVector3{Float64}
	signal
	spiketrains::AbstratVector{SpikeTrain}
end

type ChannelGroup
	name::String
	channels::AbstratVector{Channel}
end

type AnalogSignal
	name::String
	description
	channel::Uint64
	fs::Float64
	value::AbstratVector{Float64}
	startime::Float64
end

type AnalogSignalArray
  name::String
  signals::AbstratVector{AnalogSignal}
end

type Event
	name::String
	value
	time::Float64
end

type EventArray
	name::String
	events::AbstratVector{Event}
end

type Epoch
	name::String
	time::Float64
	duration::Float64
	value
end

type EpochArray
	name::String
	epochs::AbstratVector{Epoch}
end

type Cell
	name::String
	celltype
	coordinate::AbstratVector3{Float64}
	spiketrain::AbstratVector{Float64}
end

type CellAssemble
	name::String
	cells::AbstratVector{Cell}
end

type Segment
  name::String
	description
	starttime
  endtime
  duration
	settings::Dict
  eventarrays::AbstratVector{EventArray}
  epocharrays::AbstratVector{EpochArray}
  cellassembles::AbstratVector{CellAssemble}
  channelgroups::AbstratVector{ChannelGroup}
end

type Block
	name::String
	description
	source
	starttime
  endtime
	duration
	settings::Dict
  segments::AbstratVector{Segment}
end

type RecordSession
	name::String
	description
	region
	date
	experimenters::AbstratVector{String}
	blocks::AbstratVector{Block}
end

type Subject
	name::String
	description
	contact
	gender
	age
	height
	weight
	recordsessions::AbstratVector{RecordSession}
end

type Experiment
	name::String
	description
	designers::AbstratVector{String}
	subjects::AbstratVector{Subject}
end
