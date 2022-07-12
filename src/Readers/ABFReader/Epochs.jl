
"""
This object contains information about the stimulus Epoch coming from the Cmd

To initialize an Epoch use this function: 
    julia> Epoch()
Epochs will have a type associated {T} with the type of numerical data being T. 
For Julia this can be Float64, but for older files this may be Float32. 

### Epoch is a struct that contains: 
    1) 'epochLetter::String' -> The alphabetical label of the Epoch
    2) 'epochType::String' -> The type of epoch. Comes in terms of Step, Pulse, Ramp, Triangle, Saw.  
    3) 'dacNum::T' -> The Digital to Analog channel that contains the stimulus. 
    4) 'epochNumber::T' -> The number ID related to the Epoch
    5) 'level::T' -> The holding level of the Epoch. The value related to the command is found in the "dacUnits"
    6) 'levelDelta::T' -> The change in level for each sweep. 
    7) 'duration::T' -> How long the DAC will be held at a certai level
    8) 'durationDelta::T' -> From sweep to sweep the change in the duration
    9) 'digitalPattern::Vector'{T} -> If a digital channel is used, a digital bit pattern will be used (8 bits)
    10) 'pulsePeriod::T' -> The periodicity of the pulse (Sine Wave)
    11) 'pulseWidth::T' -> The width of the pulse (Sine Wave)


"""
mutable struct Epoch{T}
    epochLetter::String
    epochType::String
    dacNum::Int64
    epochNumber::Int64
    level::T
    levelDelta::T
    duration::T
    durationDelta::T
    digitalPattern::Vector{Int64}
    pulsePeriod::T
    pulseWidth::T
end

Epoch() = Epoch(" ", "Off", 0, 0, 0.0, 0.0, 0.0, 0.0, zeros(Int64, 8), 0.0, 0.0)
#EpochTable will be instantiated last

mutable struct EpochSweepWaveform
    p1s::Vector
    p2s::Vector
    levels::Vector
    types::Vector
    pulseWidths::Vector
    pulsePeriods::Vector
    digitalStates::Vector
end

mutable struct EpochTable
    sampleRateHz
    holdingLevel
    sweepPointCount
    channel
    epochs::Vector{Epoch}
    epochWaveformBySweep::Vector{EpochSweepWaveform}
end


EpochSweepWaveform() = EpochSweepWaveform([], [], [], [], [], [], [])

function addEpoch(e::EpochSweepWaveform, pt1, pt2, level, type, pulseWidth, pulsePeriod, digitalState)
    push!(e.p1s, round(Int64, pt1))
    push!(e.p2s, round(Int64, pt2))
    push!(e.levels, level)
    push!(e.types, type)
    push!(e.pulseWidths, pulseWidth)
    push!(e.pulsePeriods, pulsePeriod)
    push!(e.digitalStates, digitalState)
end


function EpochTable(abf::Dict{String,Any}, channel::Int64)
    #channel has to be a channel number in this case
    sampleRateHz = abf["dataRate"]
    holdingLevel = abf["holdingCommand"][channel]
    sweepPointCount = abf["sweepPointCount"]
    epochs = Epoch[]

    returnToHold = false
    if abf["abfVersionDict"]["major"] == 1 && abf["nWaveformEnable"][channel] == 1
        if channel > 1
            channel = 0
        end
        #not fully implemented yet
    elseif abf["abfVersionDict"]["major"] == 2 && abf["DACSection"]["nWaveformEnable"][channel] == 1

        #Create a list of Epochs
        for (i, epochDACNum) in enumerate(abf["EpochPerDACSection"]["nDACNum"])
            #only create a table for the selected channel
            if epochDACNum == (channel - 1)
                epoch = Epoch()
                epoch.dacNum = epochDACNum
                epochNumber = abf["EpochPerDACSection"]["nEpochNum"][i]
                epoch.epochNumber = epochNumber
                epoch.epochType = epoch_type[abf["EpochPerDACSection"]["nEpochType"][i]]
                epoch.level = abf["EpochPerDACSection"]["fEpochInitLevel"][i]
                epoch.levelDelta = abf["EpochPerDACSection"]["fEpochLevelInc"][i]
                epoch.duration = abf["EpochPerDACSection"]["lEpochInitDuration"][i]
                epoch.durationDelta = abf["EpochPerDACSection"]["lEpochDurationInc"][i]
                epoch.pulsePeriod = abf["EpochPerDACSection"]["lEpochPulsePeriod"][i]
                epoch.pulseWidth = abf["EpochPerDACSection"]["lEpochPulseWidth"][i]

                #Add the digital channel
                if epochDACNum == abf["ProtocolSection"]["nActiveDACChannel"]
                    digOut = abf["EpochSection"]["nEpochDigitalOutput"][i]
                    #println("Digital pattern: $digOut")
                    epoch.digitalPattern = digits(digOut, base = 2, pad = 8) #convert the digital output to a bit array
                else
                    epoch.digitalPattern = digits(0, base = 2, pad = 8) #convert the digital output to a bit array
                end

                #Add the Epoch Letter
                epochLetter = String[]
                num = epochNumber
                while num >= 0
                    push!(epochLetter, Char(num % 26 + 65) |> string)
                    num -= 26
                end
                epoch.epochLetter = join(epochLetter)
                #add the epochType
                push!(epochs, epoch)
            end
        end
        returnToHold = abf["DACSection"]["nInterEpisodeLevel"][channel] == 1
    end

    epochWaveformsBySweep = []
    #Create a list of waveform objects by sweep    
    lastSweepLastLevel = holdingLevel
    for sweep in abf["sweepList"]
        ep = EpochSweepWaveform()
        #Add pre epoch values
        preEpochEndPoint = abf["sweepPointCount"] / 64.0
        pt2 = preEpochEndPoint
        addEpoch(ep, 0.0, preEpochEndPoint, lastSweepLastLevel, "Step", 0, 0, zeros(Int64, 8))

        position = preEpochEndPoint
        level = holdingLevel
        for epoch in epochs
            duration = epoch.duration + epoch.durationDelta * sweep
            pt1, pt2 = (position, position + duration)
            level = epoch.level + epoch.levelDelta * sweep
            addEpoch(ep, pt1, pt2, level, epoch.epochType, epoch.pulseWidth, epoch.pulsePeriod, epoch.digitalPattern)
            position = pt2
        end
        if returnToHold
            lastSweepLastLevel = level |> Float64
        else
            lastSweepLastLevel = holdingLevel |> Float64
        end
        addEpoch(ep, pt2, abf["sweepPointCount"], lastSweepLastLevel, "Step", 0, 0, zeros(8))
        push!(epochWaveformsBySweep, ep)
    end

    return EpochTable(sampleRateHz, holdingLevel, sweepPointCount, channel, epochs, epochWaveformsBySweep)
end

function getAnalogWaveform(e::EpochSweepWaveform)
    sweepC = zeros(e.p2s[end])
    for i = 1:length(e.levels)
        #Easier access to epoch
        epochType = e.types[i]
        chunkSize = e.p2s[i] - e.p1s[i]
        pulsePeriod = e.pulsePeriods[i]
        pulseWidth = e.pulseWidths[i]
        level = e.levels[i]
        if i == 1
            levelBefore = level
        else
            levelBefore = e.levels[i]
        end
        levelDelta = level - levelBefore
        if e.pulsePeriods[i] > 0
            pulseCount = Int64(chunkSize / e.pulsePeriods[i])
        else
            pulseCount = 0
        end

        if epochType == "Step"
            chunk = fill(level, chunkSize)
        elseif epochType == "Ramp"
            chunk = LinRange(levelBefore, level, chunkSize)
        elseif epochType == "Pulse"
            chunk = fill(levelBefore, chunkSize)
            for pulse = 1:pulseCount
                p1 = Int64(pulsePeriod * pulse)
                p2 = Int64(p1 + pulseWidth)
                chunk[p1:p2] = level
            end
        elseif epochType == "Tri"
            println("to be implemented")
        end
        #digitalStateForChannel = digitalState[channel]
        #add the chunk to the sweep
        sweepC[(e.p1s[i]+1):(e.p2s[i])]
    end
    return sweepC
end

function getDigitalWaveform(e::EpochSweepWaveform, channel)
    sweepD = zeros(e.p2s[end])
    for i = 1:length(e.levels)-1
        digitalState = e.digitalStates[i]
        digitalStateForChannel = digitalState[channel] * 5 #for voltage output
        #println(e.p1s[i])
        #println(e.p2s[i])
        sweepD[(e.p1s[i]+1):(e.p2s[i]+1)] .= digitalStateForChannel
    end
    return sweepD
end