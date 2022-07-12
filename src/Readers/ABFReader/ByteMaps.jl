#We need a map of corresponding bytes and byte_types 
ByteDict = Dict(
    :L => Int64, :l => UInt64,
    :I => Int32, :i => UInt32,
    :H => Int16, :h => UInt16,
    :B => Int8, :b => UInt8,
    :s => String, :f => Float32,
)

digitizers = Dict(
    0 => "Unknown",
    1 => "Demo",
    2 => "MiniDigi",
    3 => "DD132X",
    4 => "OPUS",
    5 => "PATCH",
    6 => "Digidata 1440",
    7 => "MINIDIGI2",
    8 => "Digidata 1550"
)

epoch_type = Dict(
    0 => "Off",
    1 => "Step",
    2 => "Ramp",
    3 => "Pulse",
    4 => "Tri",
    5 => "Cos",
    7 => "Biphasic"
)

#This is the default ABF2 Header bytemap
header_bytemap = [
    ("sFileSignature", "4s"),
    ("fFileVersionNumber", "4b"),
    ("uFileInfoSize", "I"),
    ("lActualEpisodes", "I"),
    ("uFileStartDate", "I"),
    ("uFileStartTimeMS", "I"),
    ("uStopwatchTime", "I"),
    ("nFileType", "H"),
    ("nDataFormat", "H"),
    ("nSimultaneousScan", "H"),
    ("nCRCEnable", "H"),
    ("uFileCRC", "I"),
    ("FileGUID", "16b"),
    ("uCreatorVersion", "4b"),
    ("uCreatorNameIndex", "I"),
    ("uModifierVersion", "I"),
    ("uModifierNameIndex", "I"),
    ("uProtocolPathIndex", "I"),
    #These are the useful sections
    ("ProtocolSection", "IIL", 76),
    ("ADCSection", "IIL", 92),
    ("DACSection", "IIL", 108),
    ("EpochSection", "IIL", 124),
    ("ADCPerDACSection", "IIL", 140),
    ("EpochPerDACSection", "IIL", 156),
    ("StringsSection", "IIL", 220),
    ("DataSection", "IIL", 236),
    ("TagSection", "IIL", 252),
    # Sections I don't find useful
    ("UserListSection", "IIL", 172),
    ("StatsRegionSection", "IIL", 188),
    ("MathSection", "IIL", 204),
    ("ScopeSection", "IIL", 268),
    ("DeltaSection", "IIL", 284),
    ("VoiceTagSection", "IIL", 300),
    ("SynchArraySection", "IIL", 316),
    ("AnnotationSection", "IIL", 332),
    ("StatsSection", "IIL", 348)
]

protocol_bytemap = [
    ("nOperationMode", "h"),            # 0
    ("fADCSequenceInterval", "f"),      # 2
    ("bEnableFileCompression", "b"),    # 6
    ("_sUnused", "3b"),                 # 7
    ("uFileCompressionRatio", "I"),     # 10
    ("fSynchTimeUnit", "f"),            # 14
    ("fSecondsPerRun", "f"),            # 18
    ("lNumSamplesPerEpisode", "I"),     # 22
    ("lPreTriggerSamples", "I"),        # 26
    ("lEpisodesPerRun", "I"),           # 30
    ("lRunsPerTrial", "I"),             # 34
    ("lNumberOfTrials", "I"),           # 38
    ("nAveragingMode", "H"),            # 42
    ("nUndoRunCount", "H"),             # 44
    ("nFirstEpisodeInRun", "H"),        # 46
    ("fTriggerThreshold", "f"),         # 48
    ("nTriggerSource", "H"),            # 52
    ("nTriggerAction", "H"),            # 54
    ("nTriggerPolarity", "H"),          # 56
    ("fScopeOutputInterval", "f"),      # 58
    ("fEpisodeStartToStart", "f"),      # 62
    ("fRunStartToStart", "f"),          # 66
    ("lAverageCount", "I"),             # 70
    ("fTrialStartToStart", "f"),        # 74
    ("nAutoTriggerStrategy", "H"),      # 78
    ("fFirstRunDelayS", "f"),           # 80
    ("nChannelStatsStrategy", "H"),     # 84
    ("lSamplesPerTrace", "I"),          # 86
    ("lStartDisplayNum", "I"),          # 90
    ("lFinishDisplayNum", "I"),         # 94
    ("nShowPNRawData", "H"),            # 98
    ("fStatisticsPeriod", "f"),         # 100
    ("lStatisticsMeasurements", "I"),   # 104
    ("nStatisticsSaveStrategy", "H"),   # 108
    ("fADCRange", "f"),                 # 110
    ("fDACRange", "f"),                 # 114
    ("lADCResolution", "I"),            # 118
    ("lDACResolution", "I"),            # 122
    ("nExperimentType", "H"),           # 126
    ("nManualInfoStrategy", "H"),       # 128
    ("nCommentsEnable", "H"),           # 130
    ("lFileCommentIndex", "I"),         # 132
    ("nAutoAnalyseEnable", "H"),        # 136
    ("nSignalType", "H"),               # 138
    ("nDigitalEnable", "H"),            # 140
    ("nActiveDACChannel", "H"),         # 142
    ("nDigitalHolding", "H"),           # 144
    ("nDigitalInterEpisode", "H"),      # 146
    ("nDigitalDACChannel", "H"),        # 148
    ("nDigitalTrainActiveLogic", "H"),  # 150
    ("nStatsEnable", "H"),              # 152
    ("nStatisticsClearStrategy", "H"),  # 154
    ("nLevelHysteresis", "H"),          # 156
    ("lTimeHysteresis", "I"),           # 158
    ("nAllowExternalTags", "H"),        # 162
    ("nAverageAlgorithm", "H"),         # 164
    ("fAverageWeighting", "f"),         # 166
    ("nUndoPromptStrategy", "H"),       # 170
    ("nTrialTriggerSource", "H"),       # 172
    ("nStatisticsDisplayStrategy", "H"),# 174
    ("nExternalTagType", "H"),          # 176
    ("nScopeTriggerOut", "H"),          # 178
    ("nLTPType", "H"),                  # 180
    ("nAlternateDACOutputState", "H"),  # 182
    ("nAlternateDigitalOutputState", "H"),  # 184
    ("fCellID", "fff"),                  # 186 #This one might not work right
    ("nDigitizerADCs", "H"),            # 198
    ("nDigitizerDACs", "H"),            # 200
    ("nDigitizerTotalDigitalOuts", "H"),  # 202
    ("nDigitizerSynchDigitalOuts", "H"),  # 204
    ("nDigitizerType", "H"),            # 206
]

adc_bytemap = [
    ("nADCNum", "H"),  # 0
    ("nTelegraphEnable", "H"),  # 2
    ("nTelegraphInstrument", "H"),  # 4
    ("fTelegraphAdditGain", "f"),  # 6
    ("fTelegraphFilter", "f"),  # 10
    ("fTelegraphMembraneCap", "f"),  # 14
    ("nTelegraphMode", "H"),  # 18
    ("fTelegraphAccessResistance", "f"),  # 20
    ("nADCPtoLChannelMap", "H"),  # 24
    ("nADCSamplingSeq", "H"),  # 26
    ("fADCProgrammableGain", "f"),  # 28
    ("fADCDisplayAmplification", "f"),  # 32
    ("fADCDisplayOffset", "f"),  # 36
    ("fInstrumentScaleFactor", "f"),  # 40
    ("fInstrumentOffset", "f"),  # 44
    ("fSignalGain", "f"),  # 48
    ("fSignalOffset", "f"),  # 52
    ("fSignalLowpassFilter", "f"),  # 56
    ("fSignalHighpassFilter", "f"),  # 60
    ("nLowpassFilterType", "b"),  # 64
    ("nHighpassFilterType", "b"),  # 65
    ("fPostProcessLowpassFilter", "f"),  # 66
    ("nPostProcessLowpassFilterType", "s"),  # 70
    ("bEnabledDuringPN", "b"),  # 71
    ("nStatsChannelPolarity", "H"),  # 72
    ("lADCChannelNameIndex", "I"),  # 74
    ("lADCUnitsIndex", "I"),]


dac_bytemap = [
    ("nDACNum", "H"),  # 0
    ("nTelegraphDACScaleFactorEnable", "H"),  # 2
    ("fInstrumentHoldingLevel", "f"),  # 4
    ("fDACScaleFactor", "f"),  # 8
    ("fDACHoldingLevel", "f"),  # 12
    ("fDACCalibrationFactor", "f"),  # 16
    ("fDACCalibrationOffset", "f"),  # 20
    ("lDACChannelNameIndex", "I"),  # 24
    ("lDACChannelUnitsIndex", "I"),  # 28
    ("lDACFilePtr", "I"),  # 32
    ("lDACFileNumEpisodes", "I"),  # 36
    ("nWaveformEnable", "H"),  # 40
    ("nWaveformSource", "H"),  # 42
    ("nInterEpisodeLevel", "H"),  # 44
    ("fDACFileScale", "f"),  # 46
    ("fDACFileOffset", "f"),  # 50
    ("lDACFileEpisodeNum", "I"),  # 54
    ("nDACFileADCNum", "H"),  # 58
    ("nConditEnable", "H"),  # 60
    ("lConditNumPulses", "I"),  # 62
    ("fBaselineDuration", "f"),  # 66
    ("fBaselineLevel", "f"),  # 70
    ("fStepDuration", "f"),  # 74
    ("fStepLevel", "f"),  # 78
    ("fPostTrainPeriod", "f"),  # 82
    ("fPostTrainLevel", "f"),  # 86
    ("nMembTestEnable", "H"),  # 90
    ("nLeakSubtractType", "H"),  # 92
    ("nPNPolarity", "H"),  # 94
    ("fPNHoldingLevel", "f"),  # 96
    ("nPNNumADCChannels", "H"),  # 100
    ("nPNPosition", "H"),  # 102
    ("nPNNumPulses", "H"),  # 104
    ("fPNSettlingTime", "f"),  # 106
    ("fPNInterpulse", "f"),  # 110
    ("nLTPUsageOfDAC", "H"),  # 114
    ("nLTPPresynapticPulses", "H"),  # 116
    ("lDACFilePathIndex", "I"),  # 118
    ("fMembTestPreSettlingTimeMS", "f"),  # 122
    ("fMembTestPostSettlingTimeMS", "f"),  # 126
    ("nLeakSubtractADCIndex", "H"),  # 130
]

EpochPerDAC_bytemap = [
    ("nEpochNum", "H"),  # 0
    ("nDACNum", "H"),  # 2
    ("nEpochType", "H"),  # 4
    ("fEpochInitLevel", "f"),  # 6
    ("fEpochLevelInc", "f"),  # 10
    ("lEpochInitDuration", "I"),  # 14
    ("lEpochDurationInc", "I"),  # 18
    ("lEpochPulsePeriod", "I"),  # 22
    ("lEpochPulseWidth", "I")  # 26
]


"""
These functions handle the byte interpretations of the ABF file

"""
function readStruct(f::IO, byteType::String)
    b = UInt8[0x00] #Either this needs to be UInt8 or Int8
    if length(byteType) > 1
        first_char = tryparse(Int64, byteType[1] |> string)
        if isnothing(first_char)
            vals = map(c -> readStruct(f, string(c))[1], collect(byteType))
            return vals #Return statement vs continuing
        else
            type_conv = ByteDict[Symbol(byteType[end])]
            if type_conv == String
                n_bytes = parse(Int64, byteType[1:end-1])
            else
                n_bytes = parse(Int64, byteType[1:end-1]) * sizeof(type_conv)
            end
        end
    else
        type_conv = ByteDict[Symbol(byteType)]
        if type_conv == String
            n_bytes = 1
        else
            n_bytes = sizeof(type_conv)
        end
    end

    readbytes!(f, b, n_bytes)
    if type_conv == String
        val = b |> String
    elseif type_conv == Int16 || type_conv == UInt16
        val = Vector{type_conv}(b)[1:2:end]
    elseif type_conv == Int8 || type_conv == UInt8
        val = bytes2hex(b)
    else
        val = reinterpret(type_conv, b) |> Array
    end
    return val
end

function readStruct(f::IO, byteType::String, seekTo::Int64; repeat = 1)
    seek(f, seekTo)
    if repeat == 1
        return readStruct(f, byteType)
    else
        return map(x -> readStruct(f, byteType), 1:repeat)
    end
end

