"""
It may become useful (in the case of Pluto.jl) to read the data directly from Binary Data
"""
function readABFInfo(::Type{T}, binary_file::Vector{UInt8};
    loadData=true, data_format=[3, 2, 1]
) where {T<:Real}
    #the first 4 bytes contain the version 
    #convert the binary_file into a IOBuffer
    f = IOBuffer(binary_file)
    headerSection = readHeaderSection(f) #version is automatically determined
    dataByteStart = headerSection["dataByteStart"]
    dataPointCount = headerSection["dataPointCount"]
    nDataPoints = headerSection["nDataPoints"]
    sweepPointCount = headerSection["sweepPointCount"]
    sweepCount = headerSection["sweepCount"]
    channelCount = headerSection["channelCount"]
    dataType = headerSection["dataType"]
    dataGain = headerSection["dataGain"]
    dataOffset = headerSection["dataOffset"]

    if loadData
        seek(f, dataByteStart)
        raw = read(f, dataPointCount * sizeof(dataType)) #Read the raw data into a byte array
        raw = reinterpret(dataType, raw) #convert the byte array into the dataType
        raw = reshape(raw, (channelCount, nDataPoints)) #Reshape the raw data array
        if dataType == Int16
            raw = raw .* dataGain #Multiply the data by the gain
            raw = raw .+ dataOffset #Scale the data by the offset
        end
        #We can try to convert the data into a array of shape [sweeps, data, channels]
        raw = reshape(raw, (channelCount, sweepPointCount, sweepCount)) #Reshape the data
        raw = permutedims(raw, data_format) #permute the dims
        headerSection["data"] = Array{T}(raw)
    end
    #We want to try to read more info from the DAC

    return headerSection
end

readABFInfo(binary_file::Vector{UInt8}; kwargs...) = readABFInfo(Float64, binary_file; kwargs...)

"""
This scans the axon binary and extracts all the most useful header information

For datashape
    1 -> channels
    2 -> dataspan
    3 -> Sweeps
"""
function readABFInfo(::Type{T}, filename::String; loadData=true, data_format=[3, 2, 1]) where {T<:Real}
    #Can we go through and convert anything before loading
    file_dir = splitpath(filename)

    headerSection = Dict{String,Any}()
    open(filename, "r") do f #Do everything within this loop
        #the first 4 bytes contain the version 
        headerSection = readHeaderSection(f) #version is automatically determined
        dataByteStart = headerSection["dataByteStart"]
        dataPointCount = headerSection["dataPointCount"]
        nDataPoints = headerSection["nDataPoints"]
        sweepPointCount = headerSection["sweepPointCount"]
        sweepCount = headerSection["sweepCount"]
        channelCount = headerSection["channelCount"]
        dataType = headerSection["dataType"]
        dataGain = headerSection["dataGain"]
        dataOffset = headerSection["dataOffset"]

        if loadData
            seek(f, dataByteStart)
            raw = read(f, dataPointCount * sizeof(dataType)) #Read the raw data into a byte array
            raw = reinterpret(dataType, raw) #convert the byte array into the dataType
            raw = reshape(raw, (channelCount, nDataPoints)) #Reshape the raw data array
            if dataType == Int16
                raw = raw .* dataGain #Multiply the data by the gain
                raw = raw .+ dataOffset #Scale the data by the offset
            end
            #We can try to convert the data into a array of shape [sweeps, data, channels]
            raw = reshape(raw, (channelCount, sweepPointCount, sweepCount)) #Reshape the data
            raw = permutedims(raw, data_format) #permute the dims
            headerSection["data"] = Array{T}(raw)
        end

        #We want to try to read more info from the DAC
    end

    if length(file_dir) > 1
        headerSection["abfPath"] = filename
        headerSection["abfFolder"] = joinpath(file_dir[1:end-1]...)
    else
        headerSection["abfPath"] = filename
        headerSection["abfFolder"] = "\\"
    end

    return headerSection #Finally return the headerSection as a dictionary
end

readABFInfo(filename::String; kwargs...) = readABFInfo(Float64, filename::String; kwargs...)


"""
This function opens the ABF file in clampfit
"""
function openABF(abf_path::String)
    try
        mycmd = `explorer.exe $(abf_path)`
        run(mycmd)
    catch
        #for some reason this throws an error but still opens
    end
end

openABF(abfDict::Dict{String,Any}) = openABF(abfDict["abfPath"])