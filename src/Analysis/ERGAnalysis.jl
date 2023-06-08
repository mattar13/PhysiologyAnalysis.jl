"""
    calculate_threshold(vm_arr::AbstractArray; Z = 4, dims = -1)

Finds the threshold of a trace by calculating the average and then adding the 4x the standard deviation. 
If using a differential solution, make sure dt is set, otherwise the standard deviation will be unevenly sampled
"""
function calculate_threshold(x::Array{T, N}; Z = 4.0, dims = -1) where {T <: Real, N}
    if dims == -1
        return [sum(x)/length(x) + Z*std(x)]
    else
        n = size(x, dims)
        mean = sum(x, dims = dims)./n
        dev = Z * std(x, dims = dims)
        return mean + dev #We want these all to come out as vectors vs matrices
    end
end

"""
This is useful for finding sequences of items. Requires a bitvector
"""
function findsequential(sequence::BitVector; seq_to_find=:all)
    #first we need to do a normal findall
    sequences = Vector{Int64}[] #Save all sequences in a vector of sequences
    current = Int64[]
    for (idx, itm) in enumerate(sequence)
        if itm #If an item is true, push it to the current sequence
            push!(current, idx)
        elseif !isempty(current) && !itm #If the current sequence is not empty and the item is false
            push!(current, idx) #push the item to the current sequence
            if seq_to_find == :first #if we only want to find the first sequence, then return
                return current
            end
            push!(sequences, current) #push the current sequence to all sequences
            current = Int64[] #clear the current sequence
        end
    end
    if seq_to_find == :last
        return sequences[end]
    else
        return sequences
    end
end

function gaussian_saturation(data::Experiment{T};
    p0=[-250.0, 0.5, 80.0], ub=[0.0, 1.0, 10e3]
) where {T<:Real}
    rmaxes = zeros(size(data, 1), size(data, 3))
    f(xs, p) = map(x -> p[1] * exp((-(x - p[2])^2) / 2 * p[3]), xs)
    #model(xs, p) = map(x -> f(x, p), xs)
    t = data.t

    for swp = axes(data, 1), ch = axes(data, 3)
        ȳ = data.data_array[swp, :, ch]

        fit = curve_fit(f, t, ȳ, p0, upper=ub)
        y = f(t, fit.param) #Generate the model fit

        idx_min = argmin(y)
        resp_val = ȳ[idx_min]
        rmaxes[swp, ch] = resp_val
    end
    return rmaxes
end

function histogram_saturation(data::Experiment{T}; precision::Int64=100) where {T<:Real}

    norm_factor = minimum(data)
    rmaxes = zeros(size(data, 1), size(data, 3))
    minima = minimum(data, dims=2)[:, 1, :]

    for swp = axes(data, 1), ch = axes(data, 3)
        #Lets try to quickly zero any positive results
        #y_data = data[swp, :, ch]
        #y_data *= y_data[swp, y_data .==]
        y_data = data[swp, :, ch] ./ norm_factor
        hfit = Distributions.fit(Histogram, y_data, LinRange(0.15, 1.0, precision))
        weights = hfit.weights / maximum(hfit.weights)
        edges = collect(hfit.edges[1])[1:length(weights)]
        resp = edges[argmax(weights)]
        #println(minimum(edges))
        if resp == minimum(edges)
            #println("No Nose")
            #println(minima[swp, ch])
            rmaxes[swp, ch] = minima[swp, ch]
        else
            #println("Nose")
            #println(resp)
            rmaxes[swp, ch] = resp * norm_factor
        end
    end
    rmaxes
end

"""
This function uses a histogram method to find the saturation point. 
    - In ERG datas, a short nose component is usually present in saturated values
    - Does this same function work for the Rmax of nonsaturated responses?
    - Setting the saturated threshold to infinity will completely disregard the histogram method
"""
function saturated_response(data::Experiment{T}; mode = :Logistic, kwargs...) where {T<:Real}
    #We want to pick the region to analyze first
    if mode == :Gaussian
        #This mode is better for data that has 
        rmaxes = gaussian_saturation(data; kwargs...)
    elseif mode == :Histogram
        rmaxes = histogram_saturation(data; kwargs...)
    elseif mode == :Logistic
        rmaxes =  zeros(size(data,1), size(data,3))
        nose_peak = findNosePeak(data; kwargs...) #First use the logistic function to fit out the nose
        resp = minimum(data, dims = 2)[:, 1, :] #Then find the minimum
        for i in axes(nose_peak,1)
            vals = (resp[:, i] .>= nose_peak[i]) .* resp[:,i]
            sats = (resp[:, i] .< nose_peak[i]) .* nose_peak[i]
            rmaxes[:, i] .= (vals .+ sats)
       end
       #println(rmaxes)
       return rmaxes
    end
end

function minima_to_peak(data::Experiment; verbose=false)
    #We need to exclude the area 
    resp = zeros(size(data, 1), size(data, 3))
    for swp = axes(data, 1), ch = axes(data, 3)
        past_stim = findall(data.t .> 0.0)

        data_section = data[swp, past_stim, ch] #Isolate all the items past the stim
        # cutoff the analysis at the maximum (A-wave is before the B-wave)
        cutoff_idx = argmax(data_section)
        max_val = maximum(data_section)
        data_section = data_section[1:cutoff_idx]
        #Measure the minimum betweent the first value and the maximum
        min_val = minimum(data_section)
        if verbose
            println("Minimum: $min_val")
            println("Maximum: $max_val")
            println(max_val - min_val)
        end

        resp[swp, ch] = max_val - min_val
    end
    resp
end


function findRDIM(responses::Vector{T}, rng = (0.1, 0.4)) where T <: Real
     #This section we need to extract Rdim responses. 
     normalized_responses = abs.(responses) ./ maximum(abs.(responses))
     #println(normalized_responses)
     rdim_idxs = findall(rng[1] .< normalized_responses .< rng[2]) #Basically the rdim will now be any response under the half saturation point
     if isempty(rdim_idxs)
          rdim_idx = argmin(responses)
     else
          rdim_min = argmax(responses[rdim_idxs])
          rdim_idx = rdim_idxs[rdim_min]
     end
     rdim_idx
 end


"""
This function calculates the time to peak using the dim response properties of the concatenated file
"""
function time_to_peak(data::Experiment{T}) where {T<:Real}
    over_stim = findall(data.t .> 0.0) #We only want to extract time points after the stim
    lowest_val = map(x -> x[2], argmin(data[:, over_stim, :], dims=2))[:, 1, :]
    lowest_val .+= over_stim[1] - 1
    data.t[lowest_val] .* 1000
end

"""
This function is the amount of time that a certain trace spends in a particular bandwith. 
    I think it will be similar to the pepperburg, So this may become that function
    The "criterion" is the percent. 
    This function will measure how long a response takes to return to a specific criterion amount (iᵣ)
        -By default iᵣ is set to 0.60. 
    
    An intial problem is the tendancy for the function to pick up drift and other packets. We can eliminate non-sequential packets
    For more information on this function see 
        Pepperburg & Cornwall et al. Light-dependent delay in the falling phase of the retinal rod photoresponse

    Use: 
    >>> rmaxes = saturated_response(data1_testA)
    >>> Tᵣ = percent_recovery_interval(data1_testA, rmaxes)

"""
function percent_recovery_interval(data::Experiment{T}, rmaxes::Matrix{T}; iᵣ::T=0.50) where {T<:Real}
    #first we can normalize the data to a range
    @assert size(data,3) == size(rmaxes, 2) #rmax data matches data channels
    #Tᵣ = fill(NaN, size(data,1), size(data,3))
    Tᵣ = zeros(size(data,1), size(data,3))
    for swp in axes(data, 1), ch in axes(data, 3)
        data_percent = data.data_array[swp, :, ch] ./ rmaxes[ch]
        recovery_seqs = findsequential(data_percent .> iᵣ, seq_to_find=:all)
        #we have to eliminate all datavalues under data.t = 0.0
        after_stim = findall(map(seq -> all(data.t[seq] .> 0.0), recovery_seqs)) #this returns false anywhere where the data.t is less than 0.0
        recovery_seqs = recovery_seqs[after_stim] #this eliminates all sequences with less than 0.0 time
        if !isempty(recovery_seqs) #if there are no sequences then return 0.0
            long_seq = argmax(length.(recovery_seqs))
            recovery_idx = recovery_seqs[long_seq][end]
            Tᵣ[swp, ch] = data.t[recovery_idx]
        end
    end
    return Tᵣ
end

"""
The integration time is fit by integrating the dim flash response and dividing it by the dim flash response amplitude
- A key to note here is that the exact f(x) of the ERG data is not completely known
- The integral is therefore a defininte integral and a sum of the area under the curve
- This equation is according to 
    Baylor DA, Hodgkin AL (1973) Detection and resolution of visual stimuli
        by turtle photoreceptors. J Physiol 234:163–198.
"""
function integral(data::Experiment{T}) where {T<:Real}
    #we want this to be equal to any response after the stimuli
    data_section = data[:, data.t.>0.0, :]
    data_section = abs.(data_section)
    data_section ./= maximum(data_section, dims=2)
    return sum(data_section, dims=2) * data.dt
end
