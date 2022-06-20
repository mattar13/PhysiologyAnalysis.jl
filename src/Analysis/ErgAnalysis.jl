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

"""
This function uses a histogram method to find the saturation point. 
    - In ERG datas, a short nose component is usually present in saturated values
    - Does this same function work for the Rmax of nonsaturated responses?
    - Setting the saturated threshold to infinity will completely disregard the histogram method
"""
function saturated_response(data::Experiment{T}; precision::Int64=100) where {T<:Real}
    #We want to pick the region to analyze first
    norm_factor = minimum(data)
    rmaxes = zeros(size(data, 1), size(data, 3))
    minima = minimum(data, dims=2)[:, 1, :]

    for swp = 1:size(data, 1), ch = 1:size(data, 3)
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

function minima_to_peak(data::Experiment; verbose=false)
    #We need to exclude the area 
    resp = zeros(size(data, 1), size(data, 3))
    for swp = 1:size(data, 1), ch = 1:size(data, 3)
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

"""
This function calculates the time to peak using the dim response properties of the concatenated file
"""
function time_to_peak(data::Experiment{T}) where {T<:Real}
    over_stim = findall(data.t .> 0.0) #We only want to extract time points after the stim
    lowest_val = map(x -> x[2], argmin(data[:, over_stim, :], dims=2))[:, 1, :]
    lowest_val .+= over_stim[1] - 1
    data.t[lowest_val] .* 1000
end

#Pepperburg analysis
"""
This function conducts a Pepperburg analysis on a single data. 

    Two dispatches are available. 
    1) A rmax is provided, does not need to calculate rmaxes
    2) No rmax is provided, so one is calculated
"""
function pepperburg_analysis(data::Experiment{T}, rmaxes::Array{T,1};
    recovery_percent=0.60, kwargs...
) where {T<:Real}
    if size(data, 1) == 1
        throw(error("Pepperburg will not work on single sweeps"))
    end
    r_rec = rmaxes .* recovery_percent
    #try doing this  different way
    t_dom = zeros(T, size(data, 1), size(data, 3))
    for swp in 1:size(data, 1)
        for ch in 1:size(data, 3)
            not_recovered = findall(data[swp, :, ch] .< r_rec[ch])
            if isempty(not_recovered)
                #println("data never exceeded $(recovery_percent*100)% the Rmax")
                t_dom[swp, ch] = NaN
            elseif isempty(data.stim_protocol)
                #println("No stimulus protocol exists")
                t_dom[swp, ch] = data.t[not_recovered[end]]
            else
                t_dom[swp, ch] = data.t[not_recovered[end]] - data.t[data.stim_protocol[swp].index_range[1]]
            end
        end
    end
    t_dom
end

pepperburg_analysis(data::Experiment{T}; kwargs...) where {T<:Real} = pepperburg_analysis(data, saturated_response(data; kwargs...); kwargs...)

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
    @assert size(data,3) == size(rmaxes, 1)
    #Tᵣ = fill(NaN, size(data,1), size(data,3))
    Tᵣ = zeros(size(data,1), size(data,3))
    for swp in 1:size(data, 1), ch in 1:size(data, 3)
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

# The below functions are created by fitting a model 
"""
# Recovery Time Constant (τRec)

This function is a single exponential. 

### Function usage
[IN 1]:  Recovery(t, V⁰, τRec)

[OUT 1]: Response  
"""
REC(t, V⁰, τRec) = V⁰ * exp(-t / τRec)

"""
The dominant time constant is calculated by fitting the normalized Rdim with the response recovery equation
"""
function recovery_tau(data::Experiment{T}, resp::Union{T,Matrix{T}};
    τRec::T=1.0
) where {T<:Real}
    #Make sure the sizes are the same
    #@assert size(resp) == (size(data, 1), size(data,3))

    trec = zeros(T, size(data, 1), size(data, 3))
    gofs = zeros(T, size(data, 1), size(data, 3))
    #This function uses the recovery model and takes t as a independent variable
    model(x, p) = map(t -> REC(t, -1.0, p[2]), x)
    for swp in 1:size(data, 1), ch in 1:size(data, 3)
        # println(dim_idx[ch])
        xdata = data.t
        ydata = data[swp, :, ch]
        #Test both scenarios to ensure that
        ydata ./= minimum(ydata) #Normalize the Rdim to the minimum value
        #ydata ./= resp #Normalize the Rdim to the saturated response

        #cutoff all points below -0.5 and above -1.0
        over_1 = findall(ydata .>= 1.0)
        if !isempty(over_1)
            begin_rng = over_1[end]
        
            xdata = xdata[begin_rng:end]
            ydata = ydata[begin_rng:end]

            cutoff = findall(ydata .< 0.5)
            if isempty(cutoff)
                #println("Exception")
                end_rng = length(ydata)
            else
                end_rng = cutoff[1]
            end

            xdata = xdata[1:end_rng] .- xdata[1]
            ydata = -ydata[1:end_rng]
            p0 = [ydata[1], τRec]
            fit = curve_fit(model, xdata, ydata, p0)
            #report the goodness of fit
            SSE = sum(fit.resid .^ 2)
            ȳ = sum(model(xdata, fit.param)) / length(xdata)
            SST = sum((ydata .- ȳ) .^ 2)
            GOF = 1 - SSE / SST
            trec[swp, ch] = fit.param[2]
            gofs[swp, ch] = GOF
        end
    end
    return trec, gofs
end

#Function to fit Amplification
"""
# Amplification 

Amplification is a time series, therefore it is a function of time

## The relationship is demonstrated by
\$R = f(t)\$

\$f(t) = R_{max}(1-e^{-\\alpha(t-t_{eff})^2})\$

### Variables
- R: The response is the dependent variable
- t: Time is the independent variable.

### Parameters
- (\$t_{eff}\$): The effective time delay is a short delay between stimulus onset and response onset indicative of the biomolecuar diffusion rates
- (\$\\alpha\$): The amplification coefficient  represents the rate of the response increases from the biomolecular processes. 

### Function usage
[IN 1]:  AMP(t, α, t_eff, rmax)

[OUT 1]: Response

"""
AMP(t, α, t_eff, rmax) = t > t_eff ? rmax * (1 - exp(-α * (t - t_eff)^2)) : 0.0

function amplification(data::Experiment{T}, resp::Union{T,Matrix{T}}; #This argument should be offloaded to a single value 
    time_cutoff=0.1,
    lb::Vector{T}=[0.0, 0.001],
    p0::Vector{T}=[200.0, 0.002],
    ub::Vector{T}=[Inf, 0.040]
) where {T<:Real}

    #@assert size(resp) == (size(data, 1), size(data,3))

    amp = zeros(2, size(data, 1), size(data, 3))
    gofs = zeros(T, size(data, 1), size(data, 3))

    for swp = 1:size(data, 1), ch = 1:size(data, 3)
        if isa(resp, Matrix{T})
            resp_0 = resp[swp, ch]
        else
            resp_0 = resp
        end
        model(x, p) = map(t -> AMP(t, p[1], p[2], resp_0), x)
        idx_end = findall(data.t .>= time_cutoff)[1]
        xdata = data.t[1:idx_end]
        ydata = data[swp, 1:idx_end, ch]

        fit = curve_fit(model, xdata, ydata, p0, lower=lb, upper=ub)
        #Check Goodness of fit
        SSE = sum(fit.resid .^ 2)
        ȳ = sum(model(xdata, fit.param)) / length(xdata)
        SST = sum((ydata .- ȳ) .^ 2)
        GOF = 1 - SSE / SST
        amp[1, swp, ch] = fit.param[1] #Alpha amp value
        amp[2, swp, ch] = fit.param[2] #Effective time value
        gofs[swp, ch] = GOF
    end
    return amp, gofs
end