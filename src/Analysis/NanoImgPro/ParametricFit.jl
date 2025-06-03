#p0 = [A, τ_ON, τ_OFF, shift]

"""
    single_stim_model(t, p)

Single exponential stimulus model with rise and decay time constants.
# Parameters
- `t::T`: Time point
- `p::Vector{T}`: Parameters [A, τ_ON, τ_OFF, shift] where:
  - A: Amplitude
  - τ_ON: Rise time constant
  - τ_OFF: Decay time constant
  - shift: Time shift of the response

Returns the model value at time t.
"""
function single_stim_model(t::T, p) where T<:Real
    y_fit = p[1] * (1-exp(-(t-p[4])/ p[3]))*exp(-(t-p[4])/p[2])
    if y_fit < 0
        return 0.0
    else
        return y_fit
    end
end

"""
    single_stim_model(t_series, p)

Vectorized version of single_stim_model for time series.
# Parameters
- `t_series::Vector{T}`: Vector of time points
- `p::Vector{T}`: Parameters [A, τ_ON, τ_OFF, shift]

Returns a vector of model values for each time point.
"""
single_stim_model(t_series::Vector{T}, p) where T<:Real = map(TIME-> single_stim_model(TIME, p), t_series)

"""
    single_stim_model_drift(t, p)

Single exponential stimulus model with linear drift term.
# Parameters
- `t::T`: Time point
- `p::Vector{T}`: Parameters [A, τ_ON, τ_OFF, shift, slope, intercept] where:
  - A: Amplitude
  - τ_ON: Rise time constant
  - τ_OFF: Decay time constant
  - shift: Time shift of the response
  - slope: Linear drift slope
  - intercept: Linear drift intercept

Returns the model value at time t, including the drift term.
"""
function single_stim_model_drift(t::T, p) where T<:Real
    y_fit = p[1] * (1-exp(-(t-p[4])/ p[3]))*exp(-(t-p[4])/p[2]) 
    drift = p[5] * t + p[6]
    if y_fit < 0
        return drift
    else
        return y_fit + drift
    end
end

"""
    single_stim_model_drift(t_series, p)

Vectorized version of single_stim_model_drift for time series.
# Parameters
- `t_series::Vector{T}`: Vector of time points
- `p::Vector{T}`: Parameters [A, τ_ON, τ_OFF, shift, slope, intercept]

Returns a vector of model values for each time point.
"""
single_stim_model_drift(t_series::Vector{T}, p) where T<:Real = map(TIME-> single_stim_model_drift(TIME, p), t_series)

"""
    fit_parametric(t_data, y_data; p0=[0.50, 100.0, 10.0, 20.0], ub=[1.0, 1000.0, 1000.0, 100.0], lb=[0.0, 0.0, 0.0, 10.0])

Fits a single exponential stimulus model to the data.
# Parameters
- `t_data::Vector{T}`: Time points
- `y_data::Vector{T}`: Response values
- `p0`: Initial parameters [A, τ_ON, τ_OFF, shift]
- `ub`: Upper bounds for parameters
- `lb`: Lower bounds for parameters

Returns the fitted parameters [A, τ_ON, τ_OFF, shift] or [NaN, NaN, NaN, NaN] if fitting fails.
"""
function fit_parametric(t_data::Vector{T}, y_data::Vector{T}; 
    p0 = [0.50, 100.0, 10.0, 20.0], ub = [1.0, 1000.0, 1000.0, 100.0], lb = [0.0, 0.0, 0.0, 10.0]
) where T<:Real
    model(t_series, p) = map(TIME-> single_stim_model(TIME, p), t_series)
    try
        fit = curve_fit(model, t_data, y_data, p0, upper=ub, lower = lb)
        return fit.param
    catch
        return [NaN, NaN, NaN, NaN]
    end
end

"""
    fit_parametric_drift(t_data, y_data; p0=[0.50, 100.0, 10.0, 20.0, 0.002, 5.0], ub=[1.0, 1000.0, 1000.0, 100.0, 0.005, 10.0], lb=[0.0, 0.0, 0.0, 10.0, -0.005, -10.0])

Fits a single exponential stimulus model with linear drift to the data.
# Parameters
- `t_data::Vector{T}`: Time points
- `y_data::Vector{T}`: Response values
- `p0`: Initial parameters [A, τ_ON, τ_OFF, shift, slope, intercept]
- `ub`: Upper bounds for parameters
- `lb`: Lower bounds for parameters

Returns the fitted parameters [A, τ_ON, τ_OFF, shift, slope, intercept] or [NaN, NaN, NaN, NaN, NaN, NaN] if fitting fails.
"""
function fit_parametric_drift(t_data::Vector{T}, y_data::Vector{T}; 
    lb = [0.0,  0.0,    0.0,    0.0,   -1.0, -10.0],
    p0 = [0.0050, 100.0, 10.0, 6.0, -0.00002, 0.057], 
    ub = [1.0,  1000.0, 1000.0, 1000.0,  0.005,  10.0]
) where T<:Real
    model(t_series, p) = map(TIME-> single_stim_model_drift(TIME, p), t_series)
    try
        fit = curve_fit(model, t_data, y_data, p0, upper=ub, lower = lb)
        return fit.param
    catch error
        println(error)
        #throw(error)
        return [NaN, NaN, NaN, NaN, NaN, NaN]
    end
end

fit_parametric_drift(t_rng::StepRangeLen, y_data::Vector{T}; kwargs...) where T<:Real = fit_parametric_drift(collect(t_rng), y_data; kwargs...)

"""
    estimate_linear_drift(t_data::Vector{T}, y_data::Vector{T}) where T<:Real

Estimates linear drift by fitting a line between the start and end points of the signal.
# Parameters
- `t_data::Vector{T}`: Time points
- `y_data::Vector{T}`: Signal values

Returns a tuple containing:
- slope: The estimated drift slope
- intercept: The estimated drift intercept
"""
function estimate_linear_drift(t_data::Vector{T}, y_data::Vector{T}) where T<:Real
    # Get start and end points
    t_start, t_end = t_data[1], t_data[end]
    y_start, y_end = y_data[1], y_data[end]
    
    # Calculate slope and intercept
    slope = (y_end - y_start) / (t_end - t_start)
    intercept = y_start - slope * t_start
    
    return slope, intercept
end

#Now we can iteratively fit the exponentials to the trace
function iterative_linear_bridge(t_data::Vector{T}, y_data::Vector{T}, stim_times::Vector; 
    kwargs...
) where T<:Real

    subtracted_trace = deepcopy(y_data)
    stim_duration = stim_times[2] - stim_times[1]
    for (i, stim_time) in enumerate(stim_times)
        t_rng = (stim_time, stim_time + stim_duration)
            
        idx_rng = findall(t_data .> t_rng[1] .&& t_data .< t_rng[2])
        t_data_section = t_data[idx_rng]
        y_data_section = y_data[idx_rng]

        #Lets fit some data to get a good initial guess for the drift
        slope, intercept = estimate_linear_drift(t_data_section, y_data_section)
        println("Slope: $slope, Intercept: $intercept")
        #Replace the subtracted data with the linear bridge
        linear_bridge = slope * t_data_section .+ intercept
        subtracted_trace[idx_rng] = linear_bridge
    end
    return subtracted_trace
end


#Now we can iteratively fit the exponentials to the trace
function iterative_parametric_baseline(t_data::Vector{T}, y_data::Vector{T}, stim_times::Vector; 
    p0_drift = [0.0050, 50.0, 100.0, 6.0, -0.00002, 0.057],
    kwargs...
) where T<:Real

    subtracted_trace = deepcopy(y_data)
    param_sets = []
    idx_data_sections = []
    t_data_sections = []
    y_data_sections = []
    y_data_sections_fit = []
    y_data_sections_fit_nodrift = []

    stim_duration = stim_times[2] - stim_times[1]
    for (i, stim_time) in enumerate(stim_times)
        t_rng = (stim_time, stim_time + stim_duration)
        println("Fitting stim $i")
        
        idx_rng = findall(t_data .> t_rng[1] .&& t_data .< t_rng[2])
        push!(idx_data_sections, idx_rng)
        println(t_rng)
        t_data_section = t_data[idx_rng]
        push!(t_data_sections, t_data_section)
        
        y_data_section = y_data[idx_rng]
        push!(y_data_sections, y_data_section)
        
        #We want to pick out some good initial parameters
        dy = maximum(y_data_section) - y_data_section[1] 
        println("dy: $dy")
        p0_drift[1] = dy
        p0_drift[4] = stim_time
        #Lets fit some data to get a good initial guess for the drift
        slope, intercept = estimate_linear_drift(t_data_section, y_data_section)
        p0_drift[5] = slope
        p0_drift[6] = intercept

        println("p0_drift: $p0_drift")
        #Now we should automatically fit the data
        params = fit_parametric_drift(t_data_section, y_data_section; p0 = p0_drift, kwargs...)
        println("Params: $params")
        push!(param_sets, params)

        #Now we should automatically fit the data
        y_data_section_fit = single_stim_model_drift(t_data_section, params)
        push!(y_data_sections_fit, y_data_section_fit)

        p0_nodrift = params
        p0_nodrift[5] = 0.0
        y_data_section_fit_nodrift = single_stim_model_drift(t_data_section, p0_nodrift)
        push!(y_data_sections_fit_nodrift, y_data_section_fit_nodrift)
        subtracted_trace[idx_rng] -= y_data_section_fit_nodrift
    end
    return idx_data_sections, t_data_sections, y_data_sections, y_data_sections_fit, y_data_sections_fit_nodrift, param_sets
end