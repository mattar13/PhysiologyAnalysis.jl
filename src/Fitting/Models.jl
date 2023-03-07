"""
This function is used to calculate the photon density based on the photon energy from the calibrator

The equation for this function is as follows

E = Photon Energy
C = speed of Light
E/(C)
"""
photons(E::Float64; λ::Int64 = 525) = ( E / (6.626e-34 * 3e8/(λ*10e-9)))*10e-8


"""
This function is used to calculate the transferrance from the optical density (D)
"""
Transferrance(D) = 10^-D

#For the next equation we have 3 fit variables previously determined


"""
This function is the relationship between: 
    Independent Variables: 
        Transmittance (T) 
        LED Percent (I)
        Stimulus Time (t_stim)
    Dependent Variable
        Photons (P)
    x[1] = T
    x[2] = I
    x[3] = t_stim
"""
stimulus_model(x::Array{T,1}, p::Array{Float64,1}) where T <:Real = x[1]*(p[1]*x[2]^2 + p[2]*x[2] + p[3])*x[3]
stimulus_model(x::Array{T,2}, p::Array{Float64,1}) where T <:Real = [stimulus_model(x[i,:], p) for i in 1:size(x,1)]
stimulus_model(x::Array{T,1}) where T <:Real = stimulus_model(x, [25352.59, 43857.01, 929.56]) #Green stimuli
stimulus_model(x::Array{T,2}) where T <:Real = stimulus_model(x, [25352.59, 43857.01, 929.56]) 
f_I(ND::Float64, P::Float64, t_stim::Float64) = stimulus_model([Transferrance(ND), P, t_stim])
f_I(ND::Int64, P::Int64, t_stim::Int64) = stimulus_model([Transferrance(ND|>Float64), P|>Float64, t_stim|>Float64])
##############################These are the IR and Amplification models#############

"""
# Adult Intensity-Response models

## The relationship is demonstrated by 
\$R = f(I)\$ 

\$f(I) = R_{max}\\frac{I^n}{I^n_{1/2}+I^n}\$

if Response values are normalized to 1, then \$R_{max}\$ = 1 and can be cancelled out to form the equations

### Variables: 
- R: The response amplitude is the dependent variable
- I: The stimulus light intensity (I) is the independent variable
### Parameters: 
- R_max: Maximum saturating value(\$R_{max}\$)
- Ih: The flash strength required to elicit half of \$R_{max}\$: (\$I_{1/2}\$)
- n: The power of the equation
### Function usage
[IN 1]:  IR(I, Ih, n)

[OUT 1]: Response
"""
HILL(x, rmax, k, n) = rmax * (x^n / (k^n + x^n)) #This is the basic form of the model
HILL_MODEL(X, p) = map(x -> HILL(x, p[1], p[2], p[3]), X) #This is used for fitting a larger dataset
"""
# Developmental Intensity response (>P14)

## The relationship is demonstrated by 
\$R = f(I)\$ 
 
where 

\$f(I) =R_{max}\\left(\\alpha(1 - e^{SI}) + (1-\\alpha)\\frac{I^n}{Ih^n + S}\$

if Response values are normalized to 1, then \$R_{max}\$ = 1 and can be cancelled out to form the equations

### Variables: 
- R: The response amplitude is the dependent variable
- I: The stimulus light intensity (I) is the independent variable
### Parameters: 
- R_max: Maximum saturating value(\$R_{max}\$)
- Ih: The flash strength required to elicit half of \$R_{max}\$: (\$I_{1/2}\$)
- n: The power of the equation
- (\$\\alpha\$): The temperature-dependent weighting coefficient:  
- S: he fractional sensitivity
### Function usage
[IN 1]:  IR_dev(I, rmax, k, n, α, S)

[OUT 1]: Response_dev
"""
modHILL(x, rmax, k, n, α, S) = rmax * (α*(1-exp(S*x)) + (1-α)*(x^n / (k^n + S)))
modHILL_MODEL(X, p) = map(x -> modHILL(x, p[1], p[2], p[3], p[4], p[5]), X) #This is used for fitting a larger dataset
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
AMP(t, α, t_eff, rmax) = t > t_eff ? rmax * (1 - exp(-α*(t-t_eff)^2)) : 0.0

"""
# Recovery Time Constant (τRec)

This function is a single exponential. 

### Function usage
[IN 1]:  Recovery(t, V⁰, τRec)

[OUT 1]: Response  
"""
REC(t, V⁰, τRec) = V⁰ * exp(-t/τRec)

"""
Weber Contrast sensitivity

The
"""

WEBER(I_Feature::T, I_Background::T) where T <: Real = (I_Feature - I_Background)/I_Background

"""
Michelson Contrast


"""
MICHELSON(I_Min::T, I_Max::T) where T <: Real = (I_Max - I_Min)/(I_Max + I_Min)
MICHELSON(I::Array{T,2}) where T <: Real = (maximum(I) - minimum(I))/(maximum(I) + minimum(I))

"""
Root Mean Squared Contrast

This is the contrast of an image when the image is 
M x N in size. 

The image has i and j features (equal to 1->M-1 and 1->N-1)
"""
RMS_Contrast(I::Array{T, 2}; normalized = true) where T <: Real = (1/(size(I,1)*size(I,2))) .* sum((I.-sum(I)/length(I))^2)

#%% Lets write some fitting equations into the 
#=========================== The below functions are created by fitting a model ===========================#
"""
    function IRfit(intensity::AbstractArray{T}, response::AbstractArray{T};
        lb::AbstractArray{T} = [1.0, 1.0, 0.1], #Default rmin = 100, kmin = 0.1, nmin = 0.1 
        p0::AbstractArray{T} = [500.0, 1000.0, 2.0], #Default r = 500.0, k = 200.0, n = 2.0
        ub::AbstractArray{T} = [Inf, Inf, 10.0], #Default rmax = 2400, kmax = 800
    ) where {T<:Real}

This function takes X and Y values and fits them according to a HILL type fit
"""
function HILLfit(intensity::AbstractArray{T}, response::AbstractArray{T};
    lb::AbstractArray{T} = [1.0, 1.0, 0.1], #Default rmin = 100, kmin = 0.1, nmin = 0.1 
    p0::AbstractArray{T} = [500.0, 1000.0, 2.0], #Default r = 500.0, k = 200.0, n = 2.0
    ub::AbstractArray{T} = [Inf, Inf, 10.0], #Default rmax = 2400, kmax = 800
) where {T<:Real}
    fit = curve_fit(HILL_MODEL, intensity, response, p0, lower=lb, upper=ub)
    #Calculate the R-squared
    ss_resid = sum(fit.resid.^2)
    ss_total = sum((response .- mean(response)).^2)
    RSQ = 1 - ss_resid/ss_total
    return fit, RSQ
end

function STFfit(a_wave::AbstractArray{T}, b_wave::AbstractArray{T};
    lb::AbstractArray{T} = [0.001, 0.001, 0.1],
    p0::AbstractArray{T} = [1.0, 10.0, 2.0],
    ub::AbstractArray{T} = [Inf, Inf, 10.0],
) where {T <: Real}
    fit = curve_fit(HILL_MODEL, a_wave, b_wave, p0, lower=lb, upper=ub) #for some reason this works the best when b is in log units
    #Calculate the r squared
    ss_resid = sum(fit.resid.^2)
    ss_total = sum((b_wave .- mean(b_wave)).^2)
    RSQ = 1 - ss_resid/ss_total
    return fit, RSQ
end

"""
The recovery time constant is calculated by fitting the normalized Rdim with the response recovery equation
"""
function TAUfit(data::Experiment{T};
    τRec::T=1.0
) where {T<:Real}
    #Make sure the sizes are the same
    #@assert size(resp) == (size(data, 1), size(data,3))

    trec = zeros(T, size(data, 1), size(data, 3))
    gofs = zeros(T, size(data, 1), size(data, 3))
    #This function uses the recovery model and takes t as a independent variable
    model(x, p) = map(t -> REC(t, -1.0, p[2]), x)
    for swp in axes(data, 1), ch in axes(data, 3)
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

function AMPfit(data::Experiment{T}, resp::Union{T,Matrix{T}}; #This argument should be offloaded to a single value 
    time_cutoff=0.1,
    lb::Vector{T}=(0.0, 0.001),
    p0::Vector{T}=(200.0, 0.002),
    ub::Vector{T}=(Inf, 0.040)
) where {T<:Real}

    #@assert size(resp) == (size(data, 1), size(data,3))

    amp = zeros(2, size(data, 1), size(data, 3))
    gofs = zeros(T, size(data, 1), size(data, 3))

    for swp = axes(data, 1), ch = axes(data, 3)
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