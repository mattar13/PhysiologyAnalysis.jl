#p0 = [A, τ_ON, τ_OFF, shift]
function single_stim_model(t::T, p) where T<:Real
    y_fit = p[1] * (1-exp(-(t-p[4])/ p[3]))*exp(-(t-p[4])/p[2])
    if y_fit < 0
        return 0.0
    else
        return y_fit
    end
end

single_stim_model(t_series::Vector{T}, p) where T<:Real = map(TIME-> single_stim_model(TIME, p), t_series)

function fit_parametric(y_data, t_data; 
    p0 = [0.50, 100.0, 10.0, 100.0], ub = [1.0, 100.0, 100.0, 1000.0], lb = [0.0, 0.0, 0.0, 0.0]
)
    model(t_series, p) = map(TIME-> single_stim_model(TIME, p), t_series)
    fit = curve_fit(model, t_data, y_data, p0, upper=ub, lower = lb)
    return fit
end


function iterative_fit(exp::Experiment)


end