#p0 = [A, τ_ON, τ_OFF, shift]
function single_stim_model(t, p)
    y_fit = p[1] * (1-exp(-(t-p[4])/ p[3]))*exp(-(t-p[4])/p[2])
    if y_fit < 0
        return 0.0
    else
        return y_fit
    end
end