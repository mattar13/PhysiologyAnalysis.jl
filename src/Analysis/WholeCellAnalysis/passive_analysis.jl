function calculate_baseline(exp::Experiment{E, T}; channel = 1) where {E, T <: Real}
     baselines = Vector{Float64}()
     for (i, trial) in enumerate(eachtrial(exp))
          epoch1_idx = exp.HeaderDict["EpochTableByChannel"][channel].epochWaveformBytrial[i].p2s[2] #This monstrosity is the first point of the
          avg_base = mean(trial[1,1:epoch1_idx,2])
          push!(baselines, avg_base)
     end
     baselines
end

function calculate_peak(exp::Experiment{E, T}; channel = 1, digital_cmd = "Cmd 0") where {E, T<:Real}
     if size(exp,3) != 3
          #We need to add a new channel
          create_signal_waveform!(exp, digital_cmd)
     end
     cmd_v = exp[:, findfirst(exp[1,:,3] .!= 0.0), 3]
     I_CM = zeros(cmd_v |> size)
     for (i, tr) in enumerate(eachtrial(exp))
         if cmd_v[i]-10 < 0.0 #I don't know why we have to subtract 10. Maybe this isn't every time
             I_CM[i] = minimum(tr, dims = 2)[1, 1, channel]
         else
             I_CM[i] = maximum(tr, dims = 2)[1, 1, channel]
         end
     end
     I_CM
end

extract_timepoint(exp; timepoint = 0.5, channel = 1) = exp[:, findfirst(exp.t .>= timepoint), channel]

function calculate_resistance(exp)
     V_HOLD = extract_timepoint(exp; channel = 2)
     I_CM = calculate_peak(exp)
     I_RIN = extract_timepoint(exp)
     @. lin_model(x, p) = p[1]*x+p[2]
     Rs_fit = curve_fit(lin_model, I_CM, V_HOLD,  [1.0, 0.0])
     Rin_fit = curve_fit(lin_model, I_RIN[1:4], V_HOLD[1:4],  [1.0, 0.0]) #The higher 
     Rs = Rs_fit.param[1]*1e3
     Rin = Rin_fit.param[1]*1e3
     V_M = Rin_fit.param[2]
     V_HOLD = Rs_fit.param[2]
     return Rs, Rin, V_M, V_HOLD
end

function calculate_capacitance(exp; channel = 1)
     V_HOLD = extract_timepoint(exp; channel = 2)
     I_RIN = extract_timepoint(exp) #Find the stable state current
     epoch_idx1 = exp.HeaderDict["EpochTableByChannel"][channel].epochWaveformBytrial[1].p1s[3] #This monstrosity is the first point of the
     epoch_idx2 = exp.HeaderDict["EpochTableByChannel"][channel].epochWaveformBytrial[1].p2s[3] #This monstrosity is the first point of the
     data_section = exp[:,epoch_idx1:epoch_idx2,1].-I_RIN
     Q = sum((data_section[4,:])) * (exp.dt/1e3) * 1e-9 #This is the charge in Amp*ms
     V = V_HOLD[4,1] * 1e-3 #This is the holding voltage in V
     abs.((Q./V)*1e12) #This is the capacitance in picoFarads
end