"""
A simple stimulus
"""
function stimulus_callback(step_begin, step_end, level)
     condition_fn1(u, t, integrator) = t == step_begin
     affect_fn1!(integrator) = integrator.u[1] = level
     cb_fn1 = DiscreteCallback(condition_fn1, affect_fn1!)

     condition_fn2(u, t, integrator) = t == step_end
     affect_fn2!(integrator) = integrator.u[1] = 0.0
     cb_fn2 = DiscreteCallback(condition_fn2, affect_fn2!)
     return CallbackSet(cb_fn1, cb_fn2)
end

#u = [ui, uc] ui = stimulus, uc = artifact response
#P = [A, C] A = Amplitude, C = Capacitance 
RC(du, u, p, t) = du[2] = ((u[1] - u[2]) / (p[2]))

function RCArtifact(data::Experiment, p0)
     u0 = [0.0, 0.0]
     tstops = data.stim_protocol[1].timestamps
     tspan = (data.t[1], data.t[end])
     stim = stimulus_callback(tstops[1], tstops[2], p0[1])
     prob = ODEProblem(RC, u0, tspan, p0)
     sol = solve(prob, callback=stim, tstops=tstops)
     artifact = sol(data.t, idxs=1) - sol(data.t, idxs=2)
     #The second phase will be fitting the artifact to the data
     return artifact
end

function fitArtifact(data::Experiment; p0 = [150.0, 1.0e-3]; swp = 1, ch = 2)
     min_func(x) = MeanSquaredError(RCArtifact(data, x), data.data_array[swp,:,ch])
     results = optimize(min_func, p0)
     pOPT = Optim.minimizer(results)
     return pOPT
end


function removeArtifact(exp::Experiment)
     data = deepcopy(exp)
     for swp in 1:size(data, 1), ch in 1:size(data,3)
          pOPT = ePhys.fitArtifact(data)
          artifact = RCArtifact(data_unfilt, pOPT)
          artifact_removed = artifact .- data_unfilt.data_array[swp, :, ch]
     end
     return data
end