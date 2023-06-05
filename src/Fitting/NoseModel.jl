function DoubleLogistic(tsteps;
     R1 = 100.0, K1 = -300, 
     R2 = 5.0, K2 = -200.0
)
     @parameters r1 r2 k1 k2
     @variables t N1(t) N2(t) X(t)
     D = Differential(t)
     
     eqs = [D(N1) ~ r1*(k1 - N1), D(N2) ~ r2*(k2-k1 - N2)]
     
     @named sys = ODESystem(eqs)
     sys = structural_simplify(sys)
     
     u0 = [N1 => 0.0,N2 => 0.0,X => 0.0]
     
     p = [r1=> R1,r2=> R2,k1=> K1,k2=> K2,]
     
     prob = ODEProblem(sys, u0, (tsteps[1], tsteps[end]), p)
     sol_arr = solve(prob, Tsit5(), saveat = tsteps) |> Array
     return sol_arr[1,:]# .+ sol_arr[2,:]
end

function SingleLogistic(tsteps;
     R1 = 100.0, K1 = 30, 
)
     @parameters r1 k1
     @variables t N1(t)
     D = Differential(t)
     
     eqs = [D(N1) ~ r1*(k1 - N1)]
     
     @named sys = ODESystem(eqs)
     sys = structural_simplify(sys)
     
     u0 = [N1 => 0.0]
     
     p = [r1=> R1,k1=> K1]
     
     prob = ODEProblem(sys, u0, (tsteps[1], tsteps[end]), p)
     sol_arr = solve(prob, Tsit5(), saveat = tsteps) |> Array
     return sol_arr[1,:]# .+ sol_arr[2,:]
end

NOSEMODEL1(x, p; init_y = -300.0) = init_y .+ SingleLogistic(x; R1 = p[1], K1 = p[2])
NOSEMODEL2(x, p) = DoubleLogistic(x; R1 = p[1], K1 = p[2], R2 = p[3], K2 = p[4])

function NOSEfit(time, response;     
     p0 = [10.0, 30.0],
)
     MODEL(x, p) = NOSEMODEL1(x, p; init_y = minimum(response))
     fit = curve_fit(NOSEMODEL1, time, response, p0)
     ss_resid = sum(fit.resid.^2)
     ss_total = sum((response .- mean(response)).^2)
     RSQ = 1 - ss_resid/ss_total
     return fit, RSQ
end

function findNosePeak(data::Experiment; tmax = 0.75, kwargs...)
     t = data.t
     resp = zeros(size(data, 3))
     for ch in axes(data,3)
          ŷ = minimum(data, dims =1)[1,:,ch]
          if t[argmin(ŷ)] < tmax-data.dt
               tidxs = findall(t[argmin(ŷ)].< t .< tmax)
               fit, RSQ = NOSEfit(t[tidxs], ŷ[tidxs]; kwargs...)
               y = NOSEMODEL1(t[tidxs], fit.param; init_y = minimum(ŷ));
               if minimum(ŷ) < y[end] < 0.0 
                    resp[ch] = y[end]
               else
                    resp[ch] = minimum(ŷ)
               end
          else
               resp[ch] = minimum(ŷ)
          end
     end
     return resp
end