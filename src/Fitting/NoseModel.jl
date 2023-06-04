function runNoseModel(tspan ; r1 = 1.0, r2 = 0.1, k1 = 1.5, k2 = -1.0)
     p = @parameters r1 = r1 r2 = r2 k1 = k1 k2 = k2 
     u0 = @variables N1(t) = 0.0 N2(t) = 0.0 X = 0.0
     D = Differential(t)
     eqs = [
          D(N1) ~ r1*(k1 - N1)
          D(N2) ~ r2*(k2-k1 - N2)
          D(X) ~ (N1 + N2) - X
     ]
     @named sys = ODESystem(eqs)
     sys = structural_simplify(sys)
     prob = ODEProblem(sys, u0, tspan, p)
     sol = solve(prob, Tsit5())
     Array(sol)'
end