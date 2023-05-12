"""
This function is for computing the R-squared of a polynomial
"""
function RSQ(poly::PN.Polynomial, x, y)
     ŷ = poly.(x)
     ȳ = sum(ŷ) / length(ŷ)
     SSE = sum((y - ŷ) .^ 2)
     SST = sum((y .- ȳ) .^ 2)
     1 - SSE / SST
end

function RSQ(ŷ::Array{T}, y::Array{T}) where {T<:Real}
     ȳ = sum(ŷ) / length(ŷ)
     SSE = sum((y - ŷ) .^ 2)
     SST = sum((y .- ȳ) .^ 2)
     1 - SSE / SST
end

function sig_symbol(val)
     if val <= 0.001
          return "***"
     elseif val <= 0.005
          return "**"
     elseif val <= 0.05
          return "*"
     else
          return "-"
     end
end

function dataset_statistics(qEXP; 
     control = "WT", 
     stat_metrics = [:rmax, :rdim, :K_fit, :time_to_peak, :percent_recovery, :integration_time],
)
     #unflagged_exps = qEXP |> @filter(_.INCLUDE == true) |> DataFrame
     exps = DataFrame[]
     for stat in stat_metrics
          res_stat = unflagged_exps  |> 
               @groupby({_.Genotype, _.Age, _.Condition, _.Photoreceptor}) |> 
               @map({Genotype = key(_)[1], Age = key(_)[2], Condition = key(_)[3], Photoreceptor = key(_)[4],
                    N = length(_),
                    METRIC = string(stat),
                    AVG = 0.0, STD = 0.0, SEM = 0.0, CI = 0.0, 
                    LOWER = 0.0, UPPER = 0.0,
                    P = 0.0, SIGN = "-"
               })|> 
          DataFrame
          for (idx, info) in enumerate(eachrow(res_stat))
               #first load the control data
               ctrl_data = unflagged_exps |> 
                    @filter(_.Genotype == control && _.Age == info.Age && _.Condition == info.Condition && _.Photoreceptor == info.Photoreceptor) |> 
               DataFrame
               exp_data =  unflagged_exps |> 
                    @filter(_.Genotype == info.Genotype && _.Age == info.Age && _.Condition == info.Condition && _.Photoreceptor == info.Photoreceptor) |> 
               DataFrame
               res_stat[idx, :AVG] = mean(exp_data[:, stat])
               res_stat[idx, :STD] = std(exp_data[:, stat])
               res_stat[idx, :SEM] = sem(exp_data[:, stat])
               res_stat[idx, :CI] = CI = 1.96*sem(exp_data[:, stat])
               res_stat[idx, :LOWER] = mean(exp_data[:, stat]) - CI
               res_stat[idx, :UPPER] = mean(exp_data[:, stat]) + CI
               if size(exp_data,1) > 1 && size(ctrl_data,1) > 1 && sum(ctrl_data[:, stat]) != sum(exp_data[:, stat])
                    res_stat[idx, :P] = P = UnequalVarianceTTest(ctrl_data[:, stat], exp_data[:, stat]) |> pvalue
                    res_stat[idx, :SIGN] = sig_symbol(P)
               else
                    res_stat[idx, :P] = 1.0
                    res_stat[idx, :SIGN] = "-"
               end
          end
          push!(exps, res_stat)
     end
    
     stats = vcat(exps...)
     
     return stats
end

dataset_statistics(dataset::Dict{String, DataFrame}; kwargs...) = dataset_statistics(dataset["EXPERIMENTS"]; kwargs...)