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

function dataset_statistics(qEXP; control = "WT")
     res_rmax = qEXP |> @groupby({_.Genotype, _.Age, _.Condition, _.Photoreceptor}) |> 
          @map({Genotype = key(_)[1], Age = key(_)[2], Condition = key(_)[3], Photoreceptor = key(_)[4],
               N = length(_),
               METRIC = "",
               AVG = mean(_.rmax), STD = std(_.rmax), SEM = sem(_.rmax), CI = 1.96*sem(_.rmax), 
               LOWER = mean(_.rmax) - 1.96*sem(_.rmax), UPPER = mean(_.rmax) + 1.96*sem(_.rmax),
               P = 0.0, SIGN = "-"
          })|> @orderby(_.Age) |> @thenby(_.Genotype) |> DataFrame

     res_rdim = qEXP |> @groupby({_.Genotype, _.Age, _.Condition, _.Photoreceptor}) |> 
          @map({Genotype = key(_)[1], Age = key(_)[2], Condition = key(_)[3], Photoreceptor = key(_)[4],
               N = length(_),
               METRIC = "",
               AVG = mean(_.rdim), STD = std(_.rdim), SEM = sem(_.rdim), CI = 1.96*sem(_.rdim), 
               LOWER = mean(_.rdim) - 1.96*sem(_.rdim), UPPER = mean(_.rdim) + 1.96*sem(_.rdim),
               P = 0.0, SIGN = "-"
          })|> @orderby(_.Age) |> @thenby(_.Genotype) |> DataFrame

     res_K_fit = qEXP |> @groupby({_.Genotype, _.Age, _.Condition, _.Photoreceptor}) |> 
          @map({Genotype = key(_)[1], Age = key(_)[2], Condition = key(_)[3], Photoreceptor = key(_)[4],
               N = length(_),
               METRIC = "",
               AVG = mean(_.K_fit), STD = std(_.K_fit), SEM = sem(_.K_fit), CI = 1.96*sem(_.K_fit), 
               LOWER = mean(_.K_fit) - 1.96*sem(_.K_fit), UPPER = mean(_.K_fit) + 1.96*sem(_.K_fit),
               P = 0.0, SIGN = "-"
          })|> @orderby(_.Age) |> @thenby(_.Genotype) |> DataFrame

     res_tint = qEXP |> @groupby({_.Genotype, _.Age, _.Condition, _.Photoreceptor}) |> 
          @map({Genotype = key(_)[1], Age = key(_)[2], Condition = key(_)[3], Photoreceptor = key(_)[4],
               N = length(_),
               METRIC = "",
               AVG = mean(_.integration_time), STD = std(_.integration_time), SEM = sem(_.integration_time), CI = 1.96*sem(_.integration_time), 
               LOWER = mean(_.integration_time) - 1.96*sem(_.integration_time), UPPER = mean(_.integration_time) + 1.96*sem(_.integration_time),
               P = 0.0, SIGN = "-"
          })|> @orderby(_.Age) |> @thenby(_.Genotype) |> DataFrame
     res_tpeak = qEXP |> @groupby({_.Genotype, _.Age, _.Condition, _.Photoreceptor}) |> 
          @map({Genotype = key(_)[1], Age = key(_)[2], Condition = key(_)[3], Photoreceptor = key(_)[4],
               N = length(_),
               METRIC = "",
               AVG = mean(_.time_to_peak), STD = std(_.time_to_peak), SEM = sem(_.time_to_peak), CI = 1.96*sem(_.time_to_peak), 
               LOWER = mean(_.time_to_peak) - 1.96*sem(_.time_to_peak), UPPER = mean(_.time_to_peak) + 1.96*sem(_.time_to_peak),
               P = 0.0, SIGN = "-"
          })|> @orderby(_.Age) |> @thenby(_.Genotype) |> DataFrame
     res_rec = qEXP |> @groupby({_.Genotype, _.Age, _.Condition, _.Photoreceptor}) |> 
          @map({Genotype = key(_)[1], Age = key(_)[2], Condition = key(_)[3], Photoreceptor = key(_)[4],
               N = length(_),
               METRIC = "",
               AVG = mean(_.percent_recovery), STD = std(_.percent_recovery), SEM = sem(_.percent_recovery), CI = 1.96*sem(_.percent_recovery), 
               LOWER = mean(_.percent_recovery) - 1.96*sem(_.percent_recovery), UPPER = mean(_.percent_recovery) + 1.96*sem(_.percent_recovery),
               P = 0.0, SIGN = "-"
          })|> @orderby(_.Age) |> @thenby(_.Genotype) |> DataFrame
     
     for (idx, info) in enumerate(eachrow(res_rmax))
          #first load the control data
          ctrl_data = qEXP |> @filter(_.Age == info.Age && _.Genotype == control) |> DataFrame
          exp_data = qEXP |> @filter(_.Age == info.Age && _.Genotype == info.Genotype) |> DataFrame
          #println(size(ctrl_data))
          #println(size(exp_data))
          res_rmax[idx, :METRIC] = "RMAX"
          res_rdim[idx, :METRIC] = "RDIM"
          res_K_fit[idx, :METRIC] = "K"
          res_tint[idx, :METRIC] = "TINT"
          res_tpeak[idx, :METRIC] = "TPEAK"
          res_rec[idx, :METRIC] = "REC"
          if size(exp_data,1) > 1
               pvalue_rmax = UnequalVarianceTTest(ctrl_data.rmax, exp_data.rmax) |> pvalue
               pvalue_rdim = UnequalVarianceTTest(ctrl_data.rdim, exp_data.rdim) |> pvalue
               pvalue_K_fit = UnequalVarianceTTest(ctrl_data.K_fit, exp_data.K_fit) |> pvalue
               pvalue_integration_time = UnequalVarianceTTest(ctrl_data.integration_time, exp_data.integration_time) |> pvalue
               pvalue_time_to_peak = UnequalVarianceTTest(ctrl_data.time_to_peak, exp_data.time_to_peak) |> pvalue
               pvalue_percent_recovery = UnequalVarianceTTest(ctrl_data.percent_recovery, exp_data.percent_recovery) |> pvalue

               sig_rmax = sig_symbol(pvalue_rmax)
               sig_rdim = sig_symbol(pvalue_rdim)
               sig_K_fit = sig_symbol(pvalue_K_fit)
               sig_integration_time = sig_symbol(pvalue_integration_time)
               sig_time_to_peak = sig_symbol(pvalue_time_to_peak)
               sig_percent_recovery = sig_symbol(pvalue_percent_recovery)

               #println("Test $(idx) : $(info.Age) $(info.Genotype) = $(res_Pval) $(SIG_SYMBOL)")
               res_rmax[idx, :P] = pvalue_rmax
               res_rdim[idx, :P] = pvalue_rdim
               res_K_fit[idx, :P] = pvalue_K_fit
               res_tint[idx, :P] = pvalue_integration_time
               res_tpeak[idx, :P] = pvalue_time_to_peak
               res_rec[idx, :P] = pvalue_percent_recovery

               res_rmax[idx, :SIGN] = sig_rmax
               res_rdim[idx, :SIGN] = sig_rdim
               res_K_fit[idx, :SIGN] = sig_K_fit
               res_tint[idx, :SIGN] = sig_integration_time
               res_tpeak[idx, :SIGN] = sig_time_to_peak
               res_rec[idx, :SIGN] = sig_percent_recovery
          else
               res_rmax[idx, :P] = 1.0
               res_rdim[idx, :P] = 1.0
               res_K_fit[idx, :P] = 1.0
               res_tint[idx, :P] = 1.0
               res_tpeak[idx, :P] = 1.0
               res_rec[idx, :P] = 1.0

               res_rmax[idx, :SIGN] = "-"
               res_rdim[idx, :SIGN] = "-"
               res_K_fit[idx, :SIGN] = "-"
               res_tint[idx, :SIGN] = "-"
               res_tpeak[idx, :SIGN] = "-"
               res_rec[idx, :SIGN] = "-"
          end
     end
     stats = vcat(res_rmax, res_rdim, res_K_fit, res_tint, res_tpeak, res_rec)
     stats = stats |> 
          @orderby(_.P) |> 
          @thenby(_.Condition) |> 
          @thenby(_.Age) |> 
          @thenby(_.Genotype) |> 
     DataFrame
     
     return stats
end