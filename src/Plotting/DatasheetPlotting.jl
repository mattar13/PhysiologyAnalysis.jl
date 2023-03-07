function plot_ir_scatter(axis::PyObject, qData::DataFrame;
     x_row = :Photons, y_row = :Response, 
     color = :black, ms = 10.0, 
     xscale = "log", xbase = 10,
     yscale = "linear", ybase = 10, 
     label = nothing  
)
     X_VALS = qData[:, x_row]
     Y_VALS = qData[:, y_row]
     if isnothing(label)
          axis.scatter(X_VALS, Y_VALS, color = color, s = ms)
     else
          axis.scatter(X_VALS, Y_VALS, color = color, s = ms, label = label)
     end
     if xscale != "linear"
          axis.set_xscale(xscale, base = xbase)
     end
     
     if yscale != "linear"
          axis.set_yscale(yscale, base = ybase)
     end
end

function plot_ir_fit(axis::PyObject, fit_param, rsq::Real;
     xscale = "log", xbase = 10,
     yscale = "linear", ybase = 10,
     color_by_error = true, cmap = :RdYlGn, color = :black,
     fit_rng = 10 .^ range(-1, stop = 6, length = 1000),
     model = HILL_MODEL,
     label_rsq = true,
     lw = 2.0  
)
     xmin, xmax, ymin, ymax = plt.axis() #We should plot everything based on the 
     BEST_FIT = model(fit_rng, fit_param)

     if color_by_error
          colormap =  plt.get_cmap(cmap)
          axis.plot(fit_rng, BEST_FIT, color = colormap(rsq), linewidth = lw*2, alpha = 0.5) #This is the error bar on the outside
     end

     if label_rsq
          axis.plot(fit_rng, BEST_FIT, color = color, lw = lw, label = "R² = $(round(rsq, digits = 2))")
     else
          axis.plot(fit_rng, BEST_FIT, color = color, lw = lw)
     end
     axis.vlines(fit_param[2], ymin=ymin, ymax = fit_param[1]/2, linestyle=(0, (5, 3)), color = color, lw = lw)#10 ^ (p14RodsNR_STF.param[1] / 2), color = color, )
     axis.hlines(fit_param[1]/2, xmin=xmin, xmax = fit_param[2], linestyle=(0, (5, 3)), color = color, lw = lw)
     axis.set_xlim((xmin, xmax))
     axis.set_ylim((ymin, ymax))
     
     if xscale != "linear"
          axis.set_xscale(xscale, base = xbase)
     end

     if yscale != "linear"
          axis.set_yscale(yscale, base = ybase)
     end
end

function plot_IR(axis::PyObject, qData::DataFrame;
     x_row = :Photons, y_row = :Response, 
     xscale = "log", xbase = 10,
     yscale = "linear", ybase = 10,
     color_by_error = true, cmap = :RdYlGn, color = :black,
     fit_rng = 10 .^ range(-1, stop = 6, length = 1000),
     lw = 2.0, ms = 10.0, plot_fits = true,   
     label = nothing, label_rsq = true,
     fitting_kwargs...
)

     plot_ir_scatter(axis, qData;
          x_row = x_row, y_row = y_row, 
          color = color, ms = ms, 
          xscale = xscale, xbase = xbase,
          yscale = yscale, ybase = ybase,
          label = label 
     )
     if plot_fits
          X_VALS = qData[:, x_row]
          Y_VALS = qData[:, y_row]
          fit, rsq = HILLfit(X_VALS, Y_VALS; fitting_kwargs...)
          plot_ir_fit(axis, fit.param, rsq;
               xscale = xscale, xbase = xbase,
               yscale = yscale, ybase = ybase,
               color_by_error = color_by_error, cmap = cmap, color = color,
               fit_rng = fit_rng,
               model = HILL_MODEL,
               lw = lw, label_rsq = label_rsq
          )
     end
end


function plot_exp_fits(axis, df_EXPs::DataFrame, df_TRACEs::DataFrame;
     xlims = (-10.0, 200.0), ylims = (1.0, 3000.0)
)
     #axis.set_title("P14 Rods NR")
     axis.set_xlim(xlims)
     axis.set_ylim(ylims)
     cmap = plt.get_cmap(:RdYlGn)
     for litter in eachrow(df_EXPs)
          RET = df_TRACEs |> @filter((_.Year, _.Month, _.Date, _.Number, _.Channel) == (litter.Year, litter.Month, litter.Date, litter.Number, litter.Channel)) |> DataFrame
          #println(litter)
          axis.scatter(
               RET.Response, RET.Total, 
               alpha = 0.5, 
               color = cmap(litter.RSQ), 
               marker = "o", label = "$(litter.Year)_$(litter.Month)_$(litter.Date)_$(litter.Number)"
          )
          # Find the STF data
          #STF_RES = p14RodsDR_EXP |> @filter((_.Year, _.Month, _.Date, _.Number, _.Channel) == (litter.Year, litter.Month, litter.Date, litter.Number, litter.Channel)) |> DataFrame
          rmax = litter.RMAX; k = litter.K; n = litter.N

          FIT = model(A_rng, [rmax, k, n])
          axis.plot(A_rng, FIT, color = cmap(litter.RSQ))
          axis.vlines(k, ymin=1.0, ymax = rmax/2, linestyle=(0, (5, 3)))#10 ^ (p14RodsNR_STF.param[1] / 2), color = color, )
          axis.hlines(rmax/2, xmin=-10.0, xmax = k, linestyle=(0, (5, 3)))
     end
     axis.legend(loc = "lower right")
     #axis.set_xscale("Log")
     axis.set_xlabel("a-wave (μV)", fontsize=10.0)
     axis.set_yscale("Log", base = 10)
     axis.set_ylabel("b-wave (μV)", fontsize=10.0)
     #return fig
end

function plot_exp_fits(df_EXPs::DataFrame, df_TRACEs::DataFrame; kwargs...)
     fig, axis = plt.subplots(1)
     plot_exp_fits(axis, df_EXPs, df_TRACEs; kwargs...)
end

function plot_data_summary(qTRACE::DataFrame, qEXP::DataFrame; xlims = (-0.25, 2.0))
     channels = qEXP |> @unique(_.Channel) |> @map({_.Channel}) |> DataFrame
     conditions = qEXP |> @unique(_.Condition) |> @map({_.Condition}) |> DataFrame
     println(channels)
     fig_summary, ax_summary = plt.subplots(size(channels, 1), 4)
     ylims = (minimum(qTRACE.Response), maximum(qTRACE.Response))
     photon_lims = (10^-1, 10^4)

     if size(channels.Channel, 1) == 1
          fig_summary.subplots_adjust(wspace = 0.25, bottom=0.18, top=0.85, left = 0.1, right = 0.95)
          fig_summary.set_size_inches(10.0, 3.0)
          ch = channels.Channel[1]
          qTRACE_CH = qTRACE |> @filter(_.Channel == ch) |> DataFrame
          data_unsub = readABF(qTRACE_CH.Path)
          data_filter!(data_unsub)
          plot_experiment(ax_summary[1], data_unsub)
          ax_summary[1].set_xlim(xlims)
          if any(qTRACE_CH.SubPath .!= "NONE")
               data_sub = readABF(qTRACE_CH.SubPath)
               data_filter!(data_sub)
               plot_experiment(ax_summary[1], data_sub, color = :Red)

               data_subtracted = data_unsub - data_sub
               plot_experiment(ax_summary[2], data_subtracted, yaxes = false)
               ax_summary[2].set_xlim(xlims)
          else
               plot_experiment(ax_summary[2], data_unsub, yaxes = false)
               ax_summary[2].set_xlim(xlims)
          end
          fits = qEXP |> @filter(_.Channel == ch) |> DataFrame
          plot_ir_scatter(ax_summary[3], qTRACE_CH)
          plot_ir_fit(ax_summary[3], 
               (fits.RMAX_fit[1], fits.K_fit[1], fits.N_fit[1]),
               fits.RSQ_fit[1]
          )
          ax_summary[3].set_xlim(photon_lims)
          ax_summary[3].set_ylim(-1.0, ylims[2]*1.08)
          ax_summary[3].legend(loc = "upper right", bbox_to_anchor = (1.0, 1.2))

          qTRACE_CH.Percent_Recovery *= 1000
          rdim_traces_A = qTRACE_CH |> 
               @filter(_.Response .> qEXP.rdim[1]) |> 
          DataFrame			
          plot_IR(ax_summary[4], rdim_traces_A, y_row = :Peak_Time, 
               plot_fits = false, color = :red, label = "Peak Time"
          )
          
          tDOM_traces_A = qTRACE_CH |> 
               @filter(_.Percent_Recovery .> 0.0) |> 
          DataFrame
          plot_IR(ax_summary[4], tDOM_traces_A, y_row = :Percent_Recovery, 
               plot_fits = false, label = "Recovery"
          )
          ax_summary[4].legend(loc = "upper right", bbox_to_anchor = (1.0, 1.2))
     else
          fig_summary.subplots_adjust(wspace = 0.25, hspace = 0.1, bottom=0.15, top=0.90, left = 0.1, right = 0.95)
          fig_summary.set_size_inches(10.0, 5.0)
          for (idx, ch) in enumerate(channels.Channel)
               #println(ch)
               qTRACE_CH = qTRACE |> @filter(_.Channel == ch) |> DataFrame
               data_unsub = readABF(qTRACE_CH.Path)
               data_filter!(data_unsub)
               plot_experiment(ax_summary[idx, 1], data_unsub)
               ax_summary[idx, 1].set_xlim(xlims)
               if any(qTRACE_CH.SubPath .!= "NONE")
                    data_sub = readABF(qTRACE_CH.SubPath)
                    data_filter!(data_sub)
                    if idx == 1
                         plot_experiment(ax_summary[idx, 1], data_sub, color = :Red, xaxes = false)
                    else
                         plot_experiment(ax_summary[idx, 1], data_sub, color = :Red)
                    end

                    data_subtracted = data_unsub - data_sub
                    if idx == 1
                         plot_experiment(ax_summary[idx, 2], data_subtracted, axes = false)
                    else
                         plot_experiment(ax_summary[idx, 2], data_subtracted, yaxes = false)
                    end
                    ax_summary[idx, 2].set_xlim(xlims)
               else
                    if idx == 1
                         plot_experiment(ax_summary[idx, 2], data_unsub, axes = false)
                    else
                         plot_experiment(ax_summary[idx, 2], data_unsub, yaxes = false)
                    end
                    ax_summary[idx, 2].set_xlim(xlims)
               end
               fits = qEXP |> @filter(_.Channel == ch) |> DataFrame
               plot_ir_scatter(ax_summary[idx, 3], qTRACE_CH)
               plot_ir_fit(ax_summary[idx, 3], 
                    (fits.RMAX_fit[1], fits.K_fit[1], fits.N_fit[1]),
                    fits.RSQ_fit[1]
               )
               ax_summary[idx, 3].set_xlim(photon_lims)
               ax_summary[idx, 3].set_ylim(-1.0, ylims[2]*1.08)

               qTRACE_CH.Percent_Recovery *= 1000
               rdim_traces_A = qTRACE_CH |> 
                    @filter(_.Response .> qEXP.rdim[idx]) |> 
               DataFrame			
               plot_IR(ax_summary[idx, 4], rdim_traces_A, y_row = :Peak_Time, 
                    plot_fits = false, color = :red, label = "Peak Time"
               )

               tDOM_traces_A = qTRACE_CH |> 
                    @filter(_.Percent_Recovery .> 0.0) |> 
               DataFrame
               
               plot_IR(ax_summary[idx, 4], tDOM_traces_A, y_row = :Percent_Recovery, 
                    plot_fits = false, label = "Recovery"
               )
          end
          ax_summary[1, 3].xaxis.set_visible(false) #We want the spine to fully
          ax_summary[1, 3].legend(loc = "upper right", bbox_to_anchor = (1.0, 1.2))
          ax_summary[1, 4].xaxis.set_visible(false) #We want the spine to fully
          ax_summary[1, 4].legend(loc = "upper right", bbox_to_anchor = (1.0, 1.2))
     end
     #println(qTRACE.SubPath)
     return fig_summary
end