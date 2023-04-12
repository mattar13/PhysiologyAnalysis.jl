function plot_ir_scatter(axis::T, qData::DataFrame;
     x_row = :Photons, y_row = :Response, 
     color = :black, ms = 10.0, 
     xscale = "log", xbase = 10,
     yscale = "linear", ybase = 10, 
     label = nothing, kwargs...  
) where T
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

function plot_ir_fit(axis::T, fit_param, rsq::Real;
     xscale = "log", xbase = 10,
     yscale = "linear", ybase = 10,
     color_by_error = true, cmap = :RdYlGn, color = :black,
     fit_rng = 10 .^ range(-1, stop = 6, length = 1000),
     model = HILL_MODEL,
     label_rsq = true,
     lw = 2.0  
) where T
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
     axis.set_xlim((10^-1, 10^6))
     axis.set_ylim((ymin, ymax))
     
     if xscale != "linear"
          axis.set_xscale(xscale, base = xbase)
     end

     if yscale != "linear"
          axis.set_yscale(yscale, base = ybase)
     end
end

function plot_IR(axis::T, qData::DataFrame;
     x_row = :Photons, y_row = :Response, 
     xscale = "log", xbase = 10,
     yscale = "linear", ybase = 10,
     color_by_error = true, cmap = :RdYlGn, color = :black,
     fit_rng = 10 .^ range(-1, stop = 6, length = 1000),
     lw = 2.0, ms = 10.0, plot_fits = true,   
     label = nothing, label_rsq = true,
     fitting_kwargs...
) where T

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

function plot_IR(qData::DataFrame)
     fig, ax = plt.subplots(1)
     plot_IR(ax, qData)
     return fig
end

function plot_dataset_fits(
     qTRACES::DataFrame, qEXPS::DataFrame, qCONDS::DataFrame;
     xlims = (10^-1, 10^5), normalize = false, 
     condition = "BaCl_LAP4", photoreceptor = "Rods",
     xcats = ["P9", "P11", "P13", "Adult"], xcol = :Age,
     ycats = ["WT", "RS1KO", "R141C", "C59S"], ycol = :Genotype,
     colorcycle = [:Black, :Purple, :Orange, :Green]
)
     #Plot the fits of the experiment
     fig_fits, ax_fits = plt.subplots(length(xcats)+1, length(ycats), figsize = (10.0, 10.0))
     fig_fits.subplots_adjust(wspace = 0.6, hspace = 0.1, bottom = 0.1, top = 0.9)
     #println(categories)
     for (idxX, XCAT) in enumerate(xcats), (idxY, YCAT) in enumerate(ycats)
          #filter out experiments
          qEXP_category = qEXPS |> 
               @filter(_[xcol] == XCAT && _[ycol] == YCAT && _.Condition == condition && _.Photoreceptor == photoreceptor && _.INCLUDE) |> 
               @orderby(_.RSQ_fit) |> 
          DataFrame
          ylims = (-0.2, maximum(qEXP_category.rmax) * 1.1)
          for litter in eachrow(qEXP_category)
               #println(litter)
               RET = matchExperiment(qTRACES, litter)
               if normalize
                    RET.Response ./= maximum(RET.Response) #Normalize
                    fit_param = (1.0, litter.K_fit, litter.N_fit)
               else
                    fit_param = (litter.RMAX_fit, litter.K_fit, litter.N_fit)
               end
               plot_ir_fit(ax_fits[idxY, idxX], fit_param, litter.RSQ_fit, label_rsq = false)
               ax_fits[idxY, idxX].scatter(RET.Photons, RET.Response,
                    s = 6.0, 
                    label = "$(litter.Year)_$(litter.Month)_$(litter.Date)_$(litter.Number)_$(round(litter.RSQ_fit, digits = 2))"
               )
          end
          ax_fits[idxY, idxX].legend(loc = "upper right", bbox_to_anchor = (1.4, 1.0), fontsize = 4)
          ax_fits[idxY, idxX].set_xlim(xlims)
          ax_fits[idxY, idxX].set_ylim(ylims)
          if idxX == 1
               ax_fits[idxY, idxX].set_ylabel("$(YCAT) \n Resp. (μV)")
          end

          #if idxY == size(ycats, 1)
          #     ax_fits[idxY, idxX].set_xlabel("Intensity (μV)")
          if idxY == 1
               ax_fits[idxY, idxX].set_title(XCAT)
               ax_fits[idxY, idxX].xaxis.set_visible(false) #We want the spine to fully
          else
               ax_fits[idxY, idxX].xaxis.set_visible(false) #We want the spine to fully
          end
          qCONDS_category = qCONDS |> @filter(_[xcol] == XCAT && _[ycol] == YCAT && _.Condition == condition && _.Photoreceptor == photoreceptor) |> DataFrame
          coll_fit_params = (qCONDS_category.RMAX_COLL[1], qCONDS_category.K_COLL[1], qCONDS_category.N_COLL[1])
          coll_RSQ = qCONDS_category.RSQ_COLL[1]
          qAGE = qCONDS |> @filter(_[xcol] == XCAT && _.Condition == condition && _.Photoreceptor == photoreceptor) |> DataFrame
          ymax = maximum(qAGE.RMAX_COLL) * 1.2
          plot_ir_fit(ax_fits[5, idxX], coll_fit_params, coll_RSQ, color_by_error = false, color = colorcycle[idxY])
          ax_fits[5, idxX].set_xlabel("Intensity (μV)")
          ax_fits[5, idxX].set_xlim(xlims)
          ax_fits[5, idxX].set_ylim(-0.2, ymax)
     end
     ax_fits[5, 1].set_ylabel("Response (μV)")
     return fig_fits
end

plot_dataset_fits(dataset::Dict{String, DataFrame}; kwargs...) = plot_dataset_fits(dataset["TRACES"], dataset["EXPERIMENTS"], dataset["CONDITIONS"]; kwargs...)

function plot_dataset_vals(qTRACES::DataFrame, qEXPS::DataFrame, qCONDS::DataFrame;
     #xlims = (10^-1, 10^5), normalize = false, 
     metricX = :Rmax, metricXSEM = :Rmax_sem, photoreceptor = "Rods",
     xcats =  ["P9", "P11", "P13", "Adult"],  xcol = :Age,
     ycats = ["BaCl", "NoDrugs", "BaCl_LAP4"], ycol = :Condition,
     zcats = ["WT", "RS1KO", "R141C", "C59S"], zcol = :Genotype,
     colorcycle = [:Black, :Purple, :Orange, :Green]
)
     #Plot the fits of the experiment
     fig_vals, ax_vals = plt.subplots(length(ycats), 1, figsize = (7.5, 5.0))
     fig_vals.subplots_adjust(wspace = 0.6, hspace = 0.1, bottom = 0.1, top = 0.9)
     #println(categories)
     ax_vals[1].set_title("$(metricX)")
     for (idxY, YCAT) in enumerate(ycats), (idxZ, ZCAT) in enumerate(zcats) 
          qCOND = qCONDS |> 
               @filter(_[zcol] == ZCAT && _[ycol] == YCAT && _.Photoreceptor == photoreceptor) |> 
          DataFrame

          permute_idxs = indexin(xcats, qCOND[:, xcol])
          qCOND = qCOND[permute_idxs, :]
          ax_vals[idxY].errorbar(1:length(xcats), qCOND[:, metricX], c=colorcycle[idxZ], yerr=qCOND[:, metricXSEM], marker="o", ms=5.0, lw=3.0)
          ax_vals[idxY].set_ylabel("$(metricX) \n $(YCAT)")
          if idxY != length(ycats)
               ax_vals[idxY].xaxis.set_visible(false) #We want the spine to fully
          else
               ax_vals[idxY].set_xticks(1:length(xcats), xcats, fontsize=11.0)  
          end
     end
     return fig_vals
end

plot_dataset_vals(dataset::Dict{String, DataFrame}; kwargs...) = plot_dataset_vals(dataset["TRACES"], dataset["EXPERIMENTS"], dataset["CONDITIONS"]; kwargs...)


function plot_data_summary(qTRACE::DataFrame, qEXP::DataFrame; 
     xlims = (-0.25, 2.0)
)
     channels = qEXP |> @unique(_.Channel) |> @map({_.Channel}) |> DataFrame
     conditions = qEXP |> @unique(_.Condition) |> @map({_.Condition}) |> DataFrame
     #println(channels)
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
               if idx == 1
                    plot_experiment(ax_summary[idx, 1], data_unsub, xaxes = false)
               else
                    plot_experiment(ax_summary[idx, 1], data_unsub, xaxes = false)
               end
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
               ax_summary[idx, 4].set_xlim(photon_lims)

               ax_summary[idx, 3].legend(loc = "upper right", bbox_to_anchor = (1.0, 1.2))
               ax_summary[idx, 4].legend(loc = "upper right", bbox_to_anchor = (1.0, 1.2))
          end
          ax_summary[1, 3].xaxis.set_visible(false) #We want the spine to fully
          ax_summary[1, 4].xaxis.set_visible(false) #We want the spine to fully

     end
     #println(qTRACE.SubPath)
     return fig_summary
end

plot_data_summary(dataset::Dict{String, DataFrame}; kwargs...) = plot_data_summary(dataset["TRACES"], dataset["EXPERIMENTS"]; kwargs...)