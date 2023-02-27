#This file contains the calibration data

function extract_categories(cond_df::DataFrame) 
     categories = []
     for row in eachrow(cond_df)
          a = [row.Age, row.Genotype, row.Photoreceptor, row.Wavelength]
          push!(categories, a)
     end
     categories
end

function extractIR(trace_datafile::DataFrame, category; measure = :Response, kwargs...)

     allIR = trace_datafile |> 
          @filter(_.Age == category[1]) |>
          @filter(_.Genotype == category[2]) |>
          @filter(_.Photoreceptor == category[3]) |>
          @filter(_.Wavelength == category[4]) |>
     DataFrame

     intensity = allIR.Photons
     response = abs.(allIR[!, measure])
     #generate the fits for the equation
     r = measure == :Minima ? abs(minimum(response)) : maximum(response)
     #println(intensity)
     fit = IRfit(intensity, response; r = r, rmax = (r + 1000), kwargs...)

     model(I, p) = map(i -> p[1]*IR(i, p[2], p[3]), I)
     fitResponse = model(intensity, fit.param)

     allIR[!, :FitVariable] = fitResponse
     allIR[!, :Residual] = fit.resid

     select!(allIR, 
          [
               :Path, :Year, :Month, :Date, 
               :Age, :Number, :Genotype, :Photoreceptor, :Wavelength, :Channel, :Gain,
               :Photons, measure, :FitVariable, :Residual
          ]
     )
     return allIR, fit
end

function summarize_data(qTrace::DataFrame, qExperiment::DataFrame; kwargs...)
     qConditions = qExperiment |>
          @groupby({_.Age, _.Genotype, _.Photoreceptor, _.Wavelength}) |>
          @map({
               Age = _.Age[1], Genotype = _.Genotype[1], Photoreceptor = _.Photoreceptor[1], Wavelength = _.Wavelength[1],
               N = length(_),
               Rmax = mean(_.rmax), Rmax_sem = sem(_.rmax),
               Rdim = mean(_.rdim), Rdim_sem = sem(_.rdim),
               Integrated_Time = mean(_.integration_time), Integrated_Time_sem = sem(_.integration_time),
               Time_to_peak = mean(_.time_to_peak), Time_To_Peak_sem = sem(_.time_to_peak),
               Percent_Recovery = mean(_.percent_recovery), Percent_Recovery_sem = sem(_.percent_recovery), 
               K_IR = mean(_.K_fit), K_SEM_IR = sem(_.K_fit),
               RMAX_COLL = 0.0, K_COLL = 0.0, N_COLL = 0.0, RSQ_COLL = 0.0 
               #Recovery_Tau = mean(_.recovery_tau), Recovery_Tau_sem = sem(_.recovery_tau),
          }) |>
          DataFrame
     #iterate through each conditon and generate a collated IR curve
     for (idx, cond) in enumerate(eachrow(qConditions))
          qIND_COND = qTrace |> 
               @filter(_.Age == cond.Age) |> 
               @filter(_.Genotype == cond.Genotype) |> 
               @filter(_.Photoreceptor == cond.Photoreceptor) |> 
               @filter(_.Wavelength == cond.Wavelength) |> 
          DataFrame 
          if size(qIND_COND, 1) > 2
               fit, rsq = ePhys.IRfit(qIND_COND.Photons, qIND_COND.Response; kwargs...)
               qConditions[idx, :RMAX_COLL] = fit.param[1]
               qConditions[idx, :K_COLL] = fit.param[2]
               qConditions[idx, :N_COLL] = fit.param[3]
               qConditions[idx, :RSQ_COLL] = rsq
          end
     end
     return qConditions
end

function runTraceAnalysis(all_files::DataFrame;
     t_pre = 1.0, t_post = 2.0, #Extend these to see the end of the a-wave
     measure_minima = false,
     a_cond = "BaCl_LAP4", b_cond = "BaCl", g_cond = "NoDrugs", 
     lb = [1.0, 1.0, 0.1], #Default rmin = 100, kmin = 0.1, nmin = 0.1 
     p0 = [500.0, 1000.0, 2.0], #Default r = 500.0, k = 200.0, n = 2.0
     ub = [Inf, Inf, 10.0], #Default rmax = 2400, kmax = 800
     verbose = true
)
     dataset = Dict{String, DataFrame}(
          "ALL_FILES" => all_files
     )
     uniqueData = all_files |> @unique({_.Year, _.Month, _.Date, _.Number, _.Wavelength, _.Photoreceptor, _.Genotype, _.Condition}) |> DataFrame #Pull out all unique files
     if verbose
          println("Completed data query")
     end
     qTrace = DataFrame() #Make empty dataframes for all traces
     qExperiment = DataFrame() #Make empty dataframe for all experiments
     for (idx, i) in enumerate(eachrow(uniqueData)) #Walk through each unique data
          if verbose
               println("Completeing analysis for $idx out of $(size(uniqueData,1))")
               println("Path: $(i.Path)")
          end
          qData = all_files |> @filter(
               (_.Year, _.Month, _.Date, _.Number, _.Wavelength, _.Photoreceptor, _.Genotype) ==
               (i.Year, i.Month, i.Date, i.Number, i.Wavelength, i.Photoreceptor, i.Genotype)
          ) |> DataFrame #Pull out the data 
          #Determine whether or not a subtraction should be done
          if i.Condition == a_cond
               qTRIALa = qTRIAL = qData |> @filter(_.Condition == a_cond) |> DataFrame
               qTRIAL[!, :SubPath] .= nothing #There is no subtraction
               #pull out only A-wave files
               data = readABF(qTRIALa.Path)
               dataABF = data_filter(data, avg_swp = false, t_pre = t_pre, t_post=t_post) #This function is found in the filter pipeline 
               #println(dataABF |> size)
          elseif i.Condition == b_cond
               #println("Analysis of B-wave file")
               qTRIALb = qData |> @filter(_.Condition == b_cond) |> DataFrame
               qTRIALa = qData |> @filter(_.Condition == a_cond) |> DataFrame
               qTRIAL = SubFiles = qTRIALa |> @join(qTRIALb,
                    {_.Year, _.Month, _.Date, _.Number, _.Photons, _.Wavelength, _.Photoreceptor, _.Genotype},
                    {_.Year, _.Month, _.Date, _.Number, _.Photons, _.Wavelength, _.Photoreceptor, _.Genotype},
                    {__...,
                         SubPath = _.Path,
               }) |> DataFrame
               if isempty(qTRIAL)
                    println("\t Experiment $(i.Path) was incomplete")
                    continue
               end
               data = readABF(SubFiles.Path) #Read the AB data
               filt_data = data_filter(data, avg_swp = false, t_pre = t_pre, t_post=t_post)
               dataSUB = readABF(SubFiles.SubPath)
               filt_dataSUB = data_filter(dataSUB, avg_swp = false, t_pre = t_pre, t_post=t_post)
     
               #if we want to subtract we need to filter first
               dataABF = filt_data - filt_dataSUB
               #println(size(dataABF))
          elseif i.Condition == g_cond
               #println("Analysis of Glial files")
               qTRIALb = qData |> @filter(_.Condition == b_cond) |> DataFrame
               qTRIALg = qData |> @filter(_.Condition == g_cond) |> DataFrame
               qTRIAL = SubFiles = qTRIALb |> @join(qTRIALg,
                    {_.Year, _.Month, _.Date, _.Number, _.Photons, _.Wavelength, _.Photoreceptor, _.Genotype},
                    {_.Year, _.Month, _.Date, _.Number, _.Photons, _.Wavelength, _.Photoreceptor, _.Genotype},
                    {__...,
                         SubPath = _.Path,
               }) |> DataFrame
               if isempty(qTRIAL)
                    println("\t Experiment $(i.Path) was incomplete")
                    continue
               end
               #println(qTRIAL |> size)
               data = readABF(SubFiles.Path) #Read the AB data
               filt_data = data_filter(data, avg_swp = false, t_pre = t_pre, t_post=t_post)
               dataSUB = readABF(SubFiles.SubPath)
               filt_dataSUB = data_filter(dataSUB, avg_swp = false, t_pre = t_pre, t_post=t_post)
     
               #if we want to subtract we need to filter first
               dataABF = filt_data - filt_dataSUB
               #println(size(dataABF))
          end
          #Conduct the analysis for each channel
          for (ch_idx, data_ch) in enumerate(eachchannel(dataABF)) #walk through each row of the data iterator
               gain = data.chTelegraph[ch_idx]
               if gain == 1 #Extract the gain == 1
                    if verbose 
                         println("Gain is a different")
                    end
                    data / 100.0
               end
               minimas = minimum(data_ch, dims=2)[:, 1, :]
               maximas = maximum(data_ch, dims=2)[:, 1, :]
               if i.Condition == a_cond 
                    if measure_minima
                         responses = abs.(minimas)
                    else
                         responses = abs.(saturated_response(data_ch))
                    end
               elseif i.Condition == b_cond
                    responses = abs.(maximas)
               elseif i.Condition == g_cond
                    responses = abs.(minimas)
               end
               rmax = maximum(responses, dims = 1)
               Peak_Times = time_to_peak(data_ch)
               Integrated_Times = integral(data_ch)
               Percent_Recoveries = percent_recovery_interval(data_ch, rmax)
               
               for swp in axes(data_ch, 1) #Walk through all sweep info, since sweeps will represent individual datafiles most of the time
                    #println(swp)
                    #inside_row = qData[idx, :Path] #We can break down each individual subsection by the inside row
                    push!(qTrace,
                         (
                              Path=qTRIAL[swp, :Path], 
                              #SubPath= isnothing(SubFiles) ? SubFiles[swp, :SubPath] : nothing,
                              Year=qTRIAL[swp, :Year], Month=qTRIAL[swp, :Month], Date=qTRIAL[swp, :Date],
                              Age=qTRIAL[swp, :Age], Number=qTRIAL[swp, :Number], Genotype=qTRIAL[swp, :Genotype],
                              Condition = qTRIAL[swp, :Condition], Photoreceptor=qTRIAL[swp, :Photoreceptor], Wavelength=qTRIAL[swp, :Wavelength],
                              Photons=qTRIAL[swp, :Photons],
                              Channel=data_ch.chNames[1], Gain=gain,
                              Response=responses[swp], Minima=minimas[swp], Maxima=maximas[swp],
                              Peak_Time=Peak_Times[swp], Integrated_Time=Integrated_Times[swp],
                              Percent_Recovery=Percent_Recoveries[swp]
                              #Recovery_Tau=Recovery_Taus[swp], Tau_GOF=Tau_GOFs[swp],
                         )
                    )

               end

               #This section we need to extract Rdim responses. 
               norm_resp = responses ./ maximum(responses)
               rdim_idxs = findall(0.20 .< norm_resp .< 0.50)
               if isempty(rdim_idxs)
                    rdim_idx = argmin(responses)
               else
                    rdim_min = argmin(responses[rdim_idxs])
                    rdim_idx = rdim_idxs[rdim_min]
               end

               if size(responses,1) > 2
                    if verbose 
                         println("Fitting IR curves")
                    end
                    p0 = [maximum(responses), median(qData[:, :Photons]), 2.0]
                    fit, rsq = ePhys.IRfit(qTRIAL[:, :Photons], responses |> vec, 
                         p0 = p0
                    )
                    if verbose
                         println("Fit r-squared = $rsq")
                    end
                         
                    #Fitting each data trace to a IR curve
                    push!(qExperiment, (
                         Year=qTRIAL[1, :Year], Month=qTRIAL[1, :Month], Date=qTRIAL[1, :Date],
                         Age=qTRIAL[1, :Age], Number=qTRIAL[1, :Number], Genotype=qTRIAL[1, :Genotype],
                         Condition = qTRIAL[1, :Condition],
                         Channel = data_ch.chNames[1],
                         Photoreceptor=qTRIAL[1, :Photoreceptor], Wavelength=qTRIAL[1, :Wavelength],
                         rmax = maximum(responses),
                         RMAX_fit = fit.param[1], K_fit = fit.param[2], N_fit = fit.param[3],
                         RSQ_fit = rsq, #, MSE_fit = mse_FIT,
                         rdim=responses[rdim_idx],
                         integration_time=Integrated_Times[rdim_idx],
                         time_to_peak=Peak_Times[rdim_idx],
                         percent_recovery = mean(Percent_Recoveries) #really strange these usually are averaged
                         #recovery_tau=Recovery_Taus[rdim_idx],
                    ))
               else
                    #Fitting each data trace to a IR curve
                    println("\t ONly $(length(responses)) responses. IR curve couldn't be fit")
                    push!(qExperiment, (
                         Year=qData[1, :Year], Month=qData[1, :Month], Date=qData[1, :Date],
                         Age=qData[1, :Age], Number=qData[1, :Number], Genotype=qData[1, :Genotype],
                         Condition = qTRIAL[1, :Condition],
                         Channel = data_ch.chNames[1],
                         Photoreceptor=qData[1, :Photoreceptor], Wavelength=qData[1, :Wavelength],
                         rmax = maximum(responses),
                         RMAX_fit = 0.0, K_fit = 0.0, N_fit = 0.0,
                         RSQ_fit = 0.0, #, MSE_fit = mse_FIT,
                         rdim=responses[rdim_idx],
                         integration_time=Integrated_Times[rdim_idx],
                         time_to_peak=Peak_Times[rdim_idx],
                         percent_recovery = mean(Percent_Recoveries) #really strange these usually are averaged
                         #recovery_tau=Recovery_Taus[rdim_idx],
                    ))
               end
          end
     end
     dataset["TRACES"] = qTrace
     dataset["EXPERIMENTS"] = qExperiment
     dataset["CONDITIONS"] = summarize_data(qTrace, qExperiment; lb = lb, p0 = p0, ub = ub)
     return dataset
end

function runTraceAnalysis(dataset_dict::Dict{String, DataFrame}; kwargs...)
     return runTraceAnalysis(dataset_dict["ALL_FILES"]; kwargs...)
end

#This can be used for IR and STF, but not for Tau or LP model
function GenerateFitFrame(df_TRACE, xData, yData; 
     MODEL = HILL_MODEL, #These function
     lb = (100.0, 1.0, 0.1), #Default rmin = 100, kmin = 0.1, nmin = 0.1 
     p0 = (500.0, 200.0, 2.0), #Default r = 500.0, k = 200.0, n = 2.0
     ub = (2400, 400, 4.0), #Default rmax = 2400, kmax = 800
     verbose = true
)
     df_EXP = df_TRACE |> @unique({_.Year, _.Month, _.Date, _.Number, _.Channel}) |> 
          @map({
               _.Year, _.Month, _.Date, _.Number, _.Channel,      
               RMAX = 0.0, K = 0.0, N = 0.0, 
               RSQ = 0.0,
               rmin = lb[1], r = p0[1], rmax = ub[1], #The highest b-wave ever seen is (2400)
               kmin = lb[2], k = p0[2], kmax = ub[2], #Half of the highest a-wave ever seen (800), 
               nmin = lb[3], n = p0[3], nmax = ub[3],
               MSE = 0.0,
               SS_Resid = 0.0, SS_Total = 0.0, 

          }) |> DataFrame
     for (idx, exp) in enumerate(eachrow(df_EXP))
          #println(exp)
          YEAR = exp.Year 
          MONTH = exp.Month 
          DATE = exp.Date 
          NUMBER = exp.Number 
          CHANNEL = exp.Channel
          if verbose
               print("Fitting exp $idx: $(YEAR)_$(MONTH)_$(DATE)_$(NUMBER)_$(CHANNEL)")
          end
          exp_traces = df_TRACE |> @filter((_.Year, _.Month, _.Date, _.Number, _.Channel) == (exp.Year, exp.Month, exp.Date, exp.Number, exp.Channel)) |> DataFrame
          println(exp_traces)
          #Conduct the fitting 
          p0 = [exp.r, exp.k, exp.n]
          ub = [exp.rmax, exp.kmax, exp.nmax]
          lb = [exp.rmin, exp.kmin, exp.nmin]
          exp_STF = curve_fit(MODEL, exp_traces[:, xData], exp_traces[:, yData], p0, lower = lb, upper = ub)

          df_EXP[idx, :RMAX] = exp_STF.param[1]
          df_EXP[idx, :K] = exp_STF.param[2]
          df_EXP[idx, :N] = exp_STF.param[3]
          df_EXP[idx, :SS_Resid] = ss_resid = sum(exp_STF.resid.^2)
          df_EXP[idx, :MSE] = ss_resid/size(exp_traces,1)
          yTrue = exp_traces[:, yData]
          yMean = mean(yTrue)
          df_EXP[idx, :SS_Total] = ss_total = sum((yTrue .- yMean).^2)
          df_EXP[idx, :RSQ] = 1 - (ss_resid/ss_total)
          if verbose
               println("RSQ is at $(println(df_EXP[idx, :RSQ]))")
          end
     end
     return df_EXP
end