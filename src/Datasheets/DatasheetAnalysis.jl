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

"""

"""
function run_A_wave_analysis(all_files::DataFrame; 
          measure_minima = false, a_cond = "BaCl_LAP4", 
          t_pre = 1.0, t_post = 2.0, #Extend these to see the end of the a-wave
          peak_method = false, 
          lb = [1.0, 1.0, 0.1], #Default rmin = 100, kmin = 0.1, nmin = 0.1 
          p0 = [500.0, 1000.0, 2.0], #Default r = 500.0, k = 200.0, n = 2.0
          ub = [Inf, Inf, 10.0], #Default rmax = 2400, kmax = 800
          verbose=false, 
     )
     a_files = all_files |> @filter(_.Condition == a_cond) |> DataFrame #Extract all A-wave responses
     if isnothing(a_files)
          return nothing
     else
          a_files[!, :Path] = string.(a_files[!, :Path]) #Make sure the path is a string
          uniqueData = a_files |> @unique({_.Year, _.Month, _.Date, _.Number, _.Wavelength, _.Photoreceptor, _.Genotype}) |> DataFrame #Pull out all unique files
          if verbose
               println("Completed data query")
          end
          qTrace = DataFrame() #Make empty dataframes for all traces
          qExperiment = DataFrame() #Make empty dataframe for all experiments
          for (idx, i) in enumerate(eachrow(uniqueData)) #Walk through each unique data
               qData = a_files |> @filter(
                         (_.Year, _.Month, _.Date, _.Number, _.Wavelength, _.Photoreceptor, _.Genotype) ==
                         (i.Year, i.Month, i.Date, i.Number, i.Wavelength, i.Photoreceptor, i.Genotype)
                    ) |>
                    DataFrame #Pull out the data 
               dataFile = readABF(qData.Path)

               if verbose
                    println("Completeing A-wave analysis for $idx out of $(size(uniqueData,1))")
                    println("Path: $(i.Path)")
               end
               for data in eachchannel(dataFile) #walk through each row of the data iterator
                    #println(qData.Age)
                    if isa(qData.Age[1], String)
                         matches = match(r"P(?'Age'\d*|)", qData.Age[1])
                         isadult = match(r"Adult", qData.Age[1]) 
                         if !isnothing(matches)
                              age = parse(Int, matches[:Age]) #Extract the age
                         elseif !isnothing(isadult)
                              age = 30
                         end
                    else
                         age = qData.Age[1]
                    end
                    ch = data.chNames[1] #Extract channel information
                    gain = data.chTelegraph[1] #Extract the gain
                    if gain == 1
                         if verbose 
                              println("Gain is a different")
                         end
                         data / 100.0
                    end
                    #======================DATA ANALYSIS========================#
                    if age <= 11 #If the data is below P11 calculate the response differently
                         filt_data = data_filter(data, avg_swp = false, t_pre = t_pre, t_post=t_post) #This function is found in the filter pipeline 
                         responses = minimas = minimum(filt_data, dims=2)[:, 1, :]
                         maximas = maximum(filt_data, dims=2)[:, 1, :]
                         Resps = abs.(responses)
                         rmaxes = minimum(responses, dims=1)
                    else
                         filt_data = data_filter(data, avg_swp = false, t_pre = t_pre, t_post=t_post)
                         if peak_method
                              #println(argmin(filt_data))
                              peak_time = filt_data.t[argmin(filt_data)[1][2]]
                              #println(peak_time)
                              truncate_data!(filt_data, t_pre = -peak_time) #Cut the trace before peak
                         end
                         if measure_minima
                              responses = minimas = minimum(filt_data, dims=2)[:, 1, :]
                              #println(minimas)
                         else
                              responses = saturated_response(filt_data)
                              minimas = minimum(filt_data, dims=2)[:, 1, :]
                         end
                         maximas = maximum(filt_data, dims=2)[:, 1, :]
                         Resps = abs.(responses)
                         rmaxes = minimum(responses, dims=1)
                    end
                    Peak_Times = time_to_peak(filt_data) #Calculate the time to peak
                    Integrated_Times = abs.(integral(filt_data)) #Calculate the area under the curve
                    Percent_Recoveries = percent_recovery_interval(filt_data, rmaxes) #Calculate the 60% recovery

                    #======================GLUING TOGETHER THE QUERY========================#
                    #now we can walk through each one of the responses and add it to the qTrace 
                    for swp = 1:size(data, 1) #Walk through all sweep info, since sweeps will represent individual datafiles most of the time
                         #inside_row = qData[idx, :Path] #We can break down each individual subsection by the inside row
                         push!(qTrace,
                              (
                                   Path=qData[swp, :Path],
                                   Year=qData[swp, :Year], Month=qData[swp, :Month], Date=qData[swp, :Date],
                                   Age=qData[swp, :Age], Number=qData[swp, :Number], Genotype=qData[swp, :Genotype],
                                   Condition = qData[swp, :Condition], Photoreceptor=qData[swp, :Photoreceptor], Wavelength=qData[swp, :Wavelength],
                                   Photons=qData[swp, :Photons],
                                   Channel=ch, Gain=gain,
                                   Response=Resps[swp], Minima=minimas[swp], Maxima=maximas[swp],
                                   Peak_Time=Peak_Times[swp], Integrated_Time=Integrated_Times[swp],
                                   Percent_Recovery=Percent_Recoveries[swp]
                                   #Recovery_Tau=Recovery_Taus[swp], Tau_GOF=Tau_GOFs[swp],
                              )
                         )

                    end

                    #This section we need to extract Rdim responses. 
                    norm_resp = Resps ./ maximum(Resps)
                    rdim_idxs = findall(0.20 .< norm_resp .< 0.50)
                    if isempty(rdim_idxs)
                         rdim_idx = argmin(Resps)
                    else
                         rdim_min = argmin(Resps[rdim_idxs])
                         rdim_idx = rdim_idxs[rdim_min]
                    end
                    if size(Resps,1) > 2
                         #if Resps |> vec
                         p0 = [maximum(Resps), median(qData[:, :Photons]), 2.0]
                         fit, rsq = ePhys.IRfit(qData[:, :Photons], Resps |> vec, 
                              p0 = p0
                         )
                         #Fitting each data trace to a IR curve
                         push!(qExperiment, (
                              Year=qData[1, :Year], Month=qData[1, :Month], Date=qData[1, :Date],
                              Age=qData[1, :Age], Number=qData[1, :Number], Genotype=qData[1, :Genotype],
                              Channel = ch,
                              Photoreceptor=qData[1, :Photoreceptor], Wavelength=qData[1, :Wavelength],
                              rmax = maximum(Resps),
                              RMAX_fit = fit.param[1], K_fit = fit.param[2], N_fit = fit.param[3],
                              RSQ_fit = rsq, #, MSE_fit = mse_FIT,
                              rdim=Resps[rdim_idx],
                              integration_time=Integrated_Times[rdim_idx],
                              time_to_peak=Peak_Times[rdim_idx],
                              percent_recovery = mean(Percent_Recoveries) #really strange these usually are averaged
                              #recovery_tau=Recovery_Taus[rdim_idx],
                         ))
                    else
                         #Fitting each data trace to a IR curve
                         push!(qExperiment, (
                              Year=qData[1, :Year], Month=qData[1, :Month], Date=qData[1, :Date],
                              Age=qData[1, :Age], Number=qData[1, :Number], Genotype=qData[1, :Genotype],
                              Channel = ch,
                              Photoreceptor=qData[1, :Photoreceptor], Wavelength=qData[1, :Wavelength],
                              rmax = maximum(Resps),
                              RMAX_fit = 0.0, K_fit = 0.0, N_fit = 0.0,
                              RSQ_fit = 0.0, #, MSE_fit = mse_FIT,
                              rdim=Resps[rdim_idx],
                              integration_time=Integrated_Times[rdim_idx],
                              time_to_peak=Peak_Times[rdim_idx],
                              percent_recovery = mean(Percent_Recoveries) #really strange these usually are averaged
                              #recovery_tau=Recovery_Taus[rdim_idx],
                         ))
                    end
               end
          end
          
          qConditions = summarize_data(qTrace, qExperiment; lb = lb, p0 = p0, ub = ub)
          return qTrace, qExperiment, qConditions
     end
end

#We can update this with our updated channel analysis
function run_B_wave_analysis(all_files::DataFrame; 
     t_pre = 1.0, t_post = 2.0, #This can be very quick, the initiation of the b-wave is faster
     a_cond = "BaCl_LAP4", 
     b_cond = "BaCl",
     verbose=false, 
)
     trace_A = all_files |> @filter(_.Condition == a_cond) |> DataFrame
     trace_AB = all_files |> @filter(_.Condition == b_cond) |> DataFrame
     if isempty(trace_AB)
          return nothing
     elseif isempty(trace_A)
          #println("Here")
          #return nothing
          b_files = trace_AB |> @join(trace_AB,
               {_.Year, _.Month, _.Date, _.Number, _.Photons, _.Wavelength, _.Photoreceptor, _.Genotype},
               {_.Year, _.Month, _.Date, _.Number, _.Photons, _.Wavelength, _.Photoreceptor, _.Genotype},
               {__...,
                    A_condition = _.Condition,
                    A_Path = _.Path,
               }) |> DataFrame
          if verbose
               println("Completed data query")
          end
     else
          b_files = trace_A |> @join(trace_AB,
               {_.Year, _.Month, _.Date, _.Number, _.Photons, _.Wavelength, _.Photoreceptor, _.Genotype},
               {_.Year, _.Month, _.Date, _.Number, _.Photons, _.Wavelength, _.Photoreceptor, _.Genotype},
               {__...,
                    A_condition = _.Condition,
                    A_Path = _.Path,
               }) |> DataFrame
          if verbose
               println("Completed data query")
          end
     end
     b_files[!, :Path] = string.(b_files[!, :Path]) #XLSX.jl converts things into Vector{Any}      
     b_files[!, :A_Path] = string.(b_files[!, :A_Path]) #XLSX.jl converts things into Vector{Any}            
     uniqueData = b_files |> @unique({_.Year, _.Month, _.Date, _.Number, _.Wavelength, _.Photoreceptor, _.Genotype}) |> DataFrame
     qTrace = DataFrame()
     qExperiment = DataFrame()
     for (idx, i) in enumerate(eachrow(uniqueData)) #We ca
          if verbose
               println("Completeing B-wave analysis for $idx out of $(size(uniqueData,1))")
               println("Path: $(i.Path)")
          end
          qData = b_files |> @filter(
               (_.Year, _.Month, _.Date, _.Number, _.Wavelength, _.Photoreceptor, _.Genotype) ==
               (i.Year, i.Month, i.Date, i.Number, i.Wavelength, i.Photoreceptor, i.Genotype)
          ) |> DataFrame

          
          data_AB = readABF(qData.Path) #Read the AB data
          filt_data_AB = data_filter(data_AB, avg_swp = false, t_pre = t_pre, t_post=t_post)

          data_A = readABF(qData.A_Path)
          filt_data_A = data_filter(data_A, avg_swp = false, t_pre = t_pre, t_post=t_post)

          #if we want to subtract we need to filter first
          sub_data = filt_data_AB - filt_data_A
          for (ch_idx, data_ch) in enumerate(eachchannel(sub_data)) #walk through each row of the data iterator

               unsubtracted_data_ch = getchannel(filt_data_AB, ch_idx)
               ch = data_ch.chNames[1] #Extract channel information
               gain = data_ch.chTelegraph[1]
               #Calculate the response based on the age
               if gain == 1
                    if verbose 
                         println("Gain is a different")
                    end
                    data_ch / 100
               end
               #======================DATA ANALYSIS========================#
               #data from P11 doesn't always make sense, so we can invert it 
               #if age == 11
               #     data_ch * -1
               #end
               responses = Resps = maximas = maximum(data_ch, dims=2)[:, 1, :]
               rmaxes = maximum(responses, dims=1)
               Unsubtracted_Resp = abs.(minima_to_peak(unsubtracted_data_ch)) #This is the unsubtracted Response
               minimas = minimum(data_ch, dims=2)[:, 1, :]
               maximas = maximum(data_ch, dims=2)[:, 1, :]
               Peak_Times = time_to_peak(data_ch)
               Integrated_Times = integral(data_ch)
               Percent_Recoveries = percent_recovery_interval(data_ch, rmaxes)
               #======================GLUING TOGETHER THE QUERY========================#
               for swp = 1:size(data_ch, 1) #Walk through all sweep info, since sweeps will represent individual datafiles most of the time
                    push!(qTrace,
                         (
                              Path=qData[swp, :Path], A_Path=qData[swp, :A_Path],
                              Year=qData[swp, :Year], Month=qData[swp, :Month], Date=qData[swp, :Date],
                              Age=qData[swp, :Age], Number=qData[swp, :Number], Genotype=qData[swp, :Genotype],
                              Condition = qData[swp, :Condition], Photoreceptor=qData[swp, :Photoreceptor], Wavelength=qData[swp, :Wavelength],
                              Photons=qData[swp, :Photons],
                              Channel=ch, Gain=gain,
                              Response=Resps[swp], Unsubtracted_Response=Unsubtracted_Resp[swp],
                              Minima=minimas[swp], Maxima=maximas[swp],
                              Peak_Time=Peak_Times[swp], Integrated_Time=Integrated_Times[swp],
                              Percent_Recovery=Percent_Recoveries[swp]
                         )
                    )
               end

               norm_resp = Resps ./ maximum(Resps)
               #println(norm_resp)
               rdim_idxs = findall(0.20 .< norm_resp .< 0.50)
               if isempty(rdim_idxs)
                    rdim_idx = argmin(Resps)
               else
                    rdim_min = argmin(Resps[rdim_idxs])
                    rdim_idx = rdim_idxs[rdim_min]
               end
               if size(Resps,1) > 2
                    p0 = [maximum(Resps), median(qData[:, :Photons]), 2.0]
                    fit, rsq = ePhys.IRfit(qData[:, :Photons], Resps |> vec, 
                         p0 = p0
                    )
                    #println(fit, rsq)
                    push!(qExperiment, (
                         Year=qData[1, :Year], Month=qData[1, :Month], Date=qData[1, :Date],
                         Age=qData[1, :Age], Number=qData[1, :Number], Genotype=qData[1, :Genotype],
                         Photoreceptor=qData[1, :Photoreceptor], Wavelength=qData[1, :Wavelength],
                         Channel = ch,
                         Photons=qData[1, :Photons],
                         rmax=maximum(Resps),
                         RMAX_fit = fit.param[1], K_fit = fit.param[2], N_fit = fit.param[3],
                         RSQ_fit = rsq, #, MSE_fit = mse_FIT,
                         unsubtracted_rmax=maximum(Unsubtracted_Resp),
                         rdim=Resps[rdim_idx],
                         integration_time=Integrated_Times[rdim_idx],
                         time_to_peak=Peak_Times[rdim_idx],
                         percent_recovery=maximum(Percent_Recoveries) #sum(Percent_Recoveries) / length(Percent_Recoveries) #really strange these usually are averaged
                    ))
               else
                    push!(qExperiment, (
                         Year=qData[1, :Year], Month=qData[1, :Month], Date=qData[1, :Date],
                         Age=qData[1, :Age], Number=qData[1, :Number], Genotype=qData[1, :Genotype],
                         Photoreceptor=qData[1, :Photoreceptor], Wavelength=qData[1, :Wavelength],
                         Channel = ch,
                         Photons=qData[1, :Photons],
                         rmax=maximum(Resps),
                         RMAX_fit = 0.0, K_fit = 0.0, N_fit = 0.0,
                         RSQ_fit = 0.0, #, MSE_fit = mse_FIT,
                         unsubtracted_rmax=maximum(Unsubtracted_Resp),
                         rdim=Resps[rdim_idx],
                         integration_time=Integrated_Times[rdim_idx],
                         time_to_peak=Peak_Times[rdim_idx],
                         percent_recovery=maximum(Percent_Recoveries) #sum(Percent_Recoveries) / length(Percent_Recoveries) #really strange these usually are averaged
                    ))
               end
          end
     end
     qConditions = summarize_data(qTrace, qExperiment)
     return qTrace, qExperiment, qConditions
end

"""
There is no version of G component analysis that is not subtractive
"""
function run_G_wave_analysis(all_files::DataFrame; 
     t_pre = 1.0, t_post = 2.0, 
     verbose=false, 
)
     trace_ABG = all_files |> @filter(_.Condition == "NoDrugs") |> DataFrame
     trace_AB = all_files |> @filter(_.Condition == "BaCl") |> DataFrame
     if isempty(trace_ABG)
          return nothing
     elseif isempty(trace_AB)
          return nothing
     else
          g_files = trace_AB |> @join(trace_ABG,
               {_.Year, _.Month, _.Date, _.Number, _.Photons, _.Wavelength, _.Photoreceptor, _.Genotype},
               {_.Year, _.Month, _.Date, _.Number, _.Photons, _.Wavelength, _.Photoreceptor, _.Genotype},
               {__...,
                    AB_condition = _.Condition,
                    AB_Path = _.Path,
          }) |> DataFrame
          g_files[!, :Path] = string.(g_files[!, :Path]) #XLSX.jl converts things into Vector{Any}      
          g_files[!, :AB_Path] = string.(g_files[!, :AB_Path]) #XLSX.jl converts things into Vector{Any}            

          uniqueData = g_files |> @unique({_.Year, _.Month, _.Date, _.Number, _.Wavelength, _.Photoreceptor, _.Genotype}) |> DataFrame
          qTrace = DataFrame()
          qExperiment = DataFrame()
          for (idx, i) in enumerate(eachrow(uniqueData)) #We ca
               qData = g_files |> @filter(
                         (_.Year, _.Month, _.Date, _.Number, _.Wavelength, _.Photoreceptor, _.Genotype) ==
                         (i.Year, i.Month, i.Date, i.Number, i.Wavelength, i.Photoreceptor, i.Genotype)
               ) |> DataFrame
               data_ABG = readABF(qData.Path)
               filt_data_ABG = data_filter(data_ABG, avg_swp = false, t_pre = t_pre, t_post = t_post)
               data_AB = readABF(qData.AB_Path)
               filt_data_AB = data_filter(data_AB, avg_swp = false, t_pre = t_pre, t_post = t_post)
               #if we want to subtract we need to filter first
               sub_data = filt_data_ABG - filt_data_AB
               if verbose
                    println("Completeing Glial analysis for $idx out of $(size(uniqueData,1))")
                    println("Path: $(i.Path)")
               end
               for (ch_idx, data_ch) in enumerate(eachchannel(sub_data)) #walk through each row of the data iterator
                    ch = data_ch.chNames[1] #Extract channel information
                    gain = data_ch.chTelegraph[1]
                    if gain == 1
                         data_ch / 100.0
                    end
                    #Calculate the response based on the age
                    unsubtracted_data_ch = getchannel(filt_data_ABG, ch_idx)
                    #======================DATA ANALYSIS========================#
                    responses = minimas = minimum(data_ch, dims=2)[:, 1, :]
                    maximas = maximum(data_ch, dims=2)[:, 1, :]
                    rmaxes = minimum(responses, dims=1)
                    Resps = abs.(responses)
                    Unsubtracted_Resp = abs.(minimum(unsubtracted_data_ch, dims=2)[:, 1, :]) #This is the unsubtracted Response
                    #minimas = minimum(data_ch, dims=2)[:, 1, :]
                    Peak_Times = time_to_peak(data_ch)
                    Integrated_Times = integral(data_ch)
                    Percent_Recoveries = percent_recovery_interval(data_ch, rmaxes)

                    #======================GLUING TOGETHER THE QUERY========================#
                    #now we can walk through each one of the responses and add it to the qTrace 
                    for swp = 1:size(data_ch, 1) #Walk through all sweep info, since sweeps will represent individual datafiles most of the time
                         #inside_row = qData[idx, :Path] #We can break down each individual subsection by the inside row
                         push!(qTrace, (
                              Path=qData[swp, :Path], AB_Path=qData[swp, :AB_Path],
                              Year=qData[swp, :Year], Month=qData[swp, :Month], Date=qData[swp, :Date],
                              Age=qData[swp, :Age], Number=qData[swp, :Number], Genotype=qData[swp, :Genotype],
                              Condition = qData[swp, :Condition], Photoreceptor=qData[swp, :Photoreceptor], Wavelength=qData[swp, :Wavelength],
                              Photons=qData[swp, :Photons],
                              Channel=ch, Gain=gain,
                              Response=Resps[swp], Unsubtracted_Response=Unsubtracted_Resp[swp],
                              Minima=minimas[swp], Maxima=maximas[swp],
                              Peak_Time=Peak_Times[swp], 
                              Integrated_Time=Integrated_Times[swp],
                              Percent_Recovery=Percent_Recoveries[swp])
                         )
                    end
                    norm_resp = Resps ./ maximum(Resps)
                    #println(norm_resp)
                    rdim_idxs = findall(0.20 .< norm_resp .< 0.50)
                    if isempty(rdim_idxs)
                         rdim_idx = argmin(Resps)
                    else
                         rdim_min = argmin(Resps[rdim_idxs])
                         rdim_idx = rdim_idxs[rdim_min]
                    end
                    if size(Resps,1) > 2
                         #println(size(Resps))
                         #if Resps |> vec
                         p0 = [maximum(Resps), median(qData[:, :Photons]), 2.0]
                         fit, rsq = ePhys.IRfit(qData[:, :Photons], Resps |> vec, 
                              p0 = p0
                         )
                         push!(qExperiment, (
                              Year=qData[1, :Year], Month=qData[1, :Month], Date=qData[1, :Date],
                              Age=qData[1, :Age], Number=qData[1, :Number], Genotype=qData[1, :Genotype],
                              Channel = ch,
                              Photoreceptor=qData[1, :Photoreceptor], Wavelength=qData[1, :Wavelength],
                              Photons=qData[1, :Photons],
                              rmax=maximum(Resps),
                              RMAX_fit = fit.param[1], K_fit = fit.param[2], N_fit = fit.param[3],
                              RSQ_fit = rsq, #, MSE_fit = mse_FIT,
                              unsubtracted_rmax=maximum(Unsubtracted_Resp),
                              rdim=Resps[rdim_idx],
                              integration_time=Integrated_Times[rdim_idx],
                              time_to_peak=Peak_Times[rdim_idx],
                              percent_recovery=maximum(Percent_Recoveries) #This is maximum vs average
                         ))
                    else
                         push!(qExperiment, (
                              Year=qData[1, :Year], Month=qData[1, :Month], Date=qData[1, :Date],
                              Age=qData[1, :Age], Number=qData[1, :Number], Genotype=qData[1, :Genotype],
                              Channel = ch,
                              Photoreceptor=qData[1, :Photoreceptor], Wavelength=qData[1, :Wavelength],
                              Photons=qData[1, :Photons],
                              rmax=maximum(Resps),
                              RMAX_fit = 0.0, K_fit = 0.0, N_fit = 0.0,
                              RSQ_fit = 0.0, #, MSE_fit = mse_FIT,
                              unsubtracted_rmax=maximum(Unsubtracted_Resp),
                              rdim=Resps[rdim_idx],
                              integration_time=Integrated_Times[rdim_idx],
                              time_to_peak=Peak_Times[rdim_idx],
                              percent_recovery=maximum(Percent_Recoveries) #This is maximum vs average
                         ))
                    end
               end
          end
          qConditions = summarize_data(qTrace, qExperiment)
          return qTrace, qExperiment, qConditions
     end
end

function add_analysis_sheets(results, save_file::String; append="A")
     trace, experiments, conditions = results
     dataset = openDatasheet(save_file; sheetName = "all", typeConvert = false)
     rewrite_sheets = ["trace_$append", "experiments_$append", "conditions_$append"]
     for sheet in rewrite_sheets #We need to add any new sheets to the current sheet. Kind of a waste but it works
          if sheet == "trace_$append"
               dataset[sheet] = trace
          elseif sheet == "experiments_$append"
               dataset[sheet] = experiments
          elseif sheet == "conditions_$append"
               dataset[sheet] = conditions
          end
     end

     XLSX.openxlsx(save_file, mode = "w") do xf
          #This will always add a sheet1
          for sName in keys(dataset)
               #println("$sName")
               if sName == "All_Files"
                    XLSX.rename!(xf[1], "All_Files")
               else
                    XLSX.addsheet!(xf, sName)
               end
               XLSX.writetable!(xf[sName], dataset[sName])
          end
     end
end


function runAnalysis(datasheet::DataFrame; measure_minima = false, verbose = false)
     resA = ePhys.run_A_wave_analysis(datasheet; measure_minima = measure_minima, verbose = verbose)
     resB = ePhys.run_B_wave_analysis(datasheet, verbose = verbose)
     resG = ePhys.run_G_wave_analysis(datasheet, verbose = verbose)
     return (datasheet, resA, resB, resG)
end

function runAnalysis(datafile::String; kwargs...)
     print("Opening datafile $(datafile)... ")
     datasheet = openDatasheet(datafile, sheetName = "All_Files")
     println("complete")
     datasheet, resA, resB, resG = runAnalysis(datasheet; kwargs...)
     if !isnothing(resA)
          add_analysis_sheets(resA, datafile; append="A")
     end

     if !isnothing(resB)
          add_analysis_sheets(resB, datafile; append="B")
     end

     if !isnothing(resG)
          add_analysis_sheets(resG, datafile; append="G")
     end
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

function GenerateTimeFitFrame(df_TRACE)

end