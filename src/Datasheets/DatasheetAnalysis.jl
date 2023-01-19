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

"""

"""
function run_A_wave_analysis(all_files::DataFrame; 
          run_amp=false, verbose=true, measure_minima = false, a_cond = "BaCl_LAP4"
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
          #println(uniqueData)
          #println(size(uniqueData))
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
                    matches = match(r"P(?'Age'\d*|)", qData.Age[1])|> NamedTuple
                    println(matches.Age)
                    age = parse(Int, matches.Age) #Extract the age
                    ch = data.chNames[1] #Extract channel information
                    gain = data.chTelegraph[1] #Extract the gain
                    if gain == 1
                         data / 100.0
                    end
                    #======================DATA ANALYSIS========================#
                    if age <= 11 #If the data is below P11 calculate the response differently
                         filt_data = data_filter(data, avg_swp = false, t_post=5.0) #This function is found in the filter pipeline 
                         responses = minimas = minimum(filt_data, dims=2)[:, 1, :]
                         maximas = maximum(filt_data, dims=2)[:, 1, :]
                         Resps = abs.(responses)
                         rmaxes = minimum(responses, dims=1)
                    else
                         filt_data = data_filter(data, avg_swp = false, t_post=5.0)
                         if measure_minima
                              responses = minimas = minimum(filt_data, dims=2)[:, 1, :]
                         else
                              responses = saturated_response(filt_data)
                              minimas = minimum(filt_data, dims=2)[:, 1, :]
                         end
                         maximas = maximum(filt_data, dims=2)[:, 1, :]
                         Resps = abs.(responses)
                         rmaxes = minimum(responses, dims=1)
                    end
                    Peak_Times = time_to_peak(filt_data) #Calculate the time to peak
                    Integrated_Times = abs.(integral(filt_data))
                    rec_res = recovery_time_constant(filt_data, responses)
                    Recovery_Taus = rec_res[1] |> vec
                    Tau_GOFs = rec_res[2] |> vec

                    #The 60% Bandwith can be used to measure the percent recovery interval
                    Percent_Recoveries = percent_recovery_interval(filt_data, rmaxes)
                    #We need to program the amplification as well. But that may be longer

                    #======================GLUING TOGETHER THE QUERY========================#
                    #now we can walk through each one of the responses and add it to the qTrace 
                    for swp = 1:size(data, 1) #Walk through all sweep info, since sweeps will represent individual datafiles most of the time
                         #inside_row = qData[idx, :Path] #We can break down each individual subsection by the inside row
                         push!(qTrace,
                              (
                                   Path=qData[swp, :Path],
                                   Year=qData[swp, :Year], Month=qData[swp, :Month], Date=qData[swp, :Date],
                                   Age=qData[swp, :Age], Number=qData[swp, :Number], Genotype=qData[swp, :Genotype],
                                   Photoreceptor=qData[swp, :Photoreceptor], Wavelength=qData[swp, :Wavelength],
                                   Photons=qData[swp, :Photons],
                                   Channel=ch, Gain=gain,
                                   Response=Resps[swp], Minima=minimas[swp], Maxima=maximas[swp],
                                   Peak_Time=Peak_Times[swp], Integrated_Time=Integrated_Times[swp],
                                   Recovery_Tau=Recovery_Taus[swp], Tau_GOF=Tau_GOFs[swp],
                                   Percent_Recovery=Percent_Recoveries[swp]
                              )
                         )

                    end

                    #Each data entry will be added to the qExperiment frame
                    #This section we need to extract Rdim responses. 
                    norm_resp = Resps ./ maximum(Resps)
                    rdim_idxs = findall(0.20 .< norm_resp .< 0.50)
                    if isempty(rdim_idxs)
                         rdim_idx = argmin(Resps)
                    else
                         rdim_min = argmin(Resps[rdim_idxs])
                         rdim_idx = rdim_idxs[rdim_min]
                    end
                    #try fitting the data for the rmax and k individually
                    rmax_FIT = 0.0
                    k_FIT = 0.0
                    n_FIT = 0.0
                    try
                         fit = IRfit(qData[:, :Photons], Resps[:,1], r = maximum(Resps))
                         rmax_FIT, k_FIT, n_FIT = fit.param
                    catch
                         println("Something went wrong")
                         
                    end
                    push!(qExperiment, (
                         Year=qData[1, :Year], Month=qData[1, :Month], Date=qData[1, :Date],
                         Age=qData[1, :Age], Number=qData[1, :Number], Genotype=qData[1, :Genotype],
                         Photoreceptor=qData[1, :Photoreceptor], Wavelength=qData[1, :Wavelength],
                         #Photons=qData[1, :Photons], #We dont' need the first number of photons
                         rmax=maximum(Resps),
                         rmax_fit = rmax_FIT, k = k_FIT, n = n_FIT,
                         rdim=Resps[rdim_idx],
                         integration_time=Integrated_Times[rdim_idx],
                         time_to_peak=Peak_Times[rdim_idx],
                         recovery_tau=Recovery_Taus[rdim_idx],
                         percent_recovery=maximum(Percent_Recoveries) # sum(Percent_Recoveries) / length(Percent_Recoveries) #really strange these usually are averaged
                    ))

               end
          end

          qConditions = qExperiment |>
                    @groupby({_.Age, _.Genotype, _.Photoreceptor, _.Wavelength}) |>
                    @map({
                         Age = _.Age[1], Genotype = _.Genotype[1], Photoreceptor = _.Photoreceptor[1], Wavelength = _.Wavelength[1],
                         N = length(_),
                         Rmax = mean(_.rmax), Rmax_sem = sem(_.rmax),
                         Rdim = mean(_.rdim), Rdim_sem = sem(_.rdim),
                         Integrated_Time = mean(_.integration_time), Integrated_Time_sem = sem(_.integration_time),
                         Time_to_peak = mean(_.time_to_peak), Time_To_Peak_sem = sem(_.time_to_peak),
                         Recovery_Tau = mean(_.recovery_tau), Recovery_Tau_sem = sem(_.recovery_tau),
                         Percent_Recovery = mean(_.percent_recovery), Percent_Recovery_sem = sem(_.percent_recovery)
                    }) |>
                    DataFrame
          #retroactively go through and fit all IR data
          rmax_fit = Float64[]
          k_fit = Float64[]
          n_fit = Float64[]
          for (idx, categ) in enumerate(extract_categories(qConditions))
               allIR, fiti = extractIR(qTrace, categ, measure = :Minima)
               #if categ[1] <= 11
               #else
               #     allIR, fiti = extractIR(qTrace, categ)
               #end
               push!(rmax_fit, abs(fiti.param[1]))
               push!(k_fit, fiti.param[2])
               push!(n_fit, fiti.param[3])
          end
          qConditions[!, :RmaxFit] = rmax_fit
          qConditions[!, :k] = k_fit
          qConditions[!, :n] = n_fit
          return qTrace, qExperiment, qConditions
     end
end

#We can update this with our updated channel analysis
function run_B_wave_analysis(all_files::DataFrame; verbose=true, a_cond = "BaCl_LAP4", b_cond = "BaCl")
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
     else
          b_files = trace_A |> @join(trace_AB,
               {_.Year, _.Month, _.Date, _.Number, _.Photons, _.Wavelength, _.Photoreceptor, _.Genotype},
               {_.Year, _.Month, _.Date, _.Number, _.Photons, _.Wavelength, _.Photoreceptor, _.Genotype},
               {__...,
                    A_condition = _.Condition,
                    A_Path = _.Path,
               }) |> DataFrame
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
               ) |>
               DataFrame

          data_AB = readABF(qData.Path)
          filt_data_AB = data_filter(data_AB, avg_swp = false, t_post=5.0)

          data_A = readABF(qData.A_Path)
          filt_data_A = data_filter(data_A, avg_swp = false, t_post=5.0)
          #if we want to subtract we need to filter first
          sub_data = filt_data_AB - filt_data_A
          for (ch_idx, data_ch) in enumerate(eachchannel(sub_data)) #walk through each row of the data iterator

               unsubtracted_data_ch = getchannel(filt_data_AB, ch_idx)
               age = qData.Age[1] #Extract the age
               ch = data_ch.chNames[1] #Extract channel information
               gain = data_ch.chTelegraph[1]
               #Calculate the response based on the age
               if gain == 1
                    println("Gain is a different")
                    data_ch / 100
               end
               #======================DATA ANALYSIS========================#
               #data from P11 doesn't always make sense, so we can invert it 
               if age == 11
                    data_ch * -1
               end
               responses = Resps = maximas = maximum(data_ch, dims=2)[:, 1, :]
               rmaxes = maximum(responses, dims=1)
               Unsubtracted_Resp = abs.(minima_to_peak(unsubtracted_data_ch)) #This is the unsubtracted Response
               minimas = minimum(data_ch, dims=2)[:, 1, :]
               maximas = maximum(data_ch, dims=2)[:, 1, :]
               Peak_Times = time_to_peak(data_ch)
               Integrated_Times = integral(data_ch)
               rec_res = recovery_time_constant(data_ch, Resps)
               Recovery_Taus = rec_res[1] |> vec
               Tau_GOFs = rec_res[2] |> vec

               #We need to program the amplification as well. But that may be longer
               Percent_Recoveries = percent_recovery_interval(data_ch, rmaxes)
               #======================GLUING TOGETHER THE QUERY========================#
               for swp = 1:size(data_ch, 1) #Walk through all sweep info, since sweeps will represent individual datafiles most of the time
                    push!(qTrace,
                         (
                              Path=qData[swp, :Path], A_Path=qData[swp, :A_Path],
                              Year=qData[swp, :Year], Month=qData[swp, :Month], Date=qData[swp, :Date],
                              Age=qData[swp, :Age], Number=qData[swp, :Number], Genotype=qData[swp, :Genotype],
                              Photoreceptor=qData[swp, :Photoreceptor], Wavelength=qData[swp, :Wavelength],
                              Photons=qData[swp, :Photons],
                              Channel=ch, Gain=gain,
                              Response=Resps[swp], Unsubtracted_Response=Unsubtracted_Resp[swp],
                              Minima=minimas[swp], Maxima=maximas[swp],
                              Peak_Time=Peak_Times[swp], Integrated_Time=Integrated_Times[swp],
                              Recovery_Tau=Recovery_Taus[swp], Tau_GOF=Tau_GOFs[swp],
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

               push!(qExperiment, (
                    Year=qData[1, :Year], Month=qData[1, :Month], Date=qData[1, :Date],
                    Age=qData[1, :Age], Number=qData[1, :Number], Genotype=qData[1, :Genotype],
                    Photoreceptor=qData[1, :Photoreceptor], Wavelength=qData[1, :Wavelength],
                    Photons=qData[1, :Photons],
                    rmax=maximum(Resps),
                    unsubtracted_rmax=maximum(Unsubtracted_Resp),
                    rdim=Resps[rdim_idx],
                    integration_time=Integrated_Times[rdim_idx],
                    time_to_peak=Peak_Times[rdim_idx],
                    recovery_tau=Recovery_Taus[rdim_idx],
                    percent_recovery=maximum(Percent_Recoveries) #sum(Percent_Recoveries) / length(Percent_Recoveries) #really strange these usually are averaged
               ))
          end
     end
     println(qExperiment)
     qConditions = qExperiment |>
          @groupby({_.Age, _.Genotype, _.Photoreceptor, _.Wavelength}) |>
          @map({
               Age = _.Age[1], Genotype = _.Genotype[1], Photoreceptor = _.Photoreceptor[1], Wavelength = _.Wavelength[1],
               N = length(_),
               Rmax = mean(_.rmax), Rmax_sem = sem(_.rmax),
               Unsubtracted_Rmax = mean(_.unsubtracted_rmax), Unsubtracted_Rmax_sem = sem(_.unsubtracted_rmax),
               Rdim = mean(_.rdim), Rdim_sem = sem(_.rdim),
               Integrated_Time = mean(_.integration_time), Integrated_Time_sem = sem(_.integration_time),
               Time_to_peak = mean(_.time_to_peak), Time_To_Peak_sem = sem(_.time_to_peak),
               Recovery_Tau = mean(_.recovery_tau), Recovery_Tau_sem = sem(_.recovery_tau),
               Percent_Recovery = mean(_.percent_recovery), Percent_Recovery_sem = sem(_.percent_recovery)}) |>
          DataFrame

     return qTrace, qExperiment, qConditions
end

"""
There is no version of G component analysis that is not subtractive
"""
function run_G_wave_analysis(all_files::DataFrame; verbose=true)
     trace_ABG = all_files |> @filter(_.Condition == "NoDrugs") |> DataFrame
     trace_AB = all_files |> @filter(_.Condition == "BaCl") |> DataFrame
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
     #println(size(uniqueData))
     #println(size(trace_ABG))
     #println(size(trace_AB))
     #println(size(g_files))
     qTrace = DataFrame()
     qExperiment = DataFrame()
     for (idx, i) in enumerate(eachrow(uniqueData)) #We ca
          qData = g_files |> @filter(
                       (_.Year, _.Month, _.Date, _.Number, _.Wavelength, _.Photoreceptor, _.Genotype) ==
                       (i.Year, i.Month, i.Date, i.Number, i.Wavelength, i.Photoreceptor, i.Genotype)
                  ) |>
                  DataFrame
          data_ABG = readABF(qData.Path)
          filt_data_ABG = data_filter(data_ABG, avg_swp = false, t_post=5.0)
          data_AB = readABF(qData.AB_Path)
          filt_data_AB = data_filter(data_AB, avg_swp = false, t_post=5.0)
          #if we want to subtract we need to filter first
          #println(qData.Path)
          #println(qData.AB_Path)
          filt_data = filt_data_ABG - filt_data_AB


          if verbose
               println("Completeing Glial analysis for $idx out of $(size(uniqueData,1))")
          end
          for (ch_idx, data_ch) in enumerate(eachchannel(filt_data)) #walk through each row of the data iterator
               age = qData.Age[1] #Extract the age
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
               rec_res = recovery_time_constant(data_ch, Resps)
               Recovery_Taus = rec_res[1] |> vec
               Tau_GOFs = rec_res[2] |> vec
               Percent_Recoveries = percent_recovery_interval(data_ch, rmaxes)
               #println(Percent_Recoveries)
               #We need to program the amplification as well. But that may be longer

               #======================GLUING TOGETHER THE QUERY========================#
               #now we can walk through each one of the responses and add it to the qTrace 
               for swp = 1:size(data_ch, 1) #Walk through all sweep info, since sweeps will represent individual datafiles most of the time
                    #inside_row = qData[idx, :Path] #We can break down each individual subsection by the inside row
                    push!(qTrace, (
                         Path=qData[swp, :Path], AB_Path=qData[swp, :AB_Path],
                         Year=qData[swp, :Year], Month=qData[swp, :Month], Date=qData[swp, :Date],
                         Age=qData[swp, :Age], Number=qData[swp, :Number], Genotype=qData[swp, :Genotype],
                         Photoreceptor=qData[swp, :Photoreceptor], Wavelength=qData[swp, :Wavelength],
                         Photons=qData[swp, :Photons],
                         Channel=ch, Gain=gain,
                         Response=Resps[swp], Unsubtracted_Response=Unsubtracted_Resp[swp],
                         Minima=minimas[swp], Maxima=maximas[swp],
                         Peak_Time=Peak_Times[swp], Integrated_Time=Integrated_Times[swp],
                         Recovery_Tau=Recovery_Taus[swp], Tau_GOF=Tau_GOFs[swp],
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
               push!(qExperiment, (
                    Year=qData[1, :Year], Month=qData[1, :Month], Date=qData[1, :Date],
                    Age=qData[1, :Age], Number=qData[1, :Number], Genotype=qData[1, :Genotype],
                    Photoreceptor=qData[1, :Photoreceptor], Wavelength=qData[1, :Wavelength],
                    Photons=qData[1, :Photons],
                    rmax=maximum(Resps),
                    unsubtracted_rmax=maximum(Unsubtracted_Resp),
                    rdim=Resps[rdim_idx],
                    integration_time=Integrated_Times[rdim_idx],
                    time_to_peak=Peak_Times[rdim_idx],
                    recovery_tau=Recovery_Taus[rdim_idx],
                    #percent_recovery=sum(Percent_Recoveries) / length(Percent_Recoveries) #really strange these usually are averaged
                    percent_recovery=maximum(Percent_Recoveries) #This is maximum vs average
                    #percent_recovery= Percent_Recoveries[rdim_idx] #this is the recovery associated with the rdim
               ))
          end
     end
     qConditions = qExperiment |>
                   @groupby({_.Age, _.Genotype, _.Photoreceptor, _.Wavelength}) |>
                   @map({
                        Age = _.Age[1], Genotype = _.Genotype[1], Photoreceptor = _.Photoreceptor[1], Wavelength = _.Wavelength[1],
                        N = length(_),
                        Rmax = mean(_.rmax), Rmax_sem = sem(_.rmax),
                        Unsubtracted_Rmax = mean(_.unsubtracted_rmax), Unsubtracted_Rmax_sem = sem(_.unsubtracted_rmax),
                        Rdim = mean(_.rdim), Rdim_sem = sem(_.rdim),
                        Integrated_Time = mean(_.integration_time), Integrated_Time_sem = sem(_.integration_time),
                        Time_to_peak = mean(_.time_to_peak), Time_To_Peak_sem = sem(_.time_to_peak),
                        Recovery_Tau = mean(_.recovery_tau), Recovery_Tau_sem = sem(_.recovery_tau),
                        Percent_Recovery = mean(_.percent_recovery), Percent_Recovery_sem = sem(_.percent_recovery)}) |>
                   DataFrame

     return qTrace, qExperiment, qConditions
end

function add_analysis_sheets(results, save_file::String; append="A")
     trace, experiments, conditions = results
     XLSX.openxlsx(save_file, mode="rw") do xf
          try
               sheet = xf["trace_$(append)"] #try to open the sheet
               #clean the data from the sheet
               println("Cleaning Data")
               cleanDatasheet!(xf, "trace_$(append)")
          catch #the sheet is not made and must be created
               println("Adding trace sheets")
               XLSX.addsheet!(xf, "trace_$(append)")
          end
          XLSX.writetable!(xf["trace_$(append)"],
               collect(DataFrames.eachcol(trace)),
               DataFrames.names(trace))
     end
     #Extract experiments for A wave

     XLSX.openxlsx(save_file, mode="rw") do xf
          try
               sheet = xf["experiments_$(append)"]
               cleanDatasheet!(xf, "experiments_$(append)")
          catch #the sheet is not made and must be created
               println("Adding experiment sheets")
               XLSX.addsheet!(xf, "experiments_$(append)")
          end
          XLSX.writetable!(xf["experiments_$(append)"],
               collect(DataFrames.eachcol(experiments)),
               DataFrames.names(experiments))
     end

     XLSX.openxlsx(save_file, mode="rw") do xf
          try
               sheet = xf["conditions_$(append)"]
               cleanDatasheet!(xf, "conditions_$(append)")
          catch
               println("Adding condition sheets")
               XLSX.addsheet!(xf, "conditions_$(append)")
          end
          XLSX.writetable!(xf["conditions_$(append)"],
               collect(DataFrames.eachcol(conditions)),
               DataFrames.names(conditions))
     end
end

function runAnalysis(datafile::String; measure_minima = false)
     print("Opening datafile $(datafile)... ")
     df = openDatasheet(datafile)
     println("complete")
     #%% Test the a, b, and g wave analysis
     resA = ePhys.run_A_wave_analysis(df; measure_minima = measure_minima)
     add_analysis_sheets(resA, datafile; append="A")
     resB = ePhys.run_B_wave_analysis(df)
     add_analysis_sheets(resB, datafile; append="B")
     resG = ePhys.run_G_wave_analysis(df)
     add_analysis_sheets(resG, datafile; append="G")
     return (df, resA, resB, resG)
end

#Couldn't we potentially use this function to extract the datasheet anyways
#TODO: We want to change this to use the DataPathExtraction instead
function restructure_filesystem(prev_paths, restruc_path; 
          remove_old = false, 
          unknown_genotype = :WT, conditions_error = false, 
          verbose = false, mode = :save   
     )
     if remove_old
          rm(restruc_path, recursive=true)
     end
     if !isdir(restruc_path) #If the file does not exist
          mkdir(restruc_path) #Make the file
     end
     count = 0
     current_paths = restruc_path |> parseABF #These are all the files that exist inside of the resturctured path
     for path in prev_paths
          if verbose
               print("Extracting files from path: ")
               println(path)
          end

          date_res = findmatch(path, date_regex)
          if !isnothing(date_res) 
               YEAR = date_res.Year
               MONTH = date_res.Month
               DATE = date_res.Date
          else
               throw("Incorrect date specification") #Always throw a date error
          end

          #Find the age and 
          animal_res = findmatch(path, animal_regex)
          if !isnothing(animal_res)
               if !isnothing(animal_res.Animal)
                    ANIMAL = animal_res.Animal
               else
                    ANIMAL = "Mouse"
               end
               NUMBER = animal_res.Number
          else
               if verbose
                    println(path)
                    println("Incorrect animal specification. Setting Default to Mouse 1") #Always throw an animal error
               end
               ANIMAL = "Mouse"
               NUMBER = "1"
          end

          genotype_res = findmatch(path, genotype_regex)
          if !isnothing(genotype_res)
               GENOTYPE = genotype_res.Genotype
          elseif unknown_genotype == :UNKNOWN
               #println("Unknown Genotype")
               GENOTYPE = "UNKNOWN"
          elseif unknown_genotype == :WT
               #println("Unknown Genotype")
               GENOTYPE = "WT"
          else
               throw("Unknown genotype")
          end

          age_res = findmatch(path, age_regex)
          if !isnothing(age_res)
               AGE = "P$(age_res.Age)"
          elseif !isnothing(age_res) && parse(Int, age_res.Age) >= 30
               #println("Age over 30. Defaulting to Adult")
               AGE = "Adult"
          else
               if conditions_error
                    throw("Age information not included. Erroring")
               else
                    #println(path)
                    if verbose
                         println("Age information not included. Skipping: $path")
                    end
                    continue
               end
          end
          new_path = joinpath(restruc_path, "$(YEAR)_$(MONTH)_$(DATE)_ERG$(GENOTYPE)$(AGE)")
          new_path = joinpath(new_path, "$(ANIMAL)$(NUMBER)_$(AGE)_$(GENOTYPE)")
          
          #Find Drugs, BaCl, NoDrugs
          cond_res = findmatch(path, cond_regex)
          if !isnothing(cond_res) 
               if cond_res.Condition == "Drugs"
                    COND = "BaCl_LAP4"
               elseif cond_res.Condition == "NoDrugs" || cond_res.Condition == "No drugs" || cond_res.Condition == "No Drugs"
                    #println("No Drugs")
                    COND = "BaCl"
               else
                    COND = cond_res.Condition
               end
          else
               #println(path)
               if conditions_error
                    throw("Conditions not available")
               else
                    if verbose
                         @warn "Conditions not available"
                    end
                    COND = "UNKNOWN"
               end
          end
          new_path = joinpath(new_path, "$COND")
          #println(new_path)
          #Find Photoreceptor
          pc_res = findmatch(path, pc_regex)
          if !isnothing(pc_res)
               if pc_res.Photoreceptor == "Rods" #No further label is needed
                    PC = "Rods"
               else
                    color_res = findmatch(path, color_regex)
                    PC = "$(pc_res.Photoreceptor)_$(color_res.Color)"
               end
               new_path = joinpath(new_path, "$PC")
          else
               #first we want to find out if   
               #println("This case we are missing a photoreceptor category")
               background_res = findmatch(path, background_regex)
               if !isnothing(background_res)
                    if background_res.Background == "noback"
                         PC = "Rods"
                    elseif background_res.Background == "withback"
                         color_res = findmatch(path, color_regex)
                         if !isnothing(color_res)
                              PC = "Cones_$(color_res.Color)"
                         else
                              println(path)
                              throw("No color information")
                         end
                    else
                         #println(path)
                         #println(background_res.Background)
                         throw("Some other keyword is used")
                    end
                    new_path = joinpath(new_path, "$PC")
               else
                    #throw("Either Protocol or withback/noback is missing")
               end
          end
          
          nd_res = findmatch(path, nd_regex) #lets try to find the ND filter settings
          #println(nd_res)
          if !isnothing(nd_res)
               ND = nd_res.ND
          else
               if verbose
                    @warn ("ND information is missing. Setting default to ND0")
               end
               ND = 0
          end

          percent_res = findmatch(path, percent_regex)
          #println(percent_res)
          if !isnothing(percent_res)
               PERCENT = percent_res.Percent
          else
               if verbose
                    @warn ("Percent information is missing. Setting Default to 1%")
               end
               PERCENT = 1
          end

          flash_id = findmatch(path, r"\d") #find just a single digit
          if !isempty(flash_id) #This will be necessary whenever 
               println(flash_id)
          end
          
          corr_name = findmatch(path, avg_regex) # find the word "Average or average
          nd_file_res = findmatch(path, nd_file_regex) #or find the nd_file description (plus .abf)
          #if the script finds either a average, or a nd filter and a percent. 
          if !isnothing(corr_name) || !isnothing(nd_file_res)
               print("Checking if file exists:  ")
               count_str = string(count)
               #println(count_str)
               count += 1
               if count > 9999
                    count = 0
               end

               if length(count_str) < 4
                    addN = 4 - length(count_str)
               else
                    addN = 0
               end
               str_lid = "$(join(repeat("0", addN)))$count_str"
               copy_path = joinpath(new_path, "nd$(ND)_$(PERCENT)p_$str_lid.abf")

               if copy_path ∉ current_paths
                    #if true
                    print("NO! New file: $copy_path")
                    #end
                    if mode == :save
                         mkpath(new_path) #Disable this for the time being
                    end
               else#if verbose
                    print("YES! Ignore: ")
                    println(copy_path)
                    #println(new_path)
               end
               #println(!ispath(copy_path))
               if !ispath(copy_path)
                    println("Copying file to path: $path")
                    if mode == :save
                         cp(path, copy_path) #Disable these for the time being
                    end
               else#if verbose
                    println("Path already copied: $path")
               end
          end
     end
end