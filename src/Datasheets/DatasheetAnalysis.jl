#This file contains the calibration data

function extract_categories(cond_df::DataFrame) 
     categories = []
     for row in eachrow(cond_df)
          a = [row.Age, row.Genotype, row.Photoreceptor, row.Wavelength]
          push!(categories, a)
     end
     categories
end

function extractIR(trace_datafile::DataFrame, category;
          measure = :Response,
          k = 1000.0, n = 2.0,
          rmin = -1000.0, rmax = 2000.0, #This is the maximum ERG response we have gotten
          kmin = 1.0, kmax = 10e6, 
          nmin = 1.0, nmax = 10.0     
     )

     allIR = trace_datafile |> 
          @filter(_.Age == category[1]) |>
          @filter(_.Genotype == category[2]) |>
          @filter(_.Photoreceptor == category[3]) |>
          @filter(_.Wavelength == category[4]) |>
     DataFrame

     intensity = allIR.Photons
     response = allIR[!, measure]
     #generate the fits for the equation
     r = measure == :Minima ? maximum(response) : minimum(response)
     p0 = [r, k, n]
     ub = [rmax, kmax, nmax]
     lb = [rmin, kmin, nmin]
     println(response)
     model(I, p) = map(i -> p[1]*IR(i, p[2], p[3]), I)
     fit = curve_fit(model, intensity, response, p0, lower = lb, upper = ub)
     fitResponse = model(intensity, fit.param)
     allIR[!, :FitVariable] = fitResponse
     allIR[!, :Residual] = (intensity.-fitResponse).^2

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
function run_A_wave_analysis(all_files::DataFrame; run_amp=false, verbose=true)
     a_files = all_files |> @filter(_.Condition == "BaCl_LAP4") |> DataFrame
     a_files[!, :Path] = string.(a_files[!, :Path])
     uniqueData = a_files |> @unique({_.Year, _.Month, _.Date, _.Number, _.Wavelength, _.Photoreceptor, _.Genotype}) |> DataFrame
     if verbose
          println("Completed data query")
     end

     qTrace = DataFrame()
     qExperiment = DataFrame()
     for (idx, i) in enumerate(eachrow(uniqueData)) #We can walk through each experiment and extract the experiments based on that
          qData = a_files |> @filter(
                       (_.Year, _.Month, _.Date, _.Number, _.Wavelength, _.Photoreceptor, _.Genotype) ==
                       (i.Year, i.Month, i.Date, i.Number, i.Wavelength, i.Photoreceptor, i.Genotype)
                  ) |>
                  DataFrame
          dataFile = readABF(qData.Path)

          if verbose
               println("Completeing A-wave analysis for $idx out of $(size(uniqueData,1))")
          end
          for data in eachchannel(dataFile) #walk through each row of the data iterator
               age = qData.Age[1] #Extract the age
               ch = data.chNames[1] #Extract channel information
               gain = data.chTelegraph[1]
               #Calculate the response based on the age
               #println(data |> size)
               if gain == 1
                    data / 100.0
               end
               #======================DATA ANALYSIS========================#
               if age <= 11
                    filt_data = data_filter(data, t_post=0.5)
                    responses = minimas = minimum(filt_data, dims=2)[:, 1, :]
                    maximas = maximum(filt_data, dims=2)[:, 1, :]
                    Resps = abs.(responses)
                    rmaxes = minimum(responses, dims=1)
               else
                    filt_data = data_filter(data, t_post=1.0)
                    responses = saturated_response(filt_data)
                    minimas = minimum(filt_data, dims=2)[:, 1, :]
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
               #rdim_idxs = 
               #println(Resps)
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

     return qTrace, qExperiment, qConditions
end

#We can update this with our updated channel analysis
function run_B_wave_analysis(all_files::DataFrame; verbose=false)
     trace_A = all_files |> @filter(_.Condition == "BaCl_LAP4" || _.Condition == "LAP4_BaCl") |> DataFrame
     trace_AB = all_files |> @filter(_.Condition == "BaCl") |> DataFrame
     b_files = trace_A |> @join(trace_AB,
                    {_.Year, _.Month, _.Date, _.Number, _.Photons, _.Wavelength, _.Photoreceptor, _.Genotype},
                    {_.Year, _.Month, _.Date, _.Number, _.Photons, _.Wavelength, _.Photoreceptor, _.Genotype},
                    {__...,
                         A_condition = _.Condition,
                         A_Path = _.Path,
                    }) |> DataFrame
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
          filt_data_AB = data_filter(data_AB, t_post=0.5)

          data_A = readABF(qData.A_Path)
          filt_data_A = data_filter(data_A, t_post=0.5)
          #if we want to subtract we need to filter first
          sub_data = filt_data_AB - filt_data_A
          if verbose
               println("Subtraction complete")
               println("Checking analysis")
               print("Data Section: ")
               println(qData)
               print("Data section size:")
               println(qData |> size)
               print("Size of filtered AB data:")
               println(size(data_AB))
               print("Size of filtered A data:")
               println(size(data_A))
               print("Size of subtracted data:")
               println(sub_data |> size)
          end
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
               #println(data_ch)
               #println(size(data_ch))
               #println(maximum(data_ch))
               #println(minimum(data_ch))
               rec_res = recovery_time_constant(data_ch, Resps)
               Recovery_Taus = rec_res[1] |> vec
               Tau_GOFs = rec_res[2] |> vec

               #We need to program the amplification as well. But that may be longer
               Percent_Recoveries = percent_recovery_interval(data_ch, rmaxes)
               #println(Percent_Recoveries)
               #======================GLUING TOGETHER THE QUERY========================#
               #now we can walk through each one of the responses and add it to the qTrace 
               for swp = 1:size(data_ch, 1) #Walk through all sweep info, since sweeps will represent individual datafiles most of the time
                    #inside_row = qData[idx, :Path] #We can break down each individual subsection by the inside row
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
     println(size(uniqueData))
     println(size(trace_ABG))
     println(size(trace_AB))
     println(size(g_files))
     qTrace = DataFrame()
     qExperiment = DataFrame()
     for (idx, i) in enumerate(eachrow(uniqueData)) #We ca
          qData = g_files |> @filter(
                       (_.Year, _.Month, _.Date, _.Number, _.Wavelength, _.Photoreceptor, _.Genotype) ==
                       (i.Year, i.Month, i.Date, i.Number, i.Wavelength, i.Photoreceptor, i.Genotype)
                  ) |>
                  DataFrame
          data_ABG = readABF(qData.Path)
          filt_data_ABG = data_filter(data_ABG, t_post=0.5)
          data_AB = readABF(qData.AB_Path)
          filt_data_AB = data_filter(data_AB, t_post=0.5)
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
          catch #the sheet is not made and must be created
               println("Adding sheets")
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
          catch #the sheet is not made and must be created
               println("Adding sheets")
               XLSX.addsheet!(xf, "experiments_$(append)")
          end
          XLSX.writetable!(xf["experiments_$(append)"],
               collect(DataFrames.eachcol(experiments)),
               DataFrames.names(experiments))
     end

     XLSX.openxlsx(save_file, mode="rw") do xf
          try
               sheet = xf["conditions_$(append)"]
          catch
               println("Adding sheets")
               XLSX.addsheet!(xf, "conditions_$(append)")
          end
          XLSX.writetable!(xf["conditions_$(append)"],
               collect(DataFrames.eachcol(conditions)),
               DataFrames.names(conditions))
     end
end

function runAnalysis(datafile::String)
     print("Opening datafile $(datafile)... ")
     df = openDatasheet(datafile)
     println("complete")
     #%% Test the a, b, and g wave analysis
     resA = ePhys.run_A_wave_analysis(df)
     resB = ePhys.run_B_wave_analysis(df)
     resG = ePhys.run_G_wave_analysis(df)
     add_analysis_sheets(resA, datafile; append="A")
     add_analysis_sheets(resB, datafile; append="B")
     add_analysis_sheets(resG, datafile; append="G")
     return (df, resA, resB, resG)
end


