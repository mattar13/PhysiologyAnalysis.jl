"""
    photon_lookup(wavelength, nd, percent, stim_time, calibration_file[,sheet_name])

Uses the calibration file or datasheet to look up the photon density. The Photon datasheet should be either 
"""
function photon_lookup(wavelength::Real, nd::Real, percent::Real, calibration_file::String, sheet_name::String = "Current_Test")
    df = DataFrame(XLSX.readtable(calibration_file, sheet_name)...)
    Qi = df |>
         @filter(_.Wavelength == wavelength) |>
         @filter(_.ND == nd) |>
         @filter(_.Intensity == percent) |>
         #@filter(_.stim_time == stim_time) |>
         @map(_.Photons) |>
         DataFrame
    #%%
    if size(Qi, 1) != 0
        #Only return f an entry exists
        return Qi.value[1]
    end
end

function createDatasheet(all_files::Vector{String}; )
     dataframe = DataFrame()
     for (idx, file) in enumerate(all_files)
          println("Analyzing file $idx of $(size(all_files, 1)): $file")
          entry = DataPathExtraction(file)
          push!(dataframe, entry)
     end
     dataframe
end

function run_A_wave_analysis(all_files::DataFrame; run_amp=false, verbose=false)
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
               println(age)
               ch = data.chNames[1] #Extract channel information
               gain = data.chTelegraph[1]
               #Calculate the response based on the age

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
                    filt_data = filter_data(data, t_post=1.0)
                    responses = saturated_response(filt_data)
                    minimas = minimum(filt_data, dims=2)[:, 1, :]
                    maximas = maximum(filt_data, dims=2)[:, 1, :]
                    Resps = abs.(responses)
                    rmaxes = minimum(responses, dims=1)
               end

               Peak_Times = time_to_peak(filt_data)
               Integrated_Times = abs.(integral(filt_data))
               rec_res = recovery_tau(filt_data, responses)
               Recovery_Taus = rec_res[1] |> vec
               Tau_GOFs = rec_res[2] |> vec

               #We need to program the amplification as well. But that may be longer

               #The 60% Bandwith can be used to measure the percent recovery interval
               #println(size(rmaxes))
               Percent_Recoveries = percent_recovery_interval(filt_data, rmaxes)
               #println(Percent_Recoveries)
               #println(Percent_Recoveries |> size)
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