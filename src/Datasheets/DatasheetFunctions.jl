const calibration_file = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\Calibrations\photon_lookup.xlsx"


"""
    photon_lookup(wavelength, nd, percent, stim_time, calibration_file[,sheet_name])

Uses the calibration file or datasheet to look up the photon density. The Photon datasheet should be either 
"""
function photon_lookup(wavelength::Real, nd::Real, percent::Real, calibration_file::String, sheet_name::String="Current_Test")
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

#========================================================================================================================================

These functions extract the file information from the datafile path

========================================================================================================================================#

"""
Lets shift to using Regex to extract file information
"""
function DataPathExtraction(path::String, calibration_file::String;
     verbose=false,
     extract_photons=true
)
     #The last line of the string should be the nd information
     if verbose
          println(path)
     end
     path_array = splitpath(path) #first we can split the file path into seperate strings
     nt = (Path=path,)
     date_info = findmatch(path_array, date_regex; verbose=verbose)
     if !isnothing(date_info)
          nt = merge(nt, date_info)
     end

     animal_info = findmatch(path_array, animal_regex; verbose=verbose)
     if !isnothing(animal_info)
          #try to match the age Regex
          nt_keys = keys(animal_info)
          nt_vals = [values(animal_info)...]
          idx = findall(nt_keys .== :Age)[1]
          if nt_vals[idx] == "Adult"
               nt_vals[idx] = "30"
          else
               nt_vals[idx] = filter_string(Int64, nt_vals[idx])
          end
          #convert values like "P8" to 8, and Adult to 30

          animal_info = NamedTuple{nt_keys}(nt_vals)
          nt = merge(nt, animal_info)
     end

     #now lets look for a condition in the possible conds
     cond = find_condition(path_array; possible_conds=["BaCl", "BaCl_LAP4", "NoDrugs"])
     nt = merge(nt, (Condition=cond,))

     pc = find_condition(path_array; possible_conds=["Rods", "Cones"])
     if !isnothing(pc)
          nt = merge(nt, (Photoreceptor=pc,))
     else
          nt = merge(nt, (Photoreceptor="Rods",))
     end

     wv = find_condition(path_array; possible_conds=["365", "365UV", "520", "520Green", "525", "525Green"])
     if !isnothing(wv)
          wv_value = filter_string(Int64, wv)
          if wv_value == "525" #This is a file error
               println("Old file name made a mistake")
               nt = merge(nt, (Wavelength="520",))
          else
               nt = merge(nt, (Wavelength=wv_value,))
          end
     elseif any(path_array .== "UV")
          nt = merge(nt, (Wavelength="365",))
     elseif any(path_array .== "Green") || nt.Photoreceptor == "Rods"
          nt = merge(nt, (Wavelength="520",))
     end

     intensity_info = findmatch(path_array, nd_file_regex; verbose=verbose) #find and extract the intensity info
     if !isnothing(intensity_info)
          nt = merge(nt, intensity_info)
          if extract_photons
               #extract the stimulus from the data
               stim_timestamps = extract_stimulus(path)[1].timestamps
               stim_time = round(Int64, (stim_timestamps[2] - stim_timestamps[1]) * 1000)
               #println(stim_time)
               nt = merge(nt, (Stim_Time=stim_time,))
               #now lets extract photons
               #println(nt.Wavelength)
               #println(intensity_info.ND)
               #println(intensity_info.Percent)
               photons = photon_lookup(
                    parse(Int64, nt.Wavelength),
                    parse(Int64, intensity_info.ND),
                    parse(Int64, intensity_info.Percent),
                    calibration_file
               ) .* stim_time
               nt = merge(nt, (Photons=photons,))
          end
          #now we can look for any date info
     end

     #return nt |> clean_nt_numbers
     return nt |> parseNamedTuple
end

DataPathExtraction(path::String; kwargs...) = DataPathExtraction(path, calibration_file; kwargs...)


function createDatasheet(all_files::Vector{String};)
     dataframe = DataFrame()
     for (idx, file) in enumerate(all_files)
          println("Analyzing file $idx of $(size(all_files, 1)): $file")
          entry = DataPathExtraction(file)
          if isnothing(entry) #Throw this in the case that the entry cannot be fit
               println(entry)
          else
               push!(dataframe, entry)
          end
     end
     dataframe
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

#==========================================================================================
These functions can open data from the dataframes
==========================================================================================#

"""
This code is actually an extension of the ABFReader package and somewhat specific for my own needs.
     I will be some way of making this more general later

"""
function readABF(df::DataFrame; extra_channels=nothing, a_name="A_Path", ab_name="AB_Path", kwargs...)
     df_names = names(df)
     #Check to make sure path is in the dataframe
     #Check to make sure the dataframe contains channel info
     @assert "Channel" ∈ df_names
     if (a_name ∈ df_names) #&& ("AB_Path" ∈ df_names) #This is a B subtraction
          #println("B wave subtraction")
          A_paths = string.(df.A_Path)
          AB_paths = string.(df.Path)
          ch = (df.Channel |> unique) .|> String
          if !isnothing(extra_channels)
               ch = (vcat(ch..., extra_channels...))
          end
          A_data = readABF(A_paths; channels=ch, kwargs...)
          AB_data = readABF(AB_paths; channels=ch, kwargs...)
          return A_data, AB_data
     elseif (ab_name ∈ df_names) #&& ("ABG_Path" ∈ df_names) #This is the G-wave subtraction
          #println("G wave subtraction")
          AB_paths = string.(df.AB_Path)
          ABG_paths = string.(df.Path)
          ch = (df.Channel |> unique) .|> String
          if !isnothing(extra_channels)
               ch = (vcat(ch..., extra_channels...))
          end
          AB_data = readABF(AB_paths, channels=ch, kwargs...)
          ABG_data = readABF(ABG_paths, channels=ch, kwargs...)

          return AB_data, ABG_data
     elseif ("Path" ∈ df_names) #This is just the A-wave
          paths = string.(df.Path)
          ch = (df.Channel |> unique) .|> String
          if !isnothing(extra_channels)
               ch = (vcat(ch..., extra_channels...))
          end
          data = readABF(paths, channels=ch, kwargs...)
          return data
     else
          throw("There is no path key")
     end


end

readABF(df_row::DataFrameRow; kwargs...) = readABF(df_row |> DataFrame; kwargs...)