#This file contains the calibration data
const calibration_file = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\Calibrations\photon_lookup.xlsx"

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
          println(nt_vals)
          println(idx)
          println(nt_vals[idx])
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
     println(nt.Age)
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