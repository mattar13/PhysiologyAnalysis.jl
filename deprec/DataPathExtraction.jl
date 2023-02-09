"""
Lets shift to using Regex to extract file information
"""
function OLDDataPathExtraction(path::String, calibration_file::String;
     verbose=false,
     extract_photons=true
)
     #The last line of the string should be the nd information
     if verbose
          println(path)
     end
     path_array = splitpath(path) #first we can split the file path into seperate strings
     nt = (Path=path,)
     #println(path_array)
     date_res = findmatch(path_array, date_regex; verbose=verbose)
     if !isnothing(date_res)
          nt = merge(nt, date_res)
     end

     animal_res = findmatch(path_array, animal_regex; verbose=verbose)
     if !isnothing(animal_res)
          nt = merge(nt, animal_res)
     end

     genotype_res = findmatch(path, genotype_regex; verbose = verbose)
     if !isnothing(genotype_res)
          nt = merge(nt, genotype_res)
     end

     age_res = findmatch(path, age_regex; verbose = verbose)
     if !isnothing(age_res) || parse(Int, age_res.Age) > 30
          if age_res.Age == "Adult"
               nt = merge(nt, ("Age" => "30"))
          else
               nt = merge(nt, age_res)
          end
     end
     
     #now lets look for a condition in the possible conds
     cond_res = findmatch(path, cond_regex)
     if cond_res.Condition == "Drugs"
          cond = "BaCl_LAP4"
     elseif cond_res.Condition == "NoDrugs" || cond_res.Condition == "No drugs"
          cond = "BaCl"
     else
          cond = cond_res.Condition
     end
     nt = merge(nt, (Condition=cond,))
     
     pc_res = findmatch(path, pc_regex)
     if !isnothing(pc_res)
          nt = merge(nt, pc_res)
          if pc_res.Photoreceptors == "Rods" #No further label is needed
               nt = merge(nt, (Wavelength = "520",))
          else
               color_res = findmatch(path, color_regex)
               if color_res.Color == "Blue" || color_res.Color == "365" || color_res.Color == "365UV" 
                    nt = merge(nt, (Wavelength = "365",))
               elseif  color_res.Color == "Green" || color_res.Color == "525" || color_res.Color == "525Green" 
                    nt = merge(nt, (Wavelength="520",))
               elseif color_res.Color == "520" || color_res.Color == "520Green" 
                    nt = merge(nt, (Wavelength="520",))
               else
                    println(color_res)
               end
          end
     end

     nd_res = findmatch(path, nd_regex; verbose = verbose)
     if !isnothing(nd_res)
          nt = merge(nt, nd_res)
     end
     percent_res = findmatch(path, percent_regex; verbose = verbose)
     if !isnothing(nd_res)
          nt = merge(nt, percent_res)
     end

     if extract_photons
          #extract the stimulus from the data
          stim_timestamps = extract_stimulus(path)[1].timestamps
          stim_time = round(Int64, (stim_timestamps[2] - stim_timestamps[1]) * 1000)
          #println(stim_time)
          nt = merge(nt, (Stim_Time=stim_time,))
          #now lets extract photons
          if nd_res.ND == "0.5"
               #println(intensity_info.Percent)
               println(nd_res)
               println(nt.Wavelength)
               photons = photon_lookup(
                    parse(Int64, nt.Wavelength),
                    0,
                    parse(Int64, nd_res.Percent),
                    calibration_file
               )
               println(photons)
               photons = (photons*stim_time) / (10^0.5)
               nt = merge(nt, (Photons=photons,))
          else
               #println(intensity_info.Percent)
               photons = photon_lookup(
                    parse(Int64, nt.Wavelength),
                    parse(Int64, nd_res.ND),
                    parse(Int64, percent_res.Percent),
                    calibration_file
               ) .* stim_time
               nt = merge(nt, (Photons=photons,))
          end
     end
     #println(nt)
     #return nt |> clean_nt_numbers
     return nt |> parseNamedTuple
end
