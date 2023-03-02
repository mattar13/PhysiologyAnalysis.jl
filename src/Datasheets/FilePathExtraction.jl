#========================================================================================================================================

These functions extract the file information from the datafile path

========================================================================================================================================#
function DataPathExtraction(path::String, calibration_file::String;
     unknown_genotype = :WT, conditions_error = false, 
     verbose = false, extract_photons=true
)

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
     adult_res = findmatch(path, adult_regex) #This only checks if the word adult exists
     if !isnothing(age_res)
          AGE = "P$(age_res.Age)"
     elseif !isnothing(age_res) && parse(Int, age_res.Age) >= 30
          #println("Age over 30. Defaulting to Adult")
          AGE = "Adult"
     elseif !isnothing(adult_res)
          AGE = "Adult"
     else
          if conditions_error
               throw("Age information not included. Erroring")
          else
               #println(path)
               if verbose
                    println("Age information not included. Skipping: $path")
               end
               AGE = nothing
          end
     end

     #Find Drugs, BaCl, NoDrugs
     cond_res = findmatch(path, cond_regex)
     if !isnothing(cond_res)
          if cond_res.Condition == "Drugs"
               COND = "BaCl_LAP4"
          elseif cond_res.Condition == "NoDrugs" || cond_res.Condition == "No drugs" || cond_res.Condition == "No Drugs"
               #println("No Drugs")
               COND = "NoDrugs"
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

     #Find color of the light stim used
     color_res = findmatch(path, color_regex)
     if !isnothing(color_res) #This 
          if color_res.Color == "Blue" || color_res.Color == "365" || color_res.Color == "365UV" || color_res.Color == "blue"
               WAVE = 365
          elseif  color_res.Color == "green" || color_res.Color == "525" || color_res.Color == "525Green" || color_res.Color == "Green"
               WAVE = 520
          elseif color_res.Color == "520" || color_res.Color == "520Green" 
               WAVE = 520
          end
     end
     
     background_res = findmatch(path, background_regex)
     pc_res = findmatch(path, pc_regex)
     #println(pc_res)
     if !isnothing(pc_res)
          if pc_res.Photoreceptors == "Rods" #No further label is needed
               PC = "Rods"
               if isnothing(color_res) #only do this if wave is empty, but it is rods
                    WAVE = 520 #Should actually be 498
               end
          elseif pc_res.Photoreceptors == "Cones"
               PC = "Cones"          
          elseif !isnothing(background_res)
               if background_res.Background == "noback"
                    PC = "Rods"
                    WAVE = 520
               elseif background_res.Background == "withback" && !isnothing(color_res)
                    PC = "Cones"
               elseif isnothing(color_res)
                    println(path)
                    throw("No color information")
               end
          elseif !isnothing(color_res)
               throw("Either Protocol or withback/noback is missing")
          end
     elseif GENOTYPE == "GNAT-KO" #This is specific to my dataset. I don't like doing this
          #println("Here")
          PC = "Cones"
     end

     nd_res = findmatch(path, nd_regex) #lets try to find the ND filter settings
     #println(nd_res)
     if !isnothing(nd_res)
          ND = parse(Int, nd_res.ND)
     else
          if verbose
               @warn ("ND information is missing. Setting default to ND0")
          end
          ND = 0
     end

     percent_res = findmatch(path, percent_regex)
     #println(percent_res)
     if !isnothing(percent_res)
          PERCENT = parse(Int, percent_res.Percent)
     else
          if verbose
               @warn ("Percent information is missing. Setting Default to 1%")
          end
          PERCENT = 1
     end

     flash_id = findmatch(path, r"\d") #find just a single digit
     if !isempty(flash_id) #This will be necessary whenever 
          #println(flash_id)
     end
     corr_name = findmatch(path, avg_regex) # find the word "Average or average
     nd_file_res = findmatch(path, nd_file_regex) #or find the nd_file description (plus .abf)
     if !isnothing(corr_name) || !isnothing(nd_file_res)
          #println("Recording")
          if extract_photons
               #extract the stimulus from the data
               stim_timestamps = extractStimulus(path)[1]
               stim_time = round(Int64, (stim_timestamps[2] - stim_timestamps[1]) * 1000)
               STIM_TIME=stim_time
               if ND == "0.5"
                    PHOTONS = photon_lookup(
                         WAVE,
                         0,
                         PERCENT,
                         calibration_file
                    )
                    #println(PHOTONS)
                    PHOTONS = (PHOTONS*stim_time) / (10^0.5)
               else
                    PHOTONS = photon_lookup(
                         WAVE,
                         ND,
                         PERCENT,
                         calibration_file
                    ) .* stim_time
               end
          end
          nt = (
               Path = path, 
               Year = YEAR, Month = MONTH, Date = DATE, 
               Animal = ANIMAL, Number = NUMBER, Age = AGE, Genotype = GENOTYPE, 
               Condition = COND, Photoreceptor = PC, Wavelength = WAVE, 
               ND = ND, Percent = PERCENT, Stim_Time = STIM_TIME, 
               Photons = PHOTONS
          )
          return nt
          #corr_name = findmatch(path, avg_regex) # find the word "Average or average
          #nd_file_res = findmatch(path, nd_file_regex) #or find the nd_file description (plus .abf)
          #if the script finds either a average, or a nd filter and a percent. 
     else
          println("Not average nor nd")
     end

end

DataPathExtraction(path::String; kwargs...) = DataPathExtraction(path, calibration_file; kwargs...)

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