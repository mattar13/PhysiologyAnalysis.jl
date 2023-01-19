const calibration_file = raw"C:\Users\mtarc\The University of Akron\Renna Lab - General\Data\Calibrations\photon_lookup.xlsx"

"""
    photon_lookup(wavelength, nd, percent, stim_time, calibration_file[,sheet_name])

Uses the calibration file or datasheet to look up the photon density. The Photon datasheet should be either 
"""
function photon_lookup(wavelength::Real, nd::Real, percent::Real, calibration_file::String; sheet_name::String="Current_Test")
     df = DataFrame(XLSX.readtable(calibration_file, sheet_name))
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

     #Find photoreceptors
     pc_res = findmatch(path, pc_regex)
     if !isnothing(pc_res)
          if pc_res.Photoreceptors == "Rods" #No further label is needed
               PC = "Rods"
               WAVE = 520 #Shouls actually be 498
          else
               PC = "Cones"
               color_res = findmatch(path, color_regex)
               if color_res.Color == "Blue" || color_res.Color == "365" || color_res.Color == "365UV" || color_res.Color == "blue"
                    WAVE = 365
               elseif  color_res.Color == "green" || color_res.Color == "525" || color_res.Color == "525Green" || color_res.Color == "Green"
                    WAVE = 520
               elseif color_res.Color == "520" || color_res.Color == "520Green" 
                    WAVE = 520
               else
                    println(color_res)
               end
          end
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
          println(flash_id)
     end

     if extract_photons
          #extract the stimulus from the data
          stim_timestamps = extract_stimulus(path)[1].timestamps
          stim_time = round(Int64, (stim_timestamps[2] - stim_timestamps[1]) * 1000)
          #println(stim_time)
          STIM_TIME=stim_time
          #now lets extract photons
          if ND == "0.5"
               PHOTONS = photon_lookup(
                    WAVE,
                    0,
                    PERCENT,
                    calibration_file
               )
               println(PHOTONS)
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

end

DataPathExtraction(path::String; kwargs...) = DataPathExtraction(path, calibration_file; kwargs...)

"""
This function cleans the data out of a dataframe if the dataframe is already open
"""
function cleanDatasheet!(xf::XLSX.XLSXFile, sheetname::String)
     if sheetname ∈ XLSX.sheetnames(xf)
          sheet = xf[sheetname]
          nrows, ncols = size(sheet[:])
          names = map(col -> sheet[1, col], 1:ncols)
          eraser = []
          for name in names
               eraser_col = (name, fill("", nrows - 1)...)
               push!(eraser, eraser_col)
          end
          XLSX.writetable!(xf[sheetname],
               fill(Tuple(fill("", nrows)), ncols), names
          )
     else
          println("Sheetname not in sheets")
          println("Choose from one of these sheets:")
          for sn in XLSX.sheetnames(xf)
               println("-> $sn")
          end
     end
end

"""
This function cleans the data out of a dataframe and saves it
"""
function cleanDatasheet!(filename::String, sheetname::String)
     XLSX.openxlsx(filename, mode="rw") do xf
          cleanDataFrame!(xf, sheetname)
     end
end


"""
This function converts a dataframe of Any to one matching each row type. 
     catchNaN allows it to catch NaN errors from excel
"""
function safe_convert(dataframe::DataFrame)
     new_obj = DataFrame(dataframe)
     for (idx, col) in enumerate(eachcol(dataframe))
          #println(names(dataframe)[idx])
          typ = typeof(col[1]) #Check if there are 
          #We will try to convert each row. If it does not work, we can remove the NaN
          #println(col)
          if ("NaN" ∈ col) #Check if there exists a word NaN in the row (excel will call these strings)
               #print("Is NaN") #debugging statements
               whereNaN = findall(col .== "NaN")
               #println("At position $whereNaN")
               for idxNaN in whereNaN
                    #println(idxNaN)
                    col[idxNaN] = NaN #Instead use a NaN Floating point objects
               end
               new_obj[:, idx] = convert.(typ, col)
          elseif !all(isa.(col, typ))#if col[1] #This is for if there is a Int to Float64 error
               whereNotSame = findall(map(!, isa.(col, typ)))
               irregular_type = col[whereNotSame[1]] |> typeof
               #println("Column type: $typ")
               #println("Irregular type: $(irregular_type)")
               if irregular_type == Int64 && typ == Float64 #This can cause an error
                    new_obj[!, idx] = convert.(typ, col) #Convert all values to Float64
               else
                    new_obj[!, idx] = convert.(irregular_type, col) #Convert all values to Float64
               end
          else
               #conv = convert.(typ, col)
               new_obj[!, idx] = convert.(typ, col)
          end
     end
     return new_obj
end

"""
This function creates a new datasheet
"""
function createDatasheet(all_files::Vector{String}; filename="data_analysis.xlsx", verbose = false)
     dataframe = DataFrame()
     for (idx, file) in enumerate(all_files)
          if verbose
               println("Analyzing file $idx of $(size(all_files, 1)): $file")
          end
          entry = DataPathExtraction(file)
          if isnothing(entry) #Throw this in the case that the entry cannot be fit
               println(entry)
               #elseif length(entry) != size(dataframe, 2)
          #     println("Entry does not match dataframe size. Probably an extra category")
          else
               try
                    push!(dataframe, entry)
               catch
                    println(entry.ND)
                         println("Probably the ND filter. I don't know how to fix that")
                    end
               end
          end
     if !isnothing(filename)
          XLSX.openxlsx(filename, mode="w") do xf
               XLSX.rename!(xf[1], "All_Files") #Rename sheet 1
               try
                    sheet = xf["All_Files"] #Try opening the All_Files
               catch
                    println("Adding sheets")
                    XLSX.addsheet!(xf, "All_Files")
               end
               XLSX.writetable!(xf["All_Files"],
                    collect(DataFrames.eachcol(dataframe)),
                    DataFrames.names(dataframe))
          end
     end
     dataframe
end

"""
This function opens an old datasheet
"""

function openDatasheet(data_file::String; sheetName::String="All_Files", typeConvert=true)
     xf = readxlsx(data_file)
     if sheetName == "all"
          sheetnames = XLSX.sheetnames(xf)
          df_set = Dict()
          for sn in sheetnames
               #println(sn) #Use this to debug 
               df_set[sn] = openDatasheet(data_file, sheetName=sn)
          end
          return df_set
     else
          s = xf[sheetName]
          df = XLSX.eachtablerow(s) |> DataFrame
          #We can walk through and try to convert each row to either an integer, Float, or String
          if typeConvert
               df = safe_convert(df) #This converts the categories to a type in the first position
          end
          return df
     end
end

"""
This function will read and update the datasheet with the new files
"""
function updateDatasheet(data_file::String, all_files::Vector{String}; reset::Bool=false, savefile::Bool=true)
     if reset
          #If this is selected, completely reset the analysis
     else
          df = openDatasheet(data_file) #First, open the old datasheet
          nrows, ncols = size(df)

          println("Searching for files that need to be added and removed")
          old_files = df[:, :Path] #Pull out a list of old files 

          #This searches for files that occur in all files, but not in old files, indicating they need added to the analysis
          files_to_add = findall(isnothing, indexin(all_files, old_files))

          #This searches for files that are in old files, but not not in all files, indicating they may need deleted
          files_to_remove = findall(isnothing, indexin(old_files, all_files))

          duplicate_files = findall(nonunique(df))

          if !isempty(files_to_add) #There are no files to add
               println("Adding $(length(files_to_add)) files")
               new_files = all_files[files_to_add]
               df_new = createDatasheet(new_files)
               df = vcat(df, df_new) #Add the new datasheet to the old one
          end

          if !isempty(files_to_remove)
               println("Removing $(length(files_to_add)) files")
               deleteat!(df, files_to_remove)
          end

          if !isempty(duplicate_files)
               println("Removing $(length(duplicate_files)) duplicated files")
               deleteat!(df, duplicate_files)
          end

          #Sort the dataframe
          df = df |>
               @orderby(_.Year) |> @thenby(_.Month) |> @thenby(_.Date) |>
               @thenby(_.Animal) |> @thenby(_.Genotype) |> @thenby(_.Condition) |>
               @thenby(_.Wavelength) |> @thenby(_.Photons) |>
               DataFrame

          if savefile
               
               println("Saving file... ")
               XLSX.openxlsx(data_file, mode="rw") do xf
                    print("Erasing file")
                    cleanDatasheet!(xf, "All_Files") #This cleans all data from All_Files
                    #This re-writes it
                    XLSX.writetable!(xf["All_Files"],
                         collect(DataFrames.eachcol(df)),
                         DataFrames.names(df)
                    )

               end
               println("Complete")
          end

          return df
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


#=============================================================================================
These fuctions help extract experiments
=============================================================================================#
function matchExperiment(trace::DataFrame, date::Tuple{Int64,Int64,Int64,Int64}; pc="Rods", color=520)
     result = trace |>
              @filter((_.Year, _.Month, _.Date, _.Number) == date) |>
              DataFrame

     if !isnothing(pc)
          result = result |> @filter(_.Photoreceptor == pc) |> DataFrame
     end
     result = result |>
              @filter(_.Wavelength == color) |>
              @orderby(_.Photons) |>
              DataFrame

     return result
end

function matchExperiment(trace::DataFrame, date::Tuple{Int64,Int64,Int64,Int64,String}; pc="Rods", color=520)

     result = trace |>
              @filter((_.Year, _.Month, _.Date, _.Number, _.Channel) == date) |>
              DataFrame
     if !isnothing(pc)
          result = result |> @filter(_.Photoreceptor == pc) |> DataFrame
     end
     result = result |>
              @filter(_.Wavelength == color) |>
              @orderby(_.Photons) |>
              DataFrame
     return result
end

function matchExperiment(trace::DataFrame, date::Tuple{Int64,Int64,Int64,Int64,String,Int64}; pc="Rods")

     result = trace |>
              @filter((_.Year, _.Month, _.Date, _.Number, _.Channel, _.Wavelength) == date) |>
              DataFrame
     result = result |>
              #@filter(_.Wavelength == color) |>
              @orderby(_.Photons) |>
              DataFrame
     return result
end