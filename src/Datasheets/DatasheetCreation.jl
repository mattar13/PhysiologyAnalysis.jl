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
This function creates a new datasheet
"""
function createDataset(all_files::Vector{String}; verbose = false)
     dataframe = DataFrame()
     for (idx, file) in enumerate(all_files)
          if verbose
               print("Analyzing file $idx of $(size(all_files, 1)): $file ...")
          end
          try
               entry = DataPathExtraction(file)
               if isnothing(entry) && verbose #Throw this in the case that the entry cannot be fit
                    println("Failed")
                    #elseif length(entry) != size(dataframe, 2)
               #     println("Entry does not match dataframe size. Probably an extra category")
               else
                    if verbose
                         println("Success")
                    end
                    push!(dataframe, entry)
               end
          catch error
               println(file)
               println(error)
          end
     end
     dataframe
end

function createDataset(file_root::String; verbose = false, kwargs...)
     all_files = file_root |> parseABF
     ALL_FILES_dataset = createDataset(all_files; verbose = verbose)
     return runTraceAnalysis(ALL_FILES_dataset; kwargs...)
end
"""
This function opens an old datasheet
"""
function openDataset(datafile::String; sheetName::String="all", typeConvert=true)
     if sheetName == "all"
          sheetnames = ["ALL_FILES", "TRACES", "EXPERIMENTS", "CONDITIONS"]
          df_set = Dict{String, DataFrame}()
          for sn in sheetnames
               #println(sn) #Use this to debug 
               df_set[sn] = openDataset(datafile; sheetName=sn, typeConvert = typeConvert)
          end
          return df_set
     else
          #s = xf[sheetName]
          #df = XLSX.eachtablerow(s) |> DataFrame
          df = DataFrame(XLSX.readtable(datafile, sheetName))
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
function updateDataset(data_file::String, all_files::Vector{String}; reset::Bool=false, savefile::Bool=true)
     if reset
          #If this is selected, completely reset the analysis
     else
          df = openDatasheet(data_file; sheetName = "ALL_FILES") #First, open the old datasheet
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

function copyDataset(datafile::String, sheetname = "all", backup = true)
     root = joinpath(splitpath(datafile)[1:end-1])
     if isfile(datafile) #we should only do these things if the datafile exists
          new_fn = "$root\\temp.xlsx"
          XLSX.openxlsx(new_fn, mode = "w") do xf
               sheet1 = xf[1] #Sheet 1 should be renamed
               XLSX.rename!(sheet1, "ALL_FILES")
               XLSX.addsheet!(xf, "TRACES")
               XLSX.addsheet!(xf, "EXPERIMENTS")
               XLSX.addsheet!(xf, "CONDITIONS")
               #XLSX.addsheet!(xf, "STATISTICS")
               #XLSX.addsheet!(xf, "FITTING")
               #XLSX.addsheet!(xf, "SYNAPTIC TRANSFER FUNCTION")
          end
     end
end

"""
This function quicksaves the datafile you are working with
"""
function backupDataset(datafile::String)
     date = now()
     root = joinpath(splitpath(datafile)[1:end-1])
     filename = splitpath(datafile)[end][1:end-5]
     backup_file = "$(root)\\$(year(date))_$(month(date))_$(day(date))_$(filename)_BACKUP_$(hour(date))_$(minute(date))_$(second(date)).xlsx"
     cp(datafile, backup_file)
end

function saveDataset(dataset::Dict{String, DataFrame}, filename::String)
     XLSX.openxlsx(filename, mode = "w") do xf
          sheet_ALL = xf[1] #Sheet 1 should be renamed
          XLSX.rename!(sheet_ALL, "ALL_FILES")
          XLSX.writetable!(sheet_ALL, dataset["ALL_FILES"])

          XLSX.addsheet!(xf, "TRACES")
          sheet_TRACES = xf[2]
          XLSX.writetable!(sheet_TRACES, dataset["TRACES"])

          XLSX.addsheet!(xf, "EXPERIMENTS")
          sheet_EXPERIMENTS = xf[3]
          XLSX.writetable!(sheet_EXPERIMENTS , dataset["EXPERIMENTS"])

          XLSX.addsheet!(xf, "CONDITIONS")
          sheet_CONDITIONS = xf[4]
          XLSX.writetable!(sheet_CONDITIONS , dataset["CONDITIONS"])
     end
end