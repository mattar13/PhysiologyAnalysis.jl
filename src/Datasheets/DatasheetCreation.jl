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
function createDatasheet(all_files::Vector{String}; filename="data_analysis.xlsx", verbose = false)
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
          if verbose
               println("Save file success")
          end
     end
     dataframe
end

"""
This function opens an old datasheet
"""

function openDatasheet(data_file::String; sheetName::String="all", typeConvert=true)
     xf = readxlsx(data_file)
     if sheetName == "all"
          sheetnames = XLSX.sheetnames(xf)
          df_set = Dict()
          for sn in sheetnames
               #println(sn) #Use this to debug 
               df_set[sn] = openDatasheet(data_file; sheetName=sn, typeConvert = typeConvert)
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
          df = openDatasheet(data_file; sheetName = "All_Files") #First, open the old datasheet
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