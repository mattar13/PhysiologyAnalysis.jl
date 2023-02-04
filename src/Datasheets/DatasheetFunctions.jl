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

#This file contains things like extraction and convienance functions

"""
This function converts a dataframe of Any to one matching each row type. 
     catchNaN allows it to catch NaN errors from excel
"""
function safe_convert(dataframe::DataFrame)
     new_obj = DataFrame(dataframe)
     for (idx, col) in enumerate(eachcol(dataframe))
          println(names(dataframe)[idx])
          typ = typeof(col[1]) #Check if there are 
          println(typ)
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

function parseColumn!(T::Type, dataframe::DataFrame, col::Symbol) 
     if all(isa.(dataframe[!, col], T))
          println("Already converted")
     else
          dataframe[! , col] = parse.(T, dataframe[:, col])
     end
end

parseColumn!(dataframe::DataFrame, col::Symbol) = parseColumn!(Int64, dataframe, col)

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