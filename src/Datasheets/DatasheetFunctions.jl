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
function safe_convert(dataframe::DataFrame; verbose = false)
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
                    if verbose
                         println("Indexes where a NaN is: $idxNaN")
                    end
                    col[idxNaN] = NaN #Instead use a NaN Floating point objects
               end
               new_obj[:, idx] = convert.(typ, col)
          elseif !all(isa.(col, typ))#if col[1] #This is for if there is a Int to Float64 error
               whereNotSame = findall(map(!, isa.(col, typ)))
               irregular_type = col[whereNotSame[1]] |> typeof
               if verbose
                    println(col[whereNotSame[1]])
                    println("Column type: $typ")
                    println("Irregular type: $(irregular_type)")
               end
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
"""
If you pass either a named tuple or a dataframe row, this will pull out all the related 
"""
function matchExperiment(trace::DataFrame, info::NamedTuple)
     return_traces = copy(trace)
     if haskey(info, :Year)
          return_traces = return_traces |> @filter(_.Year == info.Year) |> DataFrame
     end

     if haskey(info, :Month)
          return_traces = return_traces |> @filter(_.Month == info.Month) |> DataFrame
     end

     if haskey(info, :Date)
          return_traces = return_traces |> @filter(_.Date == info.Date) |> DataFrame
     end

     if haskey(info, :Number)
          return_traces = return_traces |> @filter(_.Number == info.Number) |> DataFrame
     end
     if haskey(info,:Photoreceptor)
          return_traces = return_traces |> @filter(_.Photoreceptor == info.Photoreceptor) |> DataFrame
     end
     if haskey(info, :Condition)
          return_traces = return_traces |> @filter(_.Condition == info.Condition) |> DataFrame
     end

     if haskey(info, :Channel)
          return_traces = return_traces |> @filter(_.Channel == info.Channel) |> DataFrame
     end

     if haskey(info, :Genotype)
          return_traces = return_traces |> @filter(_.Genotype == info.Genotype) |> DataFrame
     end

     if haskey(info, :Age)
          return_traces = return_traces |> @filter(_.Age == info.Age) |> DataFrame
     end

     return return_traces
end

matchExperiment(trace::DataFrame, row::DataFrameRow) = matchExperiment(trace, NamedTuple(row))

function matchExperiment(trace::DataFrame, rows::DataFrame)
     return_dataset = DataFrame()
     for row in eachrow(rows)
          dataset_i = matchExperiment(trace, row)
          return_dataset = vcat(return_dataset, dataset_i)
     end
     return return_dataset
end

function excludeExperiment(trace::DataFrame, info::NamedTuple)
     return_traces = copy(trace)
     if haskey(info, :Year)
          return_traces = return_traces |> @filter(_.Year != info.Year) |> DataFrame
     end

     if haskey(info, :Month)
          return_traces = return_traces |> @filter(_.Month != info.Month) |> DataFrame
     end

     if haskey(info, :Date)
          return_traces = return_traces |> @filter(_.Date != info.Date) |> DataFrame
     end

     if haskey(info, :Number)
          return_traces = return_traces |> @filter(_.Number != info.Number) |> DataFrame
     end
     if haskey(info,:Photoreceptor)
          return_traces = return_traces |> @filter(_.Photoreceptor != info.Photoreceptor) |> DataFrame
     end
     if haskey(info, :Condition)
          return_traces = return_traces |> @filter(_.Condition != info.Condition) |> DataFrame
     end

     if haskey(info, :Channel)
          return_traces = return_traces |> @filter(_.Channel != info.Channel) |> DataFrame
     end

     return return_traces
end

excludeExperiment(trace::DataFrame, row::DataFrameRow) = excludeExperiment(trace, NamedTuple(row))

function excludeExperiment(trace::DataFrame, rows::DataFrame)
     return_dataset = DataFrame()
     for row in eachrow(rows)
          dataset_i = excludeExperiment(trace, row)
          return_dataset = vcat(return_dataset, dataset_i)
     end
     return return_dataset
end

function match_excludeExperiment(trace::DataFrame, match_rows, exclude_rows)
     matched = matchExperiment(trace, match_rows)
     excluded = excludeExperiment(matched, exclude_rows)
     return excluded
end

"""
This function matches the experiments in info and then switches the include flag to false
"""
function flagExperiment(trace::DataFrame, info)
     return_traces = copy(trace)
     matched = matchExperiment(trace, info)
     matched_indexes = indexin(eachrow(matched), eachrow(trace))
     return_traces[matched_indexes, :INCLUDE] .= false
     return return_traces
end

function flagExperiment!(trace::DataFrame, info)
     matched = matchExperiment(trace, info)
     matched_indexes = indexin(eachrow(matched), eachrow(trace))
     trace[matched_indexes, :INCLUDE] .= false
end

function unflagALL!(dataset)
     dataset["TRACES"][:, :INCLUDE] .= true
     dataset["EXPERIMENTS"][:, :INCLUDE] .= true
end