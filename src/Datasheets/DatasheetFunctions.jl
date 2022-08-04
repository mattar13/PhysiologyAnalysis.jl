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
     println("Revised")
     dataframe = DataFrame()
     for file in all_files
          println("Analyzing file: $file")
          entry = DataPathExtraction(file)
          push!(dataframe, entry)
     end
     dataframe
end