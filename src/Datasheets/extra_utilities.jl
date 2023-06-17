import ElectroPhysiology.readABF

function readABF(df::DataFrame; kwargs...)
     data = readABF(df.Path; kwargs...)
     if df.SubPath[1] != "NONE"
          data_sub = readABF(df.SubPath; kwargs...)
          return data_sub, data
     else
          return data
     end
end