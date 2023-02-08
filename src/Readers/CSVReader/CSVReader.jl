function readCSV(file::String; info_cols = 1:3, label_col = 4, data_col_start = 5)
     dlm_file = readdlm(file, ',', skipblanks=true,) #Open the dlm_files
     file_info = dlm_file[info_cols , :]
     labels = dlm_file[label_col, :] 
     #remove any empty columns
     time_column = findall(labels .== "Time (s)") 
     not_empty_columns = findall( .! isempty.(labels))
     data_columns = not_empty_columns[not_empty_columns .!= time_column]
          
     time = dlm_file[data_col_start:end, time_column] |> vec
     not_empty_rows = findall( .! isempty.(time))
     data_array = dlm_file[data_col_start:end, data_columns]
     time = Vector{Float64}(time[not_empty_rows])
     data_array = data_array[not_empty_rows, :]
     data_array[isempty.(data_array)] .= 0.0
     data_array = Matrix{Float64}(data_array)
     data_array = reshape(data_array, (size(data_array,1), size(data_array,2), 1))
     data_array = permutedims(data_array, (2,1,3))
     #return time, data_array
     return Experiment(time, data_array)
end