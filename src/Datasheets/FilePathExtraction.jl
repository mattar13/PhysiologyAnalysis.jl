function traverse_root(root_path)
     path_list = String[]
     for (root, dirs, files) in walkdir(root_path)
         for file in files
             full_path = joinpath(root, file)  # Create the full path to the file
             push!(path_list, full_path)  # Print the path or do other processing
         end
     end
     return path_list
end

#%% Lets get a filename and extract all of the characteristics