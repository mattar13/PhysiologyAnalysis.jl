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

date_age_regex = r"(?'Year'\d{2,4})_(?'Month'\d{1,2})_(?'Date'\d{1,2})_(?'Genotype'.+)_P(?'Age'\d+)" #One regex will have the age
date_regex = r"(?'Year'\d{2,4})_(?'Month'\d{1,2})_(?'Date'\d{1,2})_(?'Genotype'.+)" #One regex won't have the age
cell_patch_regex = r"Cell(?'Cell'\d+)"
ca_catch_regex = r"ca_img(?'Cell'\d)"
cell_catch_regex = r"cell_img(?'Cell'\d)"

function parse_cell_details(filename)
    date_age_match = match(date_age_regex, filename)
    date_match = match(date_regex, filename)
    if !isnothing(date_age_match) #We failed to notate the age
        genotype = date_age_match["Genotype"]
        age = parse(Int64, date_age_match["Age"])
    elseif !isnothing(date_match)
        genotype = date_match["Genotype"]
        age = 30
    else
        genotype = "WT"
        age = 30
    end

    cell_match = match(cell_patch_regex, filename)
    ca_catch_match = match(ca_catch_regex, filename)
    cell_catch_match = match(cell_catch_regex, filename)
    if !isnothing(cell_match)
        cell_n = parse(Int64, cell_match["Cell"])
    elseif !isnothing(ca_catch_match)
        cell_n = parse(Int64, ca_catch_match["Cell"])
    elseif !isnothing(cell_catch_match)
        cell_n = parse(Int64, cell_catch_match["Cell"])
    else
        cell_n = 0
    end
    println(age, genotype, cell_n)
    return age, genotype, cell_n
end