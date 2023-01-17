#using Revise
using ePhys
using DataFrames, Query, XLSX
import ePhys: findmatch, date_regex, animal_n_regex, age_regex, genotype_regex, cond_regex, pc_regex, color_regex
import ePhys: nd_regex, avg_regex
#%% New file location
#cone_paths = raw"C:\Users\mtarc\The University of Akron\Renna Lab - General\Data\ERG\Paul\Rods" |> parseABF
#cone_paths = raw"C:\Users\mtarc\The University of Akron\Renna Lab - General\Data\ERG\Paul\NotDetermined" |> parseABF
cone_paths = raw"C:\Users\mtarc\The University of Akron\Renna Lab - General\Data\ERG\Paul\Cones" |> parseABF
#cone_paths = [cone_paths1..., cone_paths2..., cone_paths3...]
#%% Restructure datafile location
rec_path = raw"C:\Users\mtarc\The University of Akron\Renna Lab - General\Data\ERG\Restructured"
rm(rec_path, recursive=true)
mkdir(rec_path)
#%% Restructure these paths
previous_paths = rec_path |> parseABF
previous_paths = map(path -> joinpath(splitpath(path)[1:end-1]), previous_paths)
previous_paths = unique(previous_paths)

count = 0
for path in cone_paths
     println(path)
     date_res = findmatch(path, date_regex)
     if parse(Int, date_res.Year) < 2020 && parse(Int, date_res.Month) >= 06
          user = "Paul"
     else
          user = "Matt"
     end
     #Find the age and 
     animal_res = findmatch(path, animal_n_regex)
     genotype_res = findmatch(path, r"_(?'Genotype'WT|DR)")
     age_res = findmatch(path, age_regex)
     if isnothing(age_res)
          AGE = "Adult"
     else
          AGE = age_res.Age
     end
     #Find Drugs, BaCl, NoDrugs
     cond_res = findmatch(path, cond_regex)
     new_path = joinpath(rec_path, "$(date_res.Year)_$(date_res.Month)_$(date_res.Date)_ERG$(genotype_res.Genotype)$(AGE)")

     if AGE == "30"
          new_path = joinpath(new_path, "Mouse$(animal_res.Number)_Adult_$(genotype_res.Genotype)")
     else
          new_path = joinpath(new_path, "Mouse$(animal_res.Number)_P$(AGE)_$(genotype_res.Genotype)")
     end

     #Find the conditions
     #println(cond_res)
     if cond_res.Condition == "Drugs"
          cond = "BaCl_LAP4"
     elseif cond_res.Condition == "NoDrugs" || cond_res.Condition == "No drugs" || cond_res.Condition == "No Drugs"
          #println("No Drugs")
          cond = "BaCl"
     else
          cond = cond_res.Condition
     end
     new_path = joinpath(new_path, "$cond")
     #Find photoreceptors
     pc_res = findmatch(path, pc_regex)
     if pc_res.Photoreceptors == "Rods" #No further label is needed
          new_path = joinpath(new_path, pc_res.Photoreceptors)
     else
          color_res = findmatch(path, color_regex)
          new_path = joinpath(new_path, "$(pc_res.Photoreceptors)_$(color_res.Color)")
     end

     if new_path ∉ previous_paths
          print("Append: ")
          println(new_path)
          mkpath(new_path)
     else
          print("Already exists: ")
          println(new_path)
          println(new_path)
     end

     #lets try to find the ND filter settings
     #finall lets try to find the word "Average or average
     corr_name = findmatch(path, avg_regex)
     nd_res = findmatch(path, nd_regex)
     flash_id = findmatch(path, r"\d") #find just a single digit
     if !isempty(flash_id)
          #println(flash_id)
     end
     if !isnothing(corr_name) && !isnothing(nd_res)
          count_str = string(count)
          #println(count_str)
          count += 1
          if count > 9999
               count = 0
          end

          if length(count_str) < 4
               addN = 4 - length(count_str)
          else
               addN = 0
          end
          str_lid = "$(join(repeat("0", addN)))$count_str"
          new_path = joinpath(new_path, "nd$(nd_res.ND)_$(nd_res.Percent)p_$str_lid.abf")
          #dat = read(path)
          #println(ispath(new_path))
          if !ispath(new_path)
               cp(path, new_path)
          end
     end
end