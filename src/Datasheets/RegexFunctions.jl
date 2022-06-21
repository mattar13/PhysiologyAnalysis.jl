date_regex = r"(?'Year'\d{2,4})_(?'Month'\d{1,2})_(?'Date'\d{1,2})_(?'Description'.+)"
animal_regex = r"(?'Animal'\D+)(?'Number'\d)_(?'Age'.+)_(?'Genotype'.+)"
nd_file_regex = r"nd(?'ND'.{1,3})_(?'Percent'\d{1,3})p_.+abf"


function findmatch(str_array::Vector{String}, reg_format; verbose = false, first = true)
    matches = map(r -> match(reg_format, r), str_array)
    if any(!isnothing(matches))
        if verbose
            println("We found a format")
        end
        if first
            return matches[findall(!isnothing, matches)][1] |> NamedTuple
        else
            return matches[findall(!isnothing, matches)] .|> NamedTuple
        end
    else
        if verbose
            println("We did not find a format")
        end
    end
end