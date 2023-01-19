import Base.NamedTuple

date_regex = r"(?'Year'\d{2,4})_(?'Month'\d{1,2})_(?'Date'\d{1,2})_(?'Description'.+)"
#animal_file_regex = r"(?'Animal'\D+)(?'Number'\d)_(?'Age'.+)_(?'Genotype'.+)"
nd_file_regex = r"nd(?'ND'.{1,3})_(?'Percent'\d{1,3})p_.+abf"

animal_regex = r"(_m|(?'Animal'Mouse|Zebrafish|Organoid)|m)(?'Number'\d)"
age_regex = r"_P(?'Age'\d*|)"
genotype_regex = r"_(?'Genotype'WT|DR|R141C|RS1KO|C59S|MelKO)"
cond_regex = r"(?'Condition'Drugs|NoDrugs|BaCl_LAP4|BaCl|No drugs|No Drugs)"
pc_regex = r"(?'Photoreceptors'Cones|Rods)"
color_regex = r"(?'Color'blue|green|Blue|Green|UV|365|365UV|520|520Green|525|525Green)"
avg_regex = r"(?'Correct'Average|average)"
background_regex = r"(?'Background'withback|noback)"
percent_regex = r"(?'Percent'\d{1,3})(%|p)"
nd_regex = r"(nd|ND)(?'ND'\d{1,3})"

#nd_regex = r"nd(?'ND'.{1,3})_(?'Percent'\d{1,3})p"
NamedTuple(m::RegexMatch) = NamedTuple{Symbol.(Tuple(keys(m)))}(values(m.captures))

function findmatch(str_array::Vector{String}, reg_format::Regex; verbose=false, first=true)
    matches = map(r -> match(reg_format, r), str_array)
    #println(matches)
    #println("Revise is working")
    #println(any(.! isnothing.(matches)))
    if any(.! isnothing.(matches))
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

function findmatch(path::String, reg_format::Regex; kwargs...)
    str_array = splitpath(path)
    return findmatch(str_array, reg_format; kwargs...)
end

function find_condition(str_array; possible_conds=["BaCl", "BaCl_LAP4", "NoDrugs"])
    for cond in possible_conds
        val = findall(str_array .== cond)
        if !isempty(val)
            cond = str_array[val]
            return cond[1]
        end
    end
end

"""
These functions keep only things parsed as the object T, which defaults to Int64
"""

filter_string(::Type{T}, str::String) where {T} = filter(!isnothing, map(c -> tryparse(T, c), split(str, ""))) |> join
filter_string(::Type{T}, str::SubString{String}) where {T} = filter_string(T, str |> string)
filter_string(str::String) = filter(!isnothing, map(c -> tryparse(Int64, c), split(str, ""))) |> join
filter_string(::Type{T}, ::Nothing) where {T} = nothing

"""
This function takes a named tuple that contains numbers and cleans those numbers
"""
function parseNamedTuple(::Type{T}, nt::NamedTuple{keys}) where {T<:Real,keys}
    nt_vals = [values(nt)...] #Fill nt_vals with the values of the NamedTuple
    new_vals = []
    for itm in nt_vals
        #println(itm |> typeof)
        if isa(itm, String) || isa(itm, SubString{String})#If the item is not nothing
            #We only want to parse the item if all of it is a number
            #println(itm)
            val = tryparse(T, itm) #Try to parse the item as a Int64
            if !isnothing(val)
                push!(new_vals, val)
            else #We should make sure that we at least convert these to a string
                push!(new_vals, String(itm))
            end
        else
            push!(new_vals, itm)
        end
    end
    NamedTuple{keys}(new_vals)
end

#If no type is provided the namedtuple is parsed as a Int64
parseNamedTuple(nt::NamedTuple{keys}) where {keys} = parseNamedTuple(Int64, nt)