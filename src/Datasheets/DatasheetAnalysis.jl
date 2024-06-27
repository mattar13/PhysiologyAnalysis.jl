function pair_experiments!(dataset::Dict{String, DataFrame};            
    vmin = 5.30,
    vmax = 6.30
)
    patch_datasheet = dataset["ABF Files"]
    img_datasheet = dataset["TIF Files"]
    paired_experiments = DataFrame(
        "ABF Experiment" => String[], 
        "TIF Experiment" => String[],
        "ABF Date" => DateTime[],
        "TIF Date" => DateTime[]
    )
    for (i, row) in enumerate(eachrow(patch_datasheet))
        if occursin("5mins", row.protocols)
            date = row.date_created
            elapsed_dates = (Dates.value.(abs.(date.- img_datasheet[:, "date_created"]))./1000)./60

            idx = findlast(vmin .< elapsed_dates .<= vmax)

            if !isnothing(idx)
                push!(paired_experiments, (patch_datasheet[i, "filename"],  img_datasheet[idx, "filename"], date, img_datasheet[idx, "date_created"]))
            end
        end
    end
    dataset["Paired Files"] = paired_experiments
    dataset
end

function IV_analysis!(dataset::Dict{String, DataFrame})
    patch_datasheet = dataset["ABF Files"]
    patch_datasheet[!, "Rs"] = fill(NaN, size(patch_datasheet,1))
    patch_datasheet[!, "Ri"] = fill(NaN, size(patch_datasheet,1))
    patch_datasheet[!, "Vh"] = fill(NaN, size(patch_datasheet,1))
    patch_datasheet[!, "Vm"] = fill(NaN, size(patch_datasheet,1))
    patch_datasheet[!, "Cm"] = fill(NaN, size(patch_datasheet,1))
    patch_datasheet[!, "τ"] = fill(NaN, size(patch_datasheet,1))
    patch_datasheet[!, "Quality"] = fill(NaN, size(patch_datasheet,1))

    for (i, row) in enumerate(eachrow(patch_datasheet))
        
        if occursin("IV", row.protocols)
            println(i)
            exp = readABF(row.filename); 
            if size(exp,3) >= 2 && size(exp,1) == 9
                println(row.protocols)
                println(row.filename)
                Rs, Rin, V_M, V_HOLD = calculate_resistance(exp)
                patch_datasheet[i, "Rs"] = Rs
                patch_datasheet[i, "Ri"] = Rin
                patch_datasheet[i, "Vh"] = V_HOLD
                patch_datasheet[i, "Vm"] = V_M
                
                patch_datasheet[i, "Cm"] = Cm = calculate_capacitance(exp)
                patch_datasheet[i, "τ"] = Rs*Cm
                patch_datasheet[i, "Quality"] = Rs/Rin * 100
            end
        end

    end
    dataset["ABF Files"] = patch_datasheet
    return dataset
end