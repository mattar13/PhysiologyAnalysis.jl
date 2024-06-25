function pair_experiments(patch_datasheet, img_datasheet;            
    vmin = 5.30,
    vmax = 6.30
)
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
    paired_experiments
end