function pair_experiments(patch_datasheet, img_datasheet; tolerance = 100.0)
    paired_experiments = DataFrame(
        "ABF Experiment" => String[], 
        "TIF Experiment" => String[],
        "ABF Date" => DateTime[],
        "TIF Date" => DateTime[]
    )
    for (i, date) in enumerate(patch_datasheet[:, "date_created"])
        #subtract the entirety of the tif_datetime
        elapsed_dates = abs.(date .- img_datasheet[:, "date_created"])
        val, idx = findmin(elapsed_dates)
        val_seconds = Dates.value(val)/1000
        if val_seconds < tolerance
            push!(paired_experiments, (patch_datasheet[i, "filename"],  img_datasheet[idx, "filename"], date, img_datasheet[idx, "date_created"]))
        end
    end
    paired_experiments
end