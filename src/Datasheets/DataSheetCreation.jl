function create2PDataSheet(img_dir, patch_dir; verbose = false)
    all_img_files = traverse_root(img_dir)
    img_files = [] #This is the list of all files actually analyzed
    img_allDates = []
    for (idx, img_file) in enumerate(all_img_files)
        if verbose
            println(img_file)
        end

        try
            img_test = readImage(img_file)
            push!(img_files, img_file)
            push!(img_allDates, img_test.HeaderDict["FileStartDateTime"])
            if verbose
                println("File $idx out of $(length(all_img_files))")
            end
        catch error
            if verbose
                println("File didn't load")
            end
        end
    end
    img_datasheet = DataFrame(
        filename = img_files, 
        date_created = img_allDates, filetype = "TIF", protocols = "Image"
    )


    patch_files = traverse_root(patch_dir)
    patch_allDates = []
    protocols = []

    for (idx, patch_file) in enumerate(patch_files)
        if verbose
            println(patch_file)
        end

        HeaderDict = ElectroPhysiology.readABFInfo(patch_file)
        push!(patch_allDates, getABF_datetime(patch_file))
        push!(protocols, HeaderDict["ProtocolPath"])
        if verbose
            println("File $idx out of $(length(patch_files))")
        end
    end
    patch_datasheet = DataFrame(filename = patch_files, 
        date_created = patch_allDates, filetype = "ABF", protocols = protocols)

    all_files = [img_datasheet; patch_datasheet] |> @orderby(_.date_created) |> DataFrame
    return Dict{String, DataFrame}(
        "All Files" => all_files,
        "TIF Files" => img_datasheet, 
        "ABF Files" => patch_datasheet
    )
end

function update2PDataSheet(img_dir, patch_dir, img_datasheet, patch_datasheet)
    all_img_files = traverse_root(img_dir)
    #Walk through each img file and check to see if it is in the datasheet
    patch_files = traverse_root(patch_dir)

end

function save2PDataSheet(filename::String, dataset::Dict{String, DataFrame})
    XLSX.openxlsx(filename, mode = "w") do xf
        for (k, v) in dataset
            if k == "All Files" #This one we want to rename sheet 1 (always first)
                sheet_K = xf[1] #Sheet 1 should be renamed
                XLSX.rename!(sheet_K, "All Files")
            else
                sheet_K = XLSX.addsheet!(xf, k)
            end
            XLSX.writetable!(sheet_K, v)
        end
    end
    println("Datasheet saved")
end

function open2PDataSheet(filename)
    xf = XLSX.readxlsx(filename)
    df_set = Dict{String, DataFrame}()
    for sn in XLSX.sheetnames(xf)
        df_set[sn] = DataFrame(XLSX.readtable(filename, sn))
    end
    return df_set
end