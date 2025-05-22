"""
This documentation should help you figure out what you have access to

- "dx" and "dy" are the x and y limit differentials
- "img_arr" is the raw data of all frames
- "red_zstack" is the red channel of the image array
- "grn_zstack" is the green channel of the image array
- "red_zproj" and "grn_zproj" are the red and green zprojection (taken with mean function)
- "red_img" and "grn_img" are those zprojections made into RGB{Float32} images
- "composite_img" is a composite of the red and grn images
- "red_trace" and "grn_trace" are the traces of the mean fluoresence by frame
- "dff_red_zstack" is the delta f/f (f0-f)/f taken of by a rolling mean
- "dff_grn_zstack" is the dff taken using just the total mean of the trace
- "dff_red_comp_zstack" and "dff_grn_comp_zstack" are the RGB{Float32} versions
- "dff_composite_zstack" is the composite
- "dff_red_zproj" and "dff_grn_zproj" are the zstack projections of the dff images
- "dff_red_img" and "dff_grn_img" are the image versions of the red and green z projections
- "dff_composite_img" is the composite of the red and green images
- "dff_red_trace" and "dff_grn_trace" are the trace versions
- "pks_vals" are the maxima peaks and valleys of the dff_green_trace
- "grn_row_sums" is the sum of all green section rows
- "grn_sect_arr" and "red_sect_arr" are the arrays formed as a result of the peakfinding
"""
function open2Pdata(filename;
        split_channel = true, main_channel = :red,
        trunc_rng = nothing, pre_event_time = 20.0, post_event_time = 60.0, 
        red_scale = 1.5, grn_scale = 3.0, 
        section_by = :green, peak_width = 40, peak_min = 0.2, 
        ic_stim_filename = nothing, #New option if we want to specify a 2P stimulus dataset
        stimulus_name = "IN 2", stimulus_threshold = 0.5,
        spike_train = false, 
        grn_lam = 1e4, red_lam = 1e4, grn_window = 5, red_window = 5
    )
    #╔═╡Do this for the k-puff ca image
    output = Dict{String, Any}()
    experiment = readImage(filename);
    if split_channel
        deinterleave!(experiment) #This seperates the movies into two seperate movies
    end

    #Make this a conditional
    if !isnothing(trunc_rng)
        t_begin = trunc_rng[1]
        t_end = trunc_rng[2]
        if isnothing(trunc_rng[2])
            t_end = experiment.t[end]
        end

        if isnothing(trunc_rng[1])
            t_begin = experiment.t[1]
        end
        truncate_data!(experiment, t_begin = t_begin, t_end = t_end)
    end
    output["experiment"] = experiment
    output["time"] = experiment.t
    output["xlims"] = xlims = experiment.HeaderDict["xrng"]
    output["ylims"] = ylims = experiment.HeaderDict["xrng"]
    output["dt"] = experiment.dt
    output["dx"] = xlims[2]-xlims[1]
    output["dy"] = ylims[2]-ylims[1]
    println("Image loaded")
   
    if split_channel
        
        #Create a 2D Gaussian filter (Kernel.gaussian(3))
        spatial_filter = Kernel.gaussian(2.0)  # Gaussian kernel size 3
        #imfilter!(experiment, spatial_filter; channel = 2) #Apply a gaussian filter

        #Average 5 frames together as a rolling mean
        mapdata!(mean, experiment, 5, channel = 2)
        println("Image filtered")

        #╔═╡Seperate out the image array 
        output["img_arr"] = img_arr = get_all_frames(experiment)
        #╔═╡Seperate the red and green channels
        output["red_zstack"] = red_zstack = img_arr[:,:,:,2]
        output["grn_zstack"] = grn_zstack = img_arr[:,:,:,1]
        output["composite_zstack"] = RGB{Float32}.(red_zstack, grn_zstack, zeros(size(red_zstack)...))
        
        #╔═╡Make the z projections
        output["red_zproj"] = red_zproj = project(experiment, dims = (3))[:,:,1,2]
        output["grn_zproj"] = grn_zproj = project(experiment, dims = (3))[:,:,1,1]

        output["red_img"] = red_img = red_zproj./maximum(red_zproj)*red_scale .* RGB{Float32}(1, 0, 0)
        output["grn_img"] = grn_img = grn_zproj./maximum(grn_zproj)*grn_scale .* RGB{Float32}(0, 1, 0)
        output["composite_img"] = grn_img + red_img
        #save("$analysis_loc\\$(savename)\\img.tif", composite_zstack)
        println("Z projection")

        output["red_trace"] = project(experiment, dims = (1,2))[1,1,:,2]
        output["grn_trace"] = project(experiment, dims = (1,2))[1,1,:,1]
        println("Z axis traces generated")
 
        # #Using the rolling mean method
        # output["dff_red_zstack"] = dff_red_zstack = deltaF_F(red_zstack; voxel_z = 200, mode = :rolling_mean, boundary_mode = "symmetric")
        # #Using the mean method
        # output["dff_grn_zstack"] = dff_grn_zstack = deltaF_F(grn_zstack; mode = :mean)

        # output["dff_red_comp_zstack"] = dff_red_comp_zstack = dff_red_zstack .* RGB{Float32}(1, 0, 0) 
        # output["dff_grn_comp_zstack"] = dff_grn_comp_zstack = dff_grn_zstack .* RGB{Float32}(0, 1, 0)
        # output["dff_composite_zstack"] = dff_grn_comp_zstack + dff_red_comp_zstack

        # output["dff_red_zproj"] = dff_red_zproj = mean(dff_red_zstack, dims = 3)[:,:,1]
        # output["dff_grn_zproj"] = dff_grn_zproj = mean(dff_grn_zstack, dims = 3)[:,:,1]

        # output["dff_red_img"] = dff_red_img = dff_red_zproj .* RGB{Float32}(1, 0, 0)
        # output["dff_grn_img"] = dff_grn_img = dff_grn_zproj .* RGB{Float32}(0, 1, 0)
        # output["dff_composite_img"] = dff_grn_img + dff_red_img #Make this a 

        # output["dff_red_trace"] = dff_red_trace = mean(dff_red_zstack, dims = (1,2))[1,1,:]
        # output["dff_grn_trace"] = dff_grn_trace = mean(dff_grn_zstack, dims = (1,2))[1,1,:]

        output["dff_grn_trace"] = dff_grn_trace = PhysiologyAnalysis.baseline_trace(output["grn_trace"], window = grn_window, lam = grn_lam, niter = 100)
        output["dff_red_trace"] = dff_red_trace = PhysiologyAnalysis.baseline_trace(output["red_trace"], window = red_window, lam = red_lam, niter = 100)
        println("delta f/f images and traces extracted")
    else
        #Create a 2D Gaussian filter (Kernel.gaussian(3))
        spatial_filter = Kernel.gaussian(2.0)  # Gaussian kernel size 3
        imfilter!(experiment, spatial_filter; channel = 1) #Apply a gaussian filter

        #Average 5 frames together as a rolling mean
        mapdata!(mean, experiment, 5, channel = 1)
        println("Image filtered")

        println("Z Stack channels extracted")

        if main_channel == :red
            #╔═╡Seperate out the image array and green and red zstacks
            output["img_arr"] = output["red_zstack"] = red_zstack = img_arr = get_all_frames(experiment)
            output["grn_zstack"] = zeros(size(red_zstack))
            output["composite_zstack"] = RGB{Float32}.(red_zstack, zeros(size(red_zstack)...), zeros(size(red_zstack)...))

            output["red_zproj"] = red_zproj = project(experiment, dims = (3))[:,:,1,1]
            output["grn_zproj"] = grn_zproj = zeros(size(red_zproj))

            output["red_img"] = output["composite_img"] = red_img = red_zproj./maximum(red_zproj)*red_scale .* RGB{Float32}(1, 0, 0)
            output["grn_img"] = grn_img = zeros(size(red_img))

            output["red_trace"] = red_trace = project(experiment, dims = (1,2))[1,1,:,1]
            output["grn_trace"] = grn_trace = zeros(size(red_trace))

            # #Using the rolling mean method
            # output["dff_red_zstack"] = dff_red_zstack = deltaF_F(red_zstack[:,:,:,1]; voxel_z = 200, mode = :rolling_mean, boundary_mode = "symmetric")
            # output["dff_grn_zstack"] = dff_grn_zstack = zeros(size(dff_red_zstack))

            # output["dff_red_comp_zstack"] = output["dff_composite_zstack"] = dff_red_comp_zstack = dff_red_zstack .* RGB{Float32}(1, 0, 0) 
            # output["dff_grn_comp_zstack"] = dff_grn_comp_zstack = zeros(size(dff_red_comp_zstack))

            # output["dff_red_zproj"] = dff_red_zproj = mean(dff_red_zstack, dims = 3)[:,:,1]
            # output["dff_grn_zproj"] = dff_grn_zproj = zeros(size(dff_red_zproj))

            # output["dff_red_img"] = output["dff_composite_img"] = dff_red_img = dff_red_zproj .* RGB{Float32}(1, 0, 0)
            # output["dff_grn_img"] = dff_grn_img = zeros(size(dff_red_img))

            output["dff_red_trace"] = dff_red_trace = PhysiologyAnalysis.baseline_trace(output["red_trace"],  window = red_window, lam = red_lam, niter = 100)
            output["dff_grn_trace"] = dff_grn_trace = zeros(size(dff_red_trace))


        elseif main_channel == :grn
             #╔═╡Seperate out the image array and green and red zstacks
             output["img_arr"] = output["grn_zstack"] = grn_zstack = img_arr = get_all_frames(experiment)
             output["red_zstack"] = zeros(size(grn_zstack))
             output["composite_zstack"] = RGB{Float32}.(grn_zstack, zeros(size(grn_zstack)...), zeros(size(grn_zstack)...))
 
             output["grn_zproj"] = grn_zproj = project(experiment, dims = (3))[:,:,1,1]
             output["red_zproj"] = red_zproj = zeros(size(grn_zproj))
 
             output["grn_img"] = output["composite_img"] = grn_img = grn_zproj./maximum(grn_zproj)*grn_scale .* RGB{Float32}(1, 0, 0)
             output["red_img"] = red_img = zeros(size(grn_img))
 
             output["grn_trace"] = grn_trace = project(experiment, dims = (1,2))[1,1,:,1]
             output["red_trace"] = red_trace = zeros(size(grn_trace))
 
            #  #Using the rolling mean method
            #  output["dff_grn_zstack"] = dff_grn_zstack = deltaF_F(grn_zstack[:,:,:,1]; voxel_z = 200, mode = :rolling_mean, boundary_mode = "symmetric")
            #  output["dff_red_zstack"] = dff_red_zstack = zeros(size(dff_grn_zstack))
 
            #  output["dff_grn_comp_zstack"] = output["dff_composite_zstack"] = dff_grn_comp_zstack = dff_grn_zstack .* RGB{Float32}(1, 0, 0) 
            #  output["dff_red_comp_zstack"] = dff_red_comp_zstack = zeros(size(dff_grn_comp_zstack))
 
            #  output["dff_grn_zproj"] = dff_grn_zproj = mean(dff_grn_zstack, dims = 3)[:,:,1]
            #  output["dff_red_zproj"] = dff_red_zproj = zeros(size(dff_grn_zproj))
 
            #  output["dff_grn_img"] = output["dff_composite_img"] = dff_grn_img = dff_grn_zproj .* RGB{Float32}(1, 0, 0)
            #  output["dff_red_img"] = dff_red_img = zeros(size(dff_grn_img))
 
             output["dff_grn_trace"] = dff_grn_trace = PhysiologyAnalysis.baseline_trace(output["grn_trace"],  window = grn_window, lam = grn_lam, niter = 100)
             output["dff_red_trace"] = dff_red_trace = zeros(size(dff_grn_trace))
        end
        println("Z axis traces generated")


        println("delta f/f images and traces extracted")
    end

    #%% Do and plot wavefinding events
    #If we specify a 
    if isnothing(ic_stim_filename)
        if section_by == :green
            println("Peak finding by green channel")
            output["pks"], output["vals"] = pks, vals = findmaxima(dff_grn_trace, peak_width)
        elseif section_by == :red
            println("Peak finding by red channel")
            output["pks"], output["vals"] = pks, vals = findmaxima(dff_red_trace, peak_width)
        end
    else
        println("Peak finding by using the digital stim in the IC stimulus instead")
        addStimulus!(experiment, ic_stim_filename, stimulus_name; flatten_episodic = true, stimulus_threshold = stimulus_threshold)
        dataIC = readABF(ic_stim_filename, flatten_episodic = true, stimulus_name = stimulus_name, stimulus_threshold = stimulus_threshold) #Open the IC data
        start2P = experiment.HeaderDict["FileStartDateTime"]-Second(3.0) #The computer clocks are off by 3 seconds
        startIC = dataIC.HeaderDict["FileStartDateTime"]
        t_offset = Millisecond(startIC - start2P).value/1000 
        time_offset!(dataIC, t_offset)
        stim_protocol = getStimulusProtocol(dataIC)
        spike_train_group!(stim_protocol, 3.0) #We only need to do this if there are spike trains
        
        output["dataIC"] = dataIC
        ElectroPhysiology.convert_channel_to_stimulus!(experiment, dataIC, "IN 3")
        
        if spike_train
            spike_train_group!(stim_protocol, 3.0) #We only need to do this if there are spike trains
            
            spike_train_protocol = getStimulusProtocol(experiment)
            spike_train_group!(spike_train_protocol, 3.0)
        end
        output["tstamps"] = t_stamps = map(sp -> sp[1][1] + t_offset, stim_protocol)
        output["pks"] = pks = round.(Int64, (t_stamps./experiment.dt))
        
        #output["pks"] = pks = round.(Int64, (t_episodes./experiment.dt))
        println(pks)
    end

    pre_event_length = floor(Int64, pre_event_time/experiment.dt)
    post_event_length = floor(Int64, post_event_time/experiment.dt)
    output["sect_time"] = new_t = LinRange(-pre_event_time, post_event_time, pre_event_length+post_event_length)
    red_sect_arr = zeros(length(new_t), length(pks))
    grn_sect_arr = zeros(length(new_t), length(pks))
    for (i, pk) in enumerate(pks)
        #println(pk-pre_event_length)
        if pk-pre_event_length < 0
            idx_start = 1
        else
            idx_start = round(Int64, pk-pre_event_length)
        end
        
        #println(pk+post_event_length)
        if pk+post_event_length > length(dff_red_trace)
            idx_end = length(dff_grn_trace) 
        else
            idx_end = round(Int64, pk+post_event_length-1)
        end
        
        idx_rng = idx_start:idx_end
        #println(idx_rng)

        grn_sect = dff_grn_trace[idx_rng]
        red_sect = dff_red_trace[idx_rng]
        red_sect_arr[1:length(red_sect), i] = red_sect
        grn_sect_arr[1:length(grn_sect), i] = grn_sect
    end
    output["grn_row_sums"] = grn_row_sums = sum(grn_sect_arr, dims = 1)[1,:]
    output["red_row_sums"] = red_row_sums = sum(red_sect_arr, dims = 1)[1,:]
    output["grn_sect_arr"] = grn_sect_arr# = grn_sect_arr[:, findall(grn_row_sums .!= 0.0)]
    output["red_sect_arr"] = red_sect_arr# = red_sect_arr[:, findall(red_row_sums .!= 0.0)]

    #We need to clean empty rows
    println("Peak finding completed")
    return output
end