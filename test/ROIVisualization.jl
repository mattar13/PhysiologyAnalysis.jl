using GLMakie
using Statistics
using ElectroPhysiology, PhysiologyAnalysis

"""
    plot_roi_analysis(data::Experiment{TWO_PHOTON}; 
        stim_idx::Int=1, channel_idx::Union{Int,Nothing}=nothing)

Create a comprehensive visualization of ROI analysis results for a single stimulus event including:
1. ROI map showing significant ROIs weighted by their amplitude, overlaid on max projection
2. Individual dF/F traces for significant ROIs
3. Mean trace with standard deviation ribbon

If channel_idx is provided, only that channel will be shown. Otherwise, all channels will be displayed.

Returns a Figure object containing all plots.
"""
function plot_roi_analysis(data::Experiment{TWO_PHOTON, T};
    stim_idx::Int=1, channel_idx::Union{Int,Nothing}=nothing) where {T <: Real}
    
    @assert haskey(data.HeaderDict, "ROI_Analysis") "Data must contain ROI analysis results in HeaderDict"

    analysis = data.HeaderDict["ROI_Analysis"]

    # Plot max projection as background
    xlims = data.HeaderDict["xrng"]
    ylims = data.HeaderDict["yrng"]

    # Determine which channels to process
    channels_to_process = isnothing(channel_idx) ? analysis.channels : [channel_idx]
    n_channels = length(channels_to_process)
    
    # Create figure with layout
    fig = Figure(size=(1200 * n_channels, 800))
    
    # Process each channel
    for (ch_idx, channel) in enumerate(channels_to_process)
        # Get significant ROIs for this channel and stimulus
        sig_rois = get_significant_rois(analysis, stim_idx, channel)
        
        # Create a 2x2 grid for this channel
        gl_channel = fig[1, ch_idx] = GridLayout()
        
        # 1. ROI amplitude map with max projection (top left)
        ax1 = Axis(gl_channel[1,1], title="Channel $channel Response Map (Stimulus $stim_idx)",
            xlabel="X Position (μm)", ylabel="Y Position (μm)",
            aspect=DataAspect())
        
        # Get max projection of original data
        max_proj = project(data, dims = (3))[:,:,1,channel]
        max_proj = rotr90(max_proj)
        norm_max_proj = (max_proj .- minimum(max_proj)) ./ (maximum(max_proj) - minimum(max_proj))
        hm1 = heatmap!(ax1, xlims, ylims, norm_max_proj)
        Colorbar(gl_channel[1,2], hm1, label="dF/F", width=15, vertical=true)
        
        # Create a sub-grid for the traces
        gl_traces = gl_channel[1,3] = GridLayout()
        
        # 2. Mean trace with std ribbon (top right, left half)
        ax2 = Axis(gl_traces[1,1], title="Mean Response ± STD (Stimulus $stim_idx)",
            xlabel="Time (s)", ylabel="ΔF/F")
        
        if !isempty(sig_rois)
            # Get all significant traces and calculate statistics
            sig_traces = [filter(t -> t.channel == channel && t.stimulus_index == stim_idx, analysis.rois[roi])[1].dfof for roi in sig_rois]
            t_series = filter(t -> t.channel == channel && t.stimulus_index == stim_idx, analysis.rois[first(sig_rois)])[1].t_series
            
            mean_trace = mean(sig_traces)
            std_trace = std(sig_traces)
            
            # Plot mean and standard deviation
            band!(ax2, t_series, 
                mean_trace .- std_trace, 
                mean_trace .+ std_trace,
                color=(:blue, 0.3))
            lines!(ax2, t_series, mean_trace, 
                color=:blue, linewidth=2,
                label="Mean (n=$(length(sig_rois)))")
        end
        
        # 3. Individual traces (top right, right half)
        ax3 = Axis(gl_traces[2,1], title="Individual ROI Traces (Stimulus $stim_idx)",
            xlabel="Time (s)", ylabel="ΔF/F")
        
        if !isempty(sig_rois)
            # Plot each significant ROI trace
            colors = cgrad(:viridis, length(sig_rois), categorical=true)
            for (i, roi_id) in enumerate(sig_rois)
                traces = filter(t -> t.channel == channel && t.stimulus_index == stim_idx, analysis.rois[roi_id])
                trace = traces[1]
                lines!(ax3, trace.t_series, trace.dfof, 
                    color=colors[i], alpha=0.5)
            end
        end
        
        # Link y-axes of the traces
        linkyaxes!(ax2, ax3)
        
        # Get ROI mask and create weighted map
        roi_mask = getROImask(data)
        roi_mask = rotr90(roi_mask)
        weighted_mask = zeros(size(roi_mask))
        tau_off_mask = zeros(size(roi_mask))
        
        for roi_id in keys(analysis.rois)
            if roi_id in sig_rois
                traces = filter(t -> t.channel == channel && t.stimulus_index == stim_idx, analysis.rois[roi_id])
                trace = traces[1]
                max_response = maximum(trace.dfof)
                weighted_mask[roi_mask .== roi_id] .= max_response
                if !isnothing(trace.fit_parameters) && length(trace.fit_parameters) >= 3
                    tau_off_mask[roi_mask .== roi_id] .= trace.fit_parameters[3]
                end
            end
        end
        
        # 4. Weighted response map (bottom left)
        ax_weighted = Axis(gl_channel[2,1], title="Channel $channel Weighted Response (Stimulus $stim_idx)",
            xlabel="X Position (μm)", ylabel="Y Position (μm)",
            aspect=DataAspect())
        
        if !isempty(sig_rois)
            hm_weighted = heatmap!(ax_weighted, xlims, ylims, weighted_mask,
                colormap=:viridis,
                colorrange=(0, maximum(weighted_mask)))
            Colorbar(gl_channel[2,2], hm_weighted, label="dF/F", width=15, vertical=true)
        else
            text!(ax_weighted, "No significant ROIs found", 
                position=(mean(xlims), mean(ylims)),
                align=(:center, :center),
                color=:red)
        end
        
        # 5. Tau Off map (bottom right)
        ax_tau = Axis(gl_channel[2,3], title="Channel $channel Tau Off (Stimulus $stim_idx)",
            xlabel="X Position (μm)", ylabel="Y Position (μm)",
            aspect=DataAspect())
        
        if !isempty(sig_rois)
            hm_tau = heatmap!(ax_tau, xlims, ylims, tau_off_mask,
                colormap=:viridis,
                colorrange=(0, maximum(filter(!iszero, tau_off_mask))))
            Colorbar(gl_channel[2,4], hm_tau, label="Tau Off (s)", width=15, vertical=true)
        else
            text!(ax_tau, "No significant ROIs found", 
                position=(mean(xlims), mean(ylims)),
                align=(:center, :center),
                color=:red)
        end
        
        # Add stimulus time indicator if available
        if haskey(analysis.analysis_parameters, :delay_time)
            delay_time = analysis.analysis_parameters[:delay_time]
            vlines!(ax2, [delay_time], color=:red, linestyle=:dash, label="Stimulus")
            vlines!(ax3, [delay_time], color=:red, linestyle=:dash, label="Stimulus")
        end
        
        # Add legend only for mean response
        axislegend(ax2, position=:rt)
    end
    
    return fig
end

"""
    plot_roi_analysis_averaged(data::Experiment{TWO_PHOTON}; 
        channel_idx::Union{Int,Nothing}=nothing)

Create a visualization showing the average responses across all stimuli for each ROI, including:
1. ROI map showing significant ROIs weighted by their average amplitude across stimuli
2. Average dF/F traces across all stimuli for each significant ROI
3. Overall mean trace with standard deviation ribbon

If channel_idx is provided, only that channel will be shown. Otherwise, all channels will be displayed.

Returns a Figure object containing all plots.
"""
function plot_roi_analysis_averaged(data::Experiment{TWO_PHOTON, T};
    channel_idx::Union{Int,Nothing}=nothing) where {T <: Real}
    
    @assert haskey(data.HeaderDict, "ROI_Analysis") "Data must contain ROI analysis results in HeaderDict"

    analysis = data.HeaderDict["ROI_Analysis"]
    # Plot max projection as background
    xlims = data.HeaderDict["xrng"]
    ylims = data.HeaderDict["yrng"]

    # Determine which channels to process
    channels_to_process = isnothing(channel_idx) ? analysis.channels : [channel_idx]
    n_channels = length(channels_to_process)
    
    # Create figure with layout
    fig = Figure(size=(1200 * n_channels, 800))
    
    # Process each channel
    for (ch_idx, channel) in enumerate(channels_to_process)
        # Get significant ROIs for this channel (significant in any stimulus)
        sig_rois = unique([id for (id, traces) in analysis.rois 
            if any(t -> t.channel == channel && t.is_significant, traces)])
        
        # Create a 2x2 grid for this channel
        gl_channel = fig[1, ch_idx] = GridLayout()
        
        # 1. ROI amplitude map with max projection (top left)
        ax1 = Axis(gl_channel[1,1], title="Channel $channel Response Map (Averaged)",
            xlabel="X Position (μm)", ylabel="Y Position (μm)",
            aspect=DataAspect())
        
        # Get max projection of original data
        max_proj = project(data, dims = (3))[:,:,1,channel]
        max_proj = rotr90(max_proj)
        norm_max_proj = (max_proj .- minimum(max_proj)) ./ (maximum(max_proj) - minimum(max_proj))
        hm1 = heatmap!(ax1, xlims, ylims, norm_max_proj)
        Colorbar(gl_channel[1,2], hm1, label="dF/F", width=15, vertical=true)
        
        # Create a sub-grid for the traces
        gl_traces = gl_channel[1,3] = GridLayout()
        
        # 2. Mean trace with std ribbon (top right, left half)
        ax2 = Axis(gl_traces[1,1], title="Mean Response ± STD (Averaged)",
            xlabel="Time (s)", ylabel="ΔF/F")
        
        if !isempty(sig_rois)
            # Calculate average traces across all stimuli for each ROI
            avg_traces = []
            t_series = nothing
            
            for roi_id in sig_rois
                traces = filter(t -> t.channel == channel, analysis.rois[roi_id])
                if !isempty(traces)
                    # Average the traces for this ROI across all stimuli
                    roi_avg = mean(t.dfof for t in traces)
                    push!(avg_traces, roi_avg)
                    if isnothing(t_series)
                        t_series = traces[1].t_series
                    end
                end
            end
            
            if !isempty(avg_traces)
                mean_trace = mean(avg_traces)
                std_trace = std(avg_traces)
                
                # Plot mean and standard deviation
                band!(ax2, t_series, 
                    mean_trace .- std_trace, 
                    mean_trace .+ std_trace,
                    color=(:blue, 0.3))
                lines!(ax2, t_series, mean_trace, 
                    color=:blue, linewidth=2,
                    label="Mean (n=$(length(sig_rois)))")
            end
        end
        
        # 3. Individual averaged traces (top right, right half)
        ax3 = Axis(gl_traces[2,1], title="Individual ROI Traces (Averaged)",
            xlabel="Time (s)", ylabel="ΔF/F")
        
        if !isempty(sig_rois)
            # Plot each ROI's average trace
            colors = cgrad(:viridis, length(sig_rois), categorical=true)
            for (i, roi_id) in enumerate(sig_rois)
                traces = filter(t -> t.channel == channel, analysis.rois[roi_id])
                if !isempty(traces)
                    roi_avg = mean(t.dfof for t in traces)
                    lines!(ax3, traces[1].t_series, roi_avg, 
                        color=colors[i], alpha=0.5)
                end
            end
        end
        
        # Link y-axes of the traces
        linkyaxes!(ax2, ax3)
        
        # Get ROI mask and create weighted map
        roi_mask = getROImask(data)
        roi_mask = rotr90(roi_mask)
        weighted_mask = zeros(size(roi_mask))
        tau_off_mask = zeros(size(roi_mask))
        
        for roi_id in keys(analysis.rois)
            if roi_id in sig_rois
                traces = filter(t -> t.channel == channel, analysis.rois[roi_id])
                if !isempty(traces)
                    # Average the maximum response across all stimuli
                    max_responses = [maximum(t.dfof) for t in traces]
                    avg_max_response = mean(max_responses)
                    weighted_mask[roi_mask .== roi_id] .= avg_max_response
                    
                    # Average tau_off across all stimuli
                    tau_offs = [t.fit_parameters[3] for t in traces 
                        if !isnothing(t.fit_parameters) && length(t.fit_parameters) >= 3]
                    if !isempty(tau_offs)
                        tau_off_mask[roi_mask .== roi_id] .= mean(tau_offs)
                    end
                end
            end
        end
        
        # 4. Weighted response map (bottom left)
        ax_weighted = Axis(gl_channel[2,1], title="Channel $channel Weighted Response (Averaged)",
            xlabel="X Position (μm)", ylabel="Y Position (μm)",
            aspect=DataAspect())
        
        if !isempty(sig_rois)
            hm_weighted = heatmap!(ax_weighted, xlims, ylims, weighted_mask,
                colormap=:viridis,
                colorrange=(0, maximum(weighted_mask)))
            Colorbar(gl_channel[2,2], hm_weighted, label="dF/F", width=15, vertical=true)
        else
            text!(ax_weighted, "No significant ROIs found", 
                position=(mean(xlims), mean(ylims)),
                align=(:center, :center),
                color=:red)
        end
        
        # 5. Tau Off map (bottom right)
        ax_tau = Axis(gl_channel[2,3], title="Channel $channel Tau Off (Averaged)",
            xlabel="X Position (μm)", ylabel="Y Position (μm)",
            aspect=DataAspect())
        
        if !isempty(sig_rois)
            hm_tau = heatmap!(ax_tau, xlims, ylims, tau_off_mask,
                colormap=:viridis,
                colorrange=(0, maximum(filter(!iszero, tau_off_mask))))
            Colorbar(gl_channel[2,4], hm_tau, label="Tau Off (s)", width=15, vertical=true)
        else
            text!(ax_tau, "No significant ROIs found", 
                position=(mean(xlims), mean(ylims)),
                align=(:center, :center),
                color=:red)
        end
        
        # Add stimulus time indicator if available
        if haskey(analysis.analysis_parameters, :delay_time)
            delay_time = analysis.analysis_parameters[:delay_time]
            vlines!(ax2, [delay_time], color=:red, linestyle=:dash, label="Stimulus")
            vlines!(ax3, [delay_time], color=:red, linestyle=:dash, label="Stimulus")
        end
        
        # Add legend only for mean response
        axislegend(ax2, position=:rt)
    end
    
    return fig
end 