using GLMakie
using Statistics
using ElectroPhysiology, PhysiologyAnalysis

"""
    plot_roi_analysis(data::Experiment{TWO_PHOTON}, analysis::ROIAnalysis; 
        stim_idx::Int=1, channel_idx::Int=1)

Create a comprehensive visualization of ROI analysis results including:
1. ROI map showing significant ROIs weighted by their amplitude, overlaid on max projection
2. Individual dF/F traces for significant ROIs
3. Mean trace with standard deviation ribbon

Returns a Figure object containing all plots.
"""
function plot_roi_analysis(data::Experiment{TWO_PHOTON, T}, analysis::ROIAnalysis;
    stim_idx::Int=1, channel_idx::Int=1) where {T <: Real}
    
    # Get significant ROIs and their traces
    sig_rois = get_significant_rois(analysis, stim_idx)
    
    # Create figure with layout
    fig = Figure(size=(1200, 600))
    
    # Create a layout with a square aspect ratio for the heatmap
    gl = fig[1, 1] = GridLayout()
    gl_traces = fig[1, 2] = GridLayout()
    colsize!(fig.layout, 1, Relative(0.4)) # Make the heatmap column slightly smaller
    
    # 1. ROI amplitude map
    ax1 = Axis(gl[1,1], title="ROI Response Map",
        xlabel="X Position (μm)", ylabel="Y Position (μm)",
        aspect=DataAspect()) # Force square aspect ratio
    
    # Get max projection of original data
    max_proj = zProject(data)
    # Normalize the max projection for overlay
    norm_max_proj = (max_proj .- minimum(max_proj)) ./ (maximum(max_proj) - minimum(max_proj))
    
    # Get ROI mask and create weighted map
    roi_mask = getROImask(data)
    weighted_mask = zeros(size(roi_mask))
    
    for roi_id in keys(analysis.rois)
        if roi_id in sig_rois
            trace = analysis.rois[roi_id][stim_idx]
            max_response = maximum(trace.dfof)
            weighted_mask[roi_mask .== roi_id] .= max_response
        end
    end
    
    # Plot max projection as background
    xlims = data.HeaderDict["xrng"]
    ylims = data.HeaderDict["yrng"]
    
    # Create alpha mask for ROIs
    alpha_mask = Float64.(weighted_mask .> 0)
    
    # Plot both layers
    image!(ax1, xlims, ylims, norm_max_proj, 
        colormap=:grays)
    hm = heatmap!(ax1, xlims, ylims, weighted_mask, 
        colormap=(:viridis, 0.7), # Make ROI overlay slightly transparent
        colorrange=(0, maximum(weighted_mask)))
    Colorbar(gl[1,2], hm, label="ΔF/F", width=15)
    
    # 2. Individual traces
    ax2 = Axis(gl_traces[1,1], title="Individual ROI Traces",
        xlabel="Time (s)", ylabel="ΔF/F")
    
    # Plot each significant ROI trace
    colors = cgrad(:Set2_8, length(sig_rois), categorical=true)
    for (i, roi_id) in enumerate(sig_rois)
        trace = analysis.rois[roi_id][stim_idx]
        lines!(ax2, trace.t_series, trace.dfof, 
            color=colors[i], alpha=0.5,
            label="ROI $roi_id")
    end
    
    # 3. Mean trace with std ribbon
    ax3 = Axis(gl_traces[2,1], title="Mean Response ± STD",
        xlabel="Time (s)", ylabel="ΔF/F")
    
    # Get all significant traces and calculate statistics
    sig_traces = [analysis.rois[roi][stim_idx].dfof for roi in sig_rois]
    t_series = analysis.rois[first(sig_rois)][stim_idx].t_series
    
    mean_trace = mean(sig_traces)
    std_trace = std(sig_traces)
    
    # Plot mean and standard deviation
    band!(ax3, t_series, 
        mean_trace .- std_trace, 
        mean_trace .+ std_trace,
        color=(:blue, 0.3))
    lines!(ax3, t_series, mean_trace, 
        color=:blue, linewidth=2,
        label="Mean (n=$(length(sig_rois)))")
    
    # Add stimulus time indicator if available
    if haskey(analysis.analysis_parameters, :delay_time)
        delay_time = analysis.analysis_parameters[:delay_time]
        vlines!(ax2, [delay_time], color=:red, linestyle=:dash, label="Stimulus")
        vlines!(ax3, [delay_time], color=:red, linestyle=:dash, label="Stimulus")
    end
    
    # Add legends
    axislegend(ax3, position=:rt)
    
    # Link x-axes of the traces
    linkyaxes!(ax2, ax3)
    linkxaxes!(ax2, ax3)
    
    return fig
end 