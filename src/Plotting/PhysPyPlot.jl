#We need an auxillary function 
function is_cmap(color)
    try
        plt.get_cmap(color)
        return true
    catch
        return false
    end
end

function plot_experiment(axis::T, exp::Experiment;
    channels=1, sweeps = :all, 
    axes=true, yaxes=true, xaxes=true, #Change this, this is confusing
    xlims = nothing, ylims = nothing,
    color = :black, cvals = nothing, clims = (0.0, 1.0), #still want to figure out how this wil work
    ylabel = nothing, xlabel = nothing,
    kwargs...
) where T
    dataX, dataY = plot_prep(exp; channels=channels, sweeps = sweeps)
    if is_cmap(color)
        cmapI = plt.get_cmap(color)
        for (idx, swp) in enumerate(axes(dataY, 2))
            axis.plot(dataX, dataY[:, swp], c = cmapI(idx/size(dataY,2)), kwargs...)
        end
    else
        axis.plot(dataX, dataY; c = color, kwargs...)
    end
    axis.spines["top"].set_visible(false)
    axis.spines["right"].set_visible(false)
    if !(yaxes) || !(axes)
        axis.spines["left"].set_visible(false)
        axis.yaxis.set_visible(false)
    end
    if !(xaxes) || !(axes)
        axis.spines["bottom"].set_visible(false) #We want the spine to fully
        axis.xaxis.set_visible(false)
    end

    if !isnothing(xlabel)
        axis.set_xlabel(xlabel)
    else
        axis.set_xlabel("Time (s)")
    end

    if !isnothing(ylabel)
        axis.set_ylabel(ylabel)
    else
        axis.set_ylabel("$(exp.chNames[channels]) ($(exp.chUnits[channels]))")
    end

    if !isnothing(xlims)
        axis.set_xlim(xlims)
    end

    if !isnothing(ylims)
        axis.set_ylim(ylims)
    end
    #end
end

function plot_experiment(axis::Vector{T}, exp::Experiment; kwargs...) where T #Is going to be a py object
    #This is for if there are multiple axes
    for (ch, axis) in enumerate(axis)
        if ch == 1
            plot_experiment(axis::T, exp::Experiment; channels=ch, include_xlabel = false, kwargs...)
        else
            plot_experiment(axis::T, exp::Experiment; channels=ch, kwargs...)
        end
    end
end

function plot_experiment(exp::Experiment; layout = nothing, channels = nothing, kwargs...)
    if !isnothing(layout)
        plot_layout = layout
    elseif !isnothing(channels)
        if isa(channels, Vector{Int64})
            plot_layout = (length(channels))
        elseif isa(channels, Int64)
            plot_layout = 1
        end
    else
        plot_layout = (size(exp, 3))
    end
    #println(plot_layout)
    fig, axis = plt.subplots(plot_layout)
    if plot_layout == 1 || plot_layout == (1)
        plot_experiment(axis, exp::Experiment; channels=1, kwargs...)
    else
        for (ch, axis) in enumerate(axis)
            plot_experiment(axis, exp::Experiment; 
                channels=ch, 
                xlabel = ch == length(exp.chNames),
                kwargs...
            )
        end
    end
    return fig
end

function add_scalebar(axis, loc::Tuple{T,T}, dloc::Tuple{T,T};
    fontsize=10.0, lw=3.0,
    xlabeldist=30.0, ylabeldist=15.0,
    xunits="ms", yunits="μV",
    xconvert=1000.0, yconvert=1.0, #this converts the units from x and y labels. x should be in ms
    xround=true, yround=true,
    kwargs...
) where {T<:Real}
    x, y = loc
    dx, dy = dloc
    x0, xmax = plt.xlim()
    data_width = (xmax - x0)
    y0, ymax = plt.ylim()
    data_height = abs(ymax - y0)

    #println(data_width / xlabeldist) #debug options 
    #println(data_height / ylabeldist) #debug options
    axis.plot([x, x + dx], [y, y], c=:black, lw=lw; kwargs...) #add a vertical line
    axis.plot([x, x], [y, y + dy], c=:black, lw=lw; kwargs...)
    if yround
        yscale = round(Int64, dy * yconvert)
    else
        yscale = dy * yconvert
    end
    if xround
        xscale = round(Int64, dx * xconvert)
    else
        xscale = dx * xconvert
    end
    axis.annotate("$yscale $yunits", (x - (data_width / xlabeldist), y + dy / 2), va="center", ha="center", rotation="vertical", fontsize=fontsize)
    axis.annotate("$xscale $xunits", (x + dx / 2, y - (data_height / ylabeldist)), va="center", ha="center", rotation="horizontal", fontsize=fontsize)
end