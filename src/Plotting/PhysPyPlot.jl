function plot_experiment(axis::PyObject, exp::Experiment;
    channel=1, axes_off=false, yaxes_off=false, xaxes_off=false,
    #cmap = nothing, cmap_direction = :sweeps,
    kwargs...
)
    #if !isnothing(cmap) && cmap_direction == :sweeps#cmap needs to be an array of colors
    #    for swp in 1:size(exp, 1)
    #        println(swp)
    #    end
    #else

    axis.plot(plot_prep(exp; channel=channel)...; kwargs...)
    if yaxes_off || axes_off
        axis.spines["left"].set_visible(false)
        axis.yaxis.set_visible(false)
    end
    if xaxes_off || axes_off
        axis.spines["bottom"].set_visible(false) #We want the spine to fully
        axis.xaxis.set_visible(false)
    end
    #end
end

function plot_experiment(axis::Vector{PyObject}, exp::Experiment; kwargs...)
    for (ch, ax) in enumerate(axis)
        plot_experiment(ax::PyObject, exp::Experiment; channel=ch, kwargs...)
    end
end

function add_scalebar(ax, loc::Tuple{T,T}, dloc::Tuple{T,T};
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
    ax.plot([x, x + dx], [y, y], c=:black, lw=lw; kwargs...) #add a vertical line
    ax.plot([x, x], [y, y + dy], c=:black, lw=lw; kwargs...)
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
    ax.annotate("$yscale $yunits", (x - (data_width / xlabeldist), y + dy / 2), va="center", ha="center", rotation="vertical", fontsize=fontsize)
    ax.annotate("$xscale $xunits", (x + dx / 2, y - (data_height / ylabeldist)), va="center", ha="center", rotation="horizontal", fontsize=fontsize)
end