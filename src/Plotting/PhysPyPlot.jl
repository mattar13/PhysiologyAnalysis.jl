import PyPlot.plt

#We need an auxillary function 
function is_cmap(color)
    try
        plt.get_cmap(color)
        return true
    catch
        return false
    end
end

function plot_experiment(axis::PyObject, exp::Experiment;
    channels=1, sweeps = :all, 
    axes_off=false, yaxes_off=false, xaxes_off=false, #Change this, this is confusing
    color = :black, clims = (0.0, 1.0), #still want to figure out how this wil work
    kwargs...
)
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
        plot_experiment(ax::PyObject, exp::Experiment; channels=ch, kwargs...)
    end
end

function plot_experiment(exp::Experiment; kwargs...)
    fig, axis = plt.subplots(size(exp, 3))
    if size(exp,3) == 1
        plot_experiment(axis::PyObject, exp::Experiment; channels=1, kwargs...)
    else
        for (ch, ax) in enumerate(axis)
            plot_experiment(ax::PyObject, exp::Experiment; channels=ch, kwargs...)
        end
    end
    return fig
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

function plot_exp_fits(ax, df_EXPs::DataFrame, df_TRACEs::DataFrame;
          xlims = (-10.0, 200.0), ylims = (1.0, 3000.0)
     )
     #ax.set_title("P14 Rods NR")
     ax.set_xlim(xlims)
     ax.set_ylim(ylims)
     cmap = plt.get_cmap(:RdYlGn)
     for litter in eachrow(df_EXPs)
          println(litter.RSQ)
          RET = df_TRACEs |> @filter((_.Year, _.Month, _.Date, _.Number, _.Channel) == (litter.Year, litter.Month, litter.Date, litter.Number, litter.Channel)) |> DataFrame
          #println(litter)
          ax.scatter(
               RET.Response, RET.Total, 
               alpha = 0.5, 
               color = cmap(litter.RSQ), 
               marker = "o", label = "$(litter.Year)_$(litter.Month)_$(litter.Date)_$(litter.Number)"
          )
          # Find the STF data
          #STF_RES = p14RodsDR_EXP |> @filter((_.Year, _.Month, _.Date, _.Number, _.Channel) == (litter.Year, litter.Month, litter.Date, litter.Number, litter.Channel)) |> DataFrame
          rmax = litter.RMAX; k = litter.K; n = litter.N

          FIT = model(A_rng, [rmax, k, n])
          ax.plot(A_rng, FIT, color = cmap(litter.RSQ))
          ax.vlines(k, ymin=1.0, ymax = rmax/2, linestyle=(0, (5, 3)))#10 ^ (p14RodsNR_STF.param[1] / 2), color = color, )
          ax.hlines(rmax/2, xmin=-10.0, xmax = k, linestyle=(0, (5, 3)))
     end
     ax.legend(loc = "lower right")
     #ax.set_xscale("Log")
     ax.set_xlabel("a-wave (μV)", fontsize=10.0)
     ax.set_yscale("Log", base = 10)
     ax.set_ylabel("b-wave (μV)", fontsize=10.0)
     #return fig
end

function plot_exp_fits(df_EXPs::DataFrame, df_TRACEs::DataFrame; kwargs...)
    fig,ax = plt.subplots(1)
    plot_exp_fits(ax, df_EXPs, df_TRACEs; kwargs...)
end