### A Pluto.jl notebook ###
# v0.19.19

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ a442e068-06ef-4d90-9228-0a03bc6d9379
#This section will use the packages loaded in the environment
begin
	using Pkg
	Pkg.activate("../../")
	using Dates, PlutoUI
	using PhysiologyAnalysis
	using DataFrames, Query
	import PhysiologyAnalysis: baseline_adjust!, truncate_data! , average_sweeps!
	import PhysiologyAnalysis: filter_data!, filter_data
	import PhysiologyAnalysis: EI_filter!
	import PhysiologyAnalysis: dwt_filter!, cwt_filter!
	#Use pyplot? for plotting
	using PyPlot
	#import PyPlot: plt
	import PhysiologyAnalysis.rcParams
	pygui(true)
end

# ╔═╡ e2fcae6f-d795-4258-a328-1aad5ea64195
md"
#### Current Date: $(Dates.now())

# Analysis of experiments:

#### 1) Point to the experiment root: 
a) Contains .abf files of a-wave, b-wave, and glial components

"

# ╔═╡ 7b14b019-7545-4441-833b-f7e660c23dc6
begin
	#enter in the file path of the file you would like to analyze here
	exp_root = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\ERG\Retinoschisis\2022_04_12_RCAdult\Mouse1_Adult_R141C\BaCl_LAP4"
	experiment_paths = exp_root |> parseABF
end

# ╔═╡ b400dd0c-5a40-4ee7-9116-7339939b7456
md"""
#### 2) Data information, channels, and Stimulus

a) Type in the channels you want to analyze here: 

"""

# ╔═╡ 971d6f11-2936-4d75-9641-36f81a94c2c4
channels = ["Vm_prime", "Vm_prime4"]

# ╔═╡ 05e38576-9650-4287-bac0-6d281db2ea9c
begin
	data = readABF(experiment_paths, channels = channels)
	data_filter!(data, 
		avg_swp = false,
	) #Use the default data filter
end

# ╔═╡ 76025c46-2977-4300-8597-de04f313c667
md"""
### Plotting:

xlims:
(
$(@bind xlim1 NumberField(-1.0:0.01:10.000; default = -0.25)),
$(@bind xlim2 NumberField(-1.0:0.01:10.000; default = 2.0))
)

ylims:
(
$(@bind ylim1 NumberField(-10000.0:0.01:10000.000; default = minimum(data))),
$(@bind ylim2 NumberField(-10000.0:0.01:10000.000; default = maximum(data)))
)

"""

# ╔═╡ 669c877b-efcf-4c6b-a70f-e14164abdbff
begin
	fig, ax = plt.subplots(size(data,3)) #This 
	
	# Plot the experiments
	if size(data,3) > 1
		for ch in 1:size(data, 3)
			ax[ch].set_xlim(xlim1, xlim2)
			ax[ch].set_ylim(ylim1, ylim2)
			plot_experiment(ax[ch], data, channels = ch, c = :black)
			ax[ch].set_ylabel("$(data.chNames[ch]) ($(data.chUnits[ch]))")
		end
		ax[size(data,3)].set_xlabel("Time (s)")
	else
		ax.set_xlim(xlim1, xlim2)
		ax.set_ylim(ylim1, ylim2)
		plot_experiment(ax, data, channels = 1, c = :black)
		ax.set_ylabel("$(data.chNames[1]) ($(data.chUnits[1]))")
		ax.set_xlabel("Time (s)")
	end	
	fig
end

# ╔═╡ 7662c0f6-da9b-448b-abde-9b20e15c53ee
plt.close("all"); clf()

# ╔═╡ 2ae9c8b5-473d-43db-90b2-1ca16f997c91
md"""
#### 3) Response Analysis
1) Response amplitudes
2) Saturated Response
3) Dim Response
4) Time to peak

#### 4) Fitting IR curves
R parameter:
RMIN $(@bind RMIN NumberField(0.1:0.01:10000; default = 0.1)) ->
INITIAL R $(@bind R0 NumberField(0.0:0.01:10000; default = 100.0)) -> 
RMAX $(@bind RMAX NumberField(0.0:0.01:10000; default = -minimum(data)))


k parameter:
kMIN $(@bind kMIN NumberField(0.0:0.01:10000; default = 0.1)) ->
INITIAL k $(@bind k0 NumberField(1.0:0.01:10000; default = 5000)) -> 
kMAX $(@bind kMAX NumberField(0.0:0.01:10000; default = 10e6))

n parameter:
nMIN $(@bind nMIN NumberField(0.1:0.01:10000; default = 0.1)) ->
INITIAL n $(@bind n0 NumberField(0.1:0.01:10000; default = 2.0)) -> 
nMAX $(@bind nMAX NumberField(0.1:0.01:10000; default = 10.0))
"""

# ╔═╡ 2b8e0811-2290-4488-8510-1c3476b42d20
begin
	dfs = DataFrame[]
	dfsSUMMARY = DataFrame[]
	for i in 1:size(data,3)
		df = DataFrame()
		dfSUMMARY = DataFrame()
		
		info = PhysiologyAnalysis.createDatasheet(experiment_paths, filename = nothing);
		df.SweepN = 1:size(data,1)
		df.channel .= data.chNames[i]
		df.Photons = info[:, :Photons]
		df.ND = info[:, :ND]
		df.Intensity = info[:,:Percent]
		df.StimTiem = info[:, :Stim_Time]
		df.Minima = vec(-minimum(data, dims = 2)[:,:,i])
		df.SaturatedResponse = vec(-PhysiologyAnalysis.saturated_response(data)[:,i])

		df.isSaturated = df.Minima .> df.SaturatedResponse
		#Calculate the dim responses
		agmin = map(id -> data.t[id[2]], argmin(data, dims = 2))
		df.TPeak = agmin[:,1,i].*1000
		df = df |> @orderby(_.Photons) |> DataFrame
		
		#calculate the dim responses
		isDim = fill(false, size(data,1))
		norm_responses = df.SaturatedResponse ./ maximum(df.SaturatedResponse)
		dim_idxs = findall(0.20 .< norm_responses .< 0.50)
		isDim[dim_idxs] .= true
		df.isDim .= isDim

		dfSUMMARY.CHANNEL = [data.chNames[i]]
		#dfSUMMARY.MINIMA = [maximum(df.Minima)]
		dfSUMMARY.RMAX = [maximum(df.SaturatedResponse)]
		if !isnothing(dim_idxs)
			dims = df[dim_idxs, :SaturatedResponse]
		else
			dims = 0.0
		end
		dfSUMMARY.RDIM = [minimum(dims)]
		dfSUMMARY.TIME_TO_PEAK = [maximum(df.TPeak)]

		#Lets fit the intensity response data
		fit = PhysiologyAnalysis.IRfit(
			df.Photons, df.SaturatedResponse,
			rmin = RMIN, r = R0, rmax = RMAX,
			kmin = kMIN, k = k0, kmax = kMAX, 
			nmin = nMIN, n = n0, nmax = nMAX		
		)
		#println(fit.param)
		dfSUMMARY.RMAX_FIT = [fit.param[1]]
		dfSUMMARY.K_FIT = [fit.param[2]]
		dfSUMMARY.N_FIT = [fit.param[3]]
		push!(dfs, df)
		push!(dfsSUMMARY, dfSUMMARY)
	end
	results = vcat(dfs...)
end

# ╔═╡ 5684843d-622f-487d-aa06-d5fd7bff1689
begin
	figALL, ax2 = plt.subplots(2, 3) #This 
	fit_I = 10 .^ LinRange(-2, 5, 1000)
	# Plot the experiments
	for ch in 1:size(data, 3)
		ax2[ch, 1].set_xlim(xlim1, xlim2)
		ax2[ch, 1].set_ylim(ylim1, ylim2)

		ax2[ch, 2].set_ylim(ylim1, ylim2)
		
		plot_experiment(ax2[ch, 1], data, channels = ch, c = :black)
		ax2[ch, 1].set_ylabel("$(data.chNames[ch]) ($(data.chUnits[ch]))")
		
		#query the channel
		PHOT = dfs[ch].Photons
		RESP = dfs[ch].SaturatedResponse
		TPEAK = dfs[ch].TPeak
		ax2[ch, 1].hlines(-RESP, xlim1, xlim2, color = :Red)
		#ax2[ch, 1].vlines(TPEAK./1000, ylim1, ylim2, color = :Blue)
		
		#Plot the intensity response curve
		ax2[ch,2].scatter(PHOT, -RESP, color = :Red)
		ax2[ch, 2].set_xscale("log")
		
		fitRMAX = dfsSUMMARY[ch].RMAX_FIT
		fitK = dfsSUMMARY[ch].K_FIT
		fitN = dfsSUMMARY[ch].N_FIT
		
		fit_R = -fitRMAX .* IR.(fit_I, fitK, fitN)
		ax2[ch, 2].plot(fit_I, fit_R, c = :Black, lw = 2.0)
		ax2[ch, 2].vlines([fitK], ylim1, -fitRMAX*0.50, 
			label = "K = $(round(fitK[1], digits = 1))",
			color = :Red
		)
		ax2[ch, 2].legend(loc = "lower right")
		#Plot the intensity time to peak
		ax2[ch, 3].scatter(PHOT, TPEAK, color = :Blue)
		ax2[ch, 3].set_xscale("log")
	end
	ax2[2, 1].set_xlabel("Time (s)")
	
	#Set the time label on the last variable
	figALL
end

# ╔═╡ 670d5073-19c5-45f7-950e-5408b1530535
resultsSUMMARY = vcat(dfsSUMMARY...)

# ╔═╡ aef0b08b-0aa3-4b49-8692-a9ea23a030db
begin
	figSUMMARY, ax3 = plt.subplots(2) #This 

	# Plot the experiments
	for ch in 1:size(data, 3)
		ax3[ch].set_xlim(xlim1, xlim2)
		ax3[ch].set_ylim(ylim1, ylim2)
		plot_experiment(ax3[ch], data, channel = ch, c = :black)
		ax3[ch].set_ylabel("$(data.chNames[ch]) ($(data.chUnits[ch]))")
	end
	ax3[2].set_xlabel("Time (s)")
	
	#Set the time label on the last variable
	figSUMMARY
end

# ╔═╡ 66a7a194-a23a-4f1d-a6ad-c2a7e8bee1e6
dominant_recovery = maximum(percent_recovery_interval(data, -saturated_response), dims = 1)

# ╔═╡ 3c110c6e-6909-4d24-95b3-5a5f092e9575
begin
	ax[1].hlines(-saturated_response[1], xmin = xlim1, xmax = xlim2, color = :green, 
		lw = 1.0
	)
	rmax = ax[2].hlines(-saturated_response[2], xmin = xlim1, xmax = xlim2, color = :green, 
		lw = 1.0
	)

	ax[1].hlines(-dim_response[1], xmin = xlim1, xmax = xlim2, color = :red, 
		lw = 1.0
	)
	rdim = ax[2].hlines(-dim_response[2], xmin = xlim1, xmax = xlim2, color = :red, 
		lw = 1.0
	)

	ax[1].vlines(time_to_peak[1], ymin = ylim1, ymax = ylim2, color = :blue, 
		lw = 1.0
	)
	tpeak = ax[2].vlines(time_to_peak[2], ymin = ylim1, ymax = ylim2, color = :blue, 
		lw = 1.0
	)
	ax[1].hlines(-saturated_response[1] .* 0.50, 
		xmin = 0.0, xmax = dominant_recovery[1], color = :cyan, lw =1.0, linestyle = :dashed
	)
	tdom = ax[2].hlines(-saturated_response[2] .* 0.50, 
		xmin = 0.0, xmax = dominant_recovery[2], color = :cyan, lw =1.0, linestyle = :dashed
	)
	
	ax[2].legend(
		handles = [rmax, rdim, tpeak, tdom], 
		labels = ["Rmax", "Rdim", "Time to peak", "Dominant Tau"], 
		loc = "lower right"
	)
	
	fig
end

# ╔═╡ 505ce94d-6a62-44cc-94ef-7a519425a250
rec_tau, gof = PhysiologyAnalysis.recovery_time_constant(data, saturated_response)

# ╔═╡ 8e7c7ea8-4ea0-4102-8a2f-37db0ecefd25


# ╔═╡ e284711a-5f0b-4204-ae20-3173d8496255
plt.close("all"); clf()

# ╔═╡ Cell order:
# ╟─a442e068-06ef-4d90-9228-0a03bc6d9379
# ╟─e2fcae6f-d795-4258-a328-1aad5ea64195
# ╠═7b14b019-7545-4441-833b-f7e660c23dc6
# ╟─b400dd0c-5a40-4ee7-9116-7339939b7456
# ╠═971d6f11-2936-4d75-9641-36f81a94c2c4
# ╟─05e38576-9650-4287-bac0-6d281db2ea9c
# ╟─76025c46-2977-4300-8597-de04f313c667
# ╟─669c877b-efcf-4c6b-a70f-e14164abdbff
# ╟─7662c0f6-da9b-448b-abde-9b20e15c53ee
# ╟─2ae9c8b5-473d-43db-90b2-1ca16f997c91
# ╟─2b8e0811-2290-4488-8510-1c3476b42d20
# ╟─5684843d-622f-487d-aa06-d5fd7bff1689
# ╟─670d5073-19c5-45f7-950e-5408b1530535
# ╟─aef0b08b-0aa3-4b49-8692-a9ea23a030db
# ╠═66a7a194-a23a-4f1d-a6ad-c2a7e8bee1e6
# ╟─3c110c6e-6909-4d24-95b3-5a5f092e9575
# ╠═505ce94d-6a62-44cc-94ef-7a519425a250
# ╠═8e7c7ea8-4ea0-4102-8a2f-37db0ecefd25
# ╠═e284711a-5f0b-4204-ae20-3173d8496255
