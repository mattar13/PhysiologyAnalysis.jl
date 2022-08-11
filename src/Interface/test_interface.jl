### A Pluto.jl notebook ###
# v0.19.9

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
	using ePhys

	import ePhys: baseline_adjust!, truncate_data! 
	import ePhys: lowpass_filter!, highpass_filter!
	import ePhys: EI_filter!
	import ePhys: dwt_filter!, cwt_filter!
	#Use pyplot? for plotting
	using PyPlot
	#import PyPlot: plt
	import ePhys.rcParams
	pygui(true)
end

# ╔═╡ e2fcae6f-d795-4258-a328-1aad5ea64195
md"
#### Current Date: $(Dates.now())

# Analysis of experiments:

#### 1) Point to the experiment root: Contains .abf files

"

# ╔═╡ 7b14b019-7545-4441-833b-f7e660c23dc6
#enter in the file path of the file you would like to analyze here
path = raw"C:\Users\mtarc\OneDrive\Documents\GithubRepositories\ePhys\test\to_analyze.abf"

# ╔═╡ b400dd0c-5a40-4ee7-9116-7339939b7456
md"""
#### 2) Data information, channels, and Stimulus

a) Type in the channels you want to analyze here: 

"""

# ╔═╡ 971d6f11-2936-4d75-9641-36f81a94c2c4
channels = ["Vm_prime", "Vm_prime4"]

# ╔═╡ 05e38576-9650-4287-bac0-6d281db2ea9c
if path != ""
	data = readABF(path, channels = channels)
end;

# ╔═╡ c8b4c855-64b8-4e18-9b2d-231260c67813
md"""
#### 3) Filter the waveform

Pre Stimulus Time   
$(@bind t_pre_pick NumberField(0.0:0.1:2.0; default = 1.0))
Post Stimulus Time
$(@bind t_post_pick NumberField(0.0:0.1:10.0; default = 4.0))

Display Stimulus
$(@bind display_stim CheckBox(default = true))
"

Highpass
$(@bind highpass_pick CheckBox(default = false))   
$(@bind highpass_val NumberField(0.001:0.001:1.0; default = 1.0))
hz
Lowpass
$(@bind lowpass_pick CheckBox(default = true))
$(@bind lowpass_val NumberField(0.0:10.0:1000.0; default = 300.0))
hz

Adaptive Interference Bandpass
$(@bind adap_pick CheckBox(default = true))
$(@bind adap_val NumberField(0.0:1.0:1000.0; default = 100.0))

##### These filters are experimental and only work on small noisy waveforms
DWT Filter
$(@bind DWT_pick CheckBox(default = false))
CWT Filter
$(@bind CWT_pick CheckBox(default = false))

Wavelet limits
$(@bind WT_low_val NumberField(1:50; default = 1)) ->
$(@bind WT_hi_val NumberField(1:50; default = 9))

"""

# ╔═╡ 5fdc0c43-9454-495d-9b8a-e47313d178b2
begin
	#This will act like the filtering pipeline. The pipeline goes as follows
	#A) Baseline, Truncate, Average
	
	baseline_adjust!(data)
	truncate_data!(data, t_pre = t_pre_pick, t_post = t_post_pick)
	if highpass_pick
		highpass_filter!(data, freq = highpass_val)
	end

	if lowpass_pick
		lowpass_filter!(data, freq = lowpass_val)
	end

	if adap_pick
		EI_filter!(data, bandpass = adap_val)
	end
	
	if CWT_pick
	  cwt_filter!(data;
		   period_window = (WT_low_val, WT_hi_val)
	  )
	end
	
	if DWT_pick
	  dwt_filter!(data;
		   period_window = (WT_low_val, WT_hi_val)
	  )
	end
	data * 1000
	"Filtering functions"
end

# ╔═╡ 76025c46-2977-4300-8597-de04f313c667
md"""
### Plotting:

xlims:
(
$(@bind xlim1 NumberField(-1.0:0.01:10.000; default = -1.0)),
$(@bind xlim2 NumberField(-1.0:0.01:10.000; default = 4.0))
)

ylims:
(
$(@bind ylim1 NumberField(-10000.0:0.01:10000.000; default = minimum(minimum(data, dims = 2)))),
$(@bind ylim2 NumberField(-10000.0:0.01:10000.000; default = maximum(maximum(data, dims = 2))))
)

"""

# ╔═╡ 669c877b-efcf-4c6b-a70f-e14164abdbff
begin
	#if we want to adjust some params adjust them here and add them to the default
	#rcParams["font.size"] = 12.0
	#C) Plot the data using Plotly (can be interactive)
	fig, ax = plt.subplots(size(data,3)) #This 
	
	# Plot the experiments
	ax[1].set_xlim(xlim1, xlim2)
	ax[2].set_xlim(xlim1, xlim2)
	ax[1].set_ylim(ylim1, ylim2)
	ax[2].set_ylim(ylim1, ylim2)
	
	for ch in 1:size(data, 3)
		plot_experiment(ax[ch], data, channel = ch, c = :black)
		ax[ch].set_ylabel("$(data.chNames[ch]) ($(data.chUnits[ch]))")
	end
	#Set the time label on the last variable
	ax[size(data,3)].set_xlabel("Time (s)")
	fig
end

# ╔═╡ f7842e98-7c16-4043-a047-268a2f611e9b
plt.close("all")

# ╔═╡ 2ae9c8b5-473d-43db-90b2-1ca16f997c91
md"""
#### 3) Response Analysis
1) Response amplitudes
2) Saturated Response
3) Dim Response

Take note that the actual values I get are in mV
If it is an a-wave, then the values are also negative
"""

# ╔═╡ 85446d00-b763-4e32-9f97-20c449341e60
responses = -minimum(data, dims = 2)[:,1,:]

# ╔═╡ b6365068-c86e-4dad-a2a3-e5f89ceb4baa
saturated_response = -minimum(ePhys.saturated_response(data), dims = 1)

# ╔═╡ 0c548bc7-a362-49c7-9ea1-c88f1af52f4b
begin
	#Used for calculating the dim response
	norm_responses = responses./minimum(saturated_response)
	dim_idxs = findall(0.20 .< norm_responses .< 0.50)
	#used for calculating the time to peak (dim peak)
	peak_idxs = argmin(data, dims = 2)[:,1,:][dim_idxs] .|> Tuple
	peak_idxs = map(x -> x[2], peak_idxs)
end

# ╔═╡ 8d7dd9ad-6eb2-4310-916b-d0ecc146543e
dim_response = responses[dim_idxs]

# ╔═╡ 60d55064-66eb-4f34-ab28-e2b89561dc43
time_to_peak = data.t[peak_idxs]

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
rec_tau, gof = ePhys.recovery_time_constant(data, saturated_response)

# ╔═╡ 8e7c7ea8-4ea0-4102-8a2f-37db0ecefd25


# ╔═╡ e284711a-5f0b-4204-ae20-3173d8496255
plt.close("all")

# ╔═╡ 694fbd4f-51c7-4d2d-983f-7a1d29c57b45
clf()

# ╔═╡ Cell order:
# ╟─a442e068-06ef-4d90-9228-0a03bc6d9379
# ╟─e2fcae6f-d795-4258-a328-1aad5ea64195
# ╠═7b14b019-7545-4441-833b-f7e660c23dc6
# ╟─b400dd0c-5a40-4ee7-9116-7339939b7456
# ╠═971d6f11-2936-4d75-9641-36f81a94c2c4
# ╟─05e38576-9650-4287-bac0-6d281db2ea9c
# ╟─c8b4c855-64b8-4e18-9b2d-231260c67813
# ╟─5fdc0c43-9454-495d-9b8a-e47313d178b2
# ╟─76025c46-2977-4300-8597-de04f313c667
# ╟─669c877b-efcf-4c6b-a70f-e14164abdbff
# ╟─f7842e98-7c16-4043-a047-268a2f611e9b
# ╟─2ae9c8b5-473d-43db-90b2-1ca16f997c91
# ╠═0c548bc7-a362-49c7-9ea1-c88f1af52f4b
# ╟─85446d00-b763-4e32-9f97-20c449341e60
# ╟─b6365068-c86e-4dad-a2a3-e5f89ceb4baa
# ╟─8d7dd9ad-6eb2-4310-916b-d0ecc146543e
# ╟─60d55064-66eb-4f34-ab28-e2b89561dc43
# ╟─66a7a194-a23a-4f1d-a6ad-c2a7e8bee1e6
# ╟─3c110c6e-6909-4d24-95b3-5a5f092e9575
# ╠═505ce94d-6a62-44cc-94ef-7a519425a250
# ╠═8e7c7ea8-4ea0-4102-8a2f-37db0ecefd25
# ╠═e284711a-5f0b-4204-ae20-3173d8496255
# ╠═694fbd4f-51c7-4d2d-983f-7a1d29c57b45
