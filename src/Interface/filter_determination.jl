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

	import ePhys: baseline_adjust!, truncate_data! , average_sweeps!
	import ePhys: filter_data!, filter_data
	import ePhys: dwt_filter!, cwt_filter!
	import ePhys: fft_spectrum
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
path = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\ERG\Organoids\2022_08_27_Organoids\Organoid1_Agar\BrainPhys_9cRAL\BF\nd0_100p_0000.abf"

# ╔═╡ b400dd0c-5a40-4ee7-9116-7339939b7456
md"""
#### 2) Data information, channels, and Stimulus

a) Type in the channels you want to analyze here: 

"""

# ╔═╡ 971d6f11-2936-4d75-9641-36f81a94c2c4
channels = ["Vm_prime4"]

# ╔═╡ 05e38576-9650-4287-bac0-6d281db2ea9c
if path != ""
	data = readABF(path, channels = channels)
	data * 1000; #Putting this here will make sure we reopen the data every time we filter
end;

# ╔═╡ c8b4c855-64b8-4e18-9b2d-231260c67813
md"""
#### 3) Filter the waveform

Pre Stimulus Time   
$(@bind t_pre_pick NumberField(0.0:0.1:2.0; default = 1.0))
Post Stimulus Time
$(@bind t_post_pick NumberField(0.0:0.1:10.0; default = 4.0))

Average Traces
$(@bind avg_swp CheckBox(default = true))

Filter Mode:
$(@bind filter_mode Select([
	"Highpass", 
	"Lowpass", 
	"Bandpass", 
	"Bandstop"
], default = "Lowpass"))
Method:
$(@bind filter_method Select([
	"Butterworth", 
	"Chebyshev1", 
	"Chebyshev2", 
	"Elliptic"
], default = "Chebyshev2"))

Freq Range: 
($(@bind freq_start NumberField(0.0:0.01:10.0; default = 1.0))
hz -> 
$(@bind freq_stop NumberField(0.0:10.0:1000.0; default = 300.0))
hz)

Ripple
$(@bind ripple_val NumberField(0.0:0.10:40.0; default = 10.0))
hz
Attenuation
$(@bind att_val NumberField(0.0:0.10:400.0; default = 100.0))
hz
Pole
$(@bind pole_val NumberField(0:1:8; default = 8))



##### These filters are experimental and only work on small noisy waveforms
DWT Filter
$(@bind DWT_pick CheckBox(default = false))
CWT Filter
$(@bind CWT_pick CheckBox(default = false))

Wavelet limits
$(@bind WT_low_val NumberField(1:50; default = 1)) ->
$(@bind WT_hi_val NumberField(1:50; default = 9))

"""

# ╔═╡ d6d59ac3-e6e8-4b49-9ac6-2a17cf72f30b
@bind go Button("Filter")

# ╔═╡ 5fdc0c43-9454-495d-9b8a-e47313d178b2
begin
	go
	#This will act like the filtering pipeline. The pipeline goes as follows
	#A) Baseline, Truncate, Average

	baseline_adjust!(data)
	truncate_data!(data, t_pre = t_pre_pick, t_post = t_post_pick)
	if avg_swp
		average_sweeps!(data)
	end

	# This data can be plotted seperately
	filtered_data = deepcopy(data)
	if filter_mode == "Highpass"
		filter_data!(filtered_data, 
			mode = :Highpass, 
			pole = pole_val,
			method = Symbol(filter_method), 
			ripple = ripple_val,
			attenuation = att_val, 
			freq = freq_start
		)
	elseif filter_mode == "Lowpass"
		filter_data!(filtered_data, 
			mode = :Lowpass, 
			pole = pole_val,
			method = Symbol(filter_method),
			ripple = ripple_val,
			attenuation = att_val, 
			freq = freq_stop
		)
	elseif filter_mode == "Bandpass"
		filter_data!(filtered_data, 
			mode = :Bandpass, 
			pole = pole_val,
			method = Symbol(filter_method),
			ripple = ripple_val,
			attenuation = att_val, 
			freq = freq_start,
			freqstop = freq_stop
		)
	elseif filter_mode == "Bandstop"
		filter_data!(filtered_data, 
			mode = :Bandstop, 
			pole = pole_val,
			method = Symbol(filter_method),
			ripple = ripple_val,
			attenuation = att_val, 
			freq = freq_start,
			freqstop = freq_stop
		)
	end
	
	if CWT_pick
		cwt = cwt_filter!(filtered_data;
		   period_window = (WT_low_val, WT_hi_val), return_cwt = true
		)
		println(cwt |> size)
	else
		cwt = nothing
	end
	
	if DWT_pick
		dwt = dwt_filter!(filtered_data;
		   period_window = (WT_low_val, WT_hi_val), return_dwt = true
		)
		println(dwt |> size)
	else
		dwt = nothing
	end
	
	"Filtering functions"
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
$(@bind ylim1 NumberField(-10000.0:0.01:10000.000; default = -10)),
$(@bind ylim2 NumberField(-10000.0:0.01:10000.000; default = 10))
)

"""

# ╔═╡ 669c877b-efcf-4c6b-a70f-e14164abdbff
begin
	#if we want to adjust some params adjust them here and add them to the default
	#rcParams["font.size"] = 12.0
	#C) Plot the data using Plotly (can be interactive)
	fig, ax = plt.subplots(size(data,3), 2) 
	freq, fft = fft_spectrum(data)
	freq_filt, fft_filt = fft_spectrum(filtered_data)
	fft_norm = abs.(fft)
	fft_norm_filt = abs.(fft_filt)
	
	# Plot the experiments
	if size(data,3) > 1
		for ch in 1:size(data, 3)
			ax[1, ch].set_xlim(xlim1, xlim2)
			ax[1, ch].set_ylim(ylim1, ylim2)
			plot_experiment(ax[1, ch], data, channel = ch, c = :black, alpha = 0.5)
			plot_experiment(ax[1, ch], filtered_data, channel = ch, c = :red)
			ax[1, ch].set_ylabel("$(data.chNames[ch]) ($(data.chUnits[ch]))")

			ax[2, ch].plot(freq, fft_norm[1, :, ch], c = :black, alpha = 0.5)
			ax[2, ch].plot(freq, fft_norm_filt[1, :, ch], c = :red)
		end
		ax[1, size(data,3)].set_xlabel("Time (s)")
	else
		ax[1, 1].set_xlim(xlim1, xlim2)
		ax[1, 1].set_ylim(ylim1, ylim2)
		plot_experiment(ax[1, 1], data, channel = 1, c = :black, alpha = 0.5)
		plot_experiment(ax[1, 1], filtered_data, channel = 1, c = :red)
		ax[1, 1].set_ylabel("$(data.chNames[1]) ($(data.chUnits[1]))")
		ax[1, 1].set_xlabel("Time (s)")


		ax[2, 1].plot(freq, fft_norm[1, :, 1], c = :black, alpha = 0.5)
		ax[2, 1].plot(freq, fft_norm_filt[1, :, 1], c = :red)
	end
	#Lets plot the frequency now
	#Set the time label on the last variable
	plt.xscale(:log)
	plt.yscale(:log)
	
	fig
end

# ╔═╡ 9d5c1286-1a13-428a-b772-67887ae1f7c8
begin
	#Plot the DWT info
	fig2, ax2 = plt.subplots(size(data,3), 3) 

	fig2
end

# ╔═╡ 7662c0f6-da9b-448b-abde-9b20e15c53ee
plt.close("all"); clf()

# ╔═╡ e284711a-5f0b-4204-ae20-3173d8496255
plt.close("all"); clf()

# ╔═╡ Cell order:
# ╟─a442e068-06ef-4d90-9228-0a03bc6d9379
# ╟─e2fcae6f-d795-4258-a328-1aad5ea64195
# ╠═7b14b019-7545-4441-833b-f7e660c23dc6
# ╟─b400dd0c-5a40-4ee7-9116-7339939b7456
# ╠═971d6f11-2936-4d75-9641-36f81a94c2c4
# ╠═05e38576-9650-4287-bac0-6d281db2ea9c
# ╟─c8b4c855-64b8-4e18-9b2d-231260c67813
# ╟─d6d59ac3-e6e8-4b49-9ac6-2a17cf72f30b
# ╠═5fdc0c43-9454-495d-9b8a-e47313d178b2
# ╟─76025c46-2977-4300-8597-de04f313c667
# ╟─669c877b-efcf-4c6b-a70f-e14164abdbff
# ╠═9d5c1286-1a13-428a-b772-67887ae1f7c8
# ╠═7662c0f6-da9b-448b-abde-9b20e15c53ee
# ╠═e284711a-5f0b-4204-ae20-3173d8496255
