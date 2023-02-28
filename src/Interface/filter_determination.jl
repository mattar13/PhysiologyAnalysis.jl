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
	import PhysiologyAnalysis: baseline_adjust!, truncate_data!, average_sweeps!
	import PhysiologyAnalysis: filter_data!, filter_data
	import PhysiologyAnalysis: dwt_filter!, cwt_filter!
	import PhysiologyAnalysis: fft_spectrum
	#Use pyplot? for plotting
	using PyPlot
	import PhysiologyAnalysis.plot_experiment
	#import PyPlot: plt
	import PhysiologyAnalysis.rcParams
	import PhysiologyAnalysis: wavelet, cDb2
	using ContinuousWavelets
	using Statistics
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
path = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\ERG\Organoids\2023_01_25_MattOrganoid\Organoid2_ILC42-3_d263\BrainPhys\BF\nd0_25p_0000.abf"

# ╔═╡ b400dd0c-5a40-4ee7-9116-7339939b7456
md"""
#### 2) Data information, channels, and Stimulus

a) Type in the channels you want to analyze here: 

"""

# ╔═╡ 971d6f11-2936-4d75-9641-36f81a94c2c4
channels = ["Vm_prime", "Vm_prime4"]

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
($(@bind freq_start NumberField(0.0:0.1:10.0; default = 0.1))
hz -> 
$(@bind freq_stop NumberField(0.0:1.0:1000.0; default = 55.0))
hz)

Pole
$(@bind pole_val NumberField(0:1:20; default = 8))
Ripple
$(@bind ripple_val NumberField(0.0:1.0:40.0; default = 15.0))
hz
Attenuation
$(@bind att_val NumberField(0.0:10.0:400.0; default = 100.0))
hz

$(@bind go Button("Filter"))
"""

# ╔═╡ 5fdc0c43-9454-495d-9b8a-e47313d178b2
begin
	go
	#This will act like the filtering pipeline. The pipeline goes as follows
	#A) Baseline, Truncate, Average
	data = readABF(path, channels=channels)
	baseline_adjust!(data)
	truncate_data!(data, t_pre=t_pre_pick, t_post=t_post_pick)
	
	
	# This data can be plotted seperately
	filtered_data = deepcopy(data)
	filter_data!(filtered_data,
	  mode=Symbol(filter_mode),
	  pole=pole_val,
	  method=Symbol(filter_method),
	  ripple=ripple_val,
	  attenuation=att_val,
	  freq_start=freq_start, freq_stop=freq_stop
	)
	dataSWEEP = deepcopy(data)
	filtered_dataSWEEP = deepcopy(filtered_data)

	 if avg_swp
	 	average_sweeps!(data)
		average_sweeps!(filtered_data)
		"Averaging Sweeps"
	 end;
	"Filtering functions"
end

# ╔═╡ 76025c46-2977-4300-8597-de04f313c667
md"""
### Plotting:

xlims:
(
$(@bind xlim1 NumberField(-1.0:0.01:10.000; default = -0.25)),
$(@bind xlim2 NumberField(-1.0:0.01:10.000; default = 5.0))
)

ylims:
(
$(@bind ylim1 NumberField(-10000.0:0.001:10000.000; default = minimum(data))),
$(@bind ylim2 NumberField(-10000.0:0.001:10000.000; default = maximum(data)))
)

"""

# ╔═╡ 2084267b-64a8-4d5b-8ce1-dd41f7feaa3a
begin #Plot individual traces
	fig1SWEEP, axSWEEP = plt.subplots(size(dataSWEEP,1), size(dataSWEEP, 3))
	fig1SWEEP.subplots_adjust(hspace=0.4, wspace=0.4, 
		left = 0.1
	)
	for swp in axes(dataSWEEP,1), ch in axes(dataSWEEP,3)
		axSWEEP[swp, ch].spines["left"].set_visible(false)
		axSWEEP[swp, ch].yaxis.set_visible(false)
		if swp != size(dataSWEEP,1)
			#remove the ylabel
			axSWEEP[swp, ch].spines["bottom"].set_visible(false) #We want the spine to fully
			axSWEEP[swp, ch].xaxis.set_visible(false)

		end
		plot_experiment(axSWEEP[swp, ch], dataSWEEP, channels = ch, 
			alpha=0.5, sweeps = swp
		)
		plot_experiment(axSWEEP[swp, ch], filtered_dataSWEEP, 
			color = :red, channels = ch, sweeps = swp
		)
				   	axSWEEP[swp, ch].set_xlim(xlim1, xlim2)
		   	axSWEEP[swp, ch].set_ylim(ylim1, ylim2)
	   axSWEEP[swp, ch].set_xlim(xlim1, xlim2)
	   #axSWEEP[swp, ch].set_ylim(ylim1, ylim2)
	end
	fig1SWEEP
end

# ╔═╡ 4daf9c1a-8b73-4049-8a13-b080e43f87cd
@bind do_plot Button("Plot")

# ╔═╡ 669c877b-efcf-4c6b-a70f-e14164abdbff
begin
	do_plot
	#if we want to adjust some params adjust them here and add them to the default
	#rcParams["font.size"] = 12.0
	#C) Plot the data using Plotly (can be interactive)
	fig1, ax = plt.subplots(size(data, 3), 2)
	fig1.subplots_adjust(hspace=0.4, wspace=0.4)
	
	freq, fft = fft_spectrum(data)
	freq_filt, fft_filt = fft_spectrum(filtered_data)
	fft_norm = abs.(fft)
	fft_norm_filt = abs.(fft_filt)
	
	# Plot the experiments
	if size(data, 3) > 1
	  for ch in 1:size(data, 3)
		   ax[1, ch].set_xlim(xlim1, xlim2)
		   ax[1, ch].set_ylim(ylim1, ylim2)
		   plot_experiment(ax[1, ch], data, channels=ch, c=:black, alpha=0.5)
		   plot_experiment(ax[1, ch], filtered_data, channels=ch, c=:red)
		   ax[1, ch].set_ylabel("$(data.chNames[ch]) ($(data.chUnits[ch]))")
		   ax[1, ch].set_xlabel("Time (s)")
		   ax[2, ch].plot(freq, fft_norm[1, :, ch], c=:black, alpha=0.5)
		   ax[2, ch].plot(freq, fft_norm_filt[1, :, ch], c=:red)
		   ax[2, ch].set_xscale("log")
		   ax[2, ch].set_yscale("log")
	  end
	else
	  ax[1, 1].set_xlim(xlim1, xlim2)
	  ax[1, 1].set_ylim(ylim1, ylim2)
	  plot_experiment(ax[1, 1], data, channels=1, c=:black, alpha=0.5)
	  plot_experiment(ax[1, 1], filtered_data, channel=1, c=:red)
	  ax[1, 1].set_ylabel("$(data.chNames[1]) ($(data.chUnits[1]))")
	  ax[1, 1].set_xlabel("Time (s)")
	
	  ax[2, 1].plot(freq, fft_norm[1, :, 1], c=:black, alpha=0.5)
	  ax[2, 1].plot(freq, fft_norm_filt[1, :, 1], c=:red)
	  ax[2, 1].set_xscale("log")
	end

	fig1
end

# ╔═╡ 9d5c1286-1a13-428a-b772-67887ae1f7c8
begin
     #Plot the DWT info
     fig2, ax2 = plt.subplots(3, size(data, 3))
     for ch in 1:size(data, 3)
          ax2[1, ch].set_xlim(xlim1, xlim2)
          ax2[1, ch].set_ylim(ylim1, ylim2)
          plot_experiment(ax2[1, ch], data, channels=ch, c=:black, alpha=0.5)
          ax2[1, ch].set_ylabel("$(data.chNames[ch]) ($(data.chUnits[ch]))")
          ax2[1, ch].set_xlabel("Time (s)")
     end
     fig2
end

# ╔═╡ 4b3cabb0-8de6-4364-b6dc-b23a1ff1af48
md"""
##### These filters are experimental and only work on small noisy waveforms
DWT Filter
$(@bind DWT_pick CheckBox(default = false))
CWT Filter
$(@bind CWT_pick CheckBox(default = false))

Period limits
$(@bind PER_lo NumberField(1:50; default = 2^0)) ->
$(@bind PER_hi NumberField(1:50; default = 2^4))

Power limits
$(@bind POW_lo NumberField(0.0:1.0; default = 0.0)) ->
$(@bind POW_hi NumberField(0.0:1.0; default = 1.0))
"""

# ╔═╡ dd65fd7a-b1c9-401d-8c37-149e2eaa3e5d
begin
	#settings
	wave = PhysiologyAnalysis.Morlet(0.50π)
	period_window = (PER_lo, PER_hi)
	power_window = (POW_lo, POW_hi)

	dataCWT, CWTi = cwt_filter(data*-1,
	  	wave=wave,
		period_window=period_window, power_window=power_window
	)

     #Plot the CWT info
    fig3, ax3 = plt.subplots(2, size(data, 3))
    for ch in 1:size(data, 3)
		CWT_PROC = PhysiologyAnalysis.CWTprocess(CWTi[1, :, :, ch])
		mu = mean(CWT_PROC[.!isinf.(CWT_PROC)])
		sig = std(CWT_PROC[.!isinf.(CWT_PROC)])
		levels = LinRange(mu-2 * sig, mu + 2 * sig, 10)
		freqs = log.(2, 1:size(CWTi, 3))
		
		plot_experiment(ax3[1, ch], data, channels=ch, c=:black, alpha=0.5)
		plot_experiment(ax3[1, ch], dataCWT*-1, channels=ch, c=:red)
		ax3[1, ch].set_xlim(xlim1, xlim2)
		ax3[1, ch].set_ylim(ylim1, ylim2)
		
		imF1 = ax3[2, ch].contourf(data.t, freqs, CWT_PROC,
		   cmap="seismic", levels=levels, extend = "both", 
		)
		cbari = fig3.colorbar(imF1, ax = ax3[2,ch], ticks=levels, aspect=5, location = "left")
		
		ax3[2, ch].set_ylabel("$(data.chNames[ch]) ($(data.chUnits[ch]))")
		ax3[2, ch].set_xlabel("Time (s)")
		ax3[2, ch].set_xlim(xlim1, xlim2)
	end
    fig3
end

# ╔═╡ e284711a-5f0b-4204-ae20-3173d8496255
begin
	fig1
	fig2
	fig3
	plt.close("all");
	clf();
end

# ╔═╡ Cell order:
# ╠═a442e068-06ef-4d90-9228-0a03bc6d9379
# ╟─e2fcae6f-d795-4258-a328-1aad5ea64195
# ╠═7b14b019-7545-4441-833b-f7e660c23dc6
# ╟─b400dd0c-5a40-4ee7-9116-7339939b7456
# ╠═971d6f11-2936-4d75-9641-36f81a94c2c4
# ╟─c8b4c855-64b8-4e18-9b2d-231260c67813
# ╟─5fdc0c43-9454-495d-9b8a-e47313d178b2
# ╟─76025c46-2977-4300-8597-de04f313c667
# ╟─2084267b-64a8-4d5b-8ce1-dd41f7feaa3a
# ╟─4daf9c1a-8b73-4049-8a13-b080e43f87cd
# ╠═669c877b-efcf-4c6b-a70f-e14164abdbff
# ╟─9d5c1286-1a13-428a-b772-67887ae1f7c8
# ╟─4b3cabb0-8de6-4364-b6dc-b23a1ff1af48
# ╟─dd65fd7a-b1c9-401d-8c37-149e2eaa3e5d
# ╠═e284711a-5f0b-4204-ae20-3173d8496255
