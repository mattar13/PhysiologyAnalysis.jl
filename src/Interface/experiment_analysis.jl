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

# ╔═╡ 8ec03fe0-866b-11ed-23a6-23f888e1717a
begin
	using Pkg
	Pkg.activate("../../")
	using Dates, PlutoUI
	using ePhys
	using DataFrames, XLSX, Query
	import ePhys: run_A_wave_analysis, run_B_wave_analysis, run_G_wave_analysis
	import ePhys: baseline_adjust!, truncate_data! , average_sweeps!
	import ePhys: filter_data!, filter_data
	import ePhys: HILL_MODEL
	#Use pyplot? for plotting
	using PyPlot
	#import PyPlot: plt
	import ePhys.rcParams
	pygui(true)
end

# ╔═╡ 9d51d150-3339-435c-ad76-59b81bdcf305
md"
#### Current Date: $(Dates.now())

# Analysis of experiments:

#### 1) Point to the experiment root: 
a) Contains .abf files of a-wave, b-wave, and glial components
#### 2) Point out the channel(s) you want to analyze

"

# ╔═╡ 1fb1daf9-59fd-42aa-b30b-bc5413a6cc67
#enter in the file path of the file you would like to analyze here
exp_root = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\ERG\Retinoschisis\2022_07_11_Adult\Mouse3_Adult_RS1KO"

# ╔═╡ bc4d42b1-6848-4c61-99ce-1e437c08497b
channels = ["Vm_prime", "Vm_prime4"]; #Specify which channels you are using

# ╔═╡ 3527f4ce-d1fb-46ac-a6a7-e9d3d05a03de
experiment_paths = exp_root |> parseABF #This pulls out all the paths

# ╔═╡ 022ec5ec-5b0f-4876-8ff4-18ed6f493f2b
md"
b) if you want to save the file, point to the savefile name, otherwise call this nothing

"

# ╔═╡ 01fab469-e7c3-4573-9e75-d2d13c7ff3de
savefile_name = nothing

# ╔═╡ 1cbc25be-bf12-43ac-a545-89aeab07ff93
all_files = createDatasheet(experiment_paths; filename = savefile_name)

# ╔═╡ 88ef46c3-7575-4a18-a8b5-15a0c2152c6a
md"""
# Analysis of the a-wave
Analysis Window end = $(@bind t_post_a NumberField(-1.0:0.01:10.000; default = 2.0))s
"""

# ╔═╡ 734c8f3e-7092-48fd-9784-391bde74fc08
begin	
	trace_A, exp_A, cond_A = run_A_wave_analysis(all_files, t_post = t_post_a)
	dataA = readABF(trace_A.Path |> unique, channels = channels)
	data_filter!(dataA, avg_swp = false, t_post = t_post_a) #Use the default data filter
end

# ╔═╡ 7ebbf522-125f-49da-b658-40e648e33d95
trace_A |> 
	@map({_.Channel, _.Photons, _.Response, _.Peak_Time, _.Percent_Recovery}) |> 
	@orderby(_.Channel) |> @thenby(_.Photons) |> 
DataFrame

# ╔═╡ 3a7ca4cb-8df2-4869-b3ef-f39774c00b58
md"""
### Plotting the A-wave:

xlims A:
(
$(@bind xlim1A NumberField(-1.0:0.01:10.000; default = -0.25)),
$(@bind xlim2A NumberField(-1.0:0.01:10.000; default = t_post_a))
)

ylims A:
	(
	$(@bind ylim1A NumberField(
		-10000.0:0.01:10000.000; default = minimum(dataA))),
	$(@bind ylim2A NumberField(
		-10000.0:0.01:10000.000; default = maximum(dataA)))
)

Photon Lims A:
(
$(@bind photon_lims1 NumberField(-1.0:0.1:10.0; default = -1.0)),
$(@bind photon_lims2 NumberField(-1.0:0.1:10.0; default = 5.0))
)

Recovery Lims A:
(
$(@bind rec_lims1A 
	NumberField(1.0:1.0:1000.0; default = 10.0))
)
$(@bind rec_lims2A 
	NumberField(0.0:1.0:2000.0; default = maximum(trace_A.Percent_Recovery)*1000))
"""

# ╔═╡ c11191b2-3c8f-4962-b5f5-505b4042a5a1
I_rng = 10 .^ range(photon_lims1, photon_lims2, length = 1000); #This is the range of points we fit over

# ╔═╡ 7e697def-11bf-4783-9e35-a8bb9df0e530
begin
	fig1_awave, ax1 = plt.subplots(length(channels), 3) #This 
	fig1_awave.subplots_adjust(
		left=0.0, right=1.0, 
		bottom=0.0, top=1.0, 
		wspace=0.3, hspace=0.2
	)
	if length(channels) > 1
		for (ch, name) in enumerate(dataA.chNames)
			#println(ch)
			ax1[ch, 1].set_ylabel("$(name) ($(dataA.chUnits[1]))")
			plot_experiment(ax1[ch, 1], dataA, channels = ch,
				xaxes_off = ch == 1
			)
			ax1[ch, 1].set_ylim((ylim1A, ylim2A))
			ax1[ch, 1].set_xlim((xlim1A, xlim2A))
			#Fit and plot the IR curves
			ch_trace_A = trace_A |> @filter(_.Channel == name) |> DataFrame
			ch_exp_A = exp_A |> @filter(_.Channel == name) |> DataFrame
			RMAX_fit_A = ch_exp_A.RMAX_fit
			K_fit_A = ch_exp_A.K_fit
			N_fit_A = ch_exp_A.N_fit
			R_fit_A = HILL_MODEL(I_rng, (RMAX_fit_A[1], K_fit_A[1], N_fit_A[1]))
			ax1[ch,2].scatter(ch_trace_A.Photons, -ch_trace_A.Response, 
				s = 5.0, c = :black
			)
			ax1[ch,2].plot(I_rng, -R_fit_A, 
				label = "R² = $(round(ch_exp_A.RSQ_fit[1],digits = 2))",
				linewidth = 2.0, c = :black
			)
			ax1[ch,2].set_xscale("log")
			ax1[ch, 2].set_xlim((10^photon_lims1, 10^photon_lims2))
			ax1[ch, 2].set_ylim((ylim1A, ylim2A))
			ax1[ch,2].legend(loc = "upper right")
			
			#Plot the time to peak and 60% recovery for all traces RDIM and higher
			rdim_idxs = findall(ch_trace_A.Response .> ch_exp_A.rdim)
			ax1[ch,3].scatter(
				ch_trace_A.Photons[rdim_idxs], ch_trace_A.Peak_Time[rdim_idxs], 
				s = 8.0, c = :black, label = "tPeak"
			)
			#We need to remove all zeros
			tDOM_idxs = findall(ch_trace_A.Percent_Recovery .> 0.0)
			ax1[ch,3].scatter(
				ch_trace_A.Photons[tDOM_idxs], ch_trace_A.Percent_Recovery[tDOM_idxs]*1000, 
				s = 8.0, c = :red, label = "tDom"
			)
			ax1[ch,3].set_yscale("log")
			ax1[ch,3].set_xlim((10^photon_lims1, 10^photon_lims2))
			ax1[ch,3].set_xscale("log")
			ax1[ch,3].set_ylim((rec_lims1A, rec_lims2A))
		end
		ax1[1,3].legend(loc = "lower left")
		ax1[1,2].spines["bottom"].set_visible(false)
        ax1[1,2].xaxis.set_visible(false)
		ax1[1,2].spines["left"].set_visible(false)
        ax1[1,2].yaxis.set_visible(false)
		ax1[1,2].set_xlabel("Time (s)")
		ax1[2,2].spines["left"].set_visible(false)
        ax1[2,2].yaxis.set_visible(false)
	end
	fig1_awave
end

# ╔═╡ 749e519e-32da-4e1b-b98e-d639d6489759
md"""
#### Report the Intensity-Response best fits
"""

# ╔═╡ 026ab1aa-2f80-4097-b6fb-1f93f45cd779
exp_A |> @map({_.Channel, _.RMAX_fit, _.K_fit, _.N_fit, _.RSQ_fit}) |> DataFrame

# ╔═╡ b43bc77b-81d1-4bc8-8a73-f41e8a8ee668
md"""
#### Report the Rmax, Rdim, time to peak, and percent_recovery
"""

# ╔═╡ eef22abd-af03-4509-a70d-9f3d7b7a79df
exp_A |> @map({_.Channel, _.rmax, _.rdim, _.time_to_peak, _.percent_recovery}) |> DataFrame

# ╔═╡ c2a9d00e-bca1-4c8b-a5fe-5602f449c95c
md"""
# Analysis of the b-wave
Analysis Window end = $(@bind t_post_b NumberField(-1.0:0.01:10.000; default = 2.0)),
"""

# ╔═╡ 2ce227e8-b518-479b-acbe-4b843118b517
begin
	trace_B, exp_B, cond_B = run_B_wave_analysis(all_files; t_post = t_post_b)
	dataAB = readABF(trace_B.Path |> unique, channels = channels)
	data_filter!(dataAB, avg_swp = false, t_post = t_post_b) #Use the default data filter
	data_ab_B = readABF(trace_B.Path |> unique, channels = channels)
	data_filter!(data_ab_B, avg_swp = false, t_post = t_post_b)
	data_a_B = readABF(trace_B.A_Path |> unique, channels = channels)
	data_filter!(data_a_B, avg_swp = false, t_post = t_post_b)
	dataB = data_ab_B - data_a_B
end;

# ╔═╡ 21dd7e8f-de62-45d6-bac6-d7f266e6487e
trace_B |> 
	@map({_.Channel, _.Photons, _.Response, _.Peak_Time, _.Percent_Recovery}) |> 
	@orderby(_.Channel) |> @thenby(_.Photons) |> 
DataFrame

# ╔═╡ 44645dde-129f-4bec-b246-07ff3941c8a3
md"""
### Plotting the B-wave:

xlims B:
(
$(@bind xlim1B NumberField(-1.0:0.01:10.000; default = -0.25)),
$(@bind xlim2B NumberField(-1.0:0.01:10.000; default = t_post_b))
)

ylims B:
	(
	$(@bind ylim1B NumberField(
		-10000.0:0.01:10000.000; default = minimum(dataAB))
	),
	$(@bind ylim2B NumberField(
		-10000.0:0.01:10000.000; default = maximum(dataB))
	)
)

Recovery Lims A:
(
$(@bind rec_lims1B 
	NumberField(1.0:1.0:1000.0; default = 10.0))
)
$(@bind rec_lims2B 
	NumberField(1.0:1.0:2000.0; default = maximum(trace_B.Percent_Recovery)*1000))
"""

# ╔═╡ 0a7f8699-e658-4592-ad38-715487e7b023
begin
	fig2_bwave, ax2 = plt.subplots(length(channels), 4) #This 
	fig2_bwave.subplots_adjust(
		left=0.0, right=1.0, 
		bottom=0.0, top=1.0, 
		wspace=0.3, hspace=0.2
	)
	if length(channels) > 1
		for (ch, name) in enumerate(dataAB.chNames)
			#println(ch)
			ax2[ch, 1].set_ylabel("$(name) ($(dataAB.chUnits[1]))")
			plot_experiment(ax2[ch, 1], data_ab_B, channels = ch,
				xaxes_off = ch == 1
			)
			plot_experiment(ax2[ch, 1], data_a_B, color = :red, channels = ch,
				xaxes_off = ch == 1
			)
			plot_experiment(ax2[ch, 2], dataB, channels = ch,
				xaxes_off = ch == 1
			)
			ax2[ch, 1].set_ylim((ylim1B, ylim2B))
			ax2[ch, 1].set_xlim((xlim1B, xlim2B))
			ax2[ch, 2].set_ylim((ylim1B, ylim2B))
			ax2[ch, 2].set_xlim((xlim1B, xlim2B))
			
			#filter the responses from trace_A for the channel
			ch_trace_B = trace_B |> @filter(_.Channel == name) |> DataFrame
			ch_exp_B = exp_B |> @filter(_.Channel == name) |> DataFrame
			RMAX_fit_B = ch_exp_B.RMAX_fit
			K_fit_B = ch_exp_B.K_fit
			N_fit_B = ch_exp_B.N_fit
			R_fit_B = HILL_MODEL(I_rng, (RMAX_fit_B[1], K_fit_B[1], N_fit_B[1]))
			ax2[ch,3].scatter(ch_trace_B.Photons, ch_trace_B.Response, 
				s = 4.0, c = :black
			)
			ax2[ch,3].plot(I_rng, R_fit_B, 
				label = "R² = $(round(ch_exp_B.RSQ_fit[1],digits = 2))",
				linewidth = 2.0, c = :black
			)
			ax2[ch,3].set_xscale("log")
			ax2[ch,3].set_ylim((ylim1B, ylim2B))
			ax2[ch,3].legend(loc = "upper right")
			ax2[ch,3].set_xlim((10^photon_lims1, 10^photon_lims2))
			
			#Plot the time to peak and 60% recovery for all traces RDIM and higher
			rdim_idxs = findall(ch_trace_B.Response .> ch_exp_B.rdim)
			ax2[ch,4].scatter(
				ch_trace_B.Photons[rdim_idxs], ch_trace_B.Peak_Time[rdim_idxs], 
				s = 8.0, c = :black, label = "tPeak"
			)
			#We need to remove all zeros
			tDOM_idxs = findall(ch_trace_B.Percent_Recovery .> 0.0)
			ax2[ch,4].scatter(
				ch_trace_B.Photons[tDOM_idxs], ch_trace_B.Percent_Recovery[tDOM_idxs]*1000, 
				s = 8.0, c = :red, label = "tDom"
			)
			ax2[ch,4].set_yscale("log")
			ax2[ch,4].set_xlim((10^photon_lims1, 10^photon_lims2))
			ax2[ch,4].set_xscale("log")	
			ax2[ch,4].set_ylim((rec_lims1B, rec_lims2B))
		end
		ax2[1,2].spines["left"].set_visible(false)
        ax2[1,2].yaxis.set_visible(false)
		ax2[1,2].spines["bottom"].set_visible(false)
        ax2[1,2].xaxis.set_visible(false)
		ax2[2,2].spines["left"].set_visible(false)
        ax2[2,2].yaxis.set_visible(false)

		ax2[1,3].spines["left"].set_visible(false)
        ax2[1,3].yaxis.set_visible(false)
		ax2[1,3].spines["bottom"].set_visible(false)
        ax2[1,3].xaxis.set_visible(false)
		ax2[2,3].spines["left"].set_visible(false)
        ax2[2,3].yaxis.set_visible(false)
	end
	fig2_bwave
end

# ╔═╡ 02138c18-d46e-4a5d-aaf7-3281ccf755b0
md"""
#### Report the Intensity-Response best fits
"""

# ╔═╡ 8fa36875-e7e8-442e-bb0c-2f24af295592
exp_B |> @map({_.Channel, _.RMAX_fit, _.K_fit, _.N_fit, _.RSQ_fit}) |> DataFrame

# ╔═╡ f571351f-992d-4c15-90b6-fd63be9c6793
md"""
#### Report the Rmax, Rdim, time to peak, and percent_recovery
"""

# ╔═╡ 25e51d82-66ef-4ec2-8362-dbb82a4dac55
exp_B |> @map({_.Channel, _.rmax, _.rdim, _.time_to_peak, _.percent_recovery}) |> DataFrame

# ╔═╡ c0993bf9-c0fe-4aee-9ebb-6e63ce4a2246
md"""
# Analysis of the Glial-component
Analysis Window end = $(@bind t_post_g NumberField(-1.0:0.01:10.000; default = 2.0)),
"""

# ╔═╡ 2d62d63c-081f-4ddf-8a8d-d1f4c92827f5
begin
	trace_G, exp_G, cond_G = run_G_wave_analysis(all_files)
	
	dataABG = readABF(trace_G.Path |> unique, channels = channels)
	data_filter!(dataABG, avg_swp = false, t_post = 5.0) #Use the default data filter
	data_abg_G = readABF(trace_G.Path |> unique, channels = channels)
	data_filter!(data_abg_G, avg_swp = false, t_post = 5.0)
	data_ab_G = readABF(trace_G.AB_Path |> unique, channels = channels)
	data_filter!(data_ab_G, avg_swp = false, t_post = 5.0)
	dataG = data_abg_G - data_ab_G
end;

# ╔═╡ 87c01068-8042-490b-a2b3-98a6574e87e8
trace_G |> 
	@map({_.Channel, _.Photons, _.Response, _.Peak_Time, _.Percent_Recovery}) |> 
	@orderby(_.Channel) |> @thenby(_.Photons) |> 
DataFrame

# ╔═╡ 1bce802f-0930-4eb7-bb98-e4df80abff26
md"""
### Plotting the G-wave:

xlims G:
(
$(@bind xlim1G NumberField(-1.0:0.01:10.000; default = -0.25)),
$(@bind xlim2G NumberField(-1.0:0.01:10.000; default = t_post_g))
)

ylims G:
	(
	$(@bind ylim1G NumberField(
		-10000.0:0.01:10000.000; default = minimum(dataG))
	),
	$(@bind ylim2G NumberField(
		-10000.0:0.01:10000.000; default = maximum(dataB))
	)
)

Recovery Lims G:
(
$(@bind rec_lims1G 
	NumberField(1.0:1.0:1000.0; default = 10.0))
)
$(@bind rec_lims2G 
	NumberField(1.0:1.0:2000.0; default = maximum(trace_G.Percent_Recovery)*1000))
"""

# ╔═╡ 3c2b3129-6d26-43d0-96e5-386bf47bddab
begin
	fig3_gwave, ax3 = plt.subplots(length(channels), 4) #This 
	fig3_gwave.subplots_adjust(
		left=0.0, right=1.0, 
		bottom=0.0, top=1.0, 
		wspace=0.1, hspace=0.2
	)
	if length(channels) > 1
		for (ch, name) in enumerate(dataABG.chNames)
			#println(ch)
			ax3[ch, 1].set_ylabel("$(name) ($(dataABG.chUnits[1]))")
			plot_experiment(ax3[ch, 1], data_abg_G, channels = ch,
				xaxes_off = ch == 1
			)
			plot_experiment(ax3[ch, 1], data_ab_G, color = :red, channels = ch,
				xaxes_off = ch == 1
			)
			ax3[ch, 1].set_ylim((ylim1G, ylim2G))
			ax3[ch, 1].set_xlim((xlim1G, xlim2G))
			
			plot_experiment(ax3[ch, 2], dataG, channels = ch,
				xaxes_off = ch == 1
			)
			ax3[ch, 2].set_ylim((ylim1G, ylim2G))
			ax3[ch, 2].set_xlim((xlim1G, xlim2G))
			
			#filter the responses from trace_A for the channel
			ch_trace_G = trace_G |> @filter(_.Channel == name) |> DataFrame
			ch_exp_G = exp_G |> @filter(_.Channel == name) |> DataFrame
			RMAX_fit_G = ch_exp_G.RMAX_fit
			K_fit_G = ch_exp_G.K_fit
			N_fit_G = ch_exp_G.N_fit
			R_fit_G = HILL_MODEL(I_rng, (RMAX_fit_G[1], K_fit_G[1], N_fit_G[1]))
			ax3[ch,3].scatter(ch_trace_G.Photons, -ch_trace_G.Response, 
				s = 4.0, c = :black
			)
			ax3[ch,3].plot(I_rng, -R_fit_G, 
				label = "R² = $(round(ch_exp_G.RSQ_fit[1],digits = 2))",
				linewidth = 2.0, c = :black
			)
			ax3[ch,3].set_xscale("log")
			ax3[ch, 3].set_xlim((10^photon_lims1, 10^photon_lims2))
			ax3[ch,3].legend(loc = "upper right")
			ax3[ch, 3].set_ylim((ylim1G, ylim2G))
			
			#Plot the time to peak and 60% recovery for all traces RDIM and higher
			rdim_idxs = findall(ch_trace_G.Response .> ch_exp_G.rdim)
			ax3[ch,4].scatter(
				ch_trace_G.Photons[rdim_idxs], ch_trace_G.Peak_Time[rdim_idxs], 
				s = 8.0, c = :black, label = "tPeak"
			)
			#We need to remove all zeros
			tDOM_idxs = findall(ch_trace_G.Percent_Recovery .> 0.0)
			ax3[ch,4].scatter(
				ch_trace_G.Photons[tDOM_idxs], ch_trace_G.Percent_Recovery[tDOM_idxs]*1000, 
				s = 8.0, c = :red, label = "tDom"
			)
			ax3[ch,4].set_yscale("log")
			ax3[ch,4].set_xlim((10^photon_lims1, 10^photon_lims2))
			ax3[ch,4].set_xscale("log")	
			ax3[ch,4].set_ylim((rec_lims1G, rec_lims2G))
		end
		ax3[1,2].spines["left"].set_visible(false)
        ax3[1,2].yaxis.set_visible(false)
		ax3[1,2].spines["bottom"].set_visible(false)
        ax3[1,2].xaxis.set_visible(false)
		ax3[2,2].spines["left"].set_visible(false)
        ax3[2,2].yaxis.set_visible(false)

		ax3[1,3].spines["left"].set_visible(false)
        ax3[1,3].yaxis.set_visible(false)
		ax3[1,3].spines["bottom"].set_visible(false)
        ax3[1,3].xaxis.set_visible(false)
		ax3[2,3].spines["left"].set_visible(false)
        ax3[2,3].yaxis.set_visible(false)
	end
	fig3_gwave
	
end

# ╔═╡ cfe8f792-84d8-4129-9c0c-061e6b79348e
md"""
#### Report the Intensity-Response best fits
"""

# ╔═╡ 1ec9990d-8c4a-4ecd-bcf2-b4407a61e9dd
exp_G |> @map({_.Channel, _.RMAX_fit, _.K_fit, _.N_fit, _.RSQ_fit}) |> DataFrame

# ╔═╡ a683e1be-b22c-4278-acfa-a2acfd6754cb
md"""
#### Report the Rmax, Rdim, time to peak, and percent_recovery
"""

# ╔═╡ b8246cf4-22f9-4f88-89d0-3fdbdaba0c1b
exp_G |> @map({_.Channel, _.rmax, _.rdim, _.time_to_peak, _.percent_recovery}) |> DataFrame

# ╔═╡ fa29d56f-93e2-4260-be26-1f985c4804f3
begin
	fig1_awave
	fig2_bwave
	fig3_gwave
	plt.close("all")
end

# ╔═╡ 5ab3b49d-3427-4971-8849-ecb614d4bf66
md"""
### Fitting the Synaptic transfer functions
"""

# ╔═╡ Cell order:
# ╠═8ec03fe0-866b-11ed-23a6-23f888e1717a
# ╟─9d51d150-3339-435c-ad76-59b81bdcf305
# ╠═1fb1daf9-59fd-42aa-b30b-bc5413a6cc67
# ╠═bc4d42b1-6848-4c61-99ce-1e437c08497b
# ╟─3527f4ce-d1fb-46ac-a6a7-e9d3d05a03de
# ╟─022ec5ec-5b0f-4876-8ff4-18ed6f493f2b
# ╟─01fab469-e7c3-4573-9e75-d2d13c7ff3de
# ╟─1cbc25be-bf12-43ac-a545-89aeab07ff93
# ╟─88ef46c3-7575-4a18-a8b5-15a0c2152c6a
# ╟─734c8f3e-7092-48fd-9784-391bde74fc08
# ╟─7ebbf522-125f-49da-b658-40e648e33d95
# ╟─3a7ca4cb-8df2-4869-b3ef-f39774c00b58
# ╠═c11191b2-3c8f-4962-b5f5-505b4042a5a1
# ╟─7e697def-11bf-4783-9e35-a8bb9df0e530
# ╟─749e519e-32da-4e1b-b98e-d639d6489759
# ╟─026ab1aa-2f80-4097-b6fb-1f93f45cd779
# ╟─b43bc77b-81d1-4bc8-8a73-f41e8a8ee668
# ╟─eef22abd-af03-4509-a70d-9f3d7b7a79df
# ╟─c2a9d00e-bca1-4c8b-a5fe-5602f449c95c
# ╟─2ce227e8-b518-479b-acbe-4b843118b517
# ╟─21dd7e8f-de62-45d6-bac6-d7f266e6487e
# ╟─44645dde-129f-4bec-b246-07ff3941c8a3
# ╟─0a7f8699-e658-4592-ad38-715487e7b023
# ╟─02138c18-d46e-4a5d-aaf7-3281ccf755b0
# ╟─8fa36875-e7e8-442e-bb0c-2f24af295592
# ╟─f571351f-992d-4c15-90b6-fd63be9c6793
# ╟─25e51d82-66ef-4ec2-8362-dbb82a4dac55
# ╟─c0993bf9-c0fe-4aee-9ebb-6e63ce4a2246
# ╟─2d62d63c-081f-4ddf-8a8d-d1f4c92827f5
# ╟─87c01068-8042-490b-a2b3-98a6574e87e8
# ╟─1bce802f-0930-4eb7-bb98-e4df80abff26
# ╟─3c2b3129-6d26-43d0-96e5-386bf47bddab
# ╟─cfe8f792-84d8-4129-9c0c-061e6b79348e
# ╟─1ec9990d-8c4a-4ecd-bcf2-b4407a61e9dd
# ╟─a683e1be-b22c-4278-acfa-a2acfd6754cb
# ╟─b8246cf4-22f9-4f88-89d0-3fdbdaba0c1b
# ╠═fa29d56f-93e2-4260-be26-1f985c4804f3
# ╟─5ab3b49d-3427-4971-8849-ecb614d4bf66
