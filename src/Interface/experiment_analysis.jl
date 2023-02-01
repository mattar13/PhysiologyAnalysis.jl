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
	using DataFrames, Query
	import ePhys: baseline_adjust!, truncate_data! , average_sweeps!
	import ePhys: filter_data!, filter_data
	import ePhys: EI_filter!
	import ePhys: dwt_filter!, cwt_filter!
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

"

# ╔═╡ bc4d42b1-6848-4c61-99ce-1e437c08497b
begin
	#enter in the file path of the file you would like to analyze here
	exp_root = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\ERG\Retinoschisis\2022_07_11_Adult\Mouse3_Adult_RS1KO"
	experiment_paths = exp_root |> parseABF
end

# ╔═╡ cc7f746a-5714-43f7-bd32-07ad350a2aef
pcs = "Rods"

# ╔═╡ bc2576c6-b878-4c39-855b-49311d992f27
begin #Split the paths into A, AB, ABG waves
	all_files = createDatasheet(experiment_paths)
	#filter out rods or cones if you have that
	all_files = all_files |> @filter(_.Photoreceptor == pcs) |> DataFrame
	@time ATrace, AExperiment, qConditions = ePhys.run_A_wave_analysis(all_files, verbose = false)
	#@time BTrace, BExperiment, qConditions = ePhys.run_B_wave_analysis(all_files, verbose = false)
	#@time GTrace, GExperiment, qConditions = ePhys.run_G_wave_analysis(all_files, verbose = false)
end

# ╔═╡ 1495d044-e661-4aee-bfd4-6e029f2cb51d
channels = ["Vm_prime", "Vm_prime4"]

# ╔═╡ 734c8f3e-7092-48fd-9784-391bde74fc08
begin
	dataABG = readABF(GTrace.Path |> unique, channels = channels)
	data_filter!(dataABG, avg_swp = false) #Use the default data filter
	data_abg_G = readABF(GTrace.Path |> unique, channels = channels)
	data_filter!(data_abg_G, avg_swp = false)
	data_ab_G = readABF(GTrace.AB_Path |> unique, channels = channels)
	data_filter!(data_ab_G, avg_swp = false)
	
	dataG = data_abg_G - data_ab_G
	
	dataAB = readABF(BTrace.Path |> unique, channels = channels)
	data_filter!(dataAB, avg_swp = false) #Use the default data filter
	data_ab_B = readABF(BTrace.Path |> unique, channels = channels)
	data_filter!(data_ab_B, avg_swp = false)
	data_a_B = readABF(BTrace.A_Path |> unique, channels = channels)
	data_filter!(data_a_B, avg_swp = false)
	dataB = data_ab_B - data_a_B
	
	dataA = readABF(ATrace.Path |> unique, channels = channels)
	data_filter!(dataA, avg_swp = false) #Use the default data filter
end

# ╔═╡ 3a7ca4cb-8df2-4869-b3ef-f39774c00b58
md"""
### Plotting:

xlims:
(
$(@bind xlim1 NumberField(-1.0:0.01:10.000; default = -0.25)),
$(@bind xlim2 NumberField(-1.0:0.01:10.000; default = 2.0))
)

ylims ABG:
	(
	$(@bind ylim1ABG NumberField(
		-10000.0:0.01:10000.000; default = minimum(dataABG))
	),
	$(@bind ylim2ABG NumberField(
		-10000.0:0.01:10000.000; default = maximum(dataB))
	)
)

"""

# ╔═╡ c11191b2-3c8f-4962-b5f5-505b4042a5a1
size(dataABG)

# ╔═╡ 248fa29c-f334-4f09-89ad-77636030bf9e
begin
	fig1, ax1 = plt.subplots(5, length(channels)) #This 
	fig1.subplots_adjust(
		left=0.0,right=1.0, bottom=0.0, top=1.0, 
		wspace=0.1, hspace=0.2
	)
	# Plot the experiments
	ax1[1, 1].set_ylabel("Voltage ($(dataABG.chUnits[1]))")
	if length(channels) > 1
		for ch in 1:size(dataABG, 3)
			ax1[1, ch].set_title(
				"$(dataABG.chNames[ch]) ($(dataABG.chUnits[ch]))"
			)
			
			ax1[1, ch].set_xlim(xlim1, xlim2)
			ax1[1, ch].set_ylim(ylim1ABG, ylim2ABG)
			if ch == length(channels)
				plot_experiment(ax1[1, ch], dataABG, axes_off = true, channels = ch, c = :black)
			else
				plot_experiment(ax1[1, ch], dataABG, xaxes_off = true, channels = ch, c = :black)
			end
			
			ax1[2, ch].set_xlim(xlim1, xlim2)
			ax1[2, ch].set_ylim(ylim1ABG, ylim2ABG)
			if ch == length(channels)
				plot_experiment(ax1[2, ch], dataG, axes_off = true, channels = ch, c = :black)
			else
				plot_experiment(ax1[2, ch], dataG, xaxes_off = true, channels = ch, c = :black)
			end
						
			ax1[3, ch].set_xlim(xlim1, xlim2)
			ax1[3, ch].set_ylim(ylim1ABG, ylim2ABG)
			if ch == length(channels)
				plot_experiment(ax1[3, ch], dataAB, axes_off = true, channels = ch, c = :black)
			else
				plot_experiment(ax1[3, ch], dataAB, xaxes_off = true, channels = ch, c = :black)
			end

			ax1[4, ch].set_xlim(xlim1, xlim2)
			ax1[4, ch].set_ylim(ylim1ABG, ylim2ABG)
			if ch == length(channels)
				plot_experiment(ax1[4, ch], dataB, axes_off = true, channels = ch, c = :black)
			else
				plot_experiment(ax1[4, ch], dataB, xaxes_off = true, channels = ch, c = :black)
			end
			
			ax1[5, ch].set_xlim(xlim1, xlim2)
			ax1[5, ch].set_ylim(ylim1ABG, ylim2ABG)
			if ch == length(channels)
				plot_experiment(ax1[5, ch], dataA, yaxes_off = true, channels = ch, c = :black)
			else
				plot_experiment(ax1[5, ch], dataA, channels = ch, c = :black)
			end
			
		end
		ax1[5, 1].set_xlabel("Time (s)")
		ax1[5, length(channels)].set_xlabel("Time (s)")
	else
		ax.set_xlim(xlim1, xlim2)
		ax.set_ylim(ylim1, ylim2)
		plot_experiment(ax, dataABG, channels = 1, c = :black)
		ax.set_ylabel("$(data.chNames[1]) ($(data.chUnits[1]))")
		ax.set_xlabel("Time (s)")
	end	
	fig1
end

# ╔═╡ 0be359d0-3802-4a60-9040-6e0e4597e8d0
md"""
#### 4) Fitting IR curves
R parameter:
RMIN $(@bind RMIN NumberField(0.1:0.01:10000; default = 0.1)) ->
INITIAL R $(@bind R0 NumberField(0.0:0.01:10000; default = 100.0)) -> 
RMAX $(@bind RMAX NumberField(0.0:0.01:10000; default = maximum(dataB)))


k parameter:
kMIN $(@bind kMIN NumberField(0.0:0.01:10000; default = 0.1)) ->
INITIAL k $(@bind k0 NumberField(1.0:0.01:10000; default = 5000)) -> 
kMAX $(@bind kMAX NumberField(0.0:0.01:10000; default = 10e6))

n parameter:
nMIN $(@bind nMIN NumberField(0.1:0.01:10000; default = 0.1)) ->
INITIAL n $(@bind n0 NumberField(0.1:0.01:10000; default = 2.0)) -> 
nMAX $(@bind nMAX NumberField(0.1:0.01:10000; default = 10.0))
"""

# ╔═╡ 5ab3b49d-3427-4971-8849-ecb614d4bf66
md"""
### Amplitude Statistics
"""

# ╔═╡ ade58738-bf2d-49c8-8419-e8d27f27497e
fitI = 10 .^ LinRange(-2, 5, 1000)

# ╔═╡ 06faab7e-e54a-4ce0-8364-14d11aaef85f
begin #Plot the amplitudes from the A-wave
	#Seperate the channels
	ylim1A = (minimum(dataA), maximum(dataA))
	AByChannel = DataFrame[]
	AirChannel = []
	for ch in unique(ATrace.Channel)
		ATrace_CH = ATrace |> @filter(_.Channel == ch) |> DataFrame
		push!(AByChannel, ATrace_CH)
	end	
	fig2A, ax2A = plt.subplots(length(AByChannel), 2)
	fig2A.subplots_adjust(wspace=0.2)
	for ch in 1:size(dataA,3)
		AChannel = AByChannel[ch] #This is the dataframe that contains the vals
		RESPa = -AChannel.Minima
		#MINa = AChannel.Minima
		PHOTa = AChannel.Photons
		
		ax2A[ch, 1].set_ylabel(
				"$(dataABG.chNames[ch]) ($(dataABG.chUnits[ch]))"
			)

		ax2A[ch, 1].set_xlim(xlim1, xlim2)
		ax2A[ch, 1].set_ylim(ylim1A)
		if ch == length(channels)
			plot_experiment(ax2A[ch, 1], dataA, channels = ch, c = :black)
			ax2A[ch, 1].set_xlabel("Time (s)")
		else
			plot_experiment(ax2A[ch, 1], dataA, xaxes_off = true, channels = ch, c = :black)
		end
		ax2A[ch, 1].hlines(-RESPa, xlim1, xlim2, color = :Red)
		
		#Plot the Intensity response
		ax2A[ch, 2].set_ylim(ylim1A)
		ax2A[ch, 2].scatter(PHOTa, -RESPa, color = :Red)
		
		#Fit the data
		fitA = ePhys.IRfit(
			PHOTa, RESPa,
			rmin = RMIN, r = R0, rmax = RMAX,
			kmin = kMIN, k = k0, kmax = kMAX, 
			nmin = nMIN, n = n0, nmax = nMAX		
		)
		push!(AirChannel, fitA)
		fitRMAXa, fitKa, fitNa = fitA.param
		fitRa = -fitRMAXa .* ePhys.IR.(fitI, fitKa, fitNa)
		ax2A[ch, 2].plot(fitI, fitRa, c = :Black, lw = 2.0)
		
		ax2A[ch, 2].vlines([fitKa], ylim1A[1], ylim1A[2], 
			label = "K = $(round(fitKa[1], digits = 1))",
			color = :Red
		)
		
		ax2A[ch, 2].set_xscale("log")
	end
	fig2A
end

# ╔═╡ 24550894-7f61-4f7f-a60a-64e4517abf59
begin
	ylim1B = (minimum(dataB), maximum(dataB)+200)
	BByChannel = DataFrame[]
	BirChannel = []
	for ch in unique(BTrace.Channel)
		BTrace_CH = BTrace |> @filter(_.Channel == ch) |> DataFrame
		push!(BByChannel, BTrace_CH)
	end
	fig2B, ax2B = plt.subplots(length(AByChannel), 2)
	fig2B.subplots_adjust(wspace=0.2)
	for ch in 1:size(dataB,3)
		BChannel = BByChannel[ch] #This is the dataframe that contains the vals
		RESPb = BChannel.Response
		MINb = BChannel.Minima
		PHOTb = BChannel.Photons
		
		ax2B[ch, 1].set_ylabel(
				"$(dataB.chNames[ch]) ($(dataB.chUnits[ch]))"
			)

		ax2B[ch, 1].set_xlim(xlim1, xlim2)
		ax2B[ch, 1].set_ylim(ylim1B)
		if ch == length(channels)
			plot_experiment(ax2B[ch, 1], dataB, channels = ch, c = :black)
			ax2B[ch, 1].set_xlabel("Time (s)")
		else
			plot_experiment(ax2B[ch, 1], dataB, xaxes_off = true, channels = ch, c = :black)
		end
		ax2B[ch, 1].hlines(RESPb, xlim1, xlim2, color = :Red)
		
		#Plot the Intensity response
		#ax2B[ch, 2].set_xlim(1e1, 1e4)
		ax2B[ch, 2].set_ylim(ylim1B)
		ax2B[ch, 2].scatter(PHOTb, RESPb, color = :Red)
		
		#Fit the data
		fitB = ePhys.IRfit(
			PHOTb, RESPb,
			rmin = RMIN, r = R0, rmax = RMAX,
			kmin = kMIN, k = k0, kmax = kMAX, 
			nmin = nMIN, n = n0, nmax = nMAX		
		)
		push!(BirChannel, fitB)
		fitRMAXb, fitKb, fitNb = fitB.param
		fitRb = fitRMAXb .* ePhys.IR.(fitI, fitKb, fitNb)
		ax2B[ch, 2].plot(fitI, fitRb, c = :Black, lw = 2.0)
		
		ax2B[ch, 2].vlines([fitKb], ylim1B[1], ylim1B[2],  
			label = "K = $(round(fitKb[1], digits = 1))",
			color = :Red
		)
		
		ax2B[ch, 2].set_xscale("log")
	end
	fig2B
	
end

# ╔═╡ 0af994e1-af6f-433d-9a49-2dedce76a95d
begin
	ylim1G = (minimum(dataG)-200, maximum(dataG))
	GByChannel = DataFrame[]
	GirChannel = []
	for ch in unique(GTrace.Channel)
		GTrace_CH = GTrace |> @filter(_.Channel == ch) |> DataFrame
		push!(GByChannel, GTrace_CH)
	end
	fig2G, ax2G = plt.subplots(length(GByChannel), 2)
	fig2G.subplots_adjust(wspace=0.2)
	for ch in 1:size(dataG,3)
		GChannel = GByChannel[ch] #This is the dataframe that contains the vals
		RESPg = GChannel.Response
		MINg = GChannel.Minima
		PHOTg = GChannel.Photons
		
		ax2G[ch, 1].set_ylabel(
				"$(dataG.chNames[ch]) ($(dataG.chUnits[ch]))"
			)

		ax2G[ch, 1].set_xlim(xlim1, xlim2)
		ax2G[ch, 1].set_ylim(ylim1G)
		if ch == length(channels)
			plot_experiment(ax2G[ch, 1], dataG, channels = ch, c = :black)
			ax2G[ch, 1].set_xlabel("Time (s)")
		else
			plot_experiment(ax2G[ch, 1], dataG, xaxes_off = true, channels = ch, c = :black)
		end
		ax2G[ch, 1].hlines(-RESPg, xlim1, xlim2, color = :Red)
		
		#Plot the Intensity response
		ax2G[ch, 2].set_ylim(ylim1G)
		ax2G[ch, 2].scatter(PHOTg, -RESPg, color = :Red)
		
		#Fit the data
		fitG = ePhys.IRfit(
			PHOTg, RESPg,
			rmin = RMIN, r = R0, rmax = RMAX,
			kmin = kMIN, k = k0, kmax = kMAX, 
			nmin = nMIN, n = n0, nmax = nMAX		
		)
		push!(GirChannel, fitG)
		fitRMAXg, fitKg, fitNg = fitG.param
		fitRg = -fitRMAXg .* ePhys.IR.(fitI, fitKg, fitNg)
		ax2G[ch, 2].plot(fitI, fitRg, c = :Black, lw = 2.0)
		
		ax2G[ch, 2].vlines([fitKg], ylim1G[1], ylim1G[2],
			label = "K = $(round(fitKg[1], digits = 1))",
			color = :Red
		)
		
		ax2G[ch, 2].set_xscale("log")
	end
	fig2G
		
end

# ╔═╡ 2501825a-093f-42cd-9c0e-0a771f1700d0
begin #overlay all the intensity response relationships
	fig3AB, ax3AB = plt.subplots(4, length(channels))
	fig3AB.subplots_adjust(hspace=1.0)
	for ch in 1:length(channels)
		plot_experiment(ax3AB[1, ch], dataA, channels = ch, c = :red)
		plot_experiment(ax3AB[1, ch], dataB, channels = ch, c = :black)
		ax3AB[1, ch].set_xlim(xlim1, xlim2)
		ax3AB[1, ch].set_ylim(ylim1A[1], ylim1B[2])
		ax3AB[1, ch].set_xlabel("Time (s)")

		fitRMAXaB, fitKaB, fitNaB = BirChannel[ch].param
		fitRMAXAb, fitKAb, fitNAb = AirChannel[ch].param
		
		#fitRaB = fitRMAXaB .* ePhys.IR.(fitI, fitKaB, fitNaB)
		#fitRAb = fitRMAXAb .* ePhys.IR.(fitI, fitKAb, fitNAb)
		#Normalize the relationships
		fitRaB = ePhys.IR.(fitI, fitKaB, fitNaB)
		fitRAb = ePhys.IR.(fitI, fitKAb, fitNAb)
		ax3AB[2, ch].plot(fitI, fitRaB, c = :Black, lw = 2.0) #plot the B IR 
		ax3AB[2, ch].plot(fitI, fitRAb, c = :red, lw = 2.0) #plot the A ir curve
		ax3AB[2, ch].set_xscale("log")

		AbChannel = AByChannel[ch] #This is the dataframe that contains the vals
		aBChannel = BByChannel[ch]
		PHOTAB = AbChannel.Photons
		RESPAb = AbChannel.Response
		RESPaB = aBChannel.Response

		ax3AB[3, ch].scatter(RESPAb, RESPaB, c = :Black, lw = 2.0) #plot the B IR 
		#ax3AB[3, ch].set_xscale("log")
		ax3AB[3, ch].set_yscale("log")
		#ax3AB[3, ch].set_xlim(abs.(ylim1A))
		ax3AB[3, ch].set_ylim(ylim1B)
		
		ax3AB[4, ch].scatter(PHOTAB, RESPaB./RESPAb, c = :Black, lw = 2.0) #plot the B IR 
		ax3AB[4, ch].set_xscale("log")
	end
	fig3AB
end

# ╔═╡ a67381ec-5725-4f31-973e-e0028df40e78
begin #overlay all the intensity response relationships
	fig3AG, ax3AG = plt.subplots(3, length(channels))
	fig3AG.subplots_adjust(hspace=0.4)
	for ch in 1:length(channels)
		plot_experiment(ax3AG[1, ch], dataA, channels = ch, c = :red)
		plot_experiment(ax3AG[1, ch], dataG, channels = ch, c = :black)
		ax3AG[1, ch].set_xlim(xlim1, xlim2)
		ax3AG[1, ch].set_ylim(ylim1G[1], ylim1G[2])

		fitRMAXaG, fitKaG, fitNaG = GirChannel[ch].param
		fitRMAXAg, fitKAg, fitNAg = AirChannel[ch].param
		
		fitRaG = fitRMAXaG .* ePhys.IR.(fitI, fitKaG, fitNaG)
		fitRAg = fitRMAXAg .* ePhys.IR.(fitI, fitKAg, fitNAg)
		
		ax3AG[2, ch].plot(fitI, fitRaG, c = :Black, lw = 2.0) #plot the G IR 
		ax3AG[2, ch].plot(fitI, fitRAg, c = :red, lw = 2.0) #plot the A ir curve
		ax3AG[2, ch].set_xscale("log")
		
		#calculate the ratio between each 
		ax3AG[1, ch].set_xlabel("Time (s)")
	end
	fig3AG
end

# ╔═╡ Cell order:
# ╠═8ec03fe0-866b-11ed-23a6-23f888e1717a
# ╟─9d51d150-3339-435c-ad76-59b81bdcf305
# ╠═bc4d42b1-6848-4c61-99ce-1e437c08497b
# ╠═cc7f746a-5714-43f7-bd32-07ad350a2aef
# ╟─bc2576c6-b878-4c39-855b-49311d992f27
# ╠═1495d044-e661-4aee-bfd4-6e029f2cb51d
# ╟─734c8f3e-7092-48fd-9784-391bde74fc08
# ╟─3a7ca4cb-8df2-4869-b3ef-f39774c00b58
# ╠═c11191b2-3c8f-4962-b5f5-505b4042a5a1
# ╟─248fa29c-f334-4f09-89ad-77636030bf9e
# ╟─0be359d0-3802-4a60-9040-6e0e4597e8d0
# ╟─5ab3b49d-3427-4971-8849-ecb614d4bf66
# ╟─ade58738-bf2d-49c8-8419-e8d27f27497e
# ╟─06faab7e-e54a-4ce0-8364-14d11aaef85f
# ╟─24550894-7f61-4f7f-a60a-64e4517abf59
# ╟─0af994e1-af6f-433d-9a49-2dedce76a95d
# ╟─2501825a-093f-42cd-9c0e-0a771f1700d0
# ╟─a67381ec-5725-4f31-973e-e0028df40e78
