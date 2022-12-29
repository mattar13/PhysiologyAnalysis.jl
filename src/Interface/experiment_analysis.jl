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
	exp_root = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\ERG\Retinoschisis\2022_04_12_RCAdult\Mouse1_Adult_R141C"
	experiment_paths = exp_root |> parseABF
end

# ╔═╡ bc2576c6-b878-4c39-855b-49311d992f27
begin #Split the paths into A, AB, ABG waves
	all_files = createDatasheet(experiment_paths)
	@time ATrace, AExperiment, qConditions = ePhys.run_A_wave_analysis(all_files, verbose = false)
	@time BTrace, BExperiment, qConditions = ePhys.run_B_wave_analysis(all_files, verbose = false)
	@time GTrace, GExperiment, qConditions = ePhys.run_G_wave_analysis(all_files, verbose = false)
end

# ╔═╡ 1787482b-cac3-4875-8c42-663cba7447d6
qTrace.Path |> unique

# ╔═╡ 1495d044-e661-4aee-bfd4-6e029f2cb51d
channels = ["Vm_prime", "Vm_prime4"]

# ╔═╡ 734c8f3e-7092-48fd-9784-391bde74fc08
begin
	dataABG = readABF(ABG_df.Path, channels = channels)
	data_filter!(dataABG, avg_swp = false) #Use the default data filter

	dataAB = readABF(AB_df.Path, channels = channels)
	data_filter!(dataAB, avg_swp = false) #Use the default data filter

	dataA = readABF(A_df.Path, channels = channels)
	data_filter!(dataA, avg_swp = false) #Use the default data filter
end

# ╔═╡ 6e25216b-56c6-4b02-a642-33d319dc2f2b
size(dataABG)

# ╔═╡ a7eaf35a-2075-471b-9950-6330c932b849
size(dataAB)

# ╔═╡ 77e9089b-4e65-4346-a0cd-a89f6a12da2b
size(dataA)

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
		-10000.0:0.01:10000.000; default = maximum(dataABG))
	)
)

"""

# ╔═╡ 248fa29c-f334-4f09-89ad-77636030bf9e
begin
	fig1, ax1 = plt.subplots(length(channels), 3) #This 
	
	# Plot the experiments
	if length(channels) > 1
		for ch in 1:length(channels)
			ax1[ch, 1].set_xlim(xlim1, xlim2)
			ax1[ch, 1].set_ylim(ylim1ABG, ylim2ABG)
			plot_experiment(ax1[ch, 1], dataABG, channels = ch, c = :black)
			ax1[ch, 1].set_ylabel(
				"$(dataABG.chNames[ch]) ($(dataABG.chUnits[ch]))"
			)
			
			ax1[ch, 2].set_xlim(xlim1, xlim2)
			ax1[ch, 2].set_ylim(ylim1ABG, ylim2ABG)
			plot_experiment(ax1[ch, 2], dataAB, channels = ch, c = :black)

			ax1[ch, 3].set_xlim(xlim1, xlim2)
			ax1[ch, 3].set_ylim(ylim1ABG, ylim2ABG)
			plot_experiment(ax1[ch, 3], dataA, channels = ch, c = :black)
			
		end
		ax1[length(channels), 1].set_xlabel("Time (s)")
		ax1[length(channels), 2].set_xlabel("Time (s)")
		ax1[length(channels), 3].set_xlabel("Time (s)")
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
RMAX $(@bind RMAX NumberField(0.0:0.01:10000; default = -minimum(dataABG)))


k parameter:
kMIN $(@bind kMIN NumberField(0.0:0.01:10000; default = 0.1)) ->
INITIAL k $(@bind k0 NumberField(1.0:0.01:10000; default = 5000)) -> 
kMAX $(@bind kMAX NumberField(0.0:0.01:10000; default = 10e6))

n parameter:
nMIN $(@bind nMIN NumberField(0.1:0.01:10000; default = 0.1)) ->
INITIAL n $(@bind n0 NumberField(0.1:0.01:10000; default = 2.0)) -> 
nMAX $(@bind nMAX NumberField(0.1:0.01:10000; default = 10.0))
"""

# ╔═╡ d1d7ce7e-22b3-42f9-b164-b6bebcfceccc


# ╔═╡ Cell order:
# ╠═8ec03fe0-866b-11ed-23a6-23f888e1717a
# ╟─9d51d150-3339-435c-ad76-59b81bdcf305
# ╠═bc4d42b1-6848-4c61-99ce-1e437c08497b
# ╠═bc2576c6-b878-4c39-855b-49311d992f27
# ╠═1787482b-cac3-4875-8c42-663cba7447d6
# ╠═1495d044-e661-4aee-bfd4-6e029f2cb51d
# ╠═734c8f3e-7092-48fd-9784-391bde74fc08
# ╠═6e25216b-56c6-4b02-a642-33d319dc2f2b
# ╠═a7eaf35a-2075-471b-9950-6330c932b849
# ╠═77e9089b-4e65-4346-a0cd-a89f6a12da2b
# ╟─3a7ca4cb-8df2-4869-b3ef-f39774c00b58
# ╟─248fa29c-f334-4f09-89ad-77636030bf9e
# ╟─0be359d0-3802-4a60-9040-6e0e4597e8d0
# ╠═d1d7ce7e-22b3-42f9-b164-b6bebcfceccc
