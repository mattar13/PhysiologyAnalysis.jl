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
	using PhysAnalysis
	#Use pyplot? for plotting
end

# ╔═╡ e2fcae6f-d795-4258-a328-1aad5ea64195
md"
#### Current Date: $(Dates.now())

# Analysis of experiments:

#### 1) Point to the experiment root: Contains .abf files

"

# ╔═╡ 7b14b019-7545-4441-833b-f7e660c23dc6
#enter in the file path of the file you would like to analyze here
path = raw"C:\Users\mtarc\OneDrive\Documents\GithubRepositories\PhysAnalysis\test\to_filter.abf"

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
end

# ╔═╡ c8b4c855-64b8-4e18-9b2d-231260c67813
md"""
#### 3) Filter the waveform

Pre Stimulus Time   
$(@bind t_pre_pick Slider(0.0:0.1:2.0; default = 1.0, show_value = true))
Post Stimulus Time
$(@bind t_post_pick Slider(0.0:0.1:10.0; default = 4.0, show_value = true))

Display Stimulus
$(@bind display_stim CheckBox(default = true))
"

Highpass
$(@bind highpass_pick CheckBox(default = false))   
$(@bind highpass_val Slider(0.001:0.001:1.0; default = 1.0, show_value = true))
hz
Lowpass
$(@bind lowpass_pick CheckBox(default = true))
$(@bind lowpass_val Slider(0.0:10.0:1000.0; default = 300.0, show_value = true))
hz

Adaptive Interference Bandpass
$(@bind adap_pick CheckBox(default = true))
$(@bind adap_val Slider(0.0:1.0:1000.0; default = 100.0, show_value = true))

##### These filters are experimental and only work on small noisy waveforms
DWT Filter
$(@bind DWT_pick CheckBox(default = false))
CWT Filter
$(@bind CWT_pick CheckBox(default = false))


$(@bind WT_low_val Slider(1:50; default = 1, show_value = true))
$(@bind WT_hi_val Slider(1:50; default = 9, show_value = true))

"""

# ╔═╡ a1ff2491-73a1-4114-a2bc-3bc6375c583a
#Now we can plot the data here

# ╔═╡ ab425ca0-65ac-4372-a4c8-3d23c8b2f4fb


# ╔═╡ Cell order:
# ╠═a442e068-06ef-4d90-9228-0a03bc6d9379
# ╟─e2fcae6f-d795-4258-a328-1aad5ea64195
# ╟─7b14b019-7545-4441-833b-f7e660c23dc6
# ╟─b400dd0c-5a40-4ee7-9116-7339939b7456
# ╟─05e38576-9650-4287-bac0-6d281db2ea9c
# ╟─971d6f11-2936-4d75-9641-36f81a94c2c4
# ╟─c8b4c855-64b8-4e18-9b2d-231260c67813
# ╠═a1ff2491-73a1-4114-a2bc-3bc6375c583a
# ╠═ab425ca0-65ac-4372-a4c8-3d23c8b2f4fb