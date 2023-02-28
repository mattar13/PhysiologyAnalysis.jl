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

# ╔═╡ c8d46060-b3fe-11ed-0887-790bb5dd2fc0
begin
	using Pkg
	Pkg.activate("../../")
	using Dates, PlutoUI
	using PhysiologyAnalysis
	using DataFrames, XLSX, Query
end

# ╔═╡ 9a5153c5-9c2d-4821-b72c-5e79aa0693b1
md"
# Datasheet analysis

### 1) Enter the path for the datasheet
"

# ╔═╡ fe369c30-a581-4b9d-a5e2-d194db4d61c6
datafile = "C:\\Users\\mtarc\\OneDrive - The University of Akron\\Projects\\Retinoschisis\\data_analysis.xlsx"

# ╔═╡ da89edc7-04ef-4bf7-923d-f48a2fa3a12c
dataset = openDataset(datafile, sheetName="all")

# ╔═╡ 1b247067-3d7d-4fb8-82a1-17328dafe734
begin

	exps = dataset["EXPERIMENTS"] |> @map({
		RSQ = round.(_.RSQ_fit, digits = 2), 
		_.Year, _.Month, _.Date, _.Number, _.Channel, _.Condition,
		}) |> 
		@orderby(_.RSQ) |>
	collect
md"

### 2) Pick a experiment from the list
We need to pick an experiment from the list of experiments loaded
$(@bind selectExps MultiSelect(exps))
"
end

# ╔═╡ bfa2b45b-fcb1-4f66-95a5-35cb4461b4cb
begin
	if !isnothing(selectExps)
		for selectExp in selectExps
			println(selectExp |> typeof)
			YEAR = selectExp[:Year]
			MONTH = selectExp[:MON]
			DAY = selectExp[:DAY]
			#NUM = 
			println(YEAR)
		end
	end
end

# ╔═╡ b1ca0549-2d4d-4c7d-8c07-29aa3e668329
haskey(selectExps[1], :Year)

# ╔═╡ a19b3fc1-6df3-4722-9ce8-76b78ecd3bda
selectExps[1].Year

# ╔═╡ Cell order:
# ╠═c8d46060-b3fe-11ed-0887-790bb5dd2fc0
# ╟─9a5153c5-9c2d-4821-b72c-5e79aa0693b1
# ╠═fe369c30-a581-4b9d-a5e2-d194db4d61c6
# ╟─da89edc7-04ef-4bf7-923d-f48a2fa3a12c
# ╟─1b247067-3d7d-4fb8-82a1-17328dafe734
# ╠═bfa2b45b-fcb1-4f66-95a5-35cb4461b4cb
# ╠═b1ca0549-2d4d-4c7d-8c07-29aa3e668329
# ╠═a19b3fc1-6df3-4722-9ce8-76b78ecd3bda
