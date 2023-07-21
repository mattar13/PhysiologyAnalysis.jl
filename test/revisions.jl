using Dates
using Revise
using PhysiologyAnalysis
using ElectroPhysiology
using DataFrames, Query, XLSX
PhysiologyAnalysis.calibration_path()
PhysiologyAnalysis.set_calibration_path(raw"C:\Users\mtarc\The University of Akron\Renna Lab - General\Data\Calibrations\photon_lookup.xlsx")

#%% We want to pull out a specific experiment or trace from the dataframe
data_file = raw"C:\Users\mtarc\OneDrive - The University of Akron\Journal Submissions\2023_05_24 Critical Timepoints in retinoschisis - Unknown Journal\data_analysis.xlsx" #The data root
dataset = openDataset(data_file , verbose = true)

exps = dataset["EXPERIMENTS"] |> 
     @filter(_.Condition == "BaCl_LAP4" &&_.Photoreceptor == "Rods") |> 
     @filter(_.Age == "P11") |> 
     @map({Date = _.Date, Number = _.Number, Channel = _.Channel, Genotype = _.Genotype}) |> 
DataFrame
RS_trace_G = dataset["TRACES"] |> @filter(_.Condition == "NoDrugs" &&_.Photoreceptor == "Rods") |> DataFrame

save_loc = raw"C:\Users\mtarc\OneDrive - The University of Akron\Journal Submissions\2023_05_24 Critical Timepoints in retinoschisis - Unknown Journal\Plot_Outputs"
for exps_sel in eachrow(exps) 
     name = "$(exps_sel.Date)_$(exps_sel.Number)_$(exps_sel.Channel)_$(exps_sel.Genotype)"
     println(name)
     #color_sel = color_dict[exps_sel.Genotype]
     #println(color_sel)
     #fig_sel = plot_timepoint(RS_trace_A, RS_trace_B, RS_trace_G, exps_sel; color = color_sel, 
     #     alims = (-100, 100), 
     #     blims = (-50, 200), 
     #     glims = (-200, 50)     
     #)
     #fig_sel.savefig("$(save_loc)\\$name.png")
     
     #qA = matchExperiment(RS_trace_A, exps_sel)
     #dataA = readABF(qA, channels = exps_sel.Channel) |> data_filter
     #writeXLSX("$(save_loc)\\$(name)_A.xlsx", dataA, :analysis)
     
     #qB = matchExperiment(RS_trace_B, exps_sel)
     #if !isempty(qB)
     #     dataB = readABF(qB, channels = exps_sel.Channel) |> data_filter
     #     writeXLSX("$(save_loc)\\$(name)_B.xlsx", dataB, :analysis)
     #end
     qG = matchExperiment(RS_trace_G, exps_sel)
     #println(qG)
     if !isempty(qG)
          dataG = readABF(qG, channels = exps_sel.Channel) |> data_filter
          writeXLSX("$(save_loc)\\$(name)_G.xlsx", dataG, :analysis)
     end
     #plt.close("all")
end

#%% Make a function that prints out a portion of a dataset
info = (Date = Date(2019, 09, 25), Number = "1")
dataset["TRACES"]
dataset

inc_dataset = matchDataset(dataset, info)
exc_dataset = excludeDataset(dataset, info)

#Now we want to pull out all files in inc dataset
reinc_dataset["ALL_FILES"].Path

#%% Make a new function for adding photon numbers to a flash stimulus
files = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\ERG\Retinoschisis\2019_03_12_AdultWT\Mouse1_Adult_WT\BaCl" |> parseABF
dataset = createDataset(files, debug = true, verbose = true)
data = readABF(files) |> data_filter
photons = dataset["ALL_FILES"].Photons ./ dataset["ALL_FILES"].Stim_Time
setIntensity(data, photons)
getIntensity(data)