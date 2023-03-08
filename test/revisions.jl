#%% Section 1. Revision of some plotting summarys
using Revise
using PhysiologyAnalysis
using PyPlot
PyPlot.pygui(true)
import PhysiologyAnalysis: readABF, parseABF
import PhysiologyAnalysis: DataFrame
using DataFrames, Query, XLSX

data_root = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\ERG\Retinoschisis" #The data root |> parseABF
datafile = "C:\\Users\\mtarc\\OneDrive - The University of Akron\\Projects\\Retinoschisis\\data_analysis.xlsx"
dataset = openDataset(datafile)
plot_dataset_fits(dataset, normalize = false)

plot_dataset_vals(dataset)

qCOND = dataset["CONDITIONS"] |> @filter(_.Age == "Adult" && _.Genotype == "WT" && _.Condition == "BaCl_LAP4" && _.Photoreceptor == "Rods") |> DataFrame
qEXP = dataset["EXPERIMENTS"] |> @filter(_.Age == "Adult" && _.Genotype == "WT" && _.Condition == "BaCl_LAP4" && _.Photoreceptor == "Rods") |> DataFrame
qTRACE = dataset["TRACES"] |> @filter(_.Age == "Adult" && _.Genotype == "WT" && _.Condition == "BaCl_LAP4" && _.Photoreceptor == "Rods") |> DataFrame

dataset["CONDITIONS"] = summarize_data(dataset)
saveDataset(dataset, datafile)

qC = dataset["CONDITIONS"] |> @filter(_.Age == "Adult" && _.Photoreceptor == "Cones" && _.Condition == "BaCl_LAP4") |> DataFrame
qC.RMAX_COLL
qC.K_COLL