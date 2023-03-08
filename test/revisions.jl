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
dataset["CONDITIONS"] = summarize_data(dataset)

plot_dataset_fits(dataset, normalize = false, condition = "NoDrugs")
plot_dataset_vals(dataset)