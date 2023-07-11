using Dates
using Revise
using PhysiologyAnalysis
using ElectroPhysiology
using DataFrames, Query, XLSX
PhysiologyAnalysis.__init__()

#%% We want to add a data analysis to writeXLSX
file_test = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\ERG\Retinoschisis\2022_12_20_WTAdult\Mouse2_WT_Adult\BaCl_LAP4\Rods" |> parseABF
data = readABF(file_test) |> data_filter
downsample!(data, 100.0)
writeXLSX("test.xlsx", data, :analysis; verbose = true)

verbose = false
filenames = file_test
dataset = createDataset(filenames, verbose = verbose)
dataset = runTraceAnalysis(dataset, verbose = verbose)
qPhotons = dataset["TRACES"] |> @unique(_.Photons) |> DataFrame
@assert length(photons.Photons) == size(data, 1)
qPhotons.Photons
setIntensity(data.stimulus_protocol, qPhotons.Photons) #Set the intensity of the flash intensities

#We want to pull out a specific experiment or trace from the dataframe
data_root = raw"C:\Users\mtarc\OneDrive - The University of Akron\Journal Submissions\2023_05_24 Critical Timepoints in retinoschisis - Unknown Journal" #The data root
data_file = joinpath(data_root, "P14_data_analysis.xlsx")
dataset = openDataset(data_file , verbose = true)

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