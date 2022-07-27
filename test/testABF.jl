#First we can test opening the files.
print("Testing opening the files... ")
test_file = "test/to_analyze.abf"
@time data = readABF(test_file);
#Test opening different test_files
println("Completed")
filtered_data = ePhys.dwt_filter(data, direction = :reverse)

data2 = readABF(test_file; channels=["Vm_prime", "Vm_prime4"]);
#=println("Testing base functionality of ABF extraction")
abf_1swp = ePhys.readABFInfo(target_path1)
abf_12swp = ePhys.readABFInfo(target_path2)

ePhys.getWaveform(abf_1swp, 1, 1; channel_type=:analog) #Test get waveform of analog 0
ePhys.getWaveform(abf_1swp, 1, 1; channel_type=:digital) #Get waveform of digital 0

ePhys.getWaveform(abf_12swp, 1, 2; channel_type=:analog) #get waveform of multisweep analog
ePhys.getWaveform(abf_12swp, 1, 2; channel_type=:digital) #get waveform of multisweep digital

ePhys.getWaveform(abf_12swp, 1; channel_type=:analog) #get waveform of multisweep analog, all sweeps
ePhys.getWaveform(abf_12swp, 1; channel_type=:digital) #get waveform of multisweep analog

#use strings to get the waveforms
ePhys.getWaveform(abf_12swp, 1, "An 0")
ePhys.getWaveform(abf_12swp, 1, "Ana 0")
ePhys.getWaveform(abf_12swp, 1, "Analog 0")
ePhys.getWaveform(abf_12swp, 1, "An 0") #Get all related sweeps to analog 0
ePhys.getWaveform(abf_12swp, 1, "Cmd 0") #Get all related sweeps to analog 0

ePhys.getWaveform(abf_12swp, 1, "D 0")
ePhys.getWaveform(abf_12swp, 1, "Dig 0")
ePhys.getWaveform(abf_12swp, 1, "Digital 0")
ePhys.getWaveform(abf_12swp, 1, "D 0") #Get all related sweeps to digital 0
=#