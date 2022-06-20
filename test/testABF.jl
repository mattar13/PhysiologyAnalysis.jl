using PhysAnalysis

#First we can test opening the files
print("Testing opening the files... ")
test_file = "test/to_analyze.abf"
@time data = readABF(test_file);
#Test opening different test_files
data2 = readABF(test_file; channels=["Vm_prime", "Vm_prime4"]);
println("Completed")