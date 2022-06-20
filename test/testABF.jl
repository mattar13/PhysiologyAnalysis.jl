using PhysAnalysis

#First we can test opening the files
print("Testing opening the files... ")
test_file = "test/to_analyze.abf"
data = readABF(test_file)
println("Completed")

