println("Testing baseline adjust... ")
print("Inplace completed: ")
@time baseline_adjust!(data);
print("Function completed: ")
@time data_baseline = baseline_adjust(data);
println("Completed! ")

print("Testing Truncate functions")
@time data_trunc = truncate_data(data);
print("Inplace: ")
@time truncate_data!(data);
println("Testing filtering functions.")

print("Lowpass filter inplace:")
@time lowpass_filter!(data);
print("Lowpass filter: ")
@time filt_data = lowpass_filter(data);
print("Lowpass filter (choosing freq): ")
@time lopass_data = lowpass_filter(data, 300.0);

print("Highpass filter inplace:")
@time highpass_filter!(data);
print("Highpass filter: ")
@time hipass_data = highpass_filter(data);
print("Highpass filter (choosing freq): ")
@time hipass_data = highpass_filter(data, 60.0);

print("Notch filter inplace:")
@time notch_filter!(data);
print("Notch filter: ")
@time hipass_data = notch_filter(data);