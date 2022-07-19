using ePhys

#finding the saturated response
print("Finding the saturated responses... ")
@time resps = saturated_response(data);
rmaxes = maximum(resps, dims = 1); 

println("Pull out the dim response")

print("Finding the minima to peak... ")
@time m2p = minima_to_peak(data);

print("Finding the time to peak")
@time t2p = time_to_peak(data);

print("Finding the dominant time constant... ")
@time t_dom = percent_recovery_interval(data, rmaxes);

println("Finding the integration time")
@time int_time = integral(data);

println("Extracting some patch clamp data")
