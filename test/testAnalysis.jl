#finding the saturated response
print("Testing the saturated response... ")
resps = saturated_response(data);
@test size(resps, 1) == size(data,1)
@test size(resps, 2) == size(data,3)

#rmaxes = maximum(resps, dims=1);

#print("Finding the minima to peak... ")
#@time m2p = minima_to_peak(data);

#print("Finding the time to peak")
#@time t2p = time_to_peak(data);

#print("Finding the dominant time constant... ")
#@time t_dom = percent_recovery_interval(data, rmaxes);

#print("Finding the integration time")
#@time int_time = integral(data);

#println("Extracting some timescale data")

