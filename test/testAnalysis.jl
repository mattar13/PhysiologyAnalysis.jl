using PhysAnalysis

#finding the saturated response
print("Finding the saturated responses... ")
@time resp = saturated_response(data);

print("Finding the minima to peak... ")
@time m2p = minima_to_peak(data);

print("Finding the time to peak")