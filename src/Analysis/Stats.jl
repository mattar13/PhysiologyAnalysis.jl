"""
This function is for computing the R-squared of a polynomial
"""
function RSQ(poly::PN.Polynomial, x, y)
     ŷ = poly.(x)
     ȳ = sum(ŷ) / length(ŷ)
     SSE = sum((y - ŷ) .^ 2)
     SST = sum((y .- ȳ) .^ 2)
     1 - SSE / SST
end

function RSQ(ŷ::Array{T}, y::Array{T}) where {T<:Real}
     ȳ = sum(ŷ) / length(ŷ)
     SSE = sum((y - ŷ) .^ 2)
     SST = sum((y .- ȳ) .^ 2)
     1 - SSE / SST
end

function sig_symbol(val)
     if val <= 0.001
          return "***"
     elseif val <= 0.005
          return "**"
     elseif val <= 0.05
          return "*"
     else
          return "-"
     end
end

using StatsBase
function cor_xy(x1, y1, x2, y2; N = 1000)
     tmin = max(minimum(x1), minimum(x2))
     tmax = min(maximum(x1), maximum(x2))
     x_common = range(tmin, tmax, length=N)
 
     # Suppose x1, y1 and x2, y2 are sorted in ascending x
     itp1 = interpolate((x1,), y1, Gridded(Linear()))
     itp2 = interpolate((x2,), y2, Gridded(Linear()))
 
     # Sample them on the common x-grid
     y1_common = itp1.(x_common)
     y2_common = itp2.(x_common)

     #corr_coeff = Statistics.cor(y2_common, y1_common)
     cc = crosscor(y1_common, y2_common)
     return cc
end