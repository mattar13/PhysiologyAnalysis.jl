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