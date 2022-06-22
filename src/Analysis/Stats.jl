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

rolling_mean(arr::AbstractArray; radius = 5) = [sum(arr[i:i+radius]) / radius for i = 1:length(arr)-radius]

ENV["R_HOME"] = raw"C:\Program Files\R\R-4.1.2"
#Have to point to the correct R path first
using RCall #used for calculating Two-Way ANOVAs

"""
This function calculates the standard error of the Mean and is safe for querying
"""
sem(x) = std(x) / sqrt(length(x))

"""
#This function calls R and computes Two-Way ANOVA
#How to call R objects
R"x = 2"
#How to retrive R objects
@rget x
#How to send objects to R
z = 10.0
@rput z
"""
function R_ANOVA(data::DataFrame)
     #install.packages("dplyr")
     r = RObject(data) #This inserts the julia DataFrame into R
     #convert Genotype into factor
     @rput r
     R"""
     r <- as.data.frame(lapply(r, unlist))#need to convert the dataframe into a new type
     #We need to convert the Age, Genotype and Photons to factor
     r$Genotype <- factor(r$Genotype)
     r$Age <- factor(r$Age) 
     #Generate a frequency table
     freq_table <- table(r$Genotype, r$Age)
     resAOV2 <- aov(Rmax ~ Genotype:Age, data = r)
     summaryAOV2 <- summary(resAOV2)
     resTUKEY <- TukeyHSD(resAOV2)
     """
     freq_table = @rget freq_table #get the frequency table back
     summaryAOV2 = @rget summaryAOV2
     resAOV2 = @rget resAOV2
     resTUKEY = @rget resTUKEY
     return (freq_table, summaryAOV2, resAOV2, resTUKEY)
end


function R_T_TEST(data::DataFrame)
     r = RObject(data) #This inserts the julia DataFrame into R
     #convert Genotype into factor
     @rput r

     R"""
         r <- as.data.frame(lapply(r, unlist))#need to convert the dataframe into a new type
         #We need to convert the Age, Genotype and Photons to factor
         r$Genotype <- factor(r$Genotype)
         freq_table <- table(r$Genotype, r$Age)
         pairTTest <- pairwise.t.test(r$Rmax, r$Genotype, p.adjust.method = "BH")
     """
     ptt = @rget pairTTest
     return ptt
end