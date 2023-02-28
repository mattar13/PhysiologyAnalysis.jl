using Documenter, PhysiologyAnalysis

makedocs(
     sitename = "PhysiologyAnalysis.jl",
     modules = [PhysiologyAnalysis],
     pages = [
          "Home" => "index.md"
          "Tutorial" => "tutorial.md"
     ]
)

#deploydocs(
     #repo = "github.com/mattar13/PhysiologyAnalysis.jl.git",
     #target = "build"
#)