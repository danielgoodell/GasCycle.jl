using Pkg; Pkg.activate(joinpath(@__DIR__, ".."))
using Test

include("test_thermo.jl")
include("test_elements.jl")
