module SolvationFreeEnergy 

export WhiteBearParameters
export solvation_free_energy, solvation_free_energy_for_hard_spheres
export solvation_free_energy_for_hard_spheres_in_bounded_container
export get_measures, get_measures_and_derivatives
export get_prefactors, get_wb_prefactors
export pressure, sigma, kappa, kappa_bar

using StaticArrays
using ..DataStructures
using ..CoordinateCalculations

include("geometric_measures.jl")
include("prefactors.jl")
include("solvation_free_energy.jl")

end