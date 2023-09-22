module Simulation

export describe
export Algorithm, algorithm_random_walk
export MetropolisParameters, Metropolis
export SAModifiedLam, EndOfSimulationInfo

export simulate
export temp_finder

export zig_zag, zig_zag_adaptive
export exponential_additive, linear_additive, quadratic_additive, geometric_cooling
export non_monotonic_adaptive_cooling_factor
export get_shift_vector, shift_centers!
export initialize_solute, initialize_atoms, initialize_two_tmv, initialize_solute_in_bounds

export plot_simulation_info

using ..DataStructures
using ..CoordinateCalculations
using ..SolvationFreeEnergy
using ..ContactGraphs
using ..Utilities

using StaticArrays
using Distributions 
using Distances
using CairoMakie

include("perturbation.jl")
include("algorithms/temperature.jl")
include("initialization/solute_initialization.jl")
include("initialization/temp_finder.jl")
include("algorithms/simulated_annealing.jl")
include("algorithms/sa_modified_lam.jl")
include("algorithms/metropolis.jl")
include("output_plots/simulation_info_plotting.jl")

struct EndOfSimulationInfo
    final_energy::Float64
    acceptance_rate::Float64
end

end #module Simulation