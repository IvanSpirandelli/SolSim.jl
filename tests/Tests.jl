module Tests

using SolSim.DataStructures
using SolSim.CoordinateCalculations
using SolSim.SolvationFreeEnergy
using SolSim.Utilities
using SolSim.Simulation

using StaticArrays
using Test 
using Mmap

export runtests

include("test_coordinate_calculations.jl")
include("test_solvation_free_energy.jl")
include("test_utilities.jl")
include("test_simulation.jl")
include("test_perturbation.jl")

function runtests()
    @testset verbose = true "SolSim Tests" begin
        run_all_coordinate_calculation_tests()
        run_all_perturbation_tests()
        run_all_solvation_free_energy_tests()
        run_all_utilities_tests()
        run_all_simulation_tests()
    end
end 

end