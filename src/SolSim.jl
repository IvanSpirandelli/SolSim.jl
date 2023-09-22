module SolSim

export DataStructures
export CoordinateCalculations
export SolvationFreeEnergy
export ContactGraphs
export Utilities
export Simulation

export Tests

# Modules
include("../modules/DataStructures/src/DataStructures.jl")
include("../modules/CoordinateCalculations/src/CoordinateCalculations.jl")
include("../modules/ContactGraphs/src/ContactGraphs.jl")
include("../modules/SolvationFreeEnergy/src/SolvationFreeEnergy.jl")
include("../modules/Utilities/src/Utilities.jl")
include("../modules/Simulation/src/Simulation.jl")

# Tests
include("../tests/Tests.jl")

end #module SolSim