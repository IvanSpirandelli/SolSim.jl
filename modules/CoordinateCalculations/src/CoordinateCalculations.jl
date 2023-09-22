module CoordinateCalculations

using ..DataStructures

using StaticArrays
using Distances

export is_intersecting_set_of_atoms, is_indexed_atom_intersecting_others
export is_center_in_bounds, are_atoms_in_bounds, is_indexed_atom_in_bounds
export are_balls_intersecting

include("coordinate_calculations.jl")

end #module CoordinateCalculations