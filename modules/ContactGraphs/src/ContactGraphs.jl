module ContactGraphs

using Graphs
using Graphs.Experimental
using StaticArrays
using Distances

using ..CoordinateCalculations

export construct_graph_from_configuration
export get_elastic_graph_with_target_edge_number
export try_to_match_edge_count_of_references

include("contact_graphs.jl")

end #module ContactGraphs
