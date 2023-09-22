# Basic variants
function solvation_free_energy(centers::Vector{SVector{3, Float64}}, radii::AbstractVector, prefactors::SVector{4, Float64})
    measures = get_measures(centers, radii)
    measures[1] * prefactors[1] + measures[2] * prefactors[2] + measures[3] * prefactors[3] + measures[4] * prefactors[4]
end

function solvation_free_energy(centers::Vector{SVector{3, Float64}}, radii::AbstractVector, rs::Float64, η::Float64)
    prefactors = get_prefactors_from_white_bear(rs, η)
    solvation_free_energy(centers, radii, prefactors)
end

function solvation_free_energy(perturb_out::Nothing, centers::Vector{SVector{3, Float64}}, radii::AbstractVector, rs::Float64, η::Float64)
    inflated_radii = get_flattened_inflated_radii(radii, rs, 1)
    prefactors = get_prefactors_from_white_bear(rs, η)    
    solvation_free_energy(centers, inflated_radii, prefactors)
end

function solvation_free_energy(index::Int, centers::Vector{SVector{3, Float64}}, radii::AbstractVector, rs::Float64, η::Float64)
    inflated_radii = get_flattened_inflated_radii(radii, rs, 1)
    prefactors = get_prefactors_from_white_bear(rs, η)
    solvation_free_energy(centers, inflated_radii, prefactors)
end

# Including Hard Sphere Check
function solvation_free_energy(perturb_out::Nothing, solute::Atoms, rs::Float64, η::Float64)
    if is_intersecting_set_of_atoms(solute) return Inf end
    inflated_radii = get_flattened_inflated_radii(solute.radii, rs, 1)
    prefactors = get_prefactors_from_white_bear(rs, η)
    solvation_free_energy(solute.centers, inflated_radii, prefactors)
end

function solvation_free_energy(index::Int, solute::Atoms, rs::Float64, η::Float64)
    if is_indexed_atom_intersecting_others(index, solute) return Inf end
    inflated_radii = get_flattened_inflated_radii(solute.radii, rs, 1)
    prefactors = get_prefactors_from_white_bear(rs, η)
    solvation_free_energy(solute.centers, inflated_radii, prefactors)
end

function solvation_free_energy_for_hard_spheres(perturb_out::Nothing, centers::Vector{SVector{3, Float64}}, radii::AbstractVector, rs::Float64, prefactors::SVector{4, Float64})
    if is_intersecting_set_of_atoms(centers, radii) return Inf end
    inflated_radii = get_flattened_inflated_radii(radii, rs, 1)
    solvation_free_energy(centers, inflated_radii, prefactors)
end

function solvation_free_energy_for_hard_spheres(index::Int, centers::Vector{SVector{3, Float64}}, radii::AbstractVector, rs::Float64, prefactors::SVector{4, Float64})
    if is_indexed_atom_intersecting_others(index, centers, radii) return Inf end
    inflated_radii = get_flattened_inflated_radii(radii, rs, 1)
    solvation_free_energy(centers, inflated_radii, prefactors)
end

function solvation_free_energy_for_hard_spheres(perturb_out::Nothing, centers::Vector{SVector{3, Float64}}, radii::AbstractVector, rs::Float64, η::Float64)
    if is_intersecting_set_of_atoms(centers, radii) return Inf end
    inflated_radii = get_flattened_inflated_radii(radii, rs, 1)
    prefactors = get_prefactors_from_white_bear(rs, η)    
    solvation_free_energy(centers, inflated_radii, prefactors)
end

function solvation_free_energy_for_hard_spheres(index::Int, centers::Vector{SVector{3, Float64}}, radii::AbstractVector, rs::Float64, η::Float64)
    if is_indexed_atom_intersecting_others(index, centers, radii) return Inf end
    inflated_radii = get_flattened_inflated_radii(radii, rs, 1)
    prefactors = get_prefactors_from_white_bear(rs, η)
    solvation_free_energy(centers, inflated_radii, prefactors)
end

function solvation_free_energy_for_hard_spheres_in_bounding_box(perturb_out::Nothing, centers::Vector{SVector{3, Float64}}, radii::AbstractVector, rs::Float64, prefactors::SVector{4, Float64}, bounds::Float64)
    if is_intersecting_set_of_atoms(centers, radii) || !are_atoms_in_bounds(centers, bounds)
        return Inf 
    end
    inflated_radii = get_flattened_inflated_radii(radii, rs, 1)
    solvation_free_energy(centers, inflated_radii, prefactors)
end

function solvation_free_energy_for_hard_spheres_in_bounding_box(index::Int, centers::Vector{SVector{3, Float64}}, radii::AbstractVector, rs::Float64, prefactors::SVector{4, Float64}, bounds::Float64)
    if is_indexed_atom_intersecting_others(index, centers, radii) || !is_indexed_atom_in_bounds(index, centers, bounds)
        return Inf 
    end
    inflated_radii = get_flattened_inflated_radii(radii, rs, 1)
    solvation_free_energy(centers, inflated_radii, prefactors)
end

function solvation_free_energy_for_hard_spheres_in_bounding_box(perturb_out::Nothing, centers::Vector{SVector{3, Float64}}, radii::AbstractVector, rs::Float64, η::Float64, bounds::Float64)
    if is_intersecting_set_of_atoms(centers, radii) || !are_atoms_in_bounds(centers, bounds)
        return Inf 
    end
    inflated_radii = get_flattened_inflated_radii(radii, rs, 1)
    prefactors = get_prefactors_from_white_bear(rs, η)    
    solvation_free_energy(centers, inflated_radii, prefactors)
end

function solvation_free_energy_for_hard_spheres_in_bounding_box(index::Int, centers::Vector{SVector{3, Float64}}, radii::AbstractVector, rs::Float64, η::Float64, bounds::Float64)
    if is_indexed_atom_intersecting_others(index, centers, radii) || !is_indexed_atom_in_bounds(index, centers, bounds)
        return Inf 
    end
    inflated_radii = get_flattened_inflated_radii(radii, rs, 1)
    prefactors = get_prefactors_from_white_bear(rs, η)
    solvation_free_energy(centers, inflated_radii, prefactors)
end

function get_flattened_inflated_radii(radii::SVector, rs::Float64, n::Int)
    flat_inflated_radii = Base.Iterators.collect(Base.Iterators.flatten([radii .+ rs for _ in 1:n]))
    SVector(ntuple(i -> flat_inflated_radii[i], length(flat_inflated_radii)))
end

function get_flattened_inflated_radii(radii::Vector, rs::Float64, n::Int)
    Base.Iterators.collect(Base.Iterators.flatten([radii .+ rs for _ in 1:n]))
end

function get_flattened_inflated_radii(radii::Vector{Vector{Float64}}, rs::Float64)
    reduce(vcat, [m_radii .+ rs for m_radii in radii])
end

