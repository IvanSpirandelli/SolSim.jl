# Basic variants
function solvation_free_energy(centers::Vector{SVector{3, Float64}}, radii::AbstractVector, rs::Float64, prefactors::SVector{4, Float64})
    measures = get_measures(centers, radii, rs)
    measures[1] * prefactors[1] + measures[2] * prefactors[2] + measures[3] * prefactors[3] + measures[4] * prefactors[4]
end

function solvation_free_energy(centers::Vector{SVector{3, Float64}}, radii::AbstractVector, rs::Float64, η::Float64)
    prefactors = get_prefactors(rs, η)
    solvation_free_energy(centers, radii, rs, prefactors)
end

function solvation_free_energy(perturb_out::Nothing, centers::Vector{SVector{3, Float64}}, radii::AbstractVector, rs::Float64, η::Float64)
    prefactors = get_prefactors(rs, η)    
    solvation_free_energy(centers, radii, rs, prefactors)
end

function solvation_free_energy(index::Int, centers::Vector{SVector{3, Float64}}, radii::AbstractVector, rs::Float64, η::Float64)
    prefactors = get_prefactors(rs, η)
    solvation_free_energy(centers, radii, rs, prefactors)
end

# Including Hard Sphere Check
function solvation_free_energy(perturb_out::Nothing, solute::Atoms, rs::Float64, η::Float64)
    if is_intersecting_set_of_atoms(solute) return Inf end
    prefactors = get_prefactors(rs, η)
    solvation_free_energy(solute.centers, solute.radii, rs, prefactors)
end

function solvation_free_energy(index::Int, solute::Atoms, rs::Float64, η::Float64)
    if is_indexed_atom_intersecting_others(index, solute) return Inf end
    prefactors = get_prefactors(rs, η)
    solvation_free_energy(solute.centers, solute.radii, rs, prefactors)
end

function solvation_free_energy_for_hard_spheres(perturb_out::Nothing, centers::Vector{SVector{3, Float64}}, radii::AbstractVector, rs::Float64, prefactors::SVector{4, Float64})
    if is_intersecting_set_of_atoms(centers, radii) return Inf end
    solvation_free_energy(centers, radii, rs, prefactors)
end

function solvation_free_energy_for_hard_spheres(index::Int, centers::Vector{SVector{3, Float64}}, radii::AbstractVector, rs::Float64, prefactors::SVector{4, Float64})
    if is_indexed_atom_intersecting_others(index, centers, radii) return Inf end
    solvation_free_energy(centers, radii, rs, prefactors)
end

function solvation_free_energy_for_hard_spheres(perturb_out::Nothing, centers::Vector{SVector{3, Float64}}, radii::AbstractVector, rs::Float64, η::Float64)
    if is_intersecting_set_of_atoms(centers, radii) return Inf end
    prefactors = get_prefactors(rs, η)    
    solvation_free_energy(centers, radii, rs, prefactors)
end

function solvation_free_energy_for_hard_spheres(index::Int, centers::Vector{SVector{3, Float64}}, radii::AbstractVector, rs::Float64, η::Float64)
    if is_indexed_atom_intersecting_others(index, centers, radii) return Inf end
    prefactors = get_prefactors(rs, η)
    solvation_free_energy(centers, radii, rs, prefactors)
end

function solvation_free_energy_for_hard_spheres_in_bounding_box(perturb_out::Nothing, centers::Vector{SVector{3, Float64}}, radii::AbstractVector, rs::Float64, prefactors::SVector{4, Float64}, bounds::Float64)
    if is_intersecting_set_of_atoms(centers, radii) || !are_atoms_in_bounds(centers, bounds)
        return Inf 
    end
    solvation_free_energy(centers, radii, rs, prefactors)
end

function solvation_free_energy_for_hard_spheres_in_bounding_box(index::Int, centers::Vector{SVector{3, Float64}}, radii::AbstractVector, rs::Float64, prefactors::SVector{4, Float64}, bounds::Float64)
    if is_indexed_atom_intersecting_others(index, centers, radii) || !is_indexed_atom_in_bounds(index, centers, bounds)
        return Inf 
    end
    solvation_free_energy(centers, radii, rs, prefactors)
end

function solvation_free_energy_for_hard_spheres_in_bounding_box(perturb_out::Nothing, centers::Vector{SVector{3, Float64}}, radii::AbstractVector, rs::Float64, η::Float64, bounds::Float64)
    if is_intersecting_set_of_atoms(centers, radii) || !are_atoms_in_bounds(centers, bounds)
        return Inf 
    end
    prefactors = (rs, η)    
    solvation_free_energy(centers, radii, rs, prefactors)
end

function solvation_free_energy_for_hard_spheres_in_bounding_box(index::Int, centers::Vector{SVector{3, Float64}}, radii::AbstractVector, rs::Float64, η::Float64, bounds::Float64)
    if is_indexed_atom_intersecting_others(index, centers, radii) || !is_indexed_atom_in_bounds(index, centers, bounds)
        return Inf 
    end
    prefactors = get_prefactors(rs, η)
    solvation_free_energy(centers, radii, rs, prefactors)
end

function get_flattened_inflated_radii(radii::SVector, rs::Float64, n::Int)
    flat_inflated_radii = collect(Base.Iterators.flatten([radii .+ rs for _ in 1:n]))
    SVector(ntuple(i -> flat_inflated_radii[i], length(flat_inflated_radii)))
end

function get_flattened_inflated_radii(radii::Vector, rs::Float64, n::Int)
    collect(Base.Iterators.flatten([radii .+ rs for _ in 1:n]))
end

function get_flattened_inflated_radii(radii::Vector{Vector{Float64}}, rs::Float64)
    reduce(vcat, [m_radii .+ rs for m_radii in radii])
end

