function is_center_in_bounds(center::SVector{3,Float64}, bounds::Float64)
    all(center .< SA_F64[bounds, bounds, bounds]) && all(center .> SA_F64[0.0, 0.0, 0.0])
end

function are_atoms_in_bounds(centers::Vector{SVector{3,Float64}}, bounds::Float64)
    all([is_center_in_bounds(center, bounds) for center in centers])
end

function is_indexed_atom_in_bounds(index::Int, centers::Vector{SVector{3,Float64}}, bounds::Float64)
    is_center_in_bounds(centers[index], bounds)
end

function are_balls_intersecting(center_one::SVector{3,Float64}, center_two::SVector{3,Float64}, radius_one::Float64, radius_two::Float64)
    return euclidean(center_one, center_two) < radius_one + radius_two
end

function is_intersecting_set_of_atoms(centers::Vector{SVector{3,Float64}}, radii::AbstractVector)
    n = length(centers)
    for i in 1:n
        for j in i+1:n
            if are_balls_intersecting(centers[i], centers[j], radii[i], radii[j])
                return true
            end
        end
    end
    false
end

function is_indexed_atom_intersecting_others(index::Int, centers::Vector{SVector{3,Float64}}, radii::AbstractVector)
    n = length(centers)
    for i in 1:n
        if i != index && are_balls_intersecting(centers[i], centers[index], radii[i], radii[index])
            return true
        end
    end
    false
end

function is_intersecting_set_of_atoms(atoms::Atoms)
    is_intersecting_set_of_atoms(atoms.centers, atoms.radii)
end

function is_indexed_atom_intersecting_others(index::Int, atoms::Atoms)
    is_indexed_atom_intersecting_others(index, atoms.centers, atoms.radii)
end