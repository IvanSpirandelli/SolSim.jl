function initialize_solute(radii::SVector)
    n = length(radii)
    mR = maximum(radii)
    initialize_atoms("", radii, n*2.0*mR)
end

function initialize_solute(radii::SVector, rs::Float64)
    n = length(radii)
    mR = maximum(radii)
    initialize_atoms("", radii, n*2.0*(mR + rs))
end

function initialize_solute_in_bounds(radii::AbstractVector, bounds::Float64)
    initialize_atoms("", radii, bounds)
end

function initialize_atoms(hot_start_filepath::String, radii::SVector, bounds::Float64; attempts_at_random_initialization::Int = 10)
    if hot_start_filepath == ""
        return initialize_atoms(radii, bounds, attempts_at_random_initialization)
    end
    return get_hard_spheres_from_hot_start_file(hot_start_filepath)
end

function get_hard_spheres_from_hot_start_file(filepath::String)
    centers, radii = ply_to_vec(filepath)
    solute = Atoms{length(centers)}(centers, SVector{length(centers)}(radii))
    @assert !is_intersecting_set_of_atoms(solute)
    solute
end

function initialize_atoms(radii::SVector, bounds::Float64, attempts_at_random_initialization::Int)
    n = length(radii)
    for i in 1:attempts_at_random_initialization
        centers = [@SVector[rand(Uniform(0.0, bounds)),rand(Uniform(0.0, bounds)),rand(Uniform(0.0, bounds))] for i in 1:n]
        solute = Atoms{length(centers)}(centers, radii)
        if !is_intersecting_set_of_atoms(solute)
            return solute
        end
    end

    elements_per_dim = ceil(cbrt(n))
    length_per_element = bounds/elements_per_dim
    offset = length_per_element / 2.0

    centers = [
        @SVector[
            i * length_per_element - offset, 
            j * length_per_element - offset, 
            k * length_per_element - offset
        ] for i in 1:elements_per_dim for j in 1:elements_per_dim for k in 1:elements_per_dim]
        
    solute = Atoms{length(centers)}(centers, radii)
    if !is_intersecting_set_of_atoms(solute)
        return solute
    end
    println("Could not initialize particles. Bounds might be to small.")
    return nothing
end
