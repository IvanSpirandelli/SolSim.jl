function ply_to_vec(filepath::String)
    coordinate_data = readlines(filepath)[12:end]
    centers = Vector{SVector{3, Float64}}([])
    radii = Vector{Float64}([])
    for line in coordinate_data
        line_array = split(line, " ")
        push!(centers, @SVector[
        parse(Float64, line_array[1]), 
        parse(Float64, line_array[2]), 
        parse(Float64, line_array[3])
        ])
        push!(radii, parse(Float64, line_array[4]))
    end
    centers, radii
end

function poly_to_vec(filepath::String)
    coordinate_data = readlines(filepath)[2:end-2]
    centers = Vector{SVector{3, Float64}}([])
    radii = Vector{Float64}([])
    # Assuming lines are formatted like "i: x y z r"
    for line in coordinate_data
        line_array = split(line, " ")
        push!(centers, @SVector[
        parse(Float64, line_array[2]), 
        parse(Float64, line_array[3]), 
        parse(Float64, line_array[4])
        ])
        push!(radii, parse(Float64, line_array[5]))
    end
    centers, radii
end

function get_references(path::String)
    lib_path = string(path)
    reference_plys = readdir(lib_path; join=true)
    ref_configurations = Vector{Tuple{Symbol, Atoms}}([])

    for ply in reference_plys
        name = split(last(split(ply, "/")), ".")[1]
        centers, radii = ply_to_vec("$(path)/$(name).ply")
        push!(ref_configurations, (Symbol(name), Atoms{length(centers)}(centers, SVector{length(centers)}(radii))))
    end

    ref_configurations
end

function get_measures_from_ply(filepath, rs)
    centers, radii = ply_to_vec(filepath)
    inflated_radii = radii .+ rs
    get_measures(centers, inflated_radii)
end

function get_energy_from_ply(filepath, rs, prefactors)
    centers, radii = ply_to_vec(filepath)
    inflated_radii = radii .+ rs
    solvation_free_energy(centers, inflated_radii, prefactors)
end

function get_energy_from_poly(filepath, rs, prefactors)
    centers, radii = poly_to_vec(filepath)
    inflated_radii = radii .+ rs
    solvation_free_energy(centers, inflated_radii, prefactors)
end

function get_maximum_contact_clusters(num_of_balls::Int)
    lib_path = string("io/input/mcc/", string(num_of_balls), "_ball_lib")
    ref_configurations = get_references(lib_path)

    ref_configurations
end

