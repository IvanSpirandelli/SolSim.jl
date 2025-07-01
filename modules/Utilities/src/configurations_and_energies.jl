function get_dispersed_conf(n::Int, interaction_radius::Float64)
    v = [@SVector[10*interaction_radius * i,0,0] for i in 1:n]
end

function get_dispersed_conf(template_particle::Vector{SVector{3, Float64}}, n::Int, interaction_radius::Float64)
    [[e + @SVector[10*interaction_radius * i, 0.0, 0.0] for e in template_particle] for i in 1:n]
end

function get_dispersed_conf(template_particle::Vector{SVector{3, Float64}}, radii::Vector{Float64}, n::Int)
    bounding_radius = get_bounding_radius(template_particle, radii)
    [[c + @SVector[10*bounding_radius * i, 0.0, 0.0] for c in template_particle] for i in 1:n]
end

function get_scaled_configuration(configuration, sf)
    [@SVector[elem[1] * sf, elem[2] * sf, elem[3] * sf] for elem in configuration]
end

function get_dispersed_energy(solute_radii::Vector{Float64}, rs::Float64, prefactors::SVector{4, Float64})
    n = length(solute_radii)
    conf = get_dispersed_conf(n, maximum(solute_radii) + rs)
    solvation_free_energy(conf, solute_radii, rs, prefactors)
end

function get_dispersed_energy(template_particle::Vector{SVector{3, Float64}}, n::Int, solute_radii::SVector, rs::Float64, prefactors::SVector{4, Float64})
    conf = get_dispersed_conf(template_particle, n, sum(solute_radii) + rs)
    solvation_free_energy(conf, solute_radii, rs, prefactors)
end

# Assumes solute radius of 1.0
function get_minimizing_scaling(conf, R, rs, prefactors; ϵ = 0.0001)
    E = Inf
    sf_min = Inf
    for sf in 1.0:ϵ:(R+rs)/R
        E_cand = solvation_free_energy(get_scaled_configuration(conf, sf), R, rs, prefactors)
        E,sf_min = E_cand < E ? (E_cand,sf) : (E,sf_min)
    end
    sf_min
end

# Assumes solute radius of 1.0
function get_minimizing_scaling_and_energy(conf, R, rs, prefactors)
    E = Inf
    sf_min = Inf
    step = rs / 10.0
    for sf in 1.0:step:(R+rs)/R
        E_cand = solvation_free_energy(get_scaled_configuration(conf, sf), R, rs, prefactors)
        E,sf_min = E_cand < E ? (E_cand,sf) : (E,sf_min)
    end
    sf_min, E
end