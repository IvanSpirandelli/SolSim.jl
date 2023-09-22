using Dates

function single_move_periodic(solute::SolSim.DataStructures.Atoms, iterations::Int, rs::Float64, η::Float64)
    σ = rs
    n = length(solute.radii)
    bounds = sum([solute.radii[i] + rs for i in 1:n]) * 2.0
    interaction_range = 2.0 * (maximum(solute.radii) + rs)
    perturb!(centers_cand, centers) = SolSim.Simulation.translate_single_modulo_bounds!(centers_cand, centers, σ, interaction_range, bounds)
    energy(perturb_out, centers) = SolSim.SolvationFreeEnergy.solvation_free_energy_for_hard_spheres(perturb_out, centers, solute.radii, rs, η)
    
    algorithm = SolSim.Simulation.SAModifiedLam(
        energy,
        perturb!,
        )

    
    start = now()
    solute, E = SolSim.Simulation.simulate(
        algorithm,
        solute,
        iterations
    ) 
    println("Simulation took  $(round((now() - start).value/60_000))min, or $(round((now() - start).value/1000))s ")

    SolSim.Utilities.vec_to_poly(solute, "io/output/modified_lam/result")
    E
end

function single_move_bounded(solute::Union{SolSim.DataStructures.Atoms, SolSim.DataStructures.EqualMolecules}, iterations::Int, rs::Float64, η::Float64, bounds::Float64; run_id::String = "")
    @assert SolSim.CoordinateCalculations.is_solute_in_bounds(solute, bounds)
    σ = rs
    perturb!, energy = get_perturbation_and_energy_for_bounded_box(solute, rs, η, bounds, σ)

    algorithm = SolSim.Simulation.SAModifiedLam(
        energy,
        perturb!,
        )

    start = now()
    solute, E = SolSim.Simulation.simulate(
        algorithm,
        solute,
        iterations
    ) 
    println("Simulation took  $(round((now() - start).value/60_000))min, or $(round((now() - start).value/1000))s ")
    
    SolSim.Utilities.vec_to_poly(solute, "io/output/modified_lam/result$(run_id)")
    E
end

function single_move_bounded_debug(solute::Union{SolSim.DataStructures.Atoms, SolSim.DataStructures.EqualMolecules}, iterations::Int, rs::Float64, η::Float64, bounds::Float64; run_id::String = "current")
    @assert SolSim.CoordinateCalculations.is_solute_in_bounds(solute, bounds)
    σ = rs
    perturb!, energy = get_perturbation_and_energy_for_bounded_box(solute, rs, η, bounds, σ)

    algorithm = SolSim.Simulation.SAModifiedLam(
        energy,
        perturb!,
        )

    info = SolSim.Simulation.DetailedSimulationInfo(Vector{Float64}([]), Vector{Float64}([]), Vector{Float64}([]), Vector{Float64}([]))
    start = now()
    solute, E = SolSim.Simulation.simulate(
        algorithm,
        solute,
        iterations,
        info
    ) 
    println("Simulation took  $(round((now() - start).value/60_000))min, or $(round((now() - start).value/1000))s ")

    path = "io/output/modified_lam/debug/$(run_id)/"
    mkpath(path)

    SolSim.Simulation.plot_simulation_info(info, path)
    SolSim.Utilities.vec_to_poly(solute, "$(path)result")
    E
end

function get_perturbation_and_energy_for_bounded_box(solute::SolSim.DataStructures.Atoms, rs, η, bounds, σ)
    perturb!(centers_cand, centers) = SolSim.Simulation.translate_single!(centers_cand, centers, σ)
    energy(perturb_out, centers) = SolSim.SolvationFreeEnergy.solvation_free_energy_for_hard_spheres_in_bounding_box(perturb_out, centers, solute.radii, rs, η, bounds)
    return perturb!, energy
end

function get_perturbation_and_energy_for_bounded_box(solute::SolSim.DataStructures.EqualMolecules, rs, η, bounds, σ)
    perturb!(centers_cand, centers) = SolSim.Simulation.translate_or_rotate_single!(centers_cand, centers, σ; mra = 0.05)
    energy(perturb_out, centers) = SolSim.SolvationFreeEnergy.solvation_free_energy_for_hard_spheres_in_bounding_box(perturb_out, centers, solute.radii, rs, η, bounds)
    return perturb!, energy
end