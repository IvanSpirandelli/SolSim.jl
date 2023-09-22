using SolSim
using StaticArrays
using Dates

function single(;
    solute_radii = @SVector[1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
    rs = 0.475,
    η = 0.475,
    iterations = 5000000,
    zig_zag = [1.0, 0.5, 0.3],
    temp_finder_iterations = 100000,
    target_acceptance_rate = 0.5,
    file_save_root = "io/output/singles/",
    hot_start_file = ""
    )

    n = length(solute_radii)
    bounds = sum([solute_radii[i] + rs for i in 1:n]) * 2.0

    solute = SolSim.Simulation.initialize_atoms(hot_start_file, solute_radii, bounds)
    prefactors = SolSim.SolvationFreeEnergy.get_prefactors_from_white_bear(rs, η)

    σ = max(0.1, rs)
    
    max_energy = SolSim.Utilities.get_dispersed_energy(solute_radii, rs, prefactors)
    T_pre = max_energy / (n * 10)
    T_pre = T_pre < 0.0 ? -1.0 * T_pre : T_pre

    T_init = SolSim.Simulation.temp_finder(
        temp_finder_iterations,
        solute,
        rs, 
        η,
        σ,
        T_pre,
        target_acceptance_rate
    )

    interaction_range = 2.0 * (maximum(solute_radii) + rs)

    perturb!(centers_cand, solute) = SolSim.Simulation.translate_single_modulo_bounds!(centers_cand, solute, σ, interaction_range, bounds)
    energy(perturb_out, centers) = SolSim.Simulation.solvation_free_energy_for_hard_spheres(perturb_out, centers, solute_radii, rs, η)
    cooling_scheme = length(zig_zag) == 0 ? _create_cooling_scheme(iterations, T_init) : _create_cooling_scheme(zig_zag, iterations, T_init)
    
    algorithm = SolSim.Simulation.SimulatedAnnealing(
        energy,
        perturb!,
        T_init,
        cooling_scheme
        )

    start = now()

    data = SolSim.Simulation.AnnealingSimulationData(
        Float64[],
        Float32[],
        Float32[],
        deepcopy(solute)
    )
    SolSim.Simulation.simulate(
        algorithm,
        solute,
        iterations,
        data;
        file_save_dir = "$(file_save_root)$(n)_$(rs)_$(η)/"
    )    
    
    println("Simulation done for: ", rs, " | ", η, " | in: ", round((now() - start).value/60_000))
end

function _create_cooling_scheme(
    iterations,
    T_init, 
)
    cooling_scheme(T, E, E_min, step) = SolSim.Simulation.quadratic_additive(
        step, iterations, T_init, 0.0) * SolSim.Simulation.non_monotonic_adaptive_cooling_factor(E, E_min)
end

function _create_cooling_scheme(
    zig_zag, 
    iterations,
    T_init
    )

    zig_cool(step, iterations, T_max, T_min) = SolSim.Simulation.quadratic_additive(step, iterations, T_max, T_min)

    cooling_scheme(temp, E, E_min, step) = SolSim.Simulation.zig_zag(
        zig_cool, 
        step, 
        iterations, 
        T_init, 
        0.0,
        zig_zag)
end

function repeated_simulations(;
    indices = 1:10,
    solute_radii = @SVector[1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
    rs = 0.475,
    η = 0.475,
    iterations = 5000000,
    zig_zag = [1.0, 0.5, 0.3],
    temp_finder_iterations = 100000,
    target_acceptance_rate = 0.5,
    file_save_root = "io/output/repeated_simulations/",
    hot_start_file = ""
    )
    for id in indices
        single(;
            solute_radii,
            rs,
            η,
            iterations,
            zig_zag,
            temp_finder_iterations,
            target_acceptance_rate,
            file_save_root = "$(file_save_root)$(id)/",
            hot_start_file
        )
    end
end

function space_scanning(;
    rss = 0.0125:0.025:0.4875,
    ηs = 0.0125:0.025:0.4875,
    solute_radii = @SVector[1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
    iterations = 5000000,
    zig_zag = [1.0, 0.5, 0.3],
    temp_finder_iterations = 100000,
    target_acceptance_rate = 0.5,
    file_save_root = "io/output/space_scanning/",
    hot_start_file = ""
    )
    for rs in rss
        for η in ηs
            single(;
            solute_radii,
            rs,
            η,
            iterations,
            zig_zag,
            temp_finder_iterations,
            target_acceptance_rate,
            file_save_root = "$(file_save_root)",
            hot_start_file
        )
        end
    end
end