# TODO: The call structure is a mess. This should call a version of METROPOLIS algortihm that generates more detailed output. Although one might also arrive at the conclusion 
# that vanilla metropolis should only be called as the special case of simulated_annealing with constant T.
function temp_finder(
    iterations::Int,
    solute::Atoms,
    rs::Float64,
    η::Float64,
    σ::Float64,
    T_pre::Float64, 
    target_acceptance_rate)

    @assert T_pre > 0.0
    Es = generate_transitions(iterations, solute, rs, η, σ, T_pre)

    transitions = []
    for i in 1:length(Es)-1
        if Es[i] > Es[i+1]
            push!(transitions, Es[i])
            push!(transitions, Es[i+1])
        end
    end

    chi_bar(T) = sum([exp(-transitions[i]/T) for i in 1:2:length(transitions)-1])/sum([exp(-transitions[i]/T) for i in 2:2:length(transitions)])
    χ_0 = target_acceptance_rate
    T_0 = T_pre
    try
        while abs(chi_bar(T_0) - χ_0) > 0.00001
            T_0 = T_0 * (log(chi_bar(T_0)) / log(χ_0 ))
        end
    catch 
        println("No energy decreasing transitions found!")
    end
    (isnan(T_0) || T_0 <= 0) ? T_pre : T_0
end

function generate_transitions(
    iterations::Int,
    solute::Atoms,
    rs::Float64,
    η::Float64,
    σ::Float64, #Normal Deviation for perturbation moves
    T_init::Float64
    )
    inflated_radii = SVector(ntuple(i -> solute.radii[i] + rs, size(solute.radii)[1]))

    interaction_range = maximum(inflated_radii) * 2.0
    bounds = interaction_range * size(solute.radii)[1]

    perturb!(x_cand, solute) = translate_single_modulo_bounds!(x_cand, solute, σ, interaction_range, bounds)
    energy(perturb_out, centers) = solvation_free_energy_for_hard_spheres(perturb_out, centers, solute.radii, rs, η)
    cooling_scheme(temp, E, E_min, iteration_step) = temp # Constant temperature for standard random walk metropolis

    algorithm = SimulatedAnnealing(
        energy,
        perturb!,
        T_init,
        cooling_scheme
        )

    data = AnnealingSimulationData(
        Float64[],
        Float32[],
        Float32[],
        deepcopy(solute)
    )

    _ = simulate(
        algorithm, 
        solute,
        iterations,
        data;
        file_save_dir = "io/output/tmp/"
    ) 

    data.Es
end