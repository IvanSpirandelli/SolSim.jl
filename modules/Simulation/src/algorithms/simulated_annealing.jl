struct SimulatedAnnealing{E,P,C}
    energy::E
    perturb!::P
    T_init::Float64
    cooling_scheme::C # Setting this to a constant temperature recovers standard Metropolis algorithm.
end

struct OutputParameters
    file_save_dir::String
end

# No detailed outputs. Therefor faster.
function simulate(
    simulated_annealing::SimulatedAnnealing, 
    solute::Atoms,
    iterations::Int
    )

    T = simulated_annealing.T_init
    E = simulated_annealing.energy(nothing, solute.centers)
    
    centers_cand = copy(solute.centers)
    accepted_steps = 0

    for iteration_step in 1:iterations
        centers_cand = copy(solute.centers)

        perturb_out = simulated_annealing.perturb!(centers_cand, solute)
        E_cand = simulated_annealing.energy(perturb_out, centers_cand)
        
        accepted = rand() <  min(1, exp(-(E_cand - E)/T))
    
        if accepted
            accepted_steps += 1
            solute.centers, centers_cand = centers_cand, solute.centers
            E = E_cand
        end

        T = simulated_annealing.cooling_scheme(T, E, results.E_min, iteration_step)
    end
    
    return solute, EndOfSimulationInfo(E, accepted_steps / iterations)
end

mutable struct AnnealingSimulationData{D}
    Es::Vector{Float64}
    Ts::Vector{Float32}
    αs::Vector{Float32}
    solute_min::D
end

# Detailed outputs. Therefor slower.
function simulate(
    simulated_annealing::SimulatedAnnealing, 
    solute::Atoms,
    iterations::Int,
    data::AnnealingSimulationData;
    file_save_dir::String = "output/last_run/",
    save_interval = Inf
    )

    T = simulated_annealing.T_init
    E = simulated_annealing.energy(nothing, solute.centers)

    initialize_output_folder(solute, file_save_dir)

    centers_cand = copy(solute.centers)
    accepted_steps = 0  
    E_min = E

    for iteration_step in 1:iterations
        centers_cand = copy(solute.centers)

        perturb_out = simulated_annealing.perturb!(centers_cand, solute)
        E_cand = simulated_annealing.energy(perturb_out, centers_cand)
        
        accepted = rand() <  min(1, exp(-(E_cand - E)/T))
    
        if accepted
            accepted_steps += 1
            solute.centers, centers_cand = centers_cand, solute.centers
            E = E_cand

            if E <= E_min
                E_min = E
                update_minimal_solute!(data, solute)
            end
            
            if save_interval != Inf && mod(save_interval, iteration_step) == 0
                vec_to_poly(solute.centers, solute.radii, string(file_save_dir,"saved/",string(accepted_steps))) 
            end

            add_data_point!(data, accepted_steps, iteration_step, E, T)
        end

        T = simulated_annealing.cooling_scheme(T, E, E_min, iteration_step)
    end
    
    vec_to_poly(data.solute_min.centers, data.solute_min.radii, string(file_save_dir,"minimal"))
    return solute, EndOfSimulationInfo(E, accepted_steps / iterations)
end

function add_data_point!(
    data::AnnealingSimulationData,
    accepted_steps::Int,
    iteration_step::Int,
    E::Float64,
    T::Float64,
    )
    push!(data.Ts, T)
    push!(data.Es, E)
    push!(data.αs, accepted_steps/iteration_step)
end

function update_minimal_solute!(
    data::AnnealingSimulationData,
    solute::Atoms,
    )
    data.solute_min.centers = solute.centers
    data.solute_min.radii = solute.radii
end

function initialize_output_folder(
    solute::Atoms,
    file_save_dir::String
    )
    try
        rm(file_save_dir, recursive=true)
    catch e 
        println(e)
    end

    mkpath(string(file_save_dir,"saved/"))
    vec_to_poly(solute.centers, solute.radii, string(file_save_dir, "initial"))
end
