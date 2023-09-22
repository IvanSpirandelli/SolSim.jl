struct SAModifiedLam{E,P}
    energy::E
    perturb!::P
end

mutable struct ModifiedLamSimulationData
    Es::Vector{Float64}
    Ts::Vector{Float16}
    αs::Vector{Float16}
    λs::Vector{Float16}
end

function simulate(
    algorithm::SAModifiedLam, 
    solute::Atoms,
    iterations::Int,
    )

    T = 0.5
    α = 0.5 #Acceptance Rate
    λ = 0.44 #Lam Rate
    
    E = algorithm.energy(nothing, solute.centers)
    
    centers_cand = copy(solute.centers)
    accepted_steps = 0

    for i in 1:iterations
        centers_cand = copy(solute.centers)

        perturb_out = algorithm.perturb!(centers_cand, solute)
        E_cand = algorithm.energy(perturb_out, centers_cand)
        
        accepted = rand() <  min(1, exp(-(E_cand - E)/T))
    
        if accepted
            accepted_steps += 1
            solute.centers, centers_cand = centers_cand, solute.centers
            E = E_cand
            α = 1.0/500.0 * (499 * α + 1)
        else
            α = 1.0/500.0 * (499 * α)
        end

        # Modfied Lam Calculations  
        if i/iterations < 0.15 
            λ = 0.44 + 0.56 * 560 ^ (-i/iterations/0.15)
        elseif 0.15 ≤ i/iterations < 0.65
            λ = 0.44
        else 
            λ = 0.44 * 440^(-(i/iterations - 0.65)/0.35)  
        end

        if α > λ
            T = 0.999 * T
        else
            T = T/0.999
        end
    end
    
    return solute, EndOfSimulationInfo(E, accepted_steps / iterations)
end

function simulate(
    algorithm::SAModifiedLam, 
    solute::Atoms,
    iterations::Int,
    info::ModifiedLamSimulationData
    )

    T = 0.5
    α = 0.5 #Acceptance Rate
    λ = 0.44 #Lam Rate
    
    E = algorithm.energy(nothing, solute.centers)
    
    add_data_point!(info, E, T, α, λ)

    centers_cand = copy(solute.centers)
    accepted_steps = 0

    for i in 1:iterations
        centers_cand = copy(solute.centers)

        perturb_out = algorithm.perturb!(centers_cand, solute)
        E_cand = algorithm.energy(perturb_out, centers_cand)
        
        accepted = rand() <  min(1, exp(-(E_cand - E)/T))
    
        if accepted
            accepted_steps += 1
            solute.centers, centers_cand = centers_cand, solute.centers
            E = E_cand
            α = 1.0/500.0 * (499 * α + 1)
        else
            α = 1.0/500.0 * (499 * α)
        end

        # Modfied Lam Calculations  
        if i/iterations < 0.15 
            λ = 0.44 + 0.56 * 560 ^ (-i/iterations/0.15)
        elseif 0.15 ≤ i/iterations < 0.65
            λ = 0.44
        else 
            λ = 0.44 * 440^(-(i/iterations - 0.65)/0.35)  
        end

        if α > λ
            T = 0.999 * T
        else
            T = T/0.999
        end


        add_data_point!(info, E, T, α, λ)
    end
    
    return solute, EndOfSimulationInfo(E, accepted_steps / iterations)
end

function add_data_point!(data::ModifiedLamSimulationData, E, T, α, λ)
    push!(data.Es, E)
    push!(data.Ts, T)
    push!(data.αs, α)
    push!(data.λs, λ)
end
    