struct Metropolis{E,P}
    energy::E
    perturb!::P
    β::Float64 # inverse temperature
end

struct MetropolisParameters
    temperature::Float64
    σ::Float64
end

function simulate(metropolis::Metropolis, solute::Atoms, iterations::Int)
    energy = metropolis.energy
    perturb! = metropolis.perturb!
    β = metropolis.β

    centers_cand = copy(solute.centers)
    E = energy(nothing, solute.centers)
    accepted_steps = 0

    for _ in 1:iterations
        centers_cand = copy(solute.centers)
        perturb_out = perturb!(centers_cand, solute)
        E_cand = energy(perturb_out, centers_cand)

        if rand() < min(1, exp(-β*(E_cand - E))) # TODO: drop the min?
            # accept step
            solute.centers, centers_cand = centers_cand, solute.centers
            E = E_cand
            accepted_steps += 1
        end
    end
    
    return solute, EndOfSimulationInfo(E, accepted_steps / iterations)
end