function zig_zag(zig_cool, step::Int, iterations::Int, T_init::Float64, T_min::Float64, level::Vector{Float64})
    iteration_level = [Int(round(i * iterations / length(level))) for i in 0:length(level)]
    current = findfirst(x -> x >= step, iteration_level)
    zig_cool(step-iteration_level[current-1], Int(round(iterations / length(level))), T_init * level[current-1], T_min)
end

# http://what-when-how.com/artificial-intelligence/a-comparison-of-cooling-schedules-for-simulated-annealing-artificial-intelligence/
function linear_additive(step::Int, iterations::Int, T_init::Float64, T_min::Float64)
    T_min + (T_init - T_min)*(iterations - step)/iterations
end

function exponential_additive(step::Int, iterations::Int, T_init::Float64, T_min::Float64)
    @assert false "faulty for T < 1. Think how to solve"
    T_min + ((T_init - T_min)/(1+exp(2.0*log(T_init - T_min)*(step-(0.5*iterations))/iterations)))
end

function quadratic_additive(step::Int, iterations::Int, T_init::Float64, T_min::Float64)
    T_min + (T_init - T_min)*((iterations - step)/iterations)^2
end

function geometric_cooling(step::Int, iterations::Int, T_init::Float64, T_min::Float64; α = 0.95)
    return T_init * α^step
end

function non_monotonic_adaptive_cooling_factor(E::Float64, E_min::Float64)
    (1+(E-E_min)/E)
end

# https://www.cicirello.org/publications/CP2007-Autonomous-Search-Workshop.pdf
function modified_lam_annealing_schedule()
    @assert false # TODO: Implement
    return
end

