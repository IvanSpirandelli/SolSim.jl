function get_normal_distributed_translation(σ::Float64)
    return SVector{3, Float64}(randn(3)) * σ
end

#Atoms
function translate_single!(
    centers_cand::Vector{SVector{3, Float64}}, 
    solute::Atoms,
    σ::Float64,
    )
    i = rand(1:length(centers_cand))
    translation = get_normal_distributed_translation(σ)
    @inbounds centers_cand[i] = solute.centers[i] + translation
    return i
end

function translate_single_modulo_bounds!(
    centers_cand::Vector{SVector{3, Float64}}, 
    solute::Atoms,
    sigma::Float64, 
    interaction_range::Float64, 
    bounds::Float64
    )
    i = rand(1:length(solute.centers))
    translation = get_normal_distributed_translation(sigma)
    @inbounds centers_cand[i] = mod.(solute.centers[i] + translation,  bounds)
    shift = get_shift_vector(centers_cand, interaction_range, bounds)
    if any(shift.>0.0) shift_centers!(centers_cand, shift, bounds) end
    return i
end


function get_shift_vector(
    centers::Vector{SVector{3, Float64}}, 
    interaction_range::Float64, 
    bounds::Float64
    )
    xs = Vector{Float64}([])
    ys = Vector{Float64}([])
    zs = Vector{Float64}([])

    shift = [0.0,0.0,0.0]

    for center in centers 
        push!(xs, center[1])
        push!(ys, center[2])
        push!(zs, center[3])
    end

    sort!(xs)
    sort!(ys)
    sort!(zs)

    if(dist_over_bound(xs[1], xs[end], bounds) < interaction_range)
        shift[1] = get_one_dimensional_shift(xs, interaction_range)
    end

    if(dist_over_bound(ys[1], ys[end], bounds) < interaction_range)
        shift[2] = get_one_dimensional_shift(ys, interaction_range)
    end

    if(dist_over_bound(zs[1], zs[end], bounds) < interaction_range)
        shift[3] = get_one_dimensional_shift(zs, interaction_range)
    end

    shift
end


function shift_centers!(centers, shift, bounds)
    for i in 1:length(centers)
        centers[i] = mod.(centers[i] - shift,  bounds)
    end
end


function dist_over_bound(x, y, b)
    lower = minimum([x,y])
    upper = maximum([x,y])
    lower + (b - upper)
end


function get_one_dimensional_shift(values, interaction_range)
    n = length(values)
    for i in 1:n-1
        if values[i+1] - values[i] > interaction_range
            return (values[i+1] + values[i]) / 2.0
        end
    end
    -1 # Can only happen if input is illegal
end