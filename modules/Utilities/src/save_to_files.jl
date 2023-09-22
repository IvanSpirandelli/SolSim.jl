function vec_to_poly(centers::Vector{SVector{3, Float64}}, radii::AbstractVector, filepath::String)
    open(string(filepath, ".poly"), "w") do io
        println(io,"POINTS")
        
        for (i, center) in enumerate(centers)
            println(io, "$(i): $(center[1]) $(center[2]) $(center[3]) $(radii[i])")
        end

        println(io,"POLYS")
        println(io,"END")
    end
end

function vec_to_poly(solute::Atoms, filepath::String)
    vec_to_poly(solute.centers, solute.radii, filepath)
end

function vec_to_ply(centers::Vector{SVector{3,Float64}}, radii::AbstractVector, filepath::String; color = (255,255,255))
    n = length(centers)

    open(string(filepath, ".ply"), "w") do io
        println(io,"ply")
        println(io,"format ascii 1.0")
        println(io,"element vertex " , n)
        println(io,
            "property float32 x\nproperty float32 y\nproperty float32 z"
        )
        println(io,"property float pscale")
        println(io,
            "property uchar red\nproperty uchar green\nproperty uchar blue"
        )
        println(io,"end_header")
        for (i,center) in enumerate(centers)
            println(io, "$(center[1]) $(center[2]) $(center[3]) $(radii[i]) $(color[1]) $(color[2]) $(color[3])")
        end
    end
end

function vec_to_ply(centers::Vector{Vector{SVector{3,Float64}}}, radii::AbstractVector, filepath::String; color = (255,255,255))
    vec_to_ply(reduce(vcat, centers), reduce(vcat, [radii for _ in 1:length(centers)]), filepath)
end

function vec_to_ply(centers::Vector{Vector{SVector{3,Float64}}}, radii::Vector{AbstractVector}, filepath::String; color = (255,255,255))
    vec_to_ply(reduce(vcat, centers), reduce(vcat, radii), filepath)
end

function vec_to_ply(solute::Atoms, filepath::String)
    vec_to_ply(solute.centers, solute.radii, filepath)
end