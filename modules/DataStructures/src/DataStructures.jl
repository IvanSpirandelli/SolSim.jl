module DataStructures

export Atoms

using StaticArrays

mutable struct Atoms{N}
    centers::Vector{SVector{3, Float64}}
    radii::SVector{N, Float64}
end

Base.deepcopy(atoms::Atoms) = Atoms{length(atoms.centers)}(deepcopy(atoms.centers), deepcopy(atoms.radii))

end #module DataStructures