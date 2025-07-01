module AlphaMolWrap
    using CxxWrap
    using AlphaMolWrapper_jll

    @wrapmodule(() -> libalphamolwrapper)

    version() = VersionNumber(unsafe_string(unsafe_load(cglobal((:__gmp_version, :libgmp), Ptr{Cchar}))))
    bits_per_limb() = Int(unsafe_load(cglobal((:__gmp_bits_per_limb, :libgmp), Cint)))
    const VERSION = version()
    const BITS_PER_LIMB = bits_per_limb()

    function __init__()
        @initcxx
        if version() != VERSION || bits_per_limb() != BITS_PER_LIMB
            msg = """The dynamically loaded GMP library (v\"$(version())\" with __gmp_bits_per_limb == $(bits_per_limb()))
                    does not correspond to the compile time version (v\"$VERSION\" with __gmp_bits_per_limb == $BITS_PER_LIMB).
                    Please rebuild Julia."""
            bits_per_limb() != BITS_PER_LIMB ? @error(msg) : @warn(msg)
        end
    end
end

function get_measures(
    atom_coordinates::Vector, 
    atom_radii::Vector, 
    probe_radius::Float64,
    delaunay_eps::Float64 = 1.0)
    outs = [0.0, 0.0, 0.0, 0.0]
    AlphaMolWrap.get_geometric_measures(
        outs,
        collect(Base.Iterators.flatten(atom_coordinates)),
        atom_radii,
        probe_radius,
        delaunay_eps
    )
    outs
end

function get_geometric_measures_with_derivatives(
    atom_coordinates::Vector, 
    atom_radii::Vector, 
    probe_radius::Float64, 
    delaunay_eps::Float64 = 1.0)
    
    measure_outs = [0.0, 0.0, 0.0, 0.0, 0.0]

    n = size(atom_coordinates)[1]

    dvol_outs = [0.0 for _ in 1:n]
    dsurf_outs = [0.0 for _ in 1:n]
    dmean_outs = [0.0 for _ in 1:n]
    dgauss_outs = [0.0 for _ in 1:n]

    AlphaMolWrap.get_geometric_measures_with_derivatives(
        measure_outs,
        dvol_outs,
        dsurf_outs, 
        dmean_outs,
        dgauss_outs,
        collect(Base.Iterators.flatten(atom_coordinates)),
        atom_radii,
        probe_radius,
        delaunay_eps
    )
    measure_outs, dvol_outs, dsurf_outs, dmean_outs, dgauss_outs
end
