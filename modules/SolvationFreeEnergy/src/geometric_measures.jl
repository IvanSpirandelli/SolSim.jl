module AlphaMolWrap
    using CxxWrap
    @wrapmodule(() -> joinpath("lib/", "libAlphaMol.so"))

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

function get_measures(coordinates::Vector{SVector{3, Float64}}, radii::Vector)
    outs = [0.0, 0.0, 0.0, 0.0]
    AlphaMolWrap.calculate_measures(
        outs,
        Base.Iterators.collect(Base.Iterators.flatten(coordinates)),
        radii,
        1.0, 1.0, 1.0, 1.0, 1, 0
    )
    outs
end

function get_measures(coordinates::Vector{SVector{3, Float64}}, radii::SVector)
    flat = Base.Iterators.collect(Base.Iterators.flatten(coordinates))
    conv_r = convert(Vector{Float64}, radii)
    outs = [0.0, 0.0, 0.0, 0.0]

    AlphaMolWrap.calculate_measures(
        outs,
        flat,
        conv_r,
        1.0, 1.0, 1.0, 1.0, 1, 0
    )
    outs
end

function get_measures_and_derivatives(coordinates::Vector{SVector{3, Float64}}, radii::SVector)
    flat = Base.Iterators.collect(Base.Iterators.flatten(coordinates))
    conv_r = convert(Vector{Float64}, radii)

    outs = [0.0, 0.0, 0.0, 0.0]

    n = size(coordinates)[1] * 3

    dvol_outs = [0.0 for _ in 1:n]
    dsurf_outs = [0.0 for _ in 1:n]
    dmean_outs = [0.0 for _ in 1:n]
    dgauss_outs = [0.0 for _ in 1:n]

    AlphaMolWrap.calculate_measures_and_derivatives(
        outs,
        dvol_outs,
        dsurf_outs,
        dmean_outs,
        dgauss_outs,
        flat,
        conv_r,
        1.0, 1.0, 1.0, 1.0, 1, 0
    )
    outs, dvol_outs, dsurf_outs, dmean_outs, dgauss_outs
end