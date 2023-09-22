include("hard_sphere_simulated_annealing_calls.jl")

attempt_100_helices() = repeated_simulations(
    indices = 1:10,
    solute_radii = @SVector[1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
    rs = 0.475,
    η = 0.475,
    iterations = 25000000,
    zig_zag = [1.0, 0.8, 0.6, 0.4, 0.2],
    temp_finder_iterations = 1000000,
    target_acceptance_rate = 0.5,
    file_save_root = "io/output/repeated_simulations/",
    hot_start_file = ""
)

four_hardsphere_scan() = space_scanning(;
    rss = 0.0125:0.025:0.4875,
    ηs = 0.0125:0.025:0.4875,
    solute_radii = @SVector[1.0, 1.0, 1.0, 1.0],
    iterations = 25000000,
    zig_zag = [1.0, 0.8, 0.6, 0.4, 0.2],
    temp_finder_iterations = 1000000,
    target_acceptance_rate = 0.5,
    file_save_root = "io/output/space_scanning/",
    hot_start_file = ""
)

five_hardsphere_scan() = space_scanning(;
    rss = 0.0125:0.025:0.4875,
    ηs = 0.0125:0.025:0.4875,
    solute_radii = @SVector[1.0, 1.0, 1.0, 1.0, 1.0],
    iterations = 25000000,
    zig_zag = [1.0, 0.8, 0.6, 0.4, 0.2],
    temp_finder_iterations = 1000000,
    target_acceptance_rate = 0.5,
    file_save_root = "io/output/space_scanning/",
    hot_start_file = ""
)

six_hardsphere_scan() = space_scanning(;
    rss = 0.0125:0.025:0.4875,
    ηs = 0.0125:0.025:0.4875,
    solute_radii = @SVector[1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
    iterations = 25000000,
    zig_zag = [1.0, 0.8, 0.6, 0.4, 0.2],
    temp_finder_iterations = 1000000,
    target_acceptance_rate = 0.5,
    file_save_root = "io/output/space_scanning/",
    hot_start_file = ""
)

seven_hardsphere_scan() = space_scanning(;
    rss = 0.0125:0.025:0.4875,
    ηs = 0.0125:0.025:0.4875,
    solute_radii = @SVector[1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
    iterations = 25000000,
    zig_zag = [1.0, 0.8, 0.6, 0.4, 0.2],
    temp_finder_iterations = 1000000,
    target_acceptance_rate = 0.5,
    file_save_root = "io/output/space_scanning/",
    hot_start_file = ""
)

eight_hardsphere_scan() = space_scanning(;
    rss = 0.0125:0.025:0.4875,
    ηs = 0.0125:0.025:0.4875,
    solute_radii = @SVector[1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
    iterations = 25000000,
    zig_zag = [1.0, 0.8, 0.6, 0.4, 0.2],
    temp_finder_iterations = 1000000,
    target_acceptance_rate = 0.5,
    file_save_root = "io/output/space_scanning/",
    hot_start_file = ""
)