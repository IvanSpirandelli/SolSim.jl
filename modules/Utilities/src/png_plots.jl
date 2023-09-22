function plot_energy_histogram(energies, th_x_tick; burn_in::Int = 1, bins::Int = 100)
    f = Figure()
    max = maximum(energies)
    ax = Axis(f[1,1], xticks = [th_x_tick, Int(round(max+1))])
    _ = hist!(ax,energies[burn_in:end];bins=bins)
    f
end

function plot_energy_histogram(energies; burn_in::Int = 1, bins::Int = 100)
    f = Figure()
    ax = Axis(f[1,1])
    _ = hist!(ax,energies[burn_in:end];bins=bins)
    f
end

function plot_temperatures(temps)
    f = Figure()
    ax = Axis(f[1,1])
    _ = lines!(ax,temps)
    f
end

function plot_energies(energies)
    f = Figure()
    ax = Axis(f[1,1])
    _ = lines!(ax,energies)
    f
end