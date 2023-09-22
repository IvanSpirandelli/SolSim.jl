function plot_simulation_info(info::ModifiedLamSimulationData, path::String)
    f = Figure(fontsize = 40)
    c = cgrad(:Dark2_4, 4, categorical = true)

    ga = f[1, 1] = GridLayout()

    ax, _ = scatter(ga[1,1], 1:length(info.Es), info.Es, color = c[1],
        axis = (
            width = 800, 
            height = 400, 
            title = ""
            )
        )

    Label(ga[1, 1, Top()], "E", valign = :bottom,
    font = :bold,
    padding = (0, 0, 5, 0))

    gb = f[1, 2] = GridLayout()
    ax, _ = scatter(gb[1,1], 1:length(info.Ts), info.Ts, color = c[2],
        axis = (
            width = 800, 
            height = 400, 
            title = ""
            )
        )

    Label(gb[1, 1, Top()], "T", valign = :bottom,
    font = :bold,
    padding = (0, 0, 5, 0))

    gc = f[2, 1] = GridLayout()
    ax, _ = scatter(gc[1,1], 1:length(info.αs), info.αs, color = c[3],
        axis = (
            width = 800, 
            height = 400, 
            title = ""
            )
        )

    Label(gc[1, 1, Top()], "α", valign = :bottom,
    font = :bold,
    padding = (0, 0, 5, 0))

    gd = f[2, 2] = GridLayout()
    ax, _ = scatter(gd[1,1], 1:length(info.λs), info.λs, color = c[4],
        axis = (
            width = 800, 
            height = 400, 
            title = ""
            )
        )

    Label(gd[1, 1, Top()], "λ", valign = :bottom,
    font = :bold,
    padding = (0, 0, 5, 0))

    resize_to_layout!(f)
    save(string("$(path)detailed_infos.png"), f, px_per_unit = 1)
end