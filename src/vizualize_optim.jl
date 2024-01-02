
function vizualize(ds::Observable{Dataset}, p; linewidth=4, markersize=10)
    f = Figure(size=(800, 600))
    nts, npd = length(ds[].timeseries), length(ds[].ponctual)
    nfigrows = max(nts, npd)
    axs = Matrix{Axis}(undef, (nfigrows, 2))
    for (i, tsname) in enumerate(name.(ds[].timeseries))
        name_sym = Symbol(tsname)
        ts = @lift $ds.timeseries[i]
        target = ts[].target
        plot_props = getproperty(p.plot, name_sym)
        axs[i, 1] = Axis(f[i, 1]; plot_props.axis...)
        draw_axis(axs[i, 1], ts, p; linewidth, markersize)
        # if (target == [:ϵ̇c]) | (target == [:ϵ̇])
        #     xlims!(axs[i,1], low = 1e-10)
        #     ylims!(axs[i,1], low = 1e-10)
        #     axs[i,1].xscale = log10
        #     axs[i,1].yscale = log10
        # end
    end
    for (i, pdname) in enumerate(name.(ds[].ponctual))
        name_sym = Symbol(pdname)
        plot_props = getproperty(p.plot, name_sym)
        axs[i, 2] = Axis(f[i, 2]; plot_props.axis...)
        pd = @lift $ds.ponctual[i]
        draw_axis(axs[i, 2], pd, p; linewidth, markersize)
    end
    assigned_axs = [isassigned(axs, i) for i in eachindex(axs)]
    on(x -> autolimits!.(axs[assigned_axs]), ds)
    display(f)
    return f
end

function vizualize(ds::Observable{<:DatasetType}, p; linewidth=4, markersize=10)
    f = Figure(size=(1000, 800))
    name_sym = name(ds[])
    plot_props = getproperty(p.plot, Symbol(name_sym))
    ax = Axis(f[1, 1]; plot_props.axis...)
    draw_axis(ax, ds, p; linewidth, markersize)
    on(x -> autolimits!(ax), ds)
    display(f)
    f, ax
end


function draw_axis(ax, ds::Observable{TimeseriesDataset}, p; linewidth=4, markersize=10)
    dg = @lift group($ds)
    name_sym = Symbol(name(ds[]))
    cols_to_plot = getproperty(p.plot, name_sym).syms
    l = s = nothing
    for i in 1:length(dg[])

        x_vec = @lift abs.($dg[i][:, cols_to_plot[1]])
        y_vec = @lift abs.($dg[i][:, cols_to_plot[2]])
        x_sim_vec = cols_to_plot[1] ∈ ds[].target ? @lift(abs.($dg[i][:, string.(cols_to_plot[1])*"_sim"])) : x_vec
        y_sim_vec = cols_to_plot[2] ∈ ds[].target ? @lift(abs.($dg[i][:, string.(cols_to_plot[2])*"_sim"])) : y_vec

        s = scatter!(ax, x_vec, y_vec;
            #marker = markers[j],
            markersize,
            color=Cycled(i)
        )

        l = lines!(ax, x_sim_vec, y_sim_vec;
            color=:black, linewidth
        )

    end
    #xlabel!(ax,x_name)
    #ylabel!(ax, string.(ds[].target)...)
    ax.xlabel = string(cols_to_plot[1])
    ax.ylabel = string(cols_to_plot[2])
    return l, s
end

function draw_axis(ax, ds::Observable{PonctualDataset}, p; group_var=:pc, linewidth=4, markersize=10)
    data = @lift $ds.data
    #colnames = names(data[][!,Not(ds[].target)])
    name_sym = Symbol(name(ds[]))
    cols_to_plot = getproperty(p.plot, name_sym).syms
    # if length(ds[].target) == 1
    #     x_name = colnames[argmax(length.(unique.(eachcol(data[][!,Not(ds[].target)]))))]
    #     y_name = string(ds[].target[1])
    # else
    #     x_name = string(ds[].target[1])
    #     y_name = string(ds[].target[2])
    # end
    dg = @lift groupby($data, group_var)
    for i in 1:length(dg[])
        g = @lift $dg[i]
        x_vec = @lift abs.($g[:, cols_to_plot[1]])
        y_vec = @lift abs.($g[:, cols_to_plot[2]])
        x_sim_vec = cols_to_plot[1] ∈ ds[].target ? @lift(abs.($g[:, string.(cols_to_plot[1])*"_sim"])) : x_vec
        y_sim_vec = cols_to_plot[2] ∈ ds[].target ? @lift(abs.($g[:, string.(cols_to_plot[2])*"_sim"])) : y_vec

        # x_vec = isnothing(xvar) ? @lift(abs.($g[:,x_name])) : @lift(abs.($g[!,xvar]))
        # if length(ds[].target) == 1
        #     x_vec_sim = isnothing(xvar) ? @lift(abs.($g[:,x_name])) : @lift(abs.($g[!,xvar]))
        # else
        #     x_vec_sim = isnothing(xvar) ? @lift(abs.($g[:,x_name*"_sim"])) : @lift(abs.($g[!,xvar]))
        # end

        # target_sim_vec = @lift abs.($g[:, y_name*"_sim"])
        # target_vec = @lift abs.($g[:, y_name])

        s = scatter!(ax, x_vec, y_vec;
            #marker = markers[j],
            markersize,
            #color = Cycled(i)
        )
        l = scatterlines!(ax, x_sim_vec, y_sim_vec;
            color=:black, linewidth, markersize
        )
    end
    ax.xlabel = string(cols_to_plot[1])
    ax.ylabel = string(cols_to_plot[2])
    #(x_name ∈ ["ϵ̇", "ϵ̇_dev"]) && (ax.xscale = log10)
    return nothing
end