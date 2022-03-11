# Automatic update plot

# function get_observable_dataset(ds::TimeseriesDataset)
#     var_to_observe = [ds.var, ds.target..., Symbol.(string.(ds.target...).*"_sim")]
#     ts = TimeseriesDataset(ds.data[:,[:t, var_to_observe...]], ds.target, ds.var)
#     return Observable(ts)
# end

# function get_observable_dataset(ds::PonctualDataset, Observe)
#     var_to_observe = [ds.var, ds.target..., Symbol.(string.(ds.target...).*"_sim")]
#     pd = PonctualDataset(ds.data[:,[:t, var_to_observe...]], ds.target)
#     return Observable(pd)
# end

# function get_observable_dataset(ds::Dataset)
#     ds_new = Dataset(TimeseriesDataset[], PonctualDataset[])
#     for tsi in ds.timeseries
#         push!(ds_new.timeseries, get_observable_dataset(tsi))
#     end
#     for pdi in ds.ponctual
#         push!(ds_new.ponctual, get_observable_dataset(pdi))
#     end
#     return Observable(ds_new)
# end

# function update_and_notify_obs!(ts_obs::Observable{TimeseriesDataset},ts::TimeseriesDataset)
#     ts_obs[].data[:, [string.(ts.target...), string.(ts.target...)*"_sim"]] .= ts.data[:, [string.(ts.target...), string.(ts.target...)*"_sim"]]
#     notify(ts_obs) 
# end

# function update_and_notify_obs!(pd_obs::Observable{PonctualDataset},pd::PonctualDataset)
#     ts_obs[].data[:, [string.(ts.target...), string.(ts.target...)*"_sim"]] .= ts.data[:, [string.(ts.target...), string.(ts.target...)*"_sim"]]
#     notify(ts_obs) 
# end

# update_and_notify_obs!(ts_obs::Observable{TimeseriesDataset},ds::Dataset) = update_and_notify_obs!(ts_obs,ds.timeseries[1])

# get_observable_dataset(ds::Union{Dataset,DatasetType}) = Observable(ds)

# function update_and_notify_obs!(ts_obs::Observable{TimeseriesDataset}, ts::TimeseriesDataset,p ; plot_transformed=true)
#     @views ts_obs[].data[:, string.(ts.target)] .= ts.data[:, string.(ts.target)]
#     @views ts_obs[].data[:, string.(ts.target).*"_sim"] .= ts.data[:, string.(ts.target).*"_sim"]
#     plot_transformed && transform_data!(ts_obs[], p)
#     notify(ts_obs)
# end

# function update_and_notify_obs!(pd_obs::Observable{PonctualDataset}, pd::PonctualDataset, p ; plot_transformed=true)
#     @views pd_obs[].data[:, string.(pd.target)] .= pd.data[:, string.(pd.target)]
#     @views pd_obs[].data[:, string.(pd.target).*"_sim"] .= pd.data[:, string.(pd.target).*"_sim"]
#     plot_transformed && transform_data!(pd_obs[], p)
#     notify(pd_obs)
# end

# function update_and_notify_obs!(ds_obs::Observable{Dataset}, ds::Dataset, p ; plot_transformed=true) 
#     for i in 1:length(ds.timeseries)
#         ts_obs = ds_obs[].timeseries[i]
#         ts = ds.timeseries[i]
#         @views ts_obs.data[:, string.(ts.target)] .= ts.data[:, string.(ts.target)]
#         @views ts_obs.data[:, string.(ts.target).*"_sim"] .= ts.data[:, string.(ts.target).*"_sim"]
#     end
#     for i in 1:length(ds.ponctual)
#         pd_obs = ds_obs[].ponctual[i]
#         pd = ds.ponctual[i]
#         @views pd_obs.data[:, string.(pd.target)] .= pd.data[:, string.(pd.target)]
#         @views pd_obs.data[:, string.(pd.target).*"_sim"] .= pd.data[:, string.(pd.target).*"_sim"]
#     end
#     plot_transformed && transform_data!(ds_obs[], p)
#     notify(ds_obs)
# end


function vizualize(ds::Observable{Dataset} ; linewidth=5, markersize=20)
    f = Figure(resolution=(1500,1200))
    max_length = max(length.([ds[].timeseries, ds[].ponctual])...)
    axs = [Axis(f[i, j], fontsize=40) for i in 1:max_length, j in 1:2]
    for i in 1:length(ds[].timeseries)
        ts = @lift $ds.timeseries[i]
        draw_axis(axs[i,1], ts ; linewidth, markersize)
    end
    for i in 1:length(ds[].ponctual)
        pd = @lift $ds.ponctual[i]
        draw_axis(axs[i,2], pd ; linewidth, markersize)
    end
    on(x->autolimits!.(axs),ds)
    display(f)
    return f
end

function vizualize(ds::Observable{<:DatasetType} ; linewidth=5, markersize=20)
    f = Figure(resolution=(1500,1200))
    ax = Axis(f[1, 1], fontsize=40)
    draw_axis(ax, ds ; linewidth, markersize)
    on(x -> autolimits!(ax), ds)
    display(f)
    f, ax
end


function draw_axis(ax, ds::Observable{TimeseriesDataset} ; xvar = nothing, linewidth=5, markersize=20)
    dg = @lift group($ds)
    l = s = nothing
    for i in 1:length(dg[])
        x_vec = isnothing(xvar) ? @lift(abs.($dg[i].t)) :  @lift(abs.($dg[i][!,xvar]))
        target_sim_vec = @lift abs.($dg[i][:, string.($ds.target...)*"_sim"])
        target_vec = @lift abs.($dg[i][:, $ds.target...])

        l = lines!(ax, x_vec, target_sim_vec;
        color = :black, linewidth
        )
    
        s = scatter!(ax, x_vec, target_vec;
            #marker = markers[j],
            markersize,
            color = Cycled(i)
        )
    end
    #xlabel!(ax,x_name)
    #ylabel!(ax, string.(ds[].target)...)
    ax.ylabel = string.(ds[].target)[1]
    return l, s
end

function draw_axis(ax, ds::Observable{PonctualDataset} ; xvar = nothing, group_var=:pc, linewidth=5, markersize=20)
    data = @lift $ds.data
    colnames = names(data[][!,Not(ds[].target)])
    if length(ds[].target) == 1
        x_name = colnames[argmax(length.(unique.(eachcol(data[][!,Not(ds[].target)]))))]
        y_name = string(ds[].target[1])
    else
        x_name = string(ds[].target[1])
        y_name = string(ds[].target[2])
    end
    dg = @lift groupby($data, group_var)
    for i in 1:length(dg[])
        g = @lift $dg[i]
        x_vec = isnothing(xvar) ? @lift(abs.($g[:,x_name])) : @lift(abs.($g[!,xvar]))
        if length(ds[].target) == 1
            x_vec_sim = isnothing(xvar) ? @lift(abs.($g[:,x_name])) : @lift(abs.($g[!,xvar]))
        else
            x_vec_sim = isnothing(xvar) ? @lift(abs.($g[:,x_name*"_sim"])) : @lift(abs.($g[!,xvar]))
        end

        target_sim_vec = @lift abs.($g[:, y_name*"_sim"])
        target_vec = @lift abs.($g[:, y_name])

        l = scatterlines!(ax, x_vec_sim, target_sim_vec;
        color = :black, linewidth, markersize
        )

        s = scatter!(ax, x_vec, target_vec;
            #marker = markers[j],
            markersize,
            #color = Cycled(i)
        )
    end
    ax.xlabel = x_name
    ax.ylabel = y_name
    (x_name ∈ ["ϵ̇", "ϵ̇_dev"]) && (ax.xscale = log10)
    return nothing
end