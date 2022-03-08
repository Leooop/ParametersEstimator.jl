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

# function update_and_notify_obs!(ts_obs::Observable{TimeseriesDataset},ts::TimeseriesDataset)
#     @views ts_obs[].data[:, string.(ts.target).*"_sim"] .= ts.data[:, string.(ts.target).*"_sim"]
#     notify(ts_obs)
# end

# function update_and_notify_obs!(pd_obs::Observable{PonctualDataset},pd::PonctualDataset)
#     @views pd_obs[].data[:, string.(pd.target).*"_sim"] .= pd.data[:, string.(pd.target).*"_sim"]
#     notify(pd_obs)
# end

# function update_and_notify_obs!(ds_obs::Observable{Dataset},ds::Dataset) 
#     for i in 1:length(ds.timeseries)
#         ts_obs = ds_obs[].timeseries[i]
#         ts = ds.timeseries[i]
#         @views ts_obs.data[:, string.(ts.target).*"_sim"] .= ts.data[:, string.(ts.target).*"_sim"]
#     end
#     for i in 1:length(ds.ponctual)
#         pd_obs = ds_obs[].ponctual[i]
#         pd = ds.ponctual[i]
#         @views pd_obs.data[:, string.(pd.target).*"_sim"] .= pd.data[:, string.(pd.target).*"_sim"]
#     end
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
    display(f)
    return f
end

function vizualize(ds::Observable{<:DatasetType} ; linewidth=5, markersize=20)
    f = Figure(resolution=(1500,1200))
    ax = Axis(f[1, 1], fontsize=40)
    draw_axis(ax, ds ; linewidth, markersize)
    display(f)
    f, ax
end


function draw_axis(ax, ds::Observable{TimeseriesDataset} ; xvar = nothing, linewidth=5, markersize=20)
    dg = @lift group($ds)
    l = s = nothing
    for i in 1:length(dg[])
        x_vec = isnothing(xvar) ? @lift($dg[i].t) :  @lift($dg[i][!,xvar])
        target_sim_vec = @lift -$dg[i][:, string.($ds.target...)*"_sim"]
        target_vec = @lift -$dg[i][:, $ds.target...]

        l = lines!(ax, x_vec, target_sim_vec;
        color = :black, linewidth
        )
    
        s = scatter!(ax, x_vec, target_vec;
            #marker = markers[j],
            markersize,
            color = Cycled(i)
        )
    end
    return l, s
end

function draw_axis(ax, ds::Observable{PonctualDataset} ; xvar = nothing, group_var=:pc, linewidth=5, markersize=20)
    data = @lift $ds.data
    colnames = names(data[][!,Not(ds[].target)])
    most_probable_x_name = colnames[argmax(length.(unique.(eachcol(data[][!,Not(ds[].target)]))))]
    dg = @lift groupby($data, group_var)
    for i in 1:length(dg[])
        g = @lift $dg[i]
        x_vec = isnothing(xvar) ? @lift($g[:,most_probable_x_name]) : @lift($g[!,xvar])
        target_sim_vec = @lift -$g[:, string.(ds[].target...)*"_sim"]
        target_vec = @lift -$g[:, ds[].target...]

        l = lines!(ax, x_vec, target_sim_vec;
        color = :black, linewidth
        )

        s = scatterlines!(ax, x_vec, target_vec;
            #marker = markers[j],
            markersize,
            linewidth
            #color = Cycled(i)
        )
    end
    return nothing
end