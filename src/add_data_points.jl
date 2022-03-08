function add_data_points!(df_target, ds::Dataset, p)
    for tsi in ds.timeseries
        add_data_points!(df_target,tsi, p)
    end
    for pdi in ds.ponctual
        add_data_points!(df_target,pdi, p)
    end
    return df_target
end

function add_data_points!(df_target, ts::TimeseriesDataset, p)
    if hasproperty(p.estimate,:extra_misfit_target)
        for extra_misfit_target_name in propertynames(p.estimate.extra_misfit_target)
            if extra_misfit_target_name ∈ ts.target
                dfg = group(ts)
                for g in dfg
                    extra_data_props = getproperty(p.estimate.extra_misfit_target,extra_misfit_target_name)
                    for i in eachindex(extra_data_props)
                        func = extra_data_props[i].f
                        colname = extra_misfit_target_name
                        val_target = func(g[!,colname])
                        val_target_sim = func(g[!,string(colname)*"_sim"])
                        for iv in eachindex(val_target)
                            row = (
                                target = val_target[iv],
                                target_sim = val_target_sim[iv],
                                std = extra_data_props[i].std,
                                target_name = extra_data_props[i].transform_target_name
                            )
                            push!(df_target,row)
                        end
                    end
                end
            end
        end
    end
    return df_target
end

function add_data_points!(df_target, pd::PonctualDataset, p)
    data = pd.data
    if hasproperty(p.estimate,:extra_misfit_target)
        for extra_misfit_target_name in propertynames(p.estimate.extra_misfit_target)
            if extra_misfit_target_name ∈ pd.target
                extra_data_props = getproperty(p.estimate.extra_misfit_target,extra_misfit_target_name)
                for i in eachindex(extra_data_props)
                    func = extra_data_props[i].f
                    colname = extra_misfit_target_name
                    val_target = func(data[!,colname])
                    val_target_sim = func(data[!,string(colname)*"_sim"])
                    for iv in eachindex(val_target)
                        row = (
                            target = val_target[iv],
                            target_sim = val_target_sim[iv],
                            std = extra_data_props[i].std,
                            target_name = extra_data_props[i].transform_target_name
                        )
                        push!(df_target,row)
                    end
                end
            end
        end
    end
    return df_target
end