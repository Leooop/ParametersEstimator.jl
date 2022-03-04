

export optimize!, get_jacobian, get_posterior_model_covariance, get_logspace_model_and_covmat

function get_jacobian(g::F, mₙ ; method=:finitediff) where F<:Function
    if method == :finitediff
        J = FiniteDiff.finite_difference_jacobian(g, mₙ)
    elseif method == :autodiff
        J = ForwardDiff.jacobian(g,mₙ)
    else
        @error "`method` kwarg must be :finitediff"
    end
    return J
end

function get_new_parameters!(g::F, mₙₑₓₜ, mₙ, mₚ, dₙ, dₒ, Cm, Cd ; μₙ=1.0, diff_method=:finitediff) where F<:Function
    G = get_jacobian(g, mₙ ; method=diff_method) 
    Cm_inv, Cd_inv = inv(Cm), inv(Cd)
    A = G'*Cd_inv*G + Cm_inv
    b = G'*Cd_inv*(dₙ - dₒ) + Cm_inv*(mₙ - mₚ)
    mₙₑₓₜ .= A \ (A*mₙ - μₙ*b)
    return mₙₑₓₜ
end

function get_new_parameters!(g::F, mₙₑₓₜ::NamedArray, mₙ::NamedArray, mₚ::NamedArray, dₙ, dₒ, Cm, Cd ; μₙ=1.0, diff_method=:finitediff) where F<:Function
    G = get_jacobian(g, mₙ ; method=diff_method)
    Cm_inv, Cd_inv = inv(Cm), inv(Cd)
    A = G'*Cd_inv*G + Cm_inv
    b = G'*Cd_inv*(dₙ - dₒ) + Cm_inv*(mₙ - mₚ).array
    mₙₑₓₜ .= A \ (A*mₙ - μₙ*b)
    return mₙₑₓₜ
end

function optimize!(dataset, p, ; 
        μₙ=1.0, reltol=1e-8, 
        maxiters=1000, 
        diff_method = :finitediff,
        allow_mismatch=false)
    reldiff = Inf
    mₙₑₓₜ, Cm = get_logspace_model_and_covmat(p)
    mₙ = copy(mₙₑₓₜ)
    mₚ = copy(mₙₑₓₜ)
    n = 0

    simulate!(dataset, p ; allow_mismatch)
    df_target = get_target_data(dataset)
    transform_data!(df_target,p)
    dₙ = df_target.target_sim
    dₒ = df_target.target
    Cd = get_data_covmat(df_target)

    function g(m)
        p.mp .= exp10.(m)
        simulate!(dataset, p ; allow_mismatch)
        df_target = get_target_data(dataset)
        transform_data!(df_target, p)
        return Vector{Float64}(df_target[!,:target_sim])
    end

    #initialize plot
    ts_obs = get_observable_timeseries(dataset)
    ts_obs = (ts_obs isa Vector) ? ts_obs[1] : ts_obs
    vizualize(ts_obs)
 

    while (reldiff > reltol) & (n < maxiters)
        n+=1 ; @show(n)
        #@show (mₙₑₓₜ.-mₚ)./mₚ
        mₙ .= mₙₑₓₜ
        #try
            get_new_parameters!(g, mₙₑₓₜ, mₙ, mₚ, dₙ, dₒ, Cm, Cd ; μₙ, diff_method)
            p.mp .= exp10.(mₙₑₓₜ)
            simulate!(dataset, p)
            df_target_next = get_target_data(dataset)
            transform_data!(df_target_next,p)
            dₙ = df_target_next.target_sim
            dₒ = df_target_next.target
            #Cd = get_data_covmat(df_target) # not needed if data always has the same length
            #reldiff = maximum(abs.((mₙₑₓₜ .- mₙ)./mₙ))
            reldiff = sum((dₙ.-dₒ).^2)
        # catch e
        #     @show mₙ mₙₑₓₜ
        #     @show maximum(abs,get_jacobian(g, mₙ ; method=:finitediff) )
        #     @show extrema(g(mₙ))
        #     throw(e)
        # end
        @show reldiff

        #plot
        update_and_notify_obs!(ts_obs,dataset)

    end

    return dataset, p
end


function get_posterior_model_covariance(g::F, m, Cm, Cd) where F<:Function
    G = get_jacobian(g, m ; method=:finitediff)
    println("G done")
    Cm*G'
    println("Cm*G' done")
    G*Cm*G' + Cd
    println("G*Cm*G' + Cd done")
    inv(G*Cm*G' + Cd)
    println("inv(G*Cm*G' + Cd) done")
    G*Cm
    println("G*Cm done")
    return Cm - Cm*G' * inv(G*Cm*G' + Cd) * G*Cm
end

function transform_data!(df_target::AbstractDataFrame, p)
    dfg = groupby(df_target,:target_name)
    if hasproperty(p.estimate,:transform_target)
        for i in 1:length(dfg)
            target_name = dfg[i][1,:target_name]
            if hasproperty(p.estimate.transform_target,target_name)
                for col in [:target,:target_sim]
                    dfg[i][!,col] = getproperty(p.estimate.transform_target, target_name)(dfg[i][!,col])
                end
            end
        end
    end
    return df_target
end

get_logspace_model_and_covmat(p) = get_logspace_model_and_covmat(p.mp, p.mp_std ; max_log_std = p.estimate.log_mp_std_max)

function get_logspace_model_and_covmat(mp, std_mp ; max_log_std=4) #TODO : change kwarg to max_log_var
    log_lb = log10.(mp.array .- min.(std_mp, mp.array .- 1e-12))
    log_ub = log10.(mp.array .+ std_mp)
    Δlog_var = min.((log_ub.-log_lb).^2, max_log_std)
    Cm_log = Diagonal(Δlog_var)
    return log10.(mp.array), Cm_log
end
# Automatic update plot

function get_observable_timeseries(ds::TimeseriesDataset)
    var_to_observe = [ds.var, ds.target..., Symbol.(string.(ds.target...).*"_sim")]
    ts = TimeseriesDataset(ds.data[:,[:t, var_to_observe...]], ds.target, ds.var)
    return Observable(ts)
end
get_observable_timeseries(ds::Dataset) = get_observable_timeseries.(ds.timeseries)


function update_and_notify_obs!(ts_obs::Observable{TimeseriesDataset},ts::TimeseriesDataset)
    ts_obs[].data[:, [string.(ts.target...), string.(ts.target...)*"_sim"]] .= ts.data[:, [string.(ts.target...), string.(ts.target...)*"_sim"]]
    notify(ts_obs) 
end
update_and_notify_obs!(ts_obs::Observable{TimeseriesDataset},ds::Dataset) = update_and_notify_obs!(ts_obs,ds.timeseries[1])

function vizualize(ds::Observable{TimeseriesDataset} )
    #ds = ds_obs[]
    dg = @lift group($ds)
    f = Figure(resolution=(1500,1200))
    ax = Axis(f[1, 1], fontsize=40)
    for i in 1:length(dg[])
            t_vec = @lift $dg[i].t
            target_sim_vec = @lift -$dg[i][:, string.($ds.target...)*"_sim"]
            target_vec = @lift -$dg[i][:, $ds.target...]

            lines!(ax, t_vec, target_sim_vec,
            color = :black, linewidth = 5
            )
        
            scatter!(ax, t_vec, target_vec,
                #marker = markers[j],
                markersize = 20,
                color = Cycled(i)
            )
    end
    display(f)
    f, ax
end