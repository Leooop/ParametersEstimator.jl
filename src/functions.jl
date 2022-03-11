
include("add_data_points.jl")
include("vizualize_optim.jl")
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
    G = get_jacobian(g, mₙ.array ; method=diff_method)
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
        allow_mismatch=false,
        plot_transformed=false)
    reldiff = Inf
    mₙₑₓₜ, Cm = get_logspace_model_and_covmat(p)
    mₙ = copy(mₙₑₓₜ)
    mₚ = copy(mₙₑₓₜ)
    n = 0
    
    simulate!(dataset, p ; allow_mismatch)
    df_target = get_target_data(dataset)
    add_data_points!(df_target, dataset, p)
    Cd = get_transformed_covmat(df_target, p) # needs to happen before data gets transformed.

    transform_data!(df_target,p)
    dₒ, dₙ = df_target.target, df_target.target_sim
    
    Cdinv, Cminv = inv(Cd), inv(Cm)

    function g(m)
        p.mp .= exp10.(m)
        simulate!(dataset, p ; allow_mismatch)
        df_target = get_target_data(dataset)
        add_data_points!(df_target, dataset, p)
        transform_data!(df_target, p)
        return Vector{Float64}(df_target[!,:target_sim])
    end

    #initialize plot
    ds_obs = Observable(dataset)#get_observable_dataset(dataset)
    #update_and_notify_obs!(ds_obs, dataset, p ; plot_transformed)
    vizualize(ds_obs)
    
    #iterate
    while (reldiff > reltol) & (n < maxiters)
        n+=1 ; println() ; @show(n)
        #@show (mₙₑₓₜ.-mₚ)./mₚ
        mₙ .= mₙₑₓₜ
        get_new_parameters!(g, mₙₑₓₜ, mₙ, mₚ, dₙ, dₒ, Cm, Cd ; μₙ, diff_method)
        p.mp .= exp10.(mₙₑₓₜ)
        simulate!(dataset, p ; allow_mismatch)
        df_target_next = get_target_data(dataset)
        add_data_points!(df_target_next, dataset, p)
        transform_data!(df_target_next,p)
        dₒ, dₙ = df_target_next.target, df_target_next.target_sim
        #Cd = get_transformed_covmat(df_target, p) # not needed if data always has the same length
        
        cost_data  = (dₙ-dₒ)'*Cdinv*(dₙ-dₒ)
        cost_model = (mₙₑₓₜ-mₚ)'*Cminv*(mₙₑₓₜ-mₚ)
        cost = cost_data + cost_model
        @show unique(Cdinv)
        @show mean((dₙ-dₒ).^2)
        @show cost_data
        @show cost_model
        @show cost

        #plot
        #update_and_notify_obs!(ds_obs, dataset, p ; plot_transformed)
        notify(ds_obs)
    end

    return dataset, p, Cm, Cd
end


function get_posterior_model_covariance(g::F, m::NamedArray, Cm, Cd) where F<:Function
    G = get_jacobian(g, m.array ; method=:finitediff)
    # println("G done")
    # Cm*G'
    # println("Cm*G' done")
    # G*Cm*G' + Cd
    # println("G*Cm*G' + Cd done")
    # inv(G*Cm*G' + Cd)
    # println("inv(G*Cm*G' + Cd) done")
    # G*Cm
    # println("G*Cm done")
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

function transform_data!(ds::Dataset, p)
    for i in eachindex(ds.timeseries)
        transform_data!(ds.timeseries[i], p)
    end
    for i in eachindex(ds.ponctual)
        transform_data!(ds.ponctual[i], p)
    end
    return ds
end

function transform_data!(ds::DatasetType, p)
    data = ds.data
    col_names = propertynames(data)
    if hasproperty(p.estimate,:transform_target)
        cols_to_transform = propertynames(p.estimate.transform_target)
        for col in cols_to_transform
            if col ∈ col_names 
                @views data[:,col] .= getproperty(p.estimate.transform_target, col)(data[:,col])
                @views data[:,string(col)*"_sim"] .= getproperty(p.estimate.transform_target, col)(data[:,string(col)*"_sim"])
            end
        end
    end
    return ds
end


function get_logspace_model_and_covmat(mp, std_mp ; max_log_std=4) #TODO : change kwarg to max_log_var
    log_lb = log10.(max.(mp.array .- std_mp, 1e-15))
    log_ub = log10.(mp.array .+ std_mp)
    Δlog_var = min.((log_ub.-log_lb).^2, max_log_std)
    Cm_log = Diagonal(Δlog_var)
    return log10.(mp.array), Cm_log
end
get_logspace_model_and_covmat(p) = get_logspace_model_and_covmat(p.mp, p.mp_std ; max_log_std = p.estimate.log_mp_std_max)

function get_transformed_covmat(dₒ_vec, std_vec, target_names, p) #TODO : change kwarg to max_log_var
    len = length(std_vec)
    Cd = Diagonal{Float64}(undef,len)
    funcs = [getproperty(p.estimate.transform_target, target_name) for target_name in target_names]
    is_positive_real_func = fill(false,length(funcs))
    println("entering try-catch block")
    for i in eachindex(funcs)
        try
            funcs[i](-1)
        catch e
            is_positive_real_func[i] = ifelse(e isa DomainError, true, false)
        end
    end
    println("exiting")

    for i in 1:len
        target_name = target_names[i]
        func_id = findfirst(target_name .== target_names)
        func = funcs[func_id]
        dₒ = dₒ_vec[i]
        std = std_vec[i]
        lb = ifelse(is_positive_real_func[i], func(max(dₒ - std, 1e-15)), func(dₒ - std))
        ub = func(dₒ + std)
        var = ifelse(is_positive_real_func[i], min((ub-lb)^2, p.estimate.log_mp_std_max), (ub-lb)^2)
        Cd[i,i] = var
    end
    return Cd
end
get_transformed_covmat(df_target, p) = get_transformed_covmat(df_target.target, df_target.std, df_target.target_name, p)

