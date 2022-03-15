
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
    Cd, pos_real_func = get_transformed_covmat(df_target, p) # needs to happen before data gets transformed.
    transform_data!(df_target,p)
    dₒ, dₙ = df_target.target, df_target.target_sim
    
    Cminv = inv(Cm)

    function g(m)
        p.mp .= exp.(m)
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
        p.mp .= exp.(mₙₑₓₜ)
        simulate!(dataset, p ; allow_mismatch)
        df_target_next = get_target_data(dataset)
        add_data_points!(df_target_next, dataset, p)
        Cd, _ = get_transformed_covmat(df_target_next, p ; is_positive_real_func=pos_real_func) # needs to happen before data gets transformed.
        #@show Cd[[1,40],[1,40]]
        Cdinv = inv(Cd)
        transform_data!(df_target_next,p)
        dₒ, dₙ = df_target_next.target, df_target_next.target_sim
        #Cd = get_transformed_covmat(df_target, p) # not needed if data always has the same length
        
        cost_data  = (dₙ-dₒ)'*Cdinv*(dₙ-dₒ)
        cost_model = (mₙₑₓₜ-mₚ)'*Cminv*(mₙₑₓₜ-mₚ)
        cost = cost_data + cost_model
        #@show unique(Cdinv)
        #@show mean((dₙ-dₒ).^2)
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
    return Cm - Cm*G' * inv(G*Cm*G' + Cd) * G*Cm
end
function get_posterior_model_covariance(g::F, m::Vector, Cm, Cd) where F<:Function
    G = get_jacobian(g, m ; method=:finitediff)
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

"log transform so that lognormal distribution σ tends toward std when m goes further away from 0, and μ is exactly equal to ln(m), i.e. the median of the lognormal"
function to_log_from_physical_m_and_std(m,std)
    σ² = log(1+std^2/m^2)
    μ = log(m)
    return μ, sqrt(σ²)
end

function get_logspace_model_and_covmat(mp, std_mp) 
    m_std_log = to_log_from_physical_m_and_std.(mp,std_mp)
    m_log = [tup[1] for tup in m_std_log]
    std_log = [tup[2] for tup in m_std_log]
    Cm_log = Diagonal(std_log.^2)
    return m_log, Cm_log
end
get_logspace_model_and_covmat(p) = get_logspace_model_and_covmat(p.mp, p.mp_std)

function get_transformed_covmat(dₒ_vec, std_vec, target_names, p ; is_positive_real_func=nothing) 
    len = length(std_vec)
    Cd = Diagonal{Float64}(undef,len)
    funcs = [getproperty(p.estimate.transform_target, target_name) for target_name in target_names]
    if isnothing(is_positive_real_func)
        is_positive_real_func = fill(false,length(funcs))
        for i in eachindex(funcs)
            try
                funcs[i](-1)
            catch e
                is_positive_real_func[i] = ifelse(e isa DomainError, true, false)
            end
        end
        for i in eachindex(funcs)
            if is_positive_real_func[i]
                if all(funcs[i].(dₒ_vec) .≈ log.(dₒ_vec))
                    continue
                else
                    throw(error("transform function associated to $(target_names[i]) is a positive real function that is not log"))
                end
            end
        end
    end

    for i in 1:len
        target_name = target_names[i]
        func_id = findfirst(target_name .== target_names)
        func = funcs[func_id]
        dₒ = dₒ_vec[i]
        std = std_vec[i]
        if is_positive_real_func[func_id]
            dlog, stdlog = to_log_from_physical_m_and_std(dₒ,std)
            var = stdlog^2
        else
            lb = func(dₒ - std)
            ub = func(dₒ + std)
            var =(ub-lb)^2
        end
        Cd[i,i] = var
    end
    return Cd, is_positive_real_func
end
get_transformed_covmat(df_target, p ; is_positive_real_func=nothing) = get_transformed_covmat(df_target.target, df_target.std, df_target.target_name, p ; is_positive_real_func)

