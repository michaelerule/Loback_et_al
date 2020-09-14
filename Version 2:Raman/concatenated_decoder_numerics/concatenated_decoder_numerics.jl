using LinearAlgebra
using Statistics
using MultivariateStats
using Random
using Distributions
using JuMP
using Ipopt
"""
Y = RX + ϵ
SSres = Var(ϵ)
SStot = Var(Y) = Var(RX) + Var(ϵ)   # assuming independence of ϵ from RX 

R^2 = 1 - SSres/SStot
    =  Var(RX)/[Var(RX) + Var(ϵ)]
Required value of Var(ϵ):
Var(ϵ) = Var(RX)[(1-R^2)/R^2]



activity spectrum stuff:
1) get covariance matrix of X
2) multiply singular values to get to desired levels 
"""


struct plotInfo
    xTicks::Vector{Float64}
    yTicks::Vector{Float64}
    name::String
end
# plotInfo(xTicks, yTicks,name::String) =  plotInfo(xTicks, yTicks, 
# [min.(xTicks), max.(xTicks)],
# [min.(yTicks), max.(yTicks)],
# name)

"""
Stuff independent of day changes
"""
abstract type generalInfo end
abstract type withDataInfo <: generalInfo end
"""
v1 was the null model of the initial submission. We could set a desired covariance spectrum when generating data X
"""
struct v1Info <: generalInfo
    name::String
    num_data_per_day::Int64
    num_regressors::Int64
    within_day_R2::Float64 #within day R squared of Y against X
    between_day_R2::Float64
    desired_cov_spectrum::Union{Nothing,Array{Float64,1}}
end
v1Info(num_data_per_day::Int64,num_regressors, within_day_R, between_day_R2) =
v1Info("unnamed", num_data_per_day,num_regressors, within_day_R, between_day_R2, nothing)
v1Info(name::String, num_data_per_day::Int64,num_regressors, within_day_R, between_day_R2) =
v1Info(name, num_data_per_day,num_regressors, within_day_R, between_day_R2, nothing)


"""
v2 allows construction of a matrix with a fixed level of sparsity
"""
struct v2Info <: generalInfo
    name::String
    num_data_per_day::Int64
    num_regressors::Int64
    within_day_R2::Float64 #within day R squared of Y against X
    between_day_R2::Float64
    sparsity_level::Float64
end

"""
v3 uses the actual mouse data, instead of a generative model.
"""
struct v3Info <: withDataInfo
    name::String
    num_data_per_day::Int64
    num_regressors::Int64
    within_day_R2::Float64 #within day R squared of Y against X
    between_day_R2::Float64
    X::Array{Float64,2} # already provided the data to use.
end


struct reducedRankInfo <: withDataInfo
    name::String
    num_data_per_day::Int64
    num_regressors::Int64
    within_day_R2::Float64 #within day R squared of Y against X
    between_day_R2::Float64
    X::Array{Float64,2} # already provided the data to use.
    num_drift_dimensions::Int64
    drift_space::Array{Float64,2} # Need to create
end

function make_drift_space(reg_d::Int64, drift_d::Int64)
    Q,_ = qr(randn(reg_d,reg_d))
    Q = Q[1:drift_d, :] # this is now an orthonormal set of drift_d vectors. the basis for the drift space
end

reducedRankInfo(name, num_data_per_day::Int64, num_regressors::Int64, within_day_R2::Float64, between_day_R2::Float64,X, num_drift_dimensions::Int64) =
reducedRankInfo(name, num_data_per_day::Int64, num_regressors::Int64, within_day_R2::Float64, between_day_R2::Float64,X, num_drift_dimensions::Int64, make_drift_space(num_regressors, num_drift_dimensions))

abstract type AbstractDailyStore end
struct dailyStore <: AbstractDailyStore
    M::Adjoint{Float64,Array{Float64,1}}  # true regressor
    X::Array{Float64,2}  # 'neural activity' matrix for day I
    Y::Adjoint{Float64, Array{Float64,1}}  # 'observed activity' matrix for day I
end
 #then make array of dailyStore
struct reducedRankDailyStore <: AbstractDailyStore
    M::Adjoint{Float64,Array{Float64,1}}  # true regressor
    X::Array{Float64,2}  # 'neural activity' matrix for day I
    Y::Adjoint{Float64, Array{Float64,1}}  # 'observed activity' matrix for day I
end



function make_regressor(inf::generalInfo)
    return rand(inf.num_regressors)'

end

function get_ΔM(M,X,Y,inf::generalInfo)
    var_ϵ = get_var_ϵ(M,X,inf::generalInfo)
    ΔM_norm2 = ((inf.num_regressors*inf.num_data_per_day)/tr(X*X'))*(
                (1-inf.between_day_R2)*var(Y) - var_ϵ
                                                                    )
    ΔM  = randn!(similar(M))
    if dot(ΔM,M) < 0
        ΔM = - ΔM
    end
    M_norm = sqrt(dot(M,M))
    proj = M*((dot(ΔM, M))/((M_norm^2)))
    ortho = ΔM - proj
    proj = proj/norm(proj)
    ortho = ortho/norm(ortho)
    c1 = -ΔM_norm2/(2*norm(M)) 
    c2 = sqrt(ΔM_norm2-c1^2) #so c1^2 + c2^2 = K
    ΔM = c1*proj + c2*ortho
    # println("c1 = $(c1), c2 = $(c2)")
    # println(c1^2 + c2^2, "vs", dot(ΔM,ΔM))
    # println(dot(proj,ortho))
    # # println("correlation is $(dot(ΔM,M)) ")
    # # println("norm2 is $(dot(ΔM, ΔM))")
    #  println("norm M vs norm M+ΔM is $(dot(M,M)) vs $(dot(M+ΔM, M + ΔM)) ")
    # # println("desired ΔM_norm2 is $(ΔM_norm2) \n")
    return ΔM
end

function get_ΔM(M,X,Y,inf::reducedRankInfo)
    ## This one is for a low dim subspace
    var_ϵ = get_var_ϵ(M,X,inf::generalInfo)
    ΔM_norm2 = ((inf.num_regressors*inf.num_data_per_day)/tr(X*X'))*(
                (1-inf.between_day_R2)*var(Y) - var_ϵ
                                                                    )
    corrs = M*inf.drift_space'
    # """
    # solve qclp: \Delta M = \sum_i c_i e_i, where e_i are the rows of drift_space
    # need
    #  \|c\|_2^2 = K
    #  \sum_i c_i corrs_i = -K/2

    #  where K = what it is in the notebook 
    # """
    model = Model(Ipopt.Optimizer)
    @variable(model, c[1:inf.num_drift_dimensions])
    @constraint(model, sum(c.^2) == ΔM_norm2)
    @constraint(model, -2*dot(c,corrs) == ΔM_norm2)
    
    for i = 1:inf.num_drift_dimensions
        set_start_value(c[i],randn())
    end

    JuMP.optimize!(model)
    cc = value.(c) 
    ΔM = cc'*inf.drift_space
    return ΔM
end



function get_var_ϵ(M,X,inf::generalInfo)
    var_ϵ = var(M*X)*((1-inf.within_day_R2)/inf.within_day_R2)
end

function random_rotation!(to_rotate, rotation_template,Q,R)
    randn!(rotation_template)
    Q,R = qr(rotation_template) #Q is now a random orthogonal matrix.
    to_rotate *= Q
    return to_rotate
end


function generate_X(inf::v1Info)
    X = randn(inf.num_regressors, inf.num_data_per_day)
    if typeof(inf.desired_cov_spectrum) != Nothing
        U,S,V = svd(X)
        X = sqrt(inf.num_data_per_day)U*Diagonal(sqrt.(inf.desired_cov_spectrum))*V' #replace singular values with desired spectrum 
    end
    return X
end

function generate_X(inf::v2Info)
    X = randn(inf.num_regressors, inf.num_data_per_day)
    b = Bernoulli(inf.sparsity_level) 
    B = 1 .- rand(b, size(X)) #more sparsity means more zeros
    return (1/(1-inf.sparsity_level)).*X.*B
end

function generate_data(inf::generalInfo, opt_regressor)
    X = generate_X(inf::generalInfo)
    ϵ = sqrt(get_var_ϵ(opt_regressor, X, inf)  
            )*randn(inf.num_data_per_day)
    Y = opt_regressor*X + ϵ'
    return X, Y
end

function generate_data(inf::withDataInfo, opt_regressor)
    X = inf.X
    ϵ = sqrt(get_var_ϵ(opt_regressor, X, inf)  
            )*randn(inf.num_data_per_day)
    Y = opt_regressor*X + ϵ'
    return X, Y
end


function generate_data(inf::generalInfo, opt_regressor, X)
    ϵ = sqrt(get_var_ϵ(opt_regressor, X, inf)  
            )*randn(inf.num_data_per_day)
    Y = opt_regressor*X + ϵ'
    return X, Y
end

function initialise_day(inf::generalInfo)
    M = make_regressor(inf)
    X, Y = generate_data(inf, M)
    return dailyStore(M, X, Y)
end

function initialise_day(inf::generalInfo, X)
    M = make_regressor(inf)
    X, Y = generate_data(inf, M, X)
    return dailyStore(M, X, Y)
end


function update_day(d::dailyStore, inf::generalInfo)
    Xnew = update_X(d.X, inf)
    Mnew = update_reg(d.M, d.X, d.Y, inf)
    Ynew = update_Y(Mnew, Xnew, inf)
    return dailyStore(Mnew, Xnew, Ynew)
end



function update_reg(M, X, Y, inf::generalInfo)
    Mnew = M + get_ΔM(M, X, Y, inf)  #randn normalised with l2 norm = inf.regressor_norm_shift
    # Rnew .*= norm(R)/norm(Rnew)
    return Mnew
end

function update_X(X, inf::generalInfo)
    # ## ISSUE: THIS WILL TEND TO INCREASE THE NORM OF X. 
    # ϵₓ = (  inf.activity_norm_shift/sqrt(length(X))     )*randn(size(X)) #randn normalised to have l2 norm = inf.activity_norm_shift
    # return Xnew = X + ϵₓ
    return X
end

function update_Y(M, X, inf::generalInfo)
    ϵ = sqrt(get_var_ϵ(M, X, inf)  
            )*randn(inf.num_data_per_day)
    Ynew = M*X + ϵ'
end

function generate_data(inf::generalInfo)
    X = randn(inf.num_regressors, inf.num_data_per_day)
end

struct concatenatedResults
    concatenated_regressors::Adjoint{Float64,Array{Float64,1}}
    R2
    MSE
end

function concatenated_decoder(arr::Array{dailyStore,1})
    X = hcat((el.X for el in arr)...)
    Y = hcat((el.Y for el in arr)...)
    concatenated_regressors = llsq(X', Y'; bias=false)'
    varMX = var(concatenated_regressors*X)
    varϵ = var(Y - concatenated_regressors*X)
    R2 = varMX/(varMX + var(Y - concatenated_regressors*X)) 
    MSE = mean((Y - concatenated_regressors*X).^2)
    return concatenatedResults(concatenated_regressors, R2, MSE)
end

function get_true_R2(d::dailyStore)
    varϵ = var(d.Y - d.M*d.X)
    varMX = var(d.M*d.X)
    R2 = 1 - (varϵ/(varϵ + varMX))
end

function get_est_reg(d::dailyStore)
     est_M = llsq(d.X',d.Y'; bias=false)'
end

function get_est_R2(d::dailyStore)
    est_M = get_est_reg(d)
    varϵ = var(d.Y - est_M*d.X)
    varMX = var(est_M*d.X)
    R2 = 1 - (varϵ/(varϵ + varMX))
end


function concat_R2_over_days(arr::Array{dailyStore,1})
    R2vec = Array{Float64,1}(undef,length(arr))
    for i in 1:length(arr)
        R2vec[i] = concatenated_decoder(arr[1:i]).R2
    end
    return R2vec
end

function concat_MSE_over_days(arr::Array{dailyStore,1})
    msevec = Array{Float64,1}(undef,length(arr))
    for i in 1:length(arr)
        msevec[i] = concatenated_decoder(arr[1:i]).MSE
    end
    return msevec
end

function fill_days!(arr::Array{dailyStore,1}, inf::generalInfo, d1)
    arr[1] = d1
    for i in 1:length(arr)
        if i == 1
            arr[i] = d1
        else
            arr[i] = update_day(arr[i-1], inf)
        end
    end
    return arr
end

function get_Y_velocity(arr::Array{dailyStore,1}, inf::generalInfo)
    velocity_vector = Array{Float64,1}(undef,length(arr))
    for (i,el) in enumerate(arr)
        if i > 1
            velocity_vector[i] = norm(el.Y - arr[i-1].Y)
        end
    end
    return velocity_vector
end

function get_est_between_day_R2(arr::Array{dailyStore,1}, inf=nothing)
    R2vec = Array{Float64,1}(undef,length(arr))
    per_day_regressor = get_est_reg.(arr)
    for (i,el) in enumerate(arr)
        if i > 1
            # R2vec[i] = 1 -   sum((el.Y- per_day_regressor[i-1]*arr[i].X).^2)/(sum((el.Y .- mean(el.Y)).^2))
            R2vec[i] = 1 - var( el.Y- per_day_regressor[i-1]*arr[i].X )/ var(el.Y)
        end
    end
    R2vec[1] = NaN
    return R2vec[2:end]
end

function get_est_mse(d::dailyStore)
    M = get_est_reg(d)
    return mean((d.Y - M*d.X).^2)
end

function get_true_mse(d::dailyStore)
    return mean((d.Y - d.M*d.X).^2)
end

function check_in_range(vec, mat)
    ## see how close vec is to the range of mat
    ext = vcat(mat,vec)
    return rank(ext) == rank(mat)
end

function check_everything(inf::generalInfo, days::Array{dailyStore,1})
    between_day_R2 = get_est_between_day_R2(days)
    within_day_R2 = get_est_R2.(days)
    @info(" Name: $(inf.name) \n
    between day R2 results: \n 
    mean: $(mean(between_day_R2))        std: $(std(between_day_R2))          
    specified: $(inf.between_day_R2)
    \n
    within day R2 results: \n 
    mean: $(mean(within_day_R2))        std: $(std(within_day_R2))          
    specified: $(inf.within_day_R2)
    \n
    " )

    if typeof(inf) <: reducedRankInfo
        Δregs = [two.M - one.M for (one, two) in zip(days[1:end-1], days[2:end])]
        in_subspace = [check_in_range(el, inf.drift_space) for el in Δregs]
        @info(" boolean: is drift constrained to subspace per day?
        \n
        $in_subspace
        ")
    end
end

function plot_data_from_everything(runtype, everything)
    readme = "each variable has 4 elements, correspondng to mouse 1,3,4,5 respectively. R2 stds are the standard deviations 
for the ribbons. they are multiplied by 1.96 for 95% confidence intervals"

    
    R2_vecs = [stuff[1] for stuff in everything]
    R2_stds = [stuff[2] for stuff in everything]
    R2_95pc_Conf = [1.96*el for el in R2_stds] 

    file = matopen("$(runtype).mat", "w")
    write(file, "R2_vecs", R2_vecs)
    write(file,"readme", readme)
    write(file, "R2_stds", R2_stds)
    write(file, "R2_95pc_Conf", R2_95pc_Conf)
    close(file)
end