"""
Graph of the rank of the subspace in which the 'true' model weights are allowed to drift, against degradation of concatenated decoder after n days.
"""


include("concatenated_decoder_numerics.jl")
using Plots
using MAT


vars = matread("../data/all_mice_day1.mat")
[vars[k] = v' for (k,v) in vars]
[vars[k] = v .- mean(v) for (k,v) in vars]






dropmean(A; dims=:) = dropdims(mean(A; dims=dims); dims=dims)
dropstd(A; dims=:) = dropdims(std(A; dims=dims); dims=dims)

const reduced_dims_to_test = [2,5,10,15,20,25,30,35,40,50,60,70,85,100]
const days_per_mouse = [10,12,10,5]
const num_repeats = 20

const true_deg =  [91.4563377325360,
93.8709428050535,
95.1427351797626,
91.5384171290284]

const mouse_numbers = [1,3,4,5]



function make_mice(data, reduced_dims)
    
    #                             #neurons, within_day_R2, between_day_R2, desired activiy covariance 
    mouse_1 = reducedRankInfo("mouse_1", size(data["mouse1_day1"])[2] , 101   , 1-0.2775  ,   1 - 0.3927, data["mouse1_day1"], reduced_dims)
    mouse_3 = reducedRankInfo("mouse_3", size(data["mouse3_day1"])[2] , 114   , 1-0.1745 ,   1 - 0.2752, data["mouse3_day1"], reduced_dims)
    mouse_4 = reducedRankInfo("mouse_4", size(data["mouse4_day1"])[2] , 146   , 1-0.1207  ,   1 - 0.1760, data["mouse4_day1"], reduced_dims)
    mouse_5 = reducedRankInfo("mouse_5", size(data["mouse5_day1"])[2] , 112   , 1-0.3208  ,   1 - 0.5722, data["mouse5_day1"], reduced_dims)

    mice = [mouse_1, mouse_3, mouse_4, mouse_5]
    return mice
end


function gen_final_mse(inf::generalInfo, num_days::Integer)
        d1 = initialise_day(inf)
        days = Array{dailyStore,1}(undef, num_days)
        days = fill_days!(days, inf, d1);
        concatenated_R2_vec = concat_R2_over_days(days)
        concat_rel_R2_vec = concatenated_R2_vec./concatenated_R2_vec[1]

        return concat_rel_R2_vec[end]
end



# mice = make_mice(vars)


function make_deg_data(ranks_to_test, num_days::Array{Int64,1}, num_repeats::Integer)
    
    
    # stores an array: for each mouse we have a matrix of which reduced rank, by repeats
     perf_storage = Array{Float64,3}(undef, length(mouse_numbers), length(ranks_to_test), num_repeats)

    for (i,r) in enumerate(ranks_to_test)
        
        for rep in 1:num_repeats
            mice = make_mice(vars, r)
            for (j,mouse) in enumerate(mice)
                perf_storage[j,i,rep] = gen_final_mse(mouse, num_days[j])
            end
        end



        # perf_storage[:,i,:] = [gen_final_mse(mouse, num_days[j]) for (j,mouse) in enumerate(mice), rep in 1:num_repeats]
    end

    
 return perf_storage
end

function make_plots(array_store)
    
    num_mice, num_ranks, reps = size(array_store)


    μs = dropmean(array_store, dims=3)
    σs = dropstd(array_store,dims=3)
    σmult = 1.96/sqrt(reps)

    subp = [plot(reduced_dims_to_test, μs[i,:], ribbon =  σmult*σs[i,:], leg=false, xlabel="rank", ylabel="normalised R2", title = "Mouse $(mouse_numbers[i])"
    ) for i in 1:num_mice]

    [hline!(p, [dg/100]) for (p,dg) in zip(subp,true_deg)]

    # plot!.(subp, ribbon=σs) 
    return subp
end

dd = make_deg_data(reduced_dims_to_test, days_per_mouse, num_repeats)
size(dd)





# 91.4563377325360 mouse 1
# 93.8709428050535 mouse 3
# 95.1427351797626 mouse 4
# 91.5384171290284 mouse 5

