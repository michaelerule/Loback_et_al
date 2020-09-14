"""
What I want to do:

similar to before. Generate data X. However, now I make X very sparse. Then I see how concatenated decoder performance drops 
over time as the weights change.

"""

include("concatenated_decoder_numerics.jl")
using Plots


# const within_day_R2 = 0.7
# const between_day_R2 = 0.6
# const num_data_per_day = 3000
# const num_regressors = 60

sparsity_level = 0.9
# const number_of_days = 14
function make_mice(sparsity_level)
    
    #                             #neurons, within_day_R2, between_day_R2, desired activiy covariance 
    mouse_1 = v2Info("mouse_1", 5000 , 101   , 1-0.2775  ,   1 - 0.3927, sparsity_level)
    mouse_3 = v2Info("mouse_3", 5000 , 114   , 1- 0.1745 ,   1 - 0.2752, sparsity_level)
    mouse_4 = v2Info("mouse_4", 5000 , 134   , 1-0.1207  ,   1 - 0.1760, sparsity_level)
    mouse_5 = v2Info("mouse_5", 5000 , 112   , 1-0.3208  ,   1 - 0.5722, sparsity_level)

    mice = [mouse_1, mouse_3, mouse_4, mouse_5]

    # 83 93 94 89
    # 14 13 13 6
    mp1 = plotInfo(collect(1:14),  [1., 0.83, 0.6, 0.5], "Mouse 1" )
    mp3 = plotInfo(collect(1:13), [1., 0.93, 0.75, 0.5], "Mouse 3"  )
    mp4 = plotInfo(collect(1:13), [1., 0.94, 0.85, 0.5], "Mouse 4")
    mp5 = plotInfo(collect(1:6), [1., 0.89, 0.7, 0.5], "Mouse 5")


    mice_plotInfo = [mp1, mp3, mp4, mp5]
    return mice, mice_plotInfo
end

function generate_graph(inf::generalInfo, pl::plotInfo)
    
    function gen_data()
        d1 = initialise_day(inf)
        days = Array{dailyStore,1}(undef, Int64(pl.xTicks[end]))
        days = fill_days!(days, inf, d1);
        concatenated_R2_vec = concat_R2_over_days(days)
        concat_rel_R2_vec = concatenated_R2_vec./concatenated_R2_vec[1]
    end
    
    num_trials = 10
    R2_vecs = [gen_data() for el in 1:num_trials]
    R2_matrix = hcat(R2_vecs...)
    R2_std = std(R2_matrix, dims=2)
    R2_mean = mean(R2_matrix, dims=2)
    p = plot(R2_mean, ribbon = 1.96*R2_std, xaxis="days", yaxis="R2 (% of day 1)",
                xticks=pl.xTicks, 
                yticks = pl.yTicks, 
                title=pl.name,
                legend=false)
    # plot!(p, xlim=[min.(pl.xTicks), max.(pl.xTicks)], ylim=[min.(pl.yTicks), max.(pl.yTicks)])
    return R2_mean, R2_std, p
end


everything_with_sparsity = [generate_graph(mouse, plot_stuff) for (mouse, plot_stuff) in zip(make_mice(sparsity_level)...)]
everything_without_sparsity = [generate_graph(mouse, plot_stuff) for (mouse, plot_stuff) in zip(make_mice(0.)...)]



plots_without_sparsity = [stuff[3] for stuff in everything_without_sparsity]
plots_with_sparsity = [stuff[3] for stuff in everything_with_sparsity]

function merged_plot()
    
    p = [plot() for el in everything_without_sparsity]
    for (pl, stuff1, stuff2) in (zip(p, everything_with_sparsity, everything_without_sparsity))
        plot!(pl, stuff1[1], ribbon=1.96*stuff1[2], xaxis="days", yaxis="R2 (% of day 1)")
    end
    return p
end


R2_vecs = [stuff[1] for stuff in everything_with_sparsity]
R2_stds = [stuff[2] for stuff in everything_with_sparsity]
the_plots = [stuff[3] for stuff in everything_with_sparsity]
dense_plots = [stuff[3] for stuff in everything_wihtout_sparsity]
# all_plots = plot(the_plots..., layout=4)
R2_95pc_Conf = [1.96*el for el in R2_stds] 


# readme = "each variable has 4 elements, correspondng to mouse 1,3,4,5 respectively. R2 stds are the standard deviations 
# for the ribbons. need to be multiplied by 1.96 for 95% confidence intervals"


using MAT

# file = matopen("with_sparsity.mat", "w")
# write(file, "R2_vecs", R2_vecs)
# write(file,"readme", readme)
# write(file, "R2_stds", R2_stds)
# write(file, "R2_95pc_Conf", R2_95pc_Conf)
# close(file)

function compare_specs(i)
R2_with = [stuff[1] for stuff in everything_with_sparsity]
R2_without = [stuff[1] for stuff in everything_without_sparsity]

p = plot(R2_with[i], label="with sparsity")
plot!(p, R2_without[i], label="without_sparsity")
end