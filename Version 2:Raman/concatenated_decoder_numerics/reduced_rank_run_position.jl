"""
Taking pseudo-neural activity as the TRUE matrix of neural activity of each mouse on day 1. Rather than something generated 


"""



include("concatenated_decoder_numerics.jl")
using Plots

using MAT

vars = matread("../data/all_mice_day1.mat")


[vars[k] = v' for (k,v) in vars]
[vars[k] = v .- mean(v) for (k,v) in vars]

reduced_dims = 50

function make_mice(data)
    
    


    #                             #neurons, within_day_R2, between_day_R2, desired activiy covariance 
    mouse_1 = reducedRankInfo("mouse_1", size(data["mouse1_day1"])[2] , 101   , 1-0.2775  ,   1 - 0.3927, data["mouse1_day1"], reduced_dims)
    mouse_3 = reducedRankInfo("mouse_3", size(data["mouse3_day1"])[2] , 114   , 1-0.1745 ,   1 - 0.2752, data["mouse3_day1"], reduced_dims)
    mouse_4 = reducedRankInfo("mouse_4", size(data["mouse4_day1"])[2] , 146   , 1-0.1207  ,   1 - 0.1760, data["mouse4_day1"], reduced_dims)
    mouse_5 = reducedRankInfo("mouse_5", size(data["mouse5_day1"])[2] , 112   , 1-0.3208  ,   1 - 0.5722, data["mouse5_day1"], reduced_dims)

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
    # plot!(p, ylim=[min.(pl.yTicks), max.(pl.yTicks)])
    return R2_mean, R2_std, p
end


everything = [generate_graph(mouse, plot_stuff) for (mouse, plot_stuff) in zip(make_mice(vars)...)]


the_plots = [stuff[3] for stuff in everything]
[plot!(el, ylim = [0.2, 1]) for el in the_plots]


function gen_days(mouse, pl)
    d1 = initialise_day(mouse)
    days = Array{dailyStore,1}(undef, Int64(pl.xTicks[end]))
    days = fill_days!(days, mouse, d1);    
end

function check_per_mouse(mice, pl)
    days_per_mouse = gen_days.(mice,pl)
    check_everything.(mice, days_per_mouse)
end


mice, pl = make_mice(vars)
check_per_mouse(mice, pl)
q = plot(the_plots...)

y = ones(3)
title = Plots.scatter(y, marker=0,markeralpha=0, annotations=(2, y[2], Plots.text("Forward position")),axis=false, leg=false,size=(200,100))

# combine the 'title' plot with your real plots
qq = Plots.plot(
    title,
    q,
    layout=grid(2,1,heights=[0.1,0.9])
)