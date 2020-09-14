include("concatenated_decoder_numerics.jl")
using Plots


# const within_day_R2 = 0.7
# const between_day_R2 = 0.6
# const num_data_per_day = 3000
# const num_regressors = 60

use_activity_spectrum = false
# const number_of_days = 14
function make_mice(use_activity_spectrum::Bool)
    
    if use_activity_spectrum == true
        m1cov = [6.0752e-01,4.7269e-01,3.4116e-01,2.4833e-01,2.0647e-01,1.7503e-01,1.5344e-01,1.4190e-01,1.2281e-01,1.0599e-01,9.5967e-02,8.5752e-02,7.9207e-02,7.3067e-02,6.8732e-02,6.4237e-02,6.0520e-02,5.6635e-02,5.2901e-02,5.0763e-02,4.8878e-02,4.5880e-02,4.3583e-02,4.1851e-02,3.9472e-02,3.8027e-02,3.6622e-02,3.5429e-02,3.4514e-02,3.3340e-02,3.2134e-02,3.1225e-02,3.0437e-02,2.9449e-02,2.8491e-02,2.7637e-02,2.6951e-02,2.6245e-02,2.5766e-02,2.5016e-02,2.4715e-02,2.4134e-02,2.3656e-02,2.3177e-02,2.2819e-02,2.2175e-02,2.1857e-02,2.1329e-02,2.0941e-02,2.0543e-02,2.0096e-02,1.9595e-02,1.9289e-02,1.8840e-02,1.8418e-02,1.8155e-02,1.7819e-02,1.7487e-02,1.7164e-02,1.6849e-02,1.6524e-02,1.6199e-02,1.5929e-02,1.5678e-02,1.5289e-02,1.4977e-02,1.4710e-02,1.4409e-02,1.4153e-02,1.3885e-02,1.3589e-02,1.3293e-02,1.3094e-02,1.2786e-02,1.2444e-02,1.2194e-02,1.1954e-02,1.1758e-02,1.1569e-02,1.1364e-02,1.1127e-02,1.0857e-02,1.0610e-02,1.0392e-02,1.0149e-02,9.9694e-03,9.7529e-03,9.3735e-03,9.1160e-03,8.7347e-03,8.5104e-03,8.2601e-03,8.0592e-03,7.7607e-03,7.3804e-03,7.0764e-03,6.6966e-03,6.4695e-03,5.9454e-03,5.1327e-03,4.5161e-03]
        m3cov = [5.1856e-01,3.5303e-01,2.8167e-01,2.2488e-01,1.8461e-01,1.5911e-01,1.3573e-01,1.2290e-01,1.1012e-01,9.9168e-02,9.0131e-02,8.3376e-02,7.6676e-02,6.9533e-02,6.5037e-02,5.9695e-02,5.6238e-02,5.1499e-02,4.7187e-02,4.4699e-02,4.1989e-02,3.9389e-02,3.6981e-02,3.4932e-02,3.3193e-02,3.1772e-02,3.0483e-02,2.9056e-02,2.7690e-02,2.6415e-02,2.5361e-02,2.4161e-02,2.3119e-02,2.1820e-02,2.1156e-02,2.0462e-02,1.9770e-02,1.9237e-02,1.8585e-02,1.7723e-02,1.7193e-02,1.6682e-02,1.6170e-02,1.5658e-02,1.5163e-02,1.4827e-02,1.4431e-02,1.3900e-02,1.3567e-02,1.3257e-02,1.2949e-02,1.2640e-02,1.2286e-02,1.2034e-02,1.1778e-02,1.1477e-02,1.1214e-02,1.0880e-02,1.0666e-02,1.0392e-02,1.0208e-02,9.9647e-03,9.7519e-03,9.4787e-03,9.2782e-03,9.1019e-03,8.9307e-03,8.7496e-03,8.5846e-03,8.4519e-03,8.2904e-03,8.1674e-03,8.0222e-03,7.8581e-03,7.7105e-03,7.5630e-03,7.4223e-03,7.2909e-03,7.1699e-03,7.0284e-03,6.9091e-03,6.7824e-03,6.6260e-03,6.5104e-03,6.3575e-03,6.2669e-03,6.1434e-03,6.0256e-03,5.9278e-03,5.8067e-03,5.6661e-03,5.5440e-03,5.4275e-03,5.3012e-03,5.1886e-03,5.0722e-03,4.9410e-03,4.8334e-03,4.7502e-03,4.5997e-03,4.4721e-03,4.4078e-03,4.2500e-03,4.1020e-03,3.9635e-03,3.8121e-03,3.6720e-03,3.5131e-03,3.3222e-03,3.1701e-03,3.0127e-03,2.8558e-03,2.6274e-03,2.3071e-03]
        m4cov = [6.8626e-01,5.4418e-01,3.6776e-01,2.9072e-01,2.5407e-01,2.1448e-01,1.8009e-01,1.6125e-01,1.4891e-01,1.3515e-01,1.2624e-01,1.1588e-01,1.0750e-01,1.0107e-01,9.4269e-02,9.0007e-02,8.5123e-02,8.0428e-02,7.6085e-02,7.2491e-02,6.7994e-02,6.5565e-02,6.2720e-02,5.9669e-02,5.7525e-02,5.4778e-02,5.2595e-02,5.0903e-02,4.8995e-02,4.7508e-02,4.5545e-02,4.4043e-02,4.2087e-02,4.0713e-02,3.9308e-02,3.7594e-02,3.6218e-02,3.5064e-02,3.3681e-02,3.2833e-02,3.1690e-02,3.0662e-02,2.9633e-02,2.8761e-02,2.8005e-02,2.7053e-02,2.6260e-02,2.5220e-02,2.4571e-02,2.3902e-02,2.3111e-02,2.2564e-02,2.1895e-02,2.1318e-02,2.0643e-02,2.0144e-02,1.9629e-02,1.9158e-02,1.8612e-02,1.8187e-02,1.7801e-02,1.7341e-02,1.6840e-02,1.6478e-02,1.6206e-02,1.5803e-02,1.5334e-02,1.4971e-02,1.4666e-02,1.4436e-02,1.4054e-02,1.3761e-02,1.3428e-02,1.3066e-02,1.2845e-02,1.2519e-02,1.2277e-02,1.2035e-02,1.1896e-02,1.1602e-02,1.1374e-02,1.1052e-02,1.0855e-02,1.0612e-02,1.0468e-02,1.0257e-02,1.0070e-02,9.8670e-03,9.7276e-03,9.5959e-03,9.3643e-03,9.1964e-03,9.0008e-03,8.8851e-03,8.7358e-03,8.5415e-03,8.3890e-03,8.2529e-03,8.0847e-03,7.9672e-03,7.7854e-03,7.6789e-03,7.5392e-03,7.4024e-03,7.2578e-03,7.1316e-03,6.9902e-03,6.8490e-03,6.7487e-03,6.6218e-03,6.5003e-03,6.3946e-03,6.2740e-03,6.1717e-03,6.0443e-03,5.9247e-03,5.8183e-03,5.7076e-03,5.5627e-03,5.4462e-03,5.3063e-03,5.1704e-03,5.0782e-03,4.9509e-03,4.7866e-03,4.6398e-03,4.4704e-03,4.3567e-03,4.0818e-03,3.9330e-03,3.7702e-03,3.5949e-03,3.2177e-03,2.8219e-03]
        m5cov = [5.2866e-01,3.4033e-01,2.7358e-01,2.0986e-01,1.7540e-01,1.5190e-01,1.3393e-01,1.1998e-01,1.1300e-01,1.0478e-01,9.2336e-02,8.7652e-02,8.3023e-02,7.8376e-02,7.2516e-02,6.9795e-02,6.6173e-02,6.2837e-02,6.0105e-02,5.7222e-02,5.4818e-02,5.1678e-02,4.9705e-02,4.7092e-02,4.5415e-02,4.3394e-02,4.1865e-02,4.0213e-02,3.8791e-02,3.7554e-02,3.5910e-02,3.5150e-02,3.3664e-02,3.2531e-02,3.1300e-02,3.0364e-02,2.9264e-02,2.8371e-02,2.7432e-02,2.6742e-02,2.5755e-02,2.4852e-02,2.4085e-02,2.3427e-02,2.2740e-02,2.2115e-02,2.1069e-02,2.0426e-02,1.9725e-02,1.9041e-02,1.8375e-02,1.7898e-02,1.7370e-02,1.6977e-02,1.6450e-02,1.5848e-02,1.5390e-02,1.4980e-02,1.4527e-02,1.4267e-02,1.3787e-02,1.3510e-02,1.3183e-02,1.2749e-02,1.2357e-02,1.2016e-02,1.1680e-02,1.1438e-02,1.1193e-02,1.0792e-02,1.0597e-02,1.0450e-02,1.0166e-02,9.9108e-03,9.5868e-03,9.2742e-03,9.0879e-03,8.8815e-03,8.6365e-03,8.4415e-03,8.2551e-03,8.0833e-03,7.9321e-03,7.7619e-03,7.5777e-03,7.3978e-03,7.2047e-03,7.0633e-03,6.9236e-03,6.7921e-03,6.6722e-03,6.4999e-03,6.3772e-03,6.2596e-03,6.1521e-03,5.9734e-03,5.8257e-03,5.6604e-03,5.5196e-03,5.3554e-03,5.1730e-03,5.0042e-03,4.8351e-03,4.6480e-03,4.4361e-03,4.1429e-03,3.9063e-03,3.7177e-03,3.3790e-03,2.9524e-03,2.5660e-03,2.1906e-03]
    else
        m1cov = nothing
        m3cov = nothing
        m4cov = nothing
        m5cov = nothing
    end


    #                             #neurons, within_day_R2, between_day_R2, desired activiy covariance 
    mouse_1 = v1Info("mouse_1", 5000 , 101   , 1-0.2775  ,   1 - 0.3927, m1cov  )
    mouse_3 = v1Info("mouse_3", 5000 , 114   , 1- 0.1745 ,   1 - 0.2752, m3cov   )
    mouse_4 = v1Info("mouse_4", 5000 , 134   , 1-0.1207  ,   1 - 0.1760, m4cov   )
    mouse_5 = v1Info("mouse_5", 5000 , 112   , 1-0.3208  ,   1 - 0.5722, m5cov  )

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
    p = plot(R2_mean, ribbon = 1.96*R2_std, xaxis="days", yaxis="R2 (% of day 1)")
    plot!(p,    xticks=pl.xTicks, 
                yticks = pl.yTicks, 
                #  xlim=[0.5, 1],
                ylim=[0.5, 1],
                title=pl.name,
                legend=false)
    # plot!(p, xlim=[min.(pl.xTicks), max.(pl.xTicks)], ylim=[min.(pl.yTicks), max.(pl.yTicks)])
    return R2_mean, R2_std, p
end


everything_with_spectrum = [generate_graph(mouse, plot_stuff) for (mouse, plot_stuff) in zip(make_mice(true)...)]
everything_without_spectrum = [generate_graph(mouse, plot_stuff) for (mouse, plot_stuff) in zip(make_mice(false)...)]



plots_without_spectrum = [stuff[3] for stuff in everything_without_spectrum]
plots_with_spectrum = [stuff[3] for stuff in everything_with_spectrum]

function merged_plot()
    
    p = [plot() for el in everything_without_spectrum]
    for (pl, stuff1, stuff2) in (zip(p, everything_with_spectrum, everything_without_spectrum))
        plot!(pl, stuff1[1], ribbon=1.96*stuff1[2], xaxis="days", yaxis="R2 (% of day 1)")
    end
    return p
end


R2_vecs = [stuff[1] for stuff in everything_with_spectrum]
R2_stds = [stuff[2] for stuff in everything_with_spectrum]
the_plots = [stuff[3] for stuff in everything_with_spectrum]
# all_plots = plot(the_plots..., layout=4)
R2_95pc_Conf = [1.96*el for el in R2_stds] 


# readme = "each variable has 4 elements, correspondng to mouse 1,3,4,5 respectively. R2 stds are the standard deviations 
# for the ribbons. need to be multiplied by 1.96 for 95% confidence intervals"


using MAT

# file = matopen("with_spectrum.mat", "w")
# write(file, "R2_vecs", R2_vecs)
# write(file,"readme", readme)
# write(file, "R2_stds", R2_stds)
# write(file, "R2_95pc_Conf", R2_95pc_Conf)
# close(file)

function compare_specs(i)
R2_with = [stuff[1] for stuff in everything_with_spectrum]
R2_without = [stuff[1] for stuff in everything_without_spectrum]

p = plot(R2_with[i], label="with spectrum")
plot!(p, R2_without[i], label="without_spectrum")
end