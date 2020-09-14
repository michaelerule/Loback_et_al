# Outline of the code:

- What we are doing mathematically is contained in concatenated_decoder_methods.pdf

- Code used to generate paper figures are:

Figure 3d data comes from running v1_results_generator.jl. 

Figure 4 supplement 1: reduced_rank_graph.jl

reduced_rank_graph.jl, which generates plots of the number of dimensions in which pseudo-neural activity can drift, against overall degradation of the concatenated decoder. 




- concatenated_decoder_numerics.jl sets up most of the code. other .jl files run scripts using functions from concatenated_decoder_numerics.jl. 

- The easiest way to see what is going on is to look up v1_results_generator.jl or v2_results_generator.jl and follow the code.

## The basic structure of the code is such:


- The data of each mouse is held in an abstract type called generalInfo. 
- There are concrete subtypes v1Info, v2Info, and v3Info <: generalInfo, that are used for different numerical experiments 

Eg we declare a mouse with:

                                #neurons, within_day_R2, between_day_R2, desired activiy covariance 
    mouse_1 = v1Info("mouse_1", 5000 ,       101   ,     1-0.2775  ,       1 - 0.3927,      m1cov  )


See concatenated_decoder_numerics.jl for descriptions of the different concrete subtypes described.

On each day of real/pseudo data, we store:
- mouse neural activity
- external variabales (eg position) to regress against
- optimal regression coefficients (linear regression
This is in the concrete type: dailystore <: AbstractDailyStore

We have a function fill_days! that populates an empty array of sequential days, given a populated initial day. On each new day, the function

- changes mouse neural activity in a manner consistent with the drift statistics specified in the mouse generalInfo container.

- updates the optimal regressors

We then have a function concat_R2_over_days() that calculates the R2 of the best static decoder trained over an array of days. 



