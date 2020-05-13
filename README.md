# Stable Task Information from an Unstable Neural Population

Source code files for the paper "Stable Task Information from an Unstable Neural Population". 

Contact LD and/or CDH to obtain the datasets necessary to reproduce these analyses.

`Version 1/Loback'` contains ARL's Matlab source code. This was used to render Figure 2, as well as in prelimenary analyses for Figures 3 and 4. There is also some plotting commands in Python. 

`Version 2/Rule/` contains MER's analysis code used to prepare figures 3, 4, 5, and the supplemental figures, in the form of iPython notebooks and Python source files. A frozen copy of the [neurotools](https://github.com/michaelerule/neurotools) Python library is included. This has a lot of dependencies, but only a limited subset of the library is actually needed to run the analyses. 

Repeating all analysis from scratch should take about a week on a high-performance laptop. Many data processing steps are shared across the notebooks. The code is set up to cache intermediate results in `.npz` files in a local directory called `PPC_cache`. You will need about 40 gigabytes of storage to hold everything. You can also contact MER to obtain a cloned copy of the `PPC_cache` directory to speed things along (it is too large to host on Github). 

This is un-optimized research code reflecting a subset of rapidly developed exploratory analyses that were included in the paper. If you are interested in reproduction, re-implementing similar or complementary analyses on the raw data may be better than attempting to run everything here from scratch.

`Version 2/Raman/` contains DVR's simulation code in Julia for sampling from the null model. 
