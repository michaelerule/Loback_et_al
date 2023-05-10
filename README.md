# Stable Task Information from an Unstable Neural Population

Rule, M. E., Loback, A. R., Raman, D. V., Driscoll, L. N., Harvey, C. D., & O'Leary, T. (2020). Stable task information from an unstable neural population. Elife, 9, e51121. [doi: https://doi.org/10.7554/eLife.51121](https://doi.org/10.7554/eLife.51121)

### Abstract

Over days and weeks, neural activity representing an animal’s position and movement in sensorimotor cortex has been found to continually reconfigure or ‘drift’ during repeated trials of learned tasks, with no obvious change in behavior. This challenges classical theories, which assume stable engrams underlie stable behavior. However, it is not known whether this drift occurs systematically, allowing downstream circuits to extract consistent information. Analyzing long-term calcium imaging recordings from posterior parietal cortex in mice (Mus musculus), we show that drift is systematically constrained far above chance, facilitating a linear weighted readout of behavioral variables. However, a significant component of drift continually degrades a fixed readout, implying that drift is not confined to a null coding space. We calculate the amount of plasticity required to compensate drift independently of any learning rule, and find that this is within physiologically achievable bounds. We demonstrate that a simple, biologically plausible local learning rule can achieve these bounds, accurately decoding behavior over many days.

### Repository Contents

This repository contains python, matlab, and julia source code files for the statistical analyses in ["Stable Task Information from an Unstable Neural Population". ](https://www.biorxiv.org/content/10.1101/796334v2) You will need the data from [Driscoll et al. (2017)](https://www.sciencedirect.com/science/article/pii/S0092867417308280) to run these analyses. Contact LND or CDH to obtain the datasets necessary to reproduce these analyses. These datasets are also available from the [Dryad repository](https://doi.org/10.5061/dryad.gqnk98sjq).

 - `Version 1/Loback` contains ARL's Matlab source code. This was used to render Figure 2, as well as in prelimenary analyses for Figures 3 and 4. There are also some plotting scipts in Python. 
 - `Version 2/Rule/` contains MER's analysis code used to prepare figures 3, 4, 5, and the supplemental figures, in the form of iPython notebooks and Python source files. A copy of the [neurotools](https://github.com/michaelerule/neurotools) Python library is included. This has a lot of dependencies, but only a limited subset of the library is needed for these analyses. Regenerating all figures from scratch should take about a week on a high-performance laptop. To avoid recomputation, intermediate results are cached in `PPC_cache`.
 - `Version 2/Raman/` contains DVR's simulation code in Julia for sampling from the null model. 
