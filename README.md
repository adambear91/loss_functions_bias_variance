# Summary

This project contains code showing how loss functions modulate the optimal bias-variance tradeoff. 

# Components

This project contains two main files: `theoretical_sims.R` and `case_study.Rmd`. The first file covers the simulation results presented in Bear & Cushman (2020), *Proceedings of the Forty-Second Annual Conference of the Cognitive Science Society*, which assume Gaussian error. The second file covers analyses done with the *k*-nearest neighbors and lasso regression algorithms for both a synthetic data set and a data set of joke ratings (see below). Note that some of the simulation code in this file takes a long time to run. For convenience, the output data files are provided in the "Data" folder. Outputted figures can be found in the "Figures" folder.

All custom functions written for this project are contained in `utilities.R`. 

# References

For the data set of joke ratings, we used the `jesterFull.Rdata` file provided in the [code](https://osf.io/gscvh/) for the following paper:

Analytis, P. P., Barkoczi, D., & Herzog, S. M. (2018). Social learning strategies for matters of taste. *Nature Human Behavior*, 2, 415–424. 

# Questions?

Please email Adam Bear ([adambear@fas.harvard.edu](mailto:adambear@fas.harvard.edu)).
