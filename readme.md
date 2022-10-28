# Fitting SPDE models with brms / Stan

This is a starting repo with some unchecked Rcode to fit SPDE models in brms and stan using the spde smooths created in the Miller paper.

## Repo orga

* figures: figure outputted by the code, right now there is only one figure comparing R2 from inla and stan models fitted via the script/spde_simu.r file
* robject: rds files saved after simulation, right now there is only the rds file generated from script/spde_simu.r
* script: some Rscript, 
    * spde_simu.r is a simulation of different data with varying spatial signal and R2 between simulated data and predicted data from inla and from a stan model with a 'plugged-in' spde smooths generated from the Miller code and plugged into the model stan code generated from brms::make_stancode
    * spde_stan.Rmd is a first try at trying different approaches with Stan and brms to fit spde models.