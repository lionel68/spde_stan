rm(list=ls())

setwd("C:/Users/joaquin/Desktop/joaquin/spde_stan")


library(rstan)
rstan_options(auto_write=TRUE)
library(INLA)
library(fields)
library(shinystan)

library(parallel)
# Calculate the number of cores
no_cores <- detectCores() - 1

#=============
# Read data
map = read.table("Leuk.map")
data(Leuk, package = "INLA")


#================
# Build the mesh 
loc = cbind(Leuk$xcoord, Leuk$ycoord)
boundary = INLA::inla.nonconvex.hull(loc)
boundary2 = INLA::inla.nonconvex.hull(loc,convex = -0.2)
mesh = INLA::inla.mesh.2d(loc=loc, boundary = list(boundary,boundary2),
                          max.edge=c(0.05, 0.15), cutoff=0.03)
A = inla.spde.make.A(mesh,loc)
spde = inla.spde2.matern(mesh, alpha=2)
spde_matrices = spde$param.inla[c("M0","M1","M2")]

#================
# Model matrix
X <- model.matrix( ~ 1 + sex, data = Leuk)



#==================
#    Stan model
#==================

# Data
data  <- list(n = length(Leuk$time),
              p = ncol(X),
              row_spar = dim(spde_matrices$M0)[1],
              col_spar = dim(spde_matrices$M0)[2],
              time=Leuk$time, 
              notcens=Leuk$cens,
              X=as.matrix(X),
              M0=as.matrix(spde_matrices$M0),
              M1=as.matrix(spde_matrices$M1),
              M2=as.matrix(spde_matrices$M2),
              A=as.matrix(A))

fit <- stan(file='stan_spde.stan', data = data, 
            control = list(max_treedepth = 5, adapt_delta = 0.6),
            iter=1000, chains = 2, warmup = 500, 
            cores = no_cores, open_progress=FALSE)

launch_shinystan(fit)

