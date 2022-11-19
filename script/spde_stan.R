rm(list=ls())

setwd("C:/Users/joaquin/Desktop/joaquin/spde_stan")


library(rstan)
rstan_options(auto_write=TRUE)
library(TMB)
library(tmbstan)
library(INLA)
library(fields)
library(shinystan)

library(parallel)
# Calculate the number of cores
no_cores <- detectCores() - 1


#Compile and load c++ code-------
compile("spde.cpp")
dyn.load(dynlib("spde"))
#--------------------------------

#Read data-----------------------
map = read.table("Leuk.map")
data(Leuk, package = "INLA")

library(dplyr)
#Leuk <- sample_n(Leuk, 200)


#--------------------------------

#Define mesh and components representing the  precision matrix----
loc = cbind(Leuk$xcoord, Leuk$ycoord)
boundary = INLA::inla.nonconvex.hull(loc)
boundary2 = INLA::inla.nonconvex.hull(loc,convex = -0.35)
mesh = INLA::inla.mesh.2d(
  loc=loc,
  boundary = list(boundary,boundary2),
  max.edge=c(0.05, 0.2),
  cutoff=0.05
)
A = inla.spde.make.A(mesh,loc)
spde = inla.spde2.matern(mesh, alpha=2)
spdeMatrices = spde$param.inla[c("M0","M1","M2")]


kappa = exp(0)
Q = kappa^4*spdeMatrices$M0 + 2.0*kappa^2*spdeMatrices$M1 + spdeMatrices$M2;



#-------------------------------------------------------------------

#Define the design matrix for the fixed effects-------
X <- model.matrix( ~ 1 + sex + age + wbc + tpi, data = Leuk)
#-----------------------------------------------------

#Define the data and parameters given to TMB----------
data <- list(time       = Leuk$time,
             notcens    = Leuk$cens,
             meshidxloc = mesh$idx$loc - 1,
             A = A,
             X          = as.matrix(X),
             #Q = Q)
             spdeMatrices = spdeMatrices)
             
parameters <- list(beta      = c(0.0,0,0,0,0),
                   log_tau   = 0,
                   log_kappa = 0,
                   log_omega = -1,
                   x         = rep(0.0, nrow(data$spdeMatrices$M0)))
#-----------------------------------------------------

#Estimating the model and extract results-------------
#data$flag = 1
startTime <- Sys.time()
obj <- MakeADFun(data, parameters, random="x", DLL="spde")
#obj <- normalize(obj, flag="flag")
opt <- nlminb(obj$par, obj$fn, obj$gr)
rep <- sdreport(obj)
endTime <- Sys.time()
timeUsed = endTime - startTime
print(timeUsed)
#-----------------------------------------------------

#==================
#  tmbstan model
#==================
fit1 <- tmbstan(obj, iter=1000, chains = 2, warmup=200, 
                control = list(max_treedepth = 10, adapt_delta = 0.8), # Stan parameters algorithm
                init='last.par.best', 
                cores=no_cores, open_progress=F)

saveRDS(fit1, 'fit1.RDS')

fit1 <- readRDS('fit1.RDS')
launch_shinystan(fit1)


#==================
#    Stan model
#==================
data2 <- list(n = length(Leuk$time),
              p = ncol(X),
              row_spar = dim(spdeMatrices$M0)[1],
              col_spar = dim(spdeMatrices$M0)[2],
              time=Leuk$time, 
              notcens=Leuk$cens,
              X=as.matrix(X),
              #Q=as.matrix(Q),
              M0=as.matrix(spdeMatrices$M0),
              M1=as.matrix(spdeMatrices$M1),
              M2=as.matrix(spdeMatrices$M2),
              A=as.matrix(A))

init.fn <- function()
  split(unname(obj$env$last.par.best),names(obj$env$last.par.best))

fit2 <- stan(file='stan_spde.stan', data = data2, 
             iter=1000, chains = 2, warmup=200, control = list(max_treedepth = 10, adapt_delta = 0.8), 
             cores = no_cores, open_progress=FALSE)

saveRDS(fit2, 'fit2.RDS')
fit2 <- readRDS('fit2.RDS')
launch_shinystan(fit2)





#=================================
#      tmbstan vs Stan 
#=================================

x1 <- as.data.frame(fit1)
x2 <- as.data.frame(fit2)

qqplot(x1[,'lp__'],x2[,'lp__'])
qqplot(x1[,'log_tau'],x2[,'log_tau'])

sum(get_elapsed_time(fit1))/60 # tmbstan
sum(get_elapsed_time(fit2))/60 # stan


kappa = 1.1

Q = kappa^4*spdeMatrices$M0 + 2.0*kappa^2*spdeMatrices$M1 + spdeMatrices$M2
class(Q)
