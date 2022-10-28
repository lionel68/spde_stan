rm(list=ls())

setwd("C:/Users/Usuario/Desktop/Projects/2021/KAUST/TMB/codes")

## Libraries used
library(TMB)
library(INLA)
library(tmbstan)
library(rstan)
library(RandomFields)
library(bayesplot)
library(tictoc)
library(ggplot2)
library(geoR)


library(parallel)
# Calculate the number of cores
no_cores <- detectCores() - 1

#============================================
#           Load the TMB model
TMB::compile('GRF1.cpp')
dyn.load( dynlib("GRF1") )
#============================================

#============================================
#         Mesh
#============================================
set.seed(201803)
inla.seed = sample.int(n=1E6, size=1)
options(width=70, digits=3)

sim_loc = matrix(c(0,0, 1, 1, 0, 1, 1, 0), nrow = 4, byrow = T)

mesh_sim = inla.mesh.2d(loc = sim_loc, max.edge=c(0.5, 1))
mesh_sim$n
plot(mesh_sim)


# Initial values for tau0 and kappa0
sigma_s0 <- 1
size <- min(c(diff(range(mesh_sim$loc[, 1])), diff(range(mesh_sim$loc[, 2]))))
rho0 <- size/5
kappa0 <- sqrt(8)/rho0
tau0 <- 1/(sqrt(4 * pi) * kappa0 * sigma_s0)

spde = inla.spde2.pcmatern(mesh_sim, alpha = 2, prior.range = c(0.5, 0.5), prior.sigma = c(0.5, 0.5))

Q = inla.spde2.precision(spde, theta=c(log(tau0), log(kappa0)))
omega_s = inla.qsample(n=1, Q = Q, seed = inla.seed)
omega_s = omega_s[ ,1]

# Simula data
n = 500
loc.data = matrix(runif(2*n), n)*1 # coordinates
plot(mesh_sim, main = "")
points(loc.data, col='red', pch = 19)
title(main = "n = 500 locations simulated",
      cex.main = 1.8, col.main= "blue")



A = inla.spde.make.A(mesh=mesh_sim, loc=loc.data)
omega_s = drop(A %*% omega_s)
# 
# 
# 
# Create the linear predictor
beta0 <- 0
beta1 <- 1
sigma_iid <- 1

x = runif(n)-0.5

#lin.pred <- beta0 + beta1*x + omega_s
lin.pred <- beta0 + omega_s

# Response variable
n <- length(lin.pred)
y <- lin.pred + rnorm(n, 0, sigma_iid)



######## ######## SPDE-based
# Build object
Data = list(y_i = as.vector(y),
#            x_1 = as.vector(x),
            M0 = spde$param.inla$M0, 
            M1 = spde$param.inla$M1, 
            M2 = spde$param.inla$M2,
            Aproj = A,
            tau0   = tau0,
            kappa0 = kappa0)

Params = list(beta0 = 0.5,
#              beta1 = 0.5,
              logsigma  = exp(0.1),  
              logtau    = exp(0.1),
              logkappa  = exp(0.1),
              omega_s = rep(0, mesh_sim$n))
              #omega_s =   rep(0, nrow(Q)))


Obj = MakeADFun(data = Data, parameters = Params, random="omega_s", DLL="GRF1" )

# Optimize
#Opt_spde = TMBhelper::Optimize( obj=Obj, newtonsteps=1, bias.correct=TRUE )
opt = nlminb(Obj$par,Obj$fn, Obj$gr)
h_spde = Obj$env$spHess(random=TRUE)
report_spde = Obj$report()

AIC = 2*opt$objective +2*length(opt$par)
AIC





#========================
#      INLA model
#========================
library(brinla)

df = data.frame(y=y, locx=loc.data[ ,1], locy=loc.data[ ,2], x = x)
summary(df)

mesh = inla.mesh.2d(loc = sim_loc, max.edge=c(0.4, 1))

A = inla.spde.make.A(mesh = mesh, loc=data.matrix(df[ , c('locx', 'locy')]))
dim(A)

plot(mesh, asp=1)
points(df[ , c('locx', 'locy')], col='red', lwd=.1)

prior.median.sd = 0.5
prior.median.range = 0.5
spde = inla.spde2.pcmatern(mesh, alpha = 2, prior.range = c(prior.median.range, .5), prior.sigma = c(prior.median.sd, .5))

stack = inla.stack(tag='est',
                   data=list(y=df$y),
                   effects=list(s=1:spde$n.spde, 
                                data.frame(intercept=1, x=df$x)),
                                A=list(A, 1))


family = "gaussian"
prior.median.sd.g = 0.5 # prior median for sigma.epsilon
control.family = list(hyper = list(prec = list(prior = "pc.prec", param = c(prior.median.sd.g, 0.5))))


formula = y ~ -1 + intercept +  x + f(s, model=spde)

initial.theta = c(-0.0658, -1.7225,  0.7257)
# - the first time you run this, set it to NULL
# - after running, set it to res$internal.summary.hyperpar$mean
# - and run the code again
# - Reason: Both faster and better inference

res = inla(formula, data=inla.stack.data(stack),
           family = family,
           control.family = control.family,
           control.predictor=list(A = inla.stack.A(stack)),
           quantiles=c(0.5, 0.025, 0.975, 0.1, 0.9, 0.25, 0.75),
           #control.compute = list(config=T, dic=T, cpo=T, waic=T), 
           # - Model comparisons
           #control.inla = list(int.strategy='grid'),
           # - More accurate integration over hyper-parameters
           control.mode = list(restart = T, theta = initial.theta),
           verbose = TRUE)

summary(res)
#res$internal.summary.hyperpar$mean


# Compare with TMB estimation
res$summary.fixed[1]
opt$par[1:2]


# sd Gaussian observations
res$summary.hyperpar
report_spde$sigma

# Range for the spatial field
report_spde$rho

# Std for the spatial field
report_spde$sigma_s




exp(opt$par[3:5])


sigmarun <- inla.tmarginal(function(x) 1/sqrt(exp(x)),res$internal.marginals.hyperpar[[2]])
sigmapos <- inla.tmarginal(function(x) 1/sqrt(exp(x)),res$internal.marginals.hyperpar[[3]])
sigmaepsilon <- inla.tmarginal(function(x) 1/sqrt(exp(x)),res$internal.marginals.hyperpar[[1]])
restab=sapply(res$marginals.fixed, function(x) inla.zmarginal(x,silent=TRUE))
restab=cbind(restab, inla.zmarginal(sigmarun,silent=TRUE))
restab=cbind(restab, inla.zmarginal(sigmapos,silent=TRUE))
restab=cbind(restab, inla.zmarginal(sigmaepsilon,silent=TRUE))
colnames(restab) = c("intercept","beta1","run","position","epsilon")
data.frame(restab)














#===========================================================================================
#                                        Fit with tmbstan
#===========================================================================================

tic("Time of estimation")
fit = tmbstan(Obj, chains = 3, open_progress = FALSE, 
              init ='last.par.best', control = list(max_treedepth = 5, adapt_delta = 0.7), 
              iter = 2000, warmup=1000, seed=483892929)
toc()


c_light <- c("#DCBCBC")
c_light_highlight <- c("#C79999")
c_mid <- c("#B97C7C")
c_mid_highlight <- c("#A25050")
c_dark <- c("#8F2727")
c_dark_highlight <- c("#7C0000")

traceplot(fit, pars=names(Obj$par), inc_warmup=TRUE)

Divergences
divergent = get_sampler_params(fit, inc_warmup=FALSE)[[1]][,'divergent__']
sum(divergent)

## Methods provided by 'rstan'
class(fit)
methods(class ="stanfit")

#launch_shinystan(fit)

## ESS and Rhat from rstan::monitor
mon = monitor(fit)
max(mon$Rhat)
min(mon$Tail_ESS)

# evalaute problem of convergence
sum(mon$Rhat > 1.01)
sum(mon$Tail_ESS < 400)

source('monitornew.R')
source('monitorplot.R')
source('stan_utility.R')

which_min_ess = which.min(mon[1:200, 'Tail_ESS'])
plot_local_ess(fit = fit, par = which_min_ess, nalpha = 10)

plot_quantile_ess(fit = fit, par = which_min_ess, nalpha = 50)

plot_change_ess(fit = fit, par = which_min_ess)

check_rhat(fit)
check_treedepth(fit, 10)
check_energy(fit)  # s?lo es una recomendaci?n muy aproximada y con estudios preliminares
check_div(fit)


#=====================================================
# rho comparison for Stan, TMb and INLA
# Stan

posterior_1 <- as.matrix(fit)

Obj$report(posterior_1[1,-ncol(posterior_1)])         # sd0 is only element
rho1 <- rep(NA, len=nrow(posterior_1))
for(i in 1:nrow(posterior_1)){
   r1 <- Obj$report(posterior_1[i,-ncol(posterior_1)])
   rho1[i] <- r1$rho
}
mean(rho1)

# TMB
Obj$report()$rho

# INLA
res$summary.hyperpar[2, 1]



#=====================================================
# sigma_s comparison for Stan, TMb and INLA
# Stan
Obj$report(posterior_1[1,-ncol(posterior_1)])         # sd0 is only element
sigma_s1 <- rep(NA, len=nrow(posterior_1))
for(i in 1:nrow(posterior_1)){
   r1 <- Obj$report(posterior_1[i,-ncol(posterior_1)])
   sigma_s1[i] <- r1$sigma_s
}
mean(sigma_s1)

# TMB
Obj$report()$sigma_s

# INLA
res$summary.hyperpar[3, 1]


