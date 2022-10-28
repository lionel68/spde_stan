# simulation SPDE Stan vs INLA

####### note 28/10/2022 ######
# the object st_model is not defined
# this should be checked but I guess that it is simply:
# brms::make_standcode that is saved in a rds object

setwd("~/Dokumente/spde/")
# load libraries
library(RandomFields)
library(brms)
library(rstan)
library(INLA)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(tidyverse)

# load helper functions
source("mgcv_spde_smooth.R") 


gen_data <- function(n = 100, range = 0.5, sigma_spat = 1, sigma_res = 1, mu = 0){
  # the location of the observations
  dat <- data.frame(x = runif(n),
                    y = runif(n))
  
  # simulate the spatial field
  dat$spat <- RFsimulate(RMmatern(1, var = sigma_spat, scale = range), x = dat$x,
                         y = dat$y, spConform = FALSE)
  
  # simulate the response
  dat$resp <- rnorm(n, mu + dat$spat, sigma_res)
  
  return(dat)
}

# function to run the simu
run_simu <- function(n = 100, range = 0.5, sigma_spat = 1, sigma_res = 1, mu = 0){
  # generate data
  dat <- gen_data(n, range, sigma_spat, sigma_res, mu)
  # create meshes and A matrix
  bnd <- inla.mesh.segment(matrix(c(0, 0, 1, 0, 1, 1, 0, 1), ncol = 2, byrow = TRUE))
  mesh <- inla.mesh.2d(boundary = bnd, max.edge = c(0.1, 10))
  sm_term <- s(x, y, bs = "spde", xt = mesh)
  sm <- mgcv::smooth.construct(sm_term, dat)
  # the projection matrix
  bs_mat <- mgcv::PredictMat(sm, data = dat[,c("x", "y")])
  
  # generate the data to pass to Stan
  stan_dat <- make_standata(resp ~ 1 + s(x, y), dat)
  stan_dat$Zs_1_1 <- bs_mat
  stan_dat$knots_1[1] <- ncol(bs_mat)
  
  # fit the Stan model
  m_stan <- sampling(st_model,
                 data = stan_dat,
                 pars = "epred",
                 open_progress = FALSE)
  
  # get the posterior average prediction
  stan_pred <- posterior_summary(rstan::extract(m_stan, "epred")$epred)
  
  # prep for the INLA model
  spde <- inla.spde2.pcmatern(mesh=mesh,
                              prior.range=c(0.1, 0.05),
                              prior.sigma=c(10, 0.05))
  # setup estimation stack
  A <- inla.spde.make.A(mesh, as.matrix(dat[,c("x", "y")]))
  intercept <- rep(1, nrow(dat))
  stk <- inla.stack(tag='est', ## tag
                    data=list(resp=dat$resp), ## response
                    A=list(A, 1), ## two projector matrix
                    effects=list(s=1:spde$n.spde, intercept = intercept))
  # model formula
  formula <- resp ~ 0 + intercept + f(s, model=spde)
  # fit with INLA
  m_inla <- inla(formula,
                 data=inla.stack.data(stk),
                 control.predictor=list(A = inla.stack.A(stk),
                                        compute=TRUE),
                 control.compute=list(config = TRUE))
  
  # extract model predictions
  ind <- inla.stack.index(stk, tag = "est")$data
  inlasamp <- inla.posterior.sample(4000, m_inla)
  inlapred <- sapply(inlasamp, FUN = function(x){x$latent[ind]})
  inla_pred <- posterior_summary(t(inlapred), probs = NA)
  
  # get correlations
  out_dat <- data.frame(type = c("data - inla",
                                 "inla - stan",
                                 "data - stan"),
                        r2 = c(cor(dat$resp, inla_pred[,1]) ** 2,
                               cor(inla_pred[,1], stan_pred[,1]) ** 2,
                               cor(dat$resp, stan_pred[,1]) ** 2),
                        range = range, sigma_spat = sigma_spat, n = n,
                        sigma_res = sigma_res, mu = mu)
  return(out_dat)
}

# run this across a range of sigma_spat and range values
param <- expand.grid(sigma_spat = c(0.1, 1, 10),
                     range = c(0.1, 1, 10),
                     sigma_res = 1, mu = 0, n = c(100, 250, 500))
# create a clean version of the model
st_model <- stan_model(file = "stan_smooth.stan")

rr_out <- NULL
n_rep <- 30
for(i in 1:nrow(param)){
  rr <- replicate(n_rep, run_simu(range = param[i, "range"],
                              sigma_spat = param[i, "sigma_spat"],
                              sigma_res = param[i, "sigma_res"],
                              mu = param[i, "mu"],
                              n = param[i, "n"]),
                  simplify = FALSE)
  rr2 <- plyr::rbind.fill(rr)
  rr2$iter <- rep(1:n_rep, each = 3)
  rr_out <- rbind(rr_out, rr2)
  print(paste0("Param set : ", i, " from ", nrow(param), " finished!\n"))
}

saveRDS(rr_out, "spde_simu.rds")

# plot the results
ss <- readRDS("spde_simu.rds")
ss$nf <-  paste0("n : ", ss$n)
gg_out <- ggplot(ss, aes(x=nf, y=r2, color = type)) +
  geom_boxplot(position = "dodge2") +
  facet_grid(range~sigma_spat, labeller = label_both) +
  labs(x = "",
       y = "R² between predicted and true data") +
  theme(axis.text.x = element_text(angle = 90))

ggsave("r2_out.png", gg_out)
