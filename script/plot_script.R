# further look at the spatial smooth simu

library(tidyverse)

simu <- readRDS("Z:/Dokumente/spde/spde_simu.rds")

# remove comp inla-stan and compute diff
simu %>%
  filter(type != "inla - stan") %>%
  group_by(iter, range, sigma_spat, n) %>%
  mutate(diff = r2 - lag(r2)) %>%
  filter(!is.na(diff)) %>%
  mutate(n = paste0("n : ", n),
         range = paste0(" spatial range : ", range),
         sigma_spat = paste0("spatial variance : ", sigma_spat)) -> dd

ggplot(dd, aes(x = n, y = diff)) +
  geom_jitter(width = 0.1) +
  stat_summary(fun.data = "mean_cl_boot", color = "red") +
  facet_grid(range ~ sigma_spat) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Sample size",
       y = "Differences",
       title = " R-square difference against the real data\nbetween Stan and INLA (R-square Stan - R-square INLA)")
