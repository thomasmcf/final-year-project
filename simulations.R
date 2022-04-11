source("simulation functions.R")

# Read the transect file and extract the points
transect_list <- read.traps(file = "transects.txt", detector = "transect")
x <- transect_list$x
y <- transect_list$y

x0 <- x[seq(1, 22, 2)]
y0 <- y[seq(1, 22, 2)]
x1 <- x[seq(2, 22, 2)]
y1 <- y[seq(2, 22, 2)]

# Make a mask for secr to use
mask <- make.mask(traps = transect_list, buffer = 0, spacing = 10)

# Times of visits at thorough search
times <- c(36500, 36530, 36560)
t0 <- 36500 - 30

# Begin parallelisation
cl <- makeCluster(20)
registerDoParallel(cl)
getDoParWorkers()

# Simulation populations under a range of abundances
N_expected <- 1:200
B <- length(N_expected)

output <- foreach(i = 1:B, .packages = c("tidyverse", "secr", "survival"), .combine = rbind) %dopar% {
  N_exp <- N_expected[i]
  pop <- prune(population(N_exp, 365 * 101), 36500 - 30)
  
  capthist <- survey_sim(x0, y0, x1, y1, pop, times = times, t0=36500-30)
  mod <- secr::secr.fit(capthist = capthist$capthist, mask = mask)
  
  data.frame(expected = N_exp,
             present = pres_det(pop, 36500)$pres,
             detectable = pres_det(pop, 36500)$det,
             estimated = region.N(mod)[2, 1]
  )
}

saveRDS(output, "exp_pres_det_est.rds")



# Simulate populations with different mean dropping life time
mean_lifes <- 7:365
B <- length(mean_lifes)

output <- foreach(i = 1:B, .packages = c("tidyverse", "secr", "survival"), .combine = rbind) %dopar% {
  mean_life <- mean_lifes[i]
  pop <- population(30, 365 * 101, mean_life = mean_life)
  
  detectable_times <- sapply(pop,
                             function(ind) max(ind$decay_times) - min(ind$drop_times))
  
  mean_detectable_time <- mean(detectable_times)
  se <- sd(detectable_times) / sqrt(length(detectable_times))
  
  slope <- mean_detectable_time / (365 + 2)
  
  pop <- prune(pop, 36500 - 30)
  capthist <- survey_sim(x0, y0, x1, y1, pop, times = times, t0=36500-30)
  mod <- secr::secr.fit(capthist = capthist$capthist, mask = mask)
  
  data.frame(estimated = region.N(mod)[2, 1],
             present = pres_det(pop, 36500)$pres,
             detectable = pres_det(pop, 36500)$det,
             mean_life = mean_life,
             mean_detectable_time = mean_detectable_time,
             se = se,
             slope = slope
  )
}

saveRDS(output, "dropping_life.rds")



# Simulate populations with different dropping rates
dropping_rate <- seq(1/7, 2, length.out = 200)
B <- length(dropping_rate)

output <- foreach(i = 1:B, .packages = c("tidyverse", "secr", "survival"), .combine = rbind) %dopar% {
  drop_rate <- dropping_rate[i]
  pop <- population(30, 365 * 101, drop_rate = drop_rate)
  
  detectable_times <- sapply(pop,
                             function(ind) max(ind$decay_times) - min(ind$drop_times))
  
  mean_detectable_time <- mean(detectable_times)
  se <- sd(detectable_times) / sqrt(length(detectable_times))
  
  slope <- mean_detectable_time / (365 + (1 / drop_rate))
  
  pop <- prune(pop, 36500 - 30)
  capthist <- survey_sim(x0, y0, x1, y1, pop, times = times, t0=36500-30)
  mod <- secr::secr.fit(capthist = capthist$capthist, mask = mask)
  
  data.frame(drop_rate = drop_rate,
             mean_detectable_time = mean_detectable_time,
             se = se,
             slope = slope,
             estimated = region.N(mod)[2, 1],
             present = pres_det(pop, 36500)$pres,
             detectable = pres_det(pop, 36500)$det
  )
}

saveRDS(output, "drop_rate.rds")



# Simulate populations with different average length of stay in population
mean_stays <- seq(365/2, 365 * 1.5, length.out = 200)
B <- length(mean_stays)

output <- foreach(i = 1:B, .packages = c("tidyverse", "secr", "survival"), .combine = rbind) %dopar% {
  mean_stay <- mean_stays[i]
  pop <- population(30, 365 * 101, mean_stay = mean_stay)
  
  detectable_times <- sapply(pop,
                             function(ind) max(ind$decay_times) - min(ind$drop_times))
  
  mean_detectable_time <- mean(detectable_times)
  se <- sd(detectable_times) / sqrt(length(detectable_times))
  
  slope <- mean_detectable_time / (mean_stay + 2)
  
  pop <- prune(pop, 36500 - 30)
  capthist <- survey_sim(x0, y0, x1, y1, pop, times = times, t0=36500-30)
  mod <- secr::secr.fit(capthist = capthist$capthist, mask = mask)
  
  data.frame(mean_stay = mean_stay,
             mean_detectable_time = mean_detectable_time,
             se = se,
             slope = slope,
             estimated = region.N(mod)[2, 1],
             present = pres_det(pop, 36500)$pres,
             detectable = pres_det(pop, 36500)$det
  )
}

saveRDS(output, "mean_stay.rds")


# Estimate mean dropping life time for different distributions
# Exponential
B <- 200
t0 <- 36500 - 30
times3 <- 36500 + 0:2 * 30
times12 <- 36500 + 0:11 * 30

output <- foreach(i = 1:B, .packages = c("tidyverse", "secr", "survival"), .combine = rbind) %dopar% {
  pop <- prune(population(30, 36500), t0)
  
  survhist3 <- survey_sim(x0, y0, x1, y1, pop, times3, t0)$survhist
  survhist12 <- survey_sim(x0, y0, x1, y1, pop, times12, t0)$survhist
  
  mod3 <- survreg(survhist3 ~ 1, dist = "exponential")
  mod12 <- survreg(survhist12 ~ 1, dist = "exponential")

  mod3.contains <- exp(confint(mod3)[1]) <= 365/2 & 365/2 <= exp(confint(mod3)[2])
  mod12.contains <- exp(confint(mod12)[1]) <= 365/2 & 365/2 <= exp(confint(mod12)[1])
  
  data.frame(mean3 = exp(coef(mod3)),
             mean12 = exp(coef(mod12)),
             contains3 = mod3.contains,
             contains12 = mod12.contains
  )
}

saveRDS(output, "survival sim exp.rds")

# Weibull
output <- foreach(i = 1:B, .packages = c("tidyverse", "secr", "survival"), .combine = rbind) %dopar% {
  pop <- prune(population(30, 36500), t0)
  
  survhist3 <- survey_sim(x0, y0, x1, y1, pop, times3, t0)$survhist
  survhist12 <- survey_sim(x0, y0, x1, y1, pop, times12, t0)$survhist
  
  mod3 <- survreg(survhist3 ~ 1, dist = "weibull")
  mod12 <- survreg(survhist12 ~ 1, dist = "weibull")
  
  mod3.confint <- predict(mod3, type = "quantile", p = c(0.025, 0.975))[1, ]
  mod12.confint <- predict(mod3, type = "quantile", p = c(0.025, 0.975))[1, ]
  
  mod3.contains <- mod3.confint[1] <= 365/2 & 365/2 <= mod3.confint[2]
  mod12.contains <- mod12.confint[1] <= 365/2 & 365/2 <= mod12.confint[1]

  data.frame(mean3 = predict(mod3)[1],
             mean12 = predict(mod12)[1],
             contains3 = mod3.contains,
             contains12 = mod12.contains
  )
}

saveRDS(output, "survival sim wei.rds")

# Log logistic
output <- foreach(i = 1:B, .packages = c("tidyverse", "secr", "survival"), .combine = rbind) %dopar% {
  pop <- prune(population(30, 36500), t0)
  
  survhist3 <- survey_sim(x0, y0, x1, y1, pop, times3, t0)$survhist
  survhist12 <- survey_sim(x0, y0, x1, y1, pop, times12, t0)$survhist
  
  mod3 <- survreg(survhist3 ~ 1, dist = "loglogistic")
  mod12 <- survreg(survhist12 ~ 1, dist = "loglogistic")
  
  mod3.confint <- predict(mod3, type = "quantile", p = c(0.025, 0.975))[1, ]
  mod12.confint <- predict(mod3, type = "quantile", p = c(0.025, 0.975))[1, ]
  
  mod3.contains <- mod3.confint[1] <= 365/2 & 365/2 <= mod3.confint[2]
  mod12.contains <- mod12.confint[1] <= 365/2 & 365/2 <= mod12.confint[1]
  
  data.frame(mean3 = predict(mod3)[1],
             mean12 = predict(mod12)[1],
             contains3 = mod3.contains,
             contains12 = mod12.contains
  )
}

saveRDS(output, "survival sim log.rds")


# Now implement the correction
B <- 200
df <- readRDS("exp_pres_det_est.rds")

 # Initial values for Markov chain
inits <- function(){
  list(beta0 = rnorm(1),
       beta1 = rnorm(1),
       sig2 = runif(1),
       m1 = runif(1))
}

# parameters to monitor
monitor <- "corrected"

nc <- 3 # Three chains
nb <- 10000 # burn in iterations
ni <- 15000 + nb # samples for inference
nt <- 1 # no thinning


output <- foreach(i = 1:B, .packages = c("tidyverse", "secr", "survival", "jagsUI"), .combine = rbind) %dopar% {
  pop <- prune(population(i, 36500), 36500 - 30)
  capthist <- survey_sim(x0, y0, x1, y1, pop, times = times, t0=36500-30)$capthist
  mod <- secr.fit(capthist = capthist, mask = mask)
  
  jagsdat <- list(N = nrow(df),
                  present = df$present,
                  estimated = df$estimated,
                  point_est = region.N(mod)[2, 1],
                  std_error = region.N(mod)[2, 2])
  
  out <- jags(data = jagsdat,
              inits = inits,
              parameters.to.save = monitor,
              model.file = "model.txt",
              n.chains = nc,
              n.iter = ni,
              n.burnin = nb,
              n.thin = nt)
  
  data.frame(uncorrected = jagsdat$point_est,
             present = pres_det(pop, 36500)$pres,
             mean = out$mean$corrected,
             median = out$q50$corrected,
             lcl = out$q2.5$corrected,
             ucl = out$q97.5$corrected
  )
}

saveRDS(output, "solution fixed.rds")
stopCluster(cl)
