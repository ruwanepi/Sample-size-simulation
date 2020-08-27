library(dplyr)
library(ggplot2)
library(plyr)
library(glmmADMB)
library(lme4)

#Derive the parameters
n <- 500
(x1_delay <- rnorm(n, mean=3, sd=1)) #derive a range of delays
(x2_population <- rnorm(n, mean=500, sd=100)) #derive a range of population sizes
(ring <- 1:500)

#Simulation: a mixed-effects regression model with random effects
#that models the effect of delay (main exposure) and effect size
#on epidemic size (main outcome)
(b0 <- 0) #intercept
(b1 <- 0.4) #x1_delay
(b2 <- 0.01) #x2_population
(b3 <- 0.4) #interaction between x1 and x2

y <- b0 + (b1*x1_delay) + (b2*x2_population) + (b3*x1_delay*x2_population) + rnorm(n, mean=0, sd=1)

model <- glmer(y~(1 + x1_delay*x2_population) + (1 | ring), family = poisson(link = "log"), nAGQ = 100)
(summary(model))

output <- summary(model)$coefficients
coefs <- output[, 1]
ps <- output[, 4]

rsq <- summary(model)$r.squared
results <- c(coefs, ps, rsq)
results

#FUNCTION
regression_sim <- function(simNum, n, b0, b1, x1_delay_mean=3, x1_delay_sd=1, err_mean=0, err_sd=1) {
  #b0(intercept),b1(delay),b2(pop),b3(delayxpop) 
  x1_delay <- rnorm(n, mean=3, sd=1) 
  x2_population <- rnorm(n, mean=500, sd=100) 
  ring <- 1:500
  Effsize <- rnorm(n, mean=0.4, sd=0.1) #derive range of effect sizes
  
  y <- b0 + (b1*x1_delay) + (b2*x2_population) + (b3*x1_delay*x2_population) + 
    rnorm(n, mean=err_mean, sd=err_sd)
  
  model <- glmer(y~(1 + x1_delay*x2_population) + (1 | ring), family = poisson(link = "log"), nAGQ = 100)
  (summary(model))
  
  output <- summary(model)$coefficients
  coefs <- output[, 1]
  ps <- output[, 4]
  
  rsq <- summary(model)$r.squared
  results <- c(coefs, ps, rsq)
  names(results) <- c('b0_coef', 'b1_coef',  
                      'b0_p','b1_p', 'rsq')
  return(results)
}

#Check that the function works
regression_sim(1, n=500, b0=0, b1=0.4)

#Run simulation 1000 times
num_sims <- 1000
(sims <- ldply(1:num_sims, regression_sim, n=500, b0=0, b1=0.4))
(power <- sum(sims$b1_p < .05) / nrow(sims))


