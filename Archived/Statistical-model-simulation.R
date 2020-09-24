library(dplyr)
library(ggplot2)
library(plyr)
library(glmmADMB)
library(lme4)


#Simulate large dataset with statistical model
n<-500 #set a sample size of rings (temporary variable)
x1_delay <- rnorm(n, mean=5, sd=2) #derive a range of delays
x2_population <- rnorm(n, mean=500, sd=50) #derive a range of population sizes
x3_time <- rnorm(n, mean=20, sd=2) #derive a range of observation times
x4_sites <- rnorm(n, mean=50, sd=5) #derive a range of numbers of sites

#Simulation: a mixed-effects regression model with random effects
#that models the effect of delay (main exposure) and effect size
#on epidemic size (main outcome)
b0 <- 0 #intercept
b1 <- 0.4 #x1_delay
b2 <- 0.01 #x2_population
b3 <- 0.01 #x3_time
b4 <- 0.01 #x4_sites
b5 <- 0.4 #interaction between x1 and x2
b6 <- 0.4 #interaction between x1 and x3
b7 <- 0.4 #interaction between x1 and x4

y <- b0 + (b1*x1_delay) + (b2*x2_population) + (b3*x3_time) + (b4*x4_sites) + 
  (b5*x1_delay*x2_population) + (b6*x1_delay*x3_time) + (b7*x1_delay*x4_sites) + 
  rnorm(n, mean=0, sd=1)

model <- lm(y~x1_delay*x2_population*x3_time*x4_sites)
(summary(model))

output <- summary(model)$coefficients
coefs <- output[, 1]
ps <- output[, 4]

rsq <- summary(model)$r.squared
results <- c(coefs, ps, rsq)
names(results) <- c('b0_coef', 'b1_coef', 'b2_coef', 'b3_coef', 
                    'b4_coef', 'b5_coef', 'b6_coef', 'b7_coef', 
                    'b8_coef', 'b9_coef', 'b10_coef','b11_coef',
                    'b12_coef','b13_coef','b14_coef','b15_coef', 
                    'b0_p','b1_p', 'b2_p', 'b3_p', 'b4_p', 'b5_p',
                    'b6_p','b7_p', 'b8_p', 'b9_p', 'b10_p', 'b11_p',
                    'b12_p','b13_p', 'b14_p', 'b15_p', 'rsq')
results

#FUNCTION
regression_sim <- function(simNum, n, b0, b1, b2, b3, b4, 
                           b5, b6, b7, x1_delay_mean=5, x1_delay_sd=2, err_mean=5, err_sd=2) {
#b0(intercept),b1(delay),b2(pop),b3(time),b4(sites),
#b5(delayxpop),b6(delayxtime),b7(popxtime),b8(delayxsites),b9(popxsites)
#b10(timexsites),b11(delayxpopxtime),b12(delayxpopxsites),b13(delayxtimesxsites)
#b14(popxtimexsites),b15(delayxpopxtimexsites)
  x1_delay <- rnorm(n, mean=5, sd=2) #derive a range of delays
  x2_population <- rnorm(n, mean=500, sd=50) #derive a range of population sizes
  x3_time <- rnorm(n, mean=20, sd=2) #derive a range of observation times
  x4_sites <- rnorm(n, mean=50, sd=5) #derive a range of numbers of sites
  Effsize <- rnorm(n, mean=0.4, sd=0.1) #derive range of effect sizes
  
  y <- b0 + (b1*x1_delay) + (b2*x2_population) + (b3*x3_time) + (b4*x4_sites) + 
    (b5*x1_delay*x2_population) + (b6*x1_delay*x3_time) + (b7*x1_delay*x4_sites) + 
    rnorm(n, mean=err_mean, sd=err_sd)
  
  model <- lm(y~x1_delay*x2_population*x3_time*x4_sites)
  summary(model)
  
  output <- summary(model)$coefficients
  coefs <- output[, 1]
  ps <- output[, 4]
  
  rsq <- summary(model)$r.squared
  
  results <- c(coefs, ps, rsq)
  names(results) <- c('b0_coef', 'b1_coef', 'b2_coef', 'b3_coef', 
                      'b4_coef', 'b5_coef', 'b6_coef', 'b7_coef', 
                      'b8_coef', 'b9_coef', 'b10_coef','b11_coef',
                      'b12_coef','b13_coef','b14_coef','b15_coef', 
                      'b0_p','b1_p', 'b2_p', 'b3_p', 'b4_p', 'b5_p',
                      'b6_p','b7_p', 'b8_p', 'b9_p', 'b10_p', 'b11_p',
                      'b12_p','b13_p', 'b14_p', 'b15_p', 'rsq')
  return(results)
}

#Check that the function works
regression_sim(1, n=1000, b0=0, b1=0.4, b2=0.01, b3=0.01, b4=0.01, b5=0.01, b6=0.01, b7=0.01)

#Run simulation 1000 times 
num_sims <- 1000
(sims <- ldply(1:num_sims, regression_sim, n=100, b0=0, b1=0.4, b2=0.01, b3=0.01, b4=0.01, 
              b5=0.01, b6=0.01, b7=0.01))
(power <- sum(sims$b5_p < .05) / nrow(sims))

#Vary the parameters
sample_sizes <- c(500, 1000, 1500, 2000)
results <- NULL
for (val in sample_sizes) {
  sims <- ldply(1:1000, regression_sim, n=100, b0=0, b1=0.4, b2=0.01, b3=0.01, b4=0.01, 
                b5=0.01, b6=0.01, b7=0.01)
  sims$n <- val  
  results <- rbind(results, sims)
}

results

power_ests <- results %>%
  group_by(n) %>%
  summarize(power=sum(b5_p < .05) / n())

ggplot(power_ests, aes(x=n, y=power)) +
  geom_point() +
  geom_line() +
  ylim(c(0, 1)) +
  theme_minimal()



#Try with just delay and population as the predictors
regression_sim2 <- function(simNum, n, b0, b1, b2, b3, x1_delay_mean=5, x1_delay_sd=2, x2_pop_mean=500, x2_pop_sd=50, err_mean=0, err_sd=1) {
  #b0(intercept),b1(delay)
  x1_delay <- rnorm(n, mean=x1_delay_mean, sd=x1_delay_mean) #derive a range of delays
  x2_population <- rnorm(n, mean=x2_pop_mean, sd=x2_pop_sd) #derive a range of population sizes
  y <- b0 + (b1*x1_delay) + (b2*x2_population) + (b3*x1_delay*x2_population) + 
    rnorm(n, mean=err_mean, sd=err_sd)
  
  model <- lm(y~x1_delay*x2_population)
  summary(model)
  
  output <- summary(model)$coefficients
  coefs <- output[, 1]
  ps <- output[, 4]
  
  rsq <- summary(model)$r.squared
  
  results <- c(coefs, ps, rsq)
  names(results) <- c('intercept_coef', 'delay_coef', 'pop_coef', 'interaction_coef', 'intercept_p','delay_p','population_p', 'interaction_p', 'Rsq')
  return(results)
}

#Check that the function works
regression_sim2(1, n=50, b0=0, b1=0.4, b2=0.1, b3=0.4)

#Run simulation 1000 times 
num_sims <- 1000
(sims <- ldply(1:num_sims, regression_sim2, n=500, b0=0, b1=0.4, b2=0.1, b3=0.4))
(power <- sum(sims$delay_p < .05) / nrow(sims))

#Vary the parameters
sample_sizes <- c(50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000)
results <- NULL
for (val in sample_sizes) {
  sims <- ldply(1:num_sims, regression_sim2, n=500, b0=0, b1=0.4, b2=0.1, b3=0.4)
  sims$n <- val  
  results <- rbind(results, sims)
}

results

power_ests <- results %>%
  group_by(n) %>%
  summarize(power=sum(delay_p < 0.05) / n())

ggplot(power_ests, aes(x=n, y=power)) +
  geom_point() +
  geom_line() +
  ylim(c(0, 1)) +
  theme_minimal()




















