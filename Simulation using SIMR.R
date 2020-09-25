#####################################################################
## Calculate power for a study using a generalised linear mixed model
## through simulation
## Author: Ruwan Ratnayake, September 2020
## Code adapted from Green and MacLeod, 2016 [SIMR package]
## <doi:10.1111/2041-210X.12504>
#####################################################################


####################
# 0 - Load libraries
####################
library(simr) 
library(truncnorm)
library(ggplot2)


#############################
# 1 - Parameters and approach
#############################

# Simulate fixed effects:
# 1.  y = incidence in 200 rings (branching model)
# 2. x1 = delay to implementation (mean, sd)
# 3. x2 = population coverage in ring (fixed %)
# 4. x3 = population size (mean, sd)
# 5. x4 = time ring is under surveillance (range of days)
# Simulate a random effect (conceptualized as )
# 1. x5 = this random effect simulates multiple undefined variables
#         that represent a grouping of the rings. This will be 
#         estimated by a score from 1 to 10, as a random normal variate)

# Approach:
# 1. Simulate data for y (from branching model)
# 2. Simulate data for all x covariates 
# 3. Refit model to simulated dependent variable, y
# 4. Calculate power. Since tested effect is known to exist, every positive test can 
#    be considered a true positive and every negative test is a 
#    Type II error.

# Input incidence from branching model
setwd("/Users/Ruwan/Desktop/Sample-size-simulation")
y.inc <- read.csv(file = 'Incidence.csv')

# Derive the mean and standard deviation by fitting a log normal distribution 
# and gamma distribution to incidence estimates
inc.logn    <- fitdist(y.inc$chol.inc.1000.30d, "lnorm")
summary(inc.logn)
inc.gamma   <- fitdist(y.inc$chol.inc.1000.30d, "gamma")
summary(inc.gamma)
# Graph the fitting graphs
par(mfrow = c(2, 2))
plot.legend <- c("lnorm", "gamma")
denscomp(list(inc.logn, inc.gamma), legendtext = plot.legend)
qqcomp(list(inc.logn, inc.gamma), legendtext = plot.legend)
cdfcomp(list(inc.logn, inc.gamma), legendtext = plot.legend)
ppcomp(list(inc.logn, inc.gamma), legendtext = plot.legend)

# Use log normal distribution (good fit) to calculate the fitted mean and SD
# = mean (20.1), sd (4.7) versus mean (20.1), sd (4.5) (very similar)
#(mean <-exp(2.976159+0.5*(0.228274^2)))
#(sd  <-exp(2.976159+0.5*(0.228274^2))*sqrt(exp(0.228274^2)-1))
# Compare with the mean and SD from the branching model incidence
mean(y.inc$chol.inc.1000.30d)
sd(y.inc$chol.inc.1000.30d)

# Simulate data for co-variates for 200 rings
n=200
x1.delay        <- rtruncnorm(n, a=1, b=7, mean=2, sd=2)
x2.coverage     <- sample(75:100, n, replace=TRUE)
x3.population   <- rnorm(n, mean=450, sd=50)
x4.time.surv    <- sample(25:30, n, replace=TRUE)
x5.random.effect<- rep(1:50, 4)

# Assemble dataset
cati.sim <- data.frame(y.inc$chol.inc.1000.30d, x1.delay, x2.coverage, x3.population, x4.time.surv, x5.random.effect)
#View(cati.sim)

#########################################
# 2 - Generate GLMM epidemic growth model
#########################################

# GLMM model using Poisson distribution

model.sim <- glmer(y.inc$chol.inc.1000.30d ~ x1.delay + x2.coverage + x3.population + x4.time.surv + 
                    (1|x5.random.effect), family="poisson", data=cati.sim)
summary(model.sim)

#evaluate power using a range of fixed effects for model from a
#10% to 50% reduction in incidence
set.seed(123)
fixef(model1)["x1_delay"] <- -0.1
powerSim(model1)
fixef(model1)["x1_delay"] <- -0.15
powerSim(model1)
fixef(model1)["x1_delay"] <- -0.2
powerSim(model1)
fixef(model1)["x1_delay"] <- -0.25
powerSim(model1)
fixef(model1)["x1_delay"] <- -0.3
powerSim(model1)
fixef(model1)["x1_delay"] <- -0.35
powerSim(model1)
fixef(model1)["x1_delay"] <- -0.4
powerSim(model1)
fixef(model1)["x1_delay"] <- -0.5
powerSim(model1)

#Power for predictor 'x1_delay'
# ****adequate power
#(95% CI): 17.10% (14.82, 19.58), effect size (-0.1)
#(95% CI): 27.70% (24.95, 30.59), effect size (-0.15)
#(95% CI): 38.30% (35.28, 41.39), effect size (-0.2)
#(95% CI): 54.90% (51.76, 58.02), effect size (-0.25)
#(95% CI): 69.80% (66.85, 72.63), effect size (-0.3)
#(95% CI): 80.50% (77.91, 82.91), effect size (-0.35)****
#(95% CI): 87.00% (84.76, 89.02), effect size (-0.4)**** (oversampling)
#(95% CI): 96.50% (95.17, 97.55), effect size (-0.5)**** (oversampling)


#High power when sample size is 100 rings and fixed effect size
#for the delay is 40% reduction in incidence. What happens at
#smaller sample sizes from 20-100 rings, holding the effect size
#at a 40% reduction in incidence?
#Power analysis of a range of sample sizes
set.seed(123)

fixef(model1)["x1_delay"] <- -0.4
pc1 <- powerCurve(model1)
print(pc1)
plot(pc1)
#Power for predictor 'x1_delay', (95% confidence interval),
#by largest value of x1_delay:
#With fixed effect size of 40%
#1.47783431771464:  0.00% ( 0.00,  0.37) - 3 rows
#1.99722517452709:  0.00% ( 0.00,  0.37) - 14 rows
#2.44107457202161: 20.80% (18.32, 23.45) - 25 rows
#2.68651460285968: 25.00% (22.34, 27.81) - 35 rows
#2.91674041350385: 40.50% (37.44, 43.62) - 46 rows
#3.25745848567784: 54.60% (51.45, 57.72) - 57 rows
#3.58784495217368: 64.40% (61.34, 67.37) - 68 rows
#3.86875328039287: 75.30% (72.50, 77.95) - 78 rows
#4.26267151868311: 82.90% (80.42, 85.18) - 89 rows****
#5.45371526042283: 89.10% (87.00, 90.96) - 100 rows**** (oversampling)
# The number of observations could be reduced to 89 

fixef(model1)["x1_delay"] <- -0.5
pc1 <- powerCurve(model1)
print(pc1)
plot(pc1)
#With fixed effect size of 50%
#1.47783431771464:  0.00% ( 0.00,  0.37) - 3 rows
#1.99722517452709:  0.00% ( 0.00,  0.37) - 14 rows
#2.44107457202161: 24.10% (21.48, 26.87) - 25 rows
#2.68651460285968: 34.80% (31.85, 37.84) - 35 rows
#2.91674041350385: 51.50% (48.35, 54.64) - 46 rows
#3.25745848567784: 64.10% (61.04, 67.08) - 57 rows
#3.58784495217368: 76.30% (73.54, 78.91) - 68 rows
#3.86875328039287: 85.90% (83.59, 88.00) - 78 rows****
#4.26267151868311: 94.20% (92.57, 95.57) - 89 rows**** (oversampling)
#5.45371526042283: 97.60% (96.45, 98.46) - 100 rows**** (oversampling)
# The number of observations could be reduced to 78 

fixef(model1)["x1_delay"] <- -0.6
pc1 <- powerCurve(model1)
print(pc1)
plot(pc1)
#1.47783431771464:  0.00% ( 0.00,  0.37) - 3 rows
#1.99722517452709:  0.00% ( 0.00,  0.37) - 14 rows
#2.44107457202161: 32.50% (29.60, 35.50) - 25 rows
#2.68651460285968: 44.70% (41.59, 47.84) - 35 rows
#2.91674041350385: 60.00% (56.89, 63.05) - 46 rows
#3.25745848567784: 74.50% (71.68, 77.18) - 57 rows
#3.58784495217368: 83.40% (80.95, 85.66) - 68 rows****
#3.86875328039287: 91.50% (89.60, 93.15) - 78 rows****
#4.26267151868311: 95.80% (94.36, 96.96) - 89 rows****
#5.45371526042283: 98.30% (97.29, 99.01) - 100 rows****
# The number of observations could be reduced to 68 


fixef(model1)["x1_delay"] <- -0.7
pc1 <- powerCurve(model1)
print(pc1)
plot(pc1)
#1.47783431771464:  0.00% ( 0.00,  0.37) - 3 rows
#1.99722517452709:  0.00% ( 0.00,  0.37) - 14 rows
#2.44107457202161: 38.20% (35.18, 41.29) - 25 rows
#2.68651460285968: 50.80% (47.65, 53.94) - 35 rows
#2.91674041350385: 67.80% (64.80, 70.69) - 46 rows
#3.25745848567784: 78.90% (76.24, 81.39) - 57 rows
#3.58784495217368: 87.80% (85.61, 89.76) - 68 rows****
#3.86875328039287: 94.70% (93.12, 96.01) - 78 rows****
#4.26267151868311: 97.70% (96.57, 98.54) - 89 rows****
#5.45371526042283: 99.40% (98.70, 99.78) - 100 rows****
# The number of observations could be reduced to 68 

fixef(model1)["x1_delay"] <- -0.8
pc1 <- powerCurve(model1)
print(pc1)
plot(pc1)
#1.47783431771464:  0.00% ( 0.00,  0.37) - 3 rows
#1.99722517452709:  0.00% ( 0.00,  0.37) - 14 rows
#2.44107457202161: 37.10% (34.10, 40.18) - 25 rows
#2.68651460285968: 57.90% (54.77, 60.98) - 35 rows
#2.91674041350385: 75.10% (72.30, 77.75) - 46 rows
#3.25745848567784: 86.50% (84.22, 88.56) - 57 rows****
#3.58784495217368: 92.20% (90.36, 93.79) - 68 rows****
#3.86875328039287: 96.40% (95.05, 97.47) - 78 rows****
#4.26267151868311: 98.50% (97.54, 99.16) - 89 rows****
#5.45371526042283: 99.80% (99.28, 99.98) - 100 rows****

fixef(model1)["x1_delay"] <- -0.85
pc1 <- powerCurve(model1)
print(pc1)
plot(pc1)
#1.47783431771464:  0.00% ( 0.00,  0.37) - 3 rows
#1.99722517452709:  0.00% ( 0.00,  0.37) - 14 rows
#2.44107457202161: 39.30% (36.26, 42.41) - 25 rows
#2.68651460285968: 58.00% (54.87, 61.08) - 35 rows
#2.91674041350385: 76.00% (73.23, 78.62) - 46 rows
#3.25745848567784: 88.20% (86.04, 90.13) - 57 rows****
#3.58784495217368: 91.90% (90.03, 93.52) - 68 rows****
#3.86875328039287: 96.10% (94.71, 97.21) - 78 rows****
#4.26267151868311: 98.20% (97.17, 98.93) - 89 rows****
#5.45371526042283: 99.30% (98.56, 99.72) - 100 rows****

fixef(model1)["x1_delay"] <- -0.9
pc1 <- powerCurve(model1)
print(pc1)
plot(pc1)
#1.47783431771464:  0.00% ( 0.00,  0.37) - 3 rows
#1.99722517452709:  0.00% ( 0.00,  0.37) - 14 rows
#2.44107457202161: 36.30% (33.31, 39.37) - 25 rows
#2.68651460285968: 56.40% (53.26, 59.50) - 35 rows
#2.91674041350385: 74.90% (72.09, 77.56) - 46 rows
#3.25745848567784: 85.40% (83.06, 87.53) - 57 rows****
#3.58784495217368: 92.40% (90.58, 93.97) - 68 rows****
#3.86875328039287: 95.90% (94.48, 97.04) - 78 rows****
#4.26267151868311: 98.10% (97.05, 98.85) - 89 rows****
#5.45371526042283: 99.30% (98.56, 99.72) - 100 rows****


#Increase the number of ring density groups ("g") from 10 to 15
#model2 <-extend(model1, along="x5_ring_density", n=15)
#pc2 <- powerCurve(model2, along="x5_ring_density")
#plot(pc2)

#Increase the sample size to 150 observations
#model2 <-extend(model1, along="x1_delay", n=15)
#fixef(model2)["x1_delay"] <- -0.1
#powerSim(model2)
#fixef(model2)["x1_delay"] <- -0.15
#powerSim(model2)
#fixef(model2)["x1_delay"] <- -0.2
#powerSim(model2)
#fixef(model2)["x1_delay"] <- -0.25
#powerSim(model2)
#fixef(model2)["x1_delay"] <- -0.3
#powerSim(model2)
#fixef(model2)["x1_delay"] <- -0.35
#powerSim(model2)
#fixef(model2)["x1_delay"] <- -0.4
#powerSim(model2)
#fixef(model2)["x1_delay"] <- -0.5
#powerSim(model2)
#(95% CI): 20.20% (17.75, 22.82), effect size (-0.1)
#(95% CI): 35.30% (32.33, 38.35), effect size (-0.15)
#(95% CI): 56.00% (52.86, 59.11), effect size (-0.2)
#(95% CI): 68.50% (65.52, 71.37), effect size (-0.25)
#(95% CI): 83.70% (81.26, 85.94), effect size (-0.3)
#(95% CI): 91.90% (90.03, 93.52), effect size (-0.35)
#(95% CI): 96.20% (94.82, 97.30), effect size (-0.4)
#(95% CI): 99.30% (98.56, 99.72), effect size (-0.5)
##adequate power at 30% reduction in incidence

