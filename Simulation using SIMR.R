#Calculate oower for generalised linear mixed model, using simulation
#Described in Green and MacLeod, 2016 <doi:10.1111/2041-210X.12504>
library(simr) 

##Assumptions:
#Fixed effects: delay (mean=3, SD=1), population (mean=450, SD=50),
#time observed (28-30d), coverage (80-100%)
#Random effect: population density (score: 1-3, random normal variate, mean 2, SD 0.25)
#Start with 100 rings (n)
##Steps:
#1. Simulate data for predictors (delay, population size per ring) and dependent variable (y, case counts)
#Simulate new values for dependent variable (y)
#3. Refit model to simulated response
#4. Since tested effect is known to exist, every positive test is a true
#positive and every negative test is a Type II error (calculate power)

n=10 #start with 5 rings to expand to 100 (rows in dataset)
(x1_delay <- rnorm(n, mean=3, sd=1)) #derive a distribution for delays
(x2_population <- rnorm(100, mean=450, sd=50)) #derive a distribution for populations of rings
(x3_time_obs <- sample(25:30, 100, replace=TRUE))
(x4_coverage <-sample(75:100, 100, replace=TRUE))
(x5_ring_density <- c(1:10)) #set rings as 10 separate groups for random effect

CATI <- expand.grid(x5_ring_density=x5_ring_density, x1_delay=x1_delay)
CATI <- data.frame(CATI)
#CATI <- cbind(CATI, x2_population, x3_time_obs, x4_coverage)
View(CATI)

(b <- c(2, -0.1)) # fixed intercept and slope, starting at 10% reduction
(V1 <- 0.5) # random intercept variance
(V2 <- matrix(c(0.5,0.05,0.05,0.1), 2)) # random intercept and slope variance-covariance matrix
(s <- 1) # residual variance

#GLMM model Poisson distribution
model1 <- makeGlmer(y ~ x1_delay + (x1_delay|x5_ring_density), family="poisson", fixef=b, VarCorr=V2, data=CATI)
print(model1)

set.seed(123)
fixef(model1)["x1_delay"] <- -0.15
fixef(model1)["x1_delay"] <- -0.2
fixef(model1)["x1_delay"] <- -0.25
fixef(model1)["x1_delay"] <- -0.3
fixef(model1)["x1_delay"] <- -0.35
fixef(model1)["x1_delay"] <- -0.4
fixef(model1)["x1_delay"] <- -0.5
powerSim(model1)
#Power for predictor 'x1_delay'
#(95% CI): 15.20% (13.03, 17.58), effect size (-0.1)
#(95% CI): 22.20% (19.66, 24.91), effect size (-0.15)
#(95% CI): 32.90% (29.99, 35.91), effect size (-0.2)
#(95% CI): 45.00% (41.89, 48.14), effect size (-0.25)
#(95% CI): 59.50% (56.38, 62.56), effect size (-0.3)
#(95% CI): 70.70% (67.77, 73.51), effect size (-0.35)
#(95% CI): 77.40% (74.68, 79.96), effect size (-0.4)
#(95% CI): 88.60% (86.47, 90.50), effect size (-0.5)

#Increase the sample size to 150 observations
model2 <-extend(model1, along="x1_delay", n=15)
fixef(model2)["x1_delay"] <- -0.1
powerSim(model2)
fixef(model2)["x1_delay"] <- -0.15
powerSim(model2)
fixef(model2)["x1_delay"] <- -0.2
powerSim(model2)
fixef(model2)["x1_delay"] <- -0.25
powerSim(model2)
fixef(model2)["x1_delay"] <- -0.3
powerSim(model2)
fixef(model2)["x1_delay"] <- -0.35
powerSim(model2)
fixef(model2)["x1_delay"] <- -0.4
powerSim(model2)
fixef(model2)["x1_delay"] <- -0.5
powerSim(model2)
#(95% CI): 20.20% (17.75, 22.82), effect size (-0.1)
#(95% CI): 35.30% (32.33, 38.35), effect size (-0.15)
#(95% CI): 56.00% (52.86, 59.11), effect size (-0.2)
#(95% CI): 68.50% (65.52, 71.37), effect size (-0.25)
#(95% CI): 83.70% (81.26, 85.94), effect size (-0.3)
#(95% CI): 91.90% (90.03, 93.52), effect size (-0.35)
#(95% CI): 96.20% (94.82, 97.30), effect size (-0.4)
#(95% CI): 99.30% (98.56, 99.72), effect size (-0.5)

#Power analysis of a range of sample sizes
pc2 <- powerCurve(model2)
print(pc2)
#Power for predictor 'x1_delay', (95% confidence interval),
#by largest value of x1_delay:
#Time elapsed: 0 h 31 m 1 s
#3: 67.80% (64.80, 70.69) - 30 rows
#4: 84.50% (82.11, 86.69) - 40 rows
#6: 94.30% (92.68, 95.65) - 60 rows
#7: 94.90% (93.35, 96.18) - 70 rows
#8: 95.30% (93.80, 96.53) - 80 rows
#10: 96.10% (94.71, 97.21) - 100 rows
#11: 96.10% (94.71, 97.21) - 110 rows
#12: 96.10% (94.71, 97.21) - 120 rows
#14: 96.10% (94.71, 97.21) - 140 rows
#15: 96.30% (94.94, 97.38) - 150 rows