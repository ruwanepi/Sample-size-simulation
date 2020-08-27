#Calculate oower for generalised linear mixed model, using simulation
#Described in Green and MacLeod, 2016 <doi:10.1111/2041-210X.12504>
library(simr) 

#Steps:
#1. Simulate data for predictors (delay, population size per ring) and dependent variable (y, case counts)
#Simulate new values for dependent variable (y)
#3. Refit model to simulated response
#4. Since tested effect is known to exist, every positive test is a true
#positive and every negative test is a Type II error (calculate power)

n=10 #start with 200 rings (rows in dataset)
(x1_delay <- rnorm(n, mean=3, sd=1)) #derive a distribution for delays
(x2_population <- rnorm(n, mean=500, sd=25)) #derive a distribution for populations of rings
(x3_time_obs <- sample(28:30, 200, replace=TRUE))
(ring <- c(1:10)) #set rings as separate groups for GLMM

CATI <- expand.grid(ring=ring, x1_delay=x1_delay)
CATI <- data.frame(CATI)
View(CATI)

(b <- c(2, -0.1)) # fixed intercept and slope
(V1 <- 0.5) # random intercept variance
(V2 <- matrix(c(0.5,0.05,0.05,0.1), 2)) # random intercept and slope variance-covariance matrix
(s <- 1) # residual variance

#GLMM model assuming Poisson distribution
model1 <- makeGlmer(z ~ x1_delay + (1|ring), family="poisson", fixef=b, VarCorr=V1, data=CATI)
print(model1)

fixef(model1)["x1_delay"]
fixef(model1)["x1_delay"] <- -0.15
set.seed(123)
powerSim(model1)

#Power if effect size for delay is 5% reduction in case count and sample size is 100 rings = 38.00% (34.98, 41.09)
#Power if effect size for delay is 10% reduction in case count and sample size is 100 rings = 84.50% (82.11, 86.69)
#Power if effect size for delay is 15% reduction in case count and sample size is 100 rings = 97.90% (96.81, 98.70)
#Power if effect size for delay is 20% reduction in case count and sample size is 100 rings = 99.70% (99.13, 99.94)
#Power if effect size for delay is 40% reduction in case count and sample size is 100 rings = 100.0% (99.63, 100.0)

pc2 <- powerCurve(model1)
print(pc2)

n2=10
(x1_delay2 <- rnorm(n2, mean=3, sd=1)) #derive a distribution for delays
(x2_population2 <- rnorm(n2, mean=500, sd=25)) #derive a distribution for populations of rings
(x3_time_obs2 <- sample(28:30, 200, replace=TRUE))
(ring2 <- c(1:5)) #set rings as separate groups for GLMM

CATI_50 <- expand.grid(ring2=ring2, x1_delay2=x1_delay2)
CATI_50 <- data.frame(CATI_50)
View(CATI_50)

#GLMM model assuming Poisson distribution
model2 <- makeGlmer(z ~ x1_delay2 + (1|ring2), family="poisson", fixef=b, VarCorr=V1, data=CATI_50)
print(model2)

fixef(model2)["x1_delay2"]
fixef(model2)["x1_delay2"] <- -0.4
set.seed(123)
powerSim(model2)

pc3 <- powerCurve(model2)
print(pc3)

#Power if effect size for delay is 5% reduction in case count and sample size is 50 rings = 21.00% (18.51, 23.66
#Power if effect size for delay is 10% reduction in case count and sample size is 50 rings = 56.40% (53.26, 59.50)
#Power if effect size for delay is 15% reduction in case count and sample size is 50 rings = 80.70% (78.11, 83.10)
#Power if effect size for delay is 20% reduction in case count and sample size is 50 rings = 94.90% (93.35, 96.18)
#Power if effect size for delay is 40% reduction in case count and sample size is 50 rings = 99.70% (99.13, 99.94)







