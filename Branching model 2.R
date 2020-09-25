#################################################################
## Branching process model for generating final epidemic size of 
## symptomatic cholera cases in a small outbreak (ring) 
## Ruwan Ratnayake, September 2020
## Code for branching process adapted from Althaus, 2015 
## <doi.org/10.1016/S1473-3099(15)70135-0>
#################################################################


##############################
# 0 - Load libraries
##############################

library(truncnorm)
library(fitdistrplus)
library(viridis)


##############################
# 1 - Parameters
##############################

# Incubation period and serial interval (gamma distribution) (Kahn, 2019, Azman, 2015) 
r.0          <- 2.5 #Azman, 2015 (early reproduction number before intervention)
ser.int      <- 5   #Kahn, 2019; Azman, 2015
ser.int.rate <- 0.1 #Kahn, 2019; Azman, 2015
ser.int.shape<- 0.5 #Kahn, 2019; Azman, 2015
k            <- 4.5 #Moore, 2014 (dispersion parameter for low potential of superspreading)

# Baseline scenario
cap.max.days <- 30 #Days under observation
cap.cases    <- 1000 #Maximum number of cases within a ring

# Estimates of reduction of R0
antibiotic.eff <- 0.66 #ACP, Reveiz 2011 (meta-analysis)
vacc.2m.eff    <- 0.87 #OCV, 2 month, Azman, 2015 (case-cohort)
water.tx.eff   <- 0.26 #POUWT, Fewtrell, 2005 (meta-analysis)
water.store.eff<- 0.21 #Safe storage, Roberts 2001 (RCT)
pop.cover      <- 0.8  #Coverage at the population level   

# Overall effectiveness of CATI
# Impact (effected on 1st day; initial). Assume no impact from vaccination.
# Estimated reduction of R0 by CATI without vaccination is 0.8
cati.efficacy      <- 1-((1-antibiotic.eff)*(1-water.tx.eff)*(1-water.store.eff))
# Estimated impact within first week is to knock down R0 to 0.998
cati.effectiveness <- pop.cover*cati.efficacy
RE.novacc      <- r.0*(1-(cati.effectiveness))

# Impact (effected on 7th day; initial). Include impact from vaccination.
# Estimated reduction of R0 by CATI without vaccination results in Re=0.97
cati.vacc.efficacy      <-  1-((1-vacc.2m.eff)*(1-antibiotic.eff)*(1-water.tx.eff)*(1-water.store.eff))
# Estimated reduction of R0 by CATI with vaccination results in Re=0.67
cati.vacc.effectiveness <- pop.cover*cati.vacc.efficacy
RE.withvacc             <- r.0*(1-(cati.vacc.effectiveness))
  

######################################
# 2 - Run the branching process model
######################################

set.seed(123)
runs          <- 200
# distribution of delays that can be randomly selected for each run
delay.distr   <- rtruncnorm(n=200, a=1, b=6, mean=3, sd=2) 
seed          <- 5

# Simulate outbreak trajectory over a 30-day period (maximum 1000 cases per ring)
plot(NA,xlim=c(0,30),ylim=c(0,50),xlab="Delay since detection of primary case (days)",
     ylab="Number of cases", cex.lab=0.6, cex.axis=0.6, frame=FALSE)
cols <- sample(viridis(runs))

# Simulate epidemic growth given varying R values as a result of R0 (day 1),
# implementation on day 1, CATI without vaccination protection (days 2-6), 
# CATI with vaccination protection (days 7-30)

total.cases <- integer(runs) # set frame with 1000 runs 
#total.cases <- data.frame(integer(runs), delay.distr)

# *** How do I include delay.distr? 
# I would like to iterate each run over the corresponding value 
# from delay.distr to force a delay on CATI implementation

for(i in 1:runs) { # Am I missing delay.distr in for loop?
  cases <- seed
  t <- rep(0,seed)
  times <- t  
  
while(cases > 0) # Am I missing delay.distr in while loop?
  
  if (length(times)<2) # Use different Re values as CATI is implemented. Here it's R0.
    {secondary <- rnbinom(cases,size=k,mu=r.0)
    t.new <- numeric()
    for(j in 1:length(secondary)) {
      t.new <- c(t.new,t[j] + rgamma(secondary[j],shape=ser.int.shape[1],
                                     scale=(1/ser.int.rate)))
    }
    cases <- length(t.new) & cases<cap.cases
    t <- t.new
    times <- c(times[times<31],t.new)
  } 
  else if (length(times)>1 & length(times)<7) { #Here it's Re with CATI (no OCV protection)
    secondary <- rnbinom(cases,size=k,mu=RE.novacc)
    t.new <- numeric()
    for(j in 1:length(secondary)) {
      t.new <- c(t.new,t[j] + rgamma(secondary[j],shape=ser.int.shape[1],
                                     scale=(1/ser.int.rate)))
    }
    cases <- length(t.new) & cases<cap.cases
    t <- t.new
    times <- c(times[times<31],t.new)
  } 
  else if (length(times)>6 & length(times)<31){ # Here it's Re with CATI (and OCV protection)
    secondary <- rnbinom(cases,size=k,mu=RE.withvacc)
    t.new <- numeric()
    for(j in 1:length(secondary)) {
      t.new <- c(t.new,t[j] + rgamma(secondary[j],shape=ser.int.shape[1],
                                     scale=(1/ser.int.rate)))
    }
    cases <- length(t.new) & cases<cap.cases
    t <- t.new
    times <- c(times[times<31],t.new)
  } 
  lines(sort(times),1:length(times),col=cols[i],lwd=0.5)
  points(max(times),length(times),col=cols[i],pch=16)
  total.cases[i] <-length(times)
}


##############################
# 3 - Generate results
##############################

# Generate final epidemic size and parameters
print(total_cases.median <- median(total.cases))
print(total_cases.sd <- sd(total.cases))
print(total_cases.range <- range(total.cases))

# Transform into log scale to estimate proportion of observations under 
# 5, 10, and 20 cases, using a log normal distribution
lsm=log(total_cases.median)-(1/2)*log((total_cases.sd/total_cases.median)^2+1)
lssd=sqrt(log((total_cases.sd/total_cases.median)^2+1))
plnorm(5, lsm, lssd)
plnorm(10, lsm, lssd)
plnorm(20, lsm, lssd)

# Create dataframe with final epidemic size and incidence across 200 runs
# Incidence = (epidemic size/500 population) * 1000
# ***Add another column for the delay used
cases.df <- data.frame("run.num"=1:200, "symp.chol.cases.30d" = total.cases,
                       "chol.inc.1000.30d" = (total.cases/500)*1000)
View(cases.df)
write.csv(cases.df,"C:\\Users\\Ruwan\\Desktop\\Sample-size-simulation\\Incidence.csv", row.names = FALSE)