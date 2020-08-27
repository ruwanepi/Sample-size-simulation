# This is code to carry out a stochastic branching model for 
# cholera infection that treats delay (days) as a continuous 
# variable and epidemic size, as an independent variable.

library(fitdistrplus)
library(viridis)

SI <- 5  
SIrate <- 0.1 
SIshape <- 0.5 
R <- 3 #(unvaccinated population)
Rrange <- seq(1, 3, length.out=5) #reproductive number range (range, 1-3)
k <- 4.5
Drange <- seq(0, 7, length.out=8) #delay range (range, 0-7 days)
cap_cases <- 1000
set.seed(645)
runs <- 1000
seed <- 3

## 6. Simulate outbreak trajectory for 30-day delay
plot(NA,xlim=c(0,14),ylim=c(0,50),xlab="Delay since primary case",
     ylab="Number of cases",frame=FALSE)
cols <- sample(viridis(runs))

total_cases <- integer(runs)

for(i in 1:runs) {
  cases <- seed
  t <- rep(0,seed)
  times <- t
  while(cases > 0) {
    secondary <- rnbinom(cases,size=k,mu=R)
    t.new <- numeric()
    for(j in 1:length(secondary)) {
      t.new <- c(t.new,t[j] + rgamma(secondary[j],shape=SIshape[1],
                                     scale=(1/SIrate)))
    }
    cases <- length(t.new) & cases<cap_cases
    t <- t.new
    times <- c(times[times<=7],t.new)
  }
  lines(sort(times),1:length(times),col=cols[i],lwd=1)
  points(max(times),length(times),col=cols[i],pch=16)
  total_cases[i] <-length(times)
}

print(total_cases.median <- median(total_cases))
print(total_cases.sd <- sd(total_cases))
print(total_cases.range <- range(total_cases))

#Transform into log scale
lsm=log(total_cases.median)-(1/2)*log((total_cases.sd/total_cases.median)^2+1)
lssd=sqrt(log((total_cases.sd/total_cases.median)^2+1))
plnorm(20, lsm, lssd)


