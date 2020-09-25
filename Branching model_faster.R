## Branching process model of epidemic spread

# This code is adapted from the study of [Riou and Althaus](https://www.biorxiv.org/content/10.1101/2020.01.23.917351v1), using the [code from GitHub](https://github.com/jriou/wcov), stripped down and rewritten for clarity.

library(data.table)
library(ggplot2)

# Set random number seed.

set.seed(1234)

# Set parameter values, using values from [Imai et al.](https://www.imperial.ac.uk/media/imperial-college/medicine/sph/ide/gida-fellowships/Imperial-2019-nCoV-transmissibility.pdf) and assuming a gamma distribution of generation times with mean `mu` 8.4 days and standard deviation `sd` of 3.8 days, taken from [Lipsitch et al. (2003](https://science.sciencemag.org/content/300/5627/1966.full).

R0 <- 2.6
k <- 0.16
mu <- 8.4
stdev <- 3.8
shape <- (mu/stdev)^2
scale <- (stdev^2)/mu

rcases <- function(n) rnbinom(n, size=k, mu=R0)
rtimes <- function(n) rgamma(n, shape = shape, scale = scale)

# Define Bellman-Harris branching process model.

bhbp <- function(
  index_cases, max_cases, max_time,
  casedistro = rcases, timedistro = rtimes
){
  cases <- index_cases
  times <- rep(0, cases)
  remaining_cases <- max_cases - cases
  while(cases & (0 < remaining_cases)) {
    secondary = casedistro(cases)
    tms <- timedistro(sum(secondary))
    t.new <- times[rep(1:cases, times = secondary)] + tms
    t.new <- t.new[t.new < max_time]
    cases <- length(t.new)
    remaining_cases <- remaining_cases - cases
    times <- c(t.new, times)
  }
  times <- sort(times)
  return(times)
}

# Plot the generation time distribution.

t <- seq(0,30,by=0.01)
g <- dgamma(t,shape=shape,scale=scale)
ggplot(
  data.table(`Generation Time`=t,Probability=g)
) + 
  aes(`Generation Time`, Probability) +
  geom_line()

# Plot the offspring distribution.

i <- seq(0,10)
d <- dnbinom(i, size=k, mu=R0)
ggplot(
  data.table(Number=i,Probability=d)
) +
  aes(x=Number,y=Probability) +
  geom_bar(stat="identity") +
  scale_x_continuous(breaks=i)

# Initial and stopping conditions.

index_cases <- 40
max_cases <- 5e4
max_time <- 90


# Set the number of simulations (note - Imai et al. used 5000).

nsims <- 500

# Run simulations.

grid <- data.table(sim=1:nsims)
compute <- grid[, {
  times <- bhbp(index_cases, max_cases, max_time)
  .(
    cases = cumsum(hist(times, breaks = 0:max_time, plot=FALSE)$counts),
    times = 1:max_time
  )
}, by=sim ]

day0 <- as.Date("2019-12-01")

summary_results <- compute[, {
  qs <- quantile(cases, probs = c(0.025, 0.5, 0.975))
  .(Lower=qs[1], Median=qs[2], Upper=qs[3])
}, by=.(Day=times)]
summary_results[, Date := day0+Day ]

# Plot trajectories over time, highlighting 4000 cases on 2020-01-18.

ntraj <- 20
ggplot(summary_results)+
  geom_line(aes(x=Date,y=Median),color="red")+
  geom_line(aes(x=Date,y=Upper),color="blue")+
  geom_line(aes(x=Date,y=Lower),color="blue")+
  geom_ribbon(aes(x=Date,ymin=Lower,ymax=Upper, linetype = NA), fill="blue", alpha = 0.1)+
  geom_line(aes(x=times+day0,y=cases,group=factor(sim)),data=compute[sim<=ntraj,],color="gray")+
  coord_cartesian(xlim=c(as.Date("2019-12-01"),as.Date("2020-01-30")),ylim=c(1,20000))+
  geom_vline(xintercept=as.Date("2020-01-18"))+
  geom_hline(yintercept=4000)+
  theme_classic()