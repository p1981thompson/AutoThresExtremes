#
# Functions to simulate from a distribution
# comprising a normal body and a Pareto tail
#
# Function originally written by Paul Thompson,
# then modified by Julian Stander
#
#--------------------------------------------------------------
#
# These functions are used in the simulation function
#

#
# Function to work out the area below u for a normal distribution
# truncated at zero
# with mean gamma and standard deviation alpha.
# A value beta is subtracted for use in get.mean.below
#
area.below <- function(gamma, alpha, u , beta){
(pnorm(u, gamma, alpha) - pnorm(0, gamma, alpha)) / (1 - pnorm(0, gamma, alpha)) - beta
}

#
# Given beta, alpha and u, the mean gamma of the above truncated normal
# is determined.  This function calculates it.
#
get.mean.below <- function(beta, alpha, u){
# beta is proportion below
# alpha is standard deviation below
# u is threshold
uniroot(area.below,c(0,u), beta = beta, alpha = alpha, u = u)$root
}

#
# Example of its use
#
get.mean.below(0.9, 0.7, 2.9)


# Function to simulated from a normal body Pareto tail distribution
# with continuous density function (see Thompson, Cai, Reeve and Stander)

# q is the proportion of the density that is below the threshold
# Above q is beta and threshold is u
# sigma.below is the standard deviation of the truncated normal body
# Above sigma.below is gamma
#
# n is the number of realizations required

r.norm.pareto <- function(n=10000,  q = 0.9, threshold, sigma.below = 0.7, xi = 0.2, do.plot = FALSE)
{
#
# Work out mean.below (alpha above)
#
mean.below <- get.mean.below(q, sigma.below, threshold)

sigma.gpd <- (1-q) *  pnorm(0, mean.below, sigma.below, lower.tail = FALSE) / dnorm(threshold, mean.below, sigma.below)

# Vector for realizations

newdist<- vector(mode = "numeric", length = n)
#
# Now fill up the vector according to the algorithm in the paper
#
for(i in 1:n)
{
  newdist[i]<-rnorm(1,mean.below,sigma.below)
  if(newdist[i]<=0)
      {
      while(newdist[i]<=0)
        {
        newdist[i]<-rnorm(1,mean.below,sigma.below)
          if(newdist[i]>threshold)
            {
            newdist[i]<-rgpd(1,xi=xi,mu=threshold,beta=sigma.gpd)
            }
        }
      }
  if(newdist[i]>threshold)
      {
       newdist[i]<-rgpd(1,xi=xi,mu=threshold,beta=sigma.gpd)
      }
}
#
# Plot
#
if(do.plot){
par(mfrow = c(2,1))
hist(newdist,nclass = 100, xlab = "Realizations", main = "Histogram of realizations from a truncated normal body and a Pareto tail")
#
title(sub = substitute(paste(u == U, ", ", beta == Beta, ", ", gamma == Gamma, ", ", alpha == Alpha, ", ", xi == Xi, ", ", sigma == Sigma),list(U = signif(threshold,3), Beta = signif(q,3), Gamma = signif(mean.below,3), Alpha = signif(sigma.below,3),  Xi = signif(xi,3), Sigma = signif(sigma.gpd,3))))
#
# Rug
#
rug(newdist)
#
# Show the threshold
#
abline(v = threshold, col = "red", lwd = 2)
#
plot(newdist)
abline(h = threshold, col = "red", lwd = 2)

}
#
#
if(do.plot){
print(paste("Proportion below threshold", length(newdist[newdist < threshold])/length(newdist)))
}
#
# Return
#
return(invisible(list(data = newdist, mean.below = mean.below, sigma.gpd = sigma.gpd)))
}
#--------------------------------------------------------------------------
#
#
#--------------------------------------------------------------------------
#
# Function to estimate a 95% return level confidence interval from a GPD fit
#
gpd.ci <- function (z,  return.year = 100)
{
a <- z$mle
u <- z$threshold
la <- z$rate
n <- z$n
npy <- z$npy
mat <- z$cov
#
# a is the mle of sigma and xi
# u is the threshold
# la is the  mle estimate of the exceedance probability
# n is the number of data points
# npy is the number of observations per year, set by default to 365 in gpd.fit
# mat is the variance-covariance matrix for sigma and xi
# return.year is the return period in years
#
# m is the number of observations, i.e. the return period in years, multiplied by
# the number of observations per year
#
    m <- npy * return.year
#
    a <- c(la, a)
    eps <- 1e-06
    a1 <- a
    a2 <- a
    a3 <- a
    a1[1] <- a[1] + eps
    a2[2] <- a[2] + eps
    a3[3] <- a[3] + eps
    q <- gpdq2(a[2:3], u, la, m)
    d1 <- (gpdq2(a1[2:3], u, la, m) - q)/eps
    d2 <- (gpdq2(a2[2:3], u, la, m) - q)/eps
    d3 <- (gpdq2(a3[2:3], u, la, m) - q)/eps
    d <- cbind(d1, d2, d3)
    mat <- matrix(c((la * (1 - la))/n, 0, 0, 0, mat[1, 1], mat[1,
        2], 0, mat[2, 1], mat[2, 2]), nc = 3)
    v <- apply(d, 1, q.form, m = mat)
    return(list(return.level = q, lower.limit = q - 1.96 * sqrt(v), upper.limit = q + 1.96 * sqrt(v)))
}
#--------------------------------------------------------------------------


# We need the library evir for generating from the Pareto tail

library(evir)
library(ismev)

# n, Proportion, threshold, sigma.below, shape

# Paper's notation

beta <- 0.9
u <- 2.9
alpha <- 0.7
xi <- 0.2

sim.data <- r.norm.pareto(365*20, beta, u,  alpha, xi, do.plot = TRUE)

# r.norm.pareto also calculate the paper's sigma (and gamma)

sigma <- sim.data$sigma.gpd
#
ny <- 365
#
# zeta is probability of being in the tail
#
zeta <- 1 - beta
#
# *** 100 year return level using formula in Coles ***
#
N <- 100
true.return.level <- u + (sigma / xi)*((N*ny*zeta)^xi - 1)
true.return.level

#
#
# Check coverage of confidence intervals
#
# Set count to zero, generate n.rep data sets and work out and check confidence interval for each one
#
i.cover <- 0
n.rep <- 1000
#
for(i in 1:n.rep){

if(i %% 10 == 0){print(paste("Iteration", i))}

#
sim.data <- r.norm.pareto(365*20, beta, u,  alpha, xi, do.plot = FALSE)
#
#
# ***** PAUL:  THESE LINES OF CODE COULD BE REPLACED WITH YOUR BOOTSTRAP CONFIDENCE INTERVAL GENERATING CODE
#
# Fit model
#
my.fit <- gpd.fit(sim.data$data, u, show = FALSE)
#
# Estimate the confidence interval
#
ci.hat <- gpd.ci(my.fit, N)
#
# ***** PAUL:  TO HERE
#
#
# Check whether it contains the true value
#
i.cover <- i.cover + ifelse(ci.hat$lower.limit < true.return.level & true.return.level < ci.hat$upper.limit, 1, 0)
}
#
# Proportion covering
#
print(paste("Proportion of intervals covering", i.cover / n.rep))


