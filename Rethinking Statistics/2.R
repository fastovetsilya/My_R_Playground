library(rethinking)

### Grid approximation: 2.3-2.5

# Define grid
p_grid <- seq(from = 0, to = 1, length.out = 20)

# Define prior
prior <- rep(1, 20)
# prior <- ifelse(p_grid < 0.5, 0, 1)
# prior <- exp(-5 * abs(p_grid - 0.5))

# Compute likelihood at each value in grid
likelihood <- dbinom(6, size = 9, prob = p_grid)

# Compute product of likelihood and prior
unstd.posterior <- likelihood * prior

# Standardize posterior, so it sums to 1
posterior <- unstd.posterior / sum(unstd.posterior)

# Display posterior distribution
plot(p_grid, posterior, type="b", 
      xlab="Probability of water", 
      ylab="Posterior probability", 
      main="20 points")


### Quadratic approximation: 2.6

# Compute quadratic approximation
globe.qa <- quap(
  alist(
    W ~ dbinom(W+L, p), # binomial likelihood
    p ~ dunif(0, 1) # uniform prior
  ), 
  data=list(W=6, L=3)
)

# Display summary of quadratic approximation
precis(globe.qa)


### Compare analytical calculation with quardatic appx: 2.7

# Analytical calculation
W <- 6
L <- 3
curve(dbeta(x, W+1, L+1), from = 0, to = 1)

# Quadratic approximation
curve(dnorm(x, 0.67, 0.16), lty = 2, add = TRUE)


### MCMC estimate (Metropolis algorithm): 2.8-2.9

# Compute MCMC estimate
n_samples <- 1000
p <- rep(NA, n_samples)
p[1] <- 0.5
W <- 6
L <- 3
for(i in 2:n_samples) {
  p_new <- rnorm(1, p[i-1], 0.1)
  if(p_new < 0) p_new <- abs(p_new)
  if(p_new > 1) p_new <- 2 - p_new
  q0 <- dbinom(W, W + L, p[i-1])
  q1 <- dbinom(W, W + L, p_new)
  p[i] <- ifelse(runif(1) < q1/q0, p_new, p[i-1])
}

# Plot the results to compare with the analytical solution
dens(p, xlim = c(0, 1))
curve(dbeta(x, W + 1, L + 1), lty = 2, add = TRUE)
