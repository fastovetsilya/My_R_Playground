library(rethinking)

### Inverting the probability: 3.1

Pr_Positive_Vampire <- 0.95
Pr_Positive_Mortal <- 0.01
Pr_Vampire <- 0.001
Pr_Positive <- Pr_Positive_Vampire * Pr_Vampire + 
               Pr_Positive_Mortal * (1 - Pr_Vampire)
(Pr_Vampire_Positive <- Pr_Positive_Vampire * Pr_Vampire / Pr_Positive)


### Sampling from a grid-approximate posterior:3.2-3.19

# Compute the grid approximation
p_grid <- seq(from = 0, to = 1, length.out = 1000)
prob_p <- rep(1, 1000)
prob_data <- dbinom(6, size = 9, prob = p_grid)
posterior <- prob_data * prob_p
posterior <- posterior / sum(posterior)

# Draw samples from the posterior
samples <- sample(p_grid, prob = posterior, size = 1e4, replace = TRUE)

# Plot the drawn samples
plot(samples)
dens(samples)

# Intervals of defined boundaries
# Compute posterior probabilities
sum(posterior[p_grid < 0.5])
sum(samples < 0.5) / 1e4
sum(samples > 0.5 & samples < 0.75) / 1e4
quantile(samples, 0.8)
quantile(samples, c(0.1, 0.9))

# Three waters in three tosses and uniform (flat) prior
p_grid <- seq(from = 0, to = 1, length.out = 1000)
prior <- seq(1, 1000)
likelihood <- dbinom(3, size = 3, prob = p_grid)
posterior <- likelihood * prior
posterior <- posterior / sum(posterior)
samples <- sample(p_grid, size = 1e4, replace = TRUE, prob = posterior)

# Confidence intervals
PI(samples, prob = 0.5) # Percentile interval
HPDI(samples, prob = 0.5) # Highest probability density interval

# Maximum a posteriori estimate (MAP)
p_grid[which.max(posterior)]

# Mode, mean, median
chainmode(samples, adj = 0.1)
mean(samples)
median(samples)

# Compute the loss

# Expected loss considering p=0.5 decision
sum(posterior * abs(0.5 - p_grid))

# Compute the loss
loss <- sapply(p_grid, function(d) sum(posterior * abs(d - p_grid)))

# Compute parameter value that minimizes loss 
p_grid[which.min(loss)]


### Sampling to simulate prediction: 3.20-3.26

# Compute analytically
dbinom(0:2, size = 2, prob = 0.7)

# Sample from distributions
rbinom(1, size = 2, prob = 0.7)
rbinom(10, size = 2, prob = 0.7)

# Generate dummy observations
dummy_w <- rbinom(10e5, size = 2, prob = 0.7)
table(dummy_w) / 10e5

# Simulate 9 tosses
dummy_w <- rbinom(1e5, size = 9, prob = 0.7)
simplehist(dummy_w, xlab = "Dummy water count")

# Draw from Binomial distribution
w <- rbinom(10e4, size = 9, prob = 0.6)
w <- rbinom(10e4, size = 9, prob = samples)
