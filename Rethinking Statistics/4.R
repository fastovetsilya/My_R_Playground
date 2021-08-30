library(rethinking)
library(MASS)

### Normal by addition: 4.1

# Simulate random walk
pos <- replicate(1000, sum(runif(16, -1, 1)))
hist(pos, freq = FALSE)
lines(density(pos))

### Normal by multiplication: 4.2-4.5

# Sample random growth rate
prod(1 + runif(12, 0, 0.1))

# Repeat multiple times
growth <- replicate(10000, prod(1 + runif(12, 0, 0.1)))
dens(growth, norm.comp = TRUE)

# Compare multiplication with big and small numbers
big <- replicate(10000, prod(1 + runif(12, 0, 0.5)))
small <- replicate(10000, prod(1 + runif(12, 0, 0.01)))
dens(big, norm.comp = TRUE)
dens(small, norm.comp = TRUE)

# Logarithm of big number multiplication
log.big <- replicate(10000, log(prod(1 + runif(12, 0, 0.5))))
dens(log.big, norm.comp = TRUE)


### Probability of water from model definition: 4.6

w <- 6; n <- 9;
p_grid <- seq(from = 0, to = 1, length.out = 100)
posterior <- dbinom(w, n, p_grid) * dunif(p_grid, 0, 1)
posterior <- posterior / sum(posterior)
dens(posterior)


### Gaussian model of height: 4.7-4.25

# Load Howell data
data(Howell1)
d <- Howell1

# Inspect dataframe
str(d)
precis(d)
d$height

# Filter the data by age
d2 <- d[d$age >= 18, ]

# Plot the priors
curve(dnorm(x, 178, 20), from = 100, to = 250)
curve(dunif(x, 0, 50), from = -10, to = 60)

# Sample heights from prior
sample_mu <- rnorm(1e4, 178, 20)
sample_sigma <- runif(1e4, 0, 50)
prior_h <- rnorm(1e4, sample_mu, sample_sigma)
dens(prior_h)

# Grid approximation of the posterior distribution

mu.list <- seq(from = 150, to = 160, length.out = 100)
sigma_list <- seq(from = 7, to = 9, length.out = 100)
post <- expand.grid(mu = mu.list, sigma = sigma_list)
post$LL <- sapply(1:nrow(post), function(i) sum(
  dnorm(d2$height, post$mu[i], post$sigma[i], log = TRUE)))
post$prod <- post$LL + dnorm(post$mu, 178, 200, TRUE) +
  dunif(post$sigma, 0, 50, TRUE)
post$prob <- exp(post$prod - max(post$prod))

contour_xyz(post$mu, post$sigma, post$prob)
image_xyz(post$mu, post$sigma, post$prob)


# Sampling from the posterior
sample.rows <- sample(1:nrow(post), size = 1e4, replace = TRUE, 
                      prob = post$prob)
sample.mu <- post$mu[sample.rows]
sample.sigma <- post$sigma[sample.rows]
plot(sample.mu, sample.sigma, cex = 0.5, pch = 16, col = col.alpha(rangi2, 0.1))

# Describe the parameters
dens(sample.mu)
dens(sample.sigma)
PI(sample.mu)
PI(sample.sigma)

# Sample size and the normality of sigma posterior
d3 <- sample(d2$height, size = 20)

mu.list <- seq(from = 150, to = 170, length.out = 200)
sigma.list <- seq(from = 4, to = 20, length.out = 200)
post2 <- expand.grid(mu = mu.list, sigma = sigma.list)
post2$LL <- sapply(1:nrow(post2), function(i)
  sum(dnorm(d3, mean = post2$mu[i], sd = post2$sigma[i], 
  log = TRUE)))
post2$prod <- post2$LL + dnorm(post2$mu, 178, 20, TRUE) + 
  dunif(post2$sigma, 0, 50, TRUE)
post2$prob <- exp(post2$prod - max(post2$prod))
sample2.rows <- sample(1:nrow(post2), size = 1e4, replace = TRUE, 
                       prob = post2$prob)
sample2.mu <- post2$mu[sample2.rows]
sample2.sigma <- post2$sigma[sample2.rows]
plot(sample2.mu, sample2.sigma, cex = 0.5, 
     col = col.alpha(rangi2, 0.1), 
     xlab = "mu", ylab = "sigma", pch = 16)
dens(sample2.sigma, norm.comp = TRUE)


### Quadratic approximation: 4.26-4.36

# Load the data
data(Howell1)
d <- Howell1
d2 <- d[d$age >= 18, ]

# Define the model in alist
flist <- alist(
  height ~ dnorm(mu, sigma), 
  mu ~ dnorm(178, 20), 
  sigma ~ dunif(0, 50)
)

# (optional) define start values for the optimizer
start <- list(
  mu = mean(d2$height), 
  sigma = sd(d2$height)
)

# Fit the model to the data
m4.1 <- quap(flist, data = d2, start = start)
precis(m4.1)

# Change prior and refit the model
m4.2 <- quap(
  alist(
    height ~ dnorm(mu, sigma), 
    mu ~ dnorm(178, 0.1), 
    sigma ~ dunif(0, 50) 
  ), data = d2)
precis(m4.2)

# Sample from the posterior
# Get variance-covariance matrix of the Gaussian distribution
vcov(m4.1)

# Decompose variance-covariance matrix
diag(vcov(m4.1))
cov2cor(vcov(m4.1))

# Sample from the quadratic approximation of the posterior
post <- extract.samples(m4.1, n = 1e4)
head(post)
precis(post)

# Closer look at extract.samples functionality
post <- mvrnorm(n = 1e4, mu = coef(m4.1), Sigma = vcov(m4.1))


### Building linear model 4.37-4.63

# Plot adult height and weight against one another
data(Howell1)
d <- Howell1
d2 <- d[d$age >= 18, ]
plot(d2$height, d2$weight)

# Simulate heights from the model using only the priors
set.seed(2971)
N <- 100
a <- rnorm(N, 178, 20)
b <- rnorm(N, 0, 10)

plot(NULL, xlim = range(d2$weight), ylim = c(-100, 400), 
     xlab = "weight", ylab = "height")
abline(h = 0, lty = 2)
abline(h = 272, lty = 1, lwd = 0.5)
mtext("b ~ dnorm(0, 10)")
xbar <- mean(d2$weight)
for(i in 1:N) curve(a[i] + b[i] * (x - xbar), 
   from = min(d2$weight), to = max(d2$weight),  
   add = TRUE, col = col.alpha("black", 0.2))

# Take a look at log-normal distribution
b <- rlnorm(1e4, 0, 1)
dens(b, xlim = c(0, 5), adj = 0.1)

# Repeat the simulation with log-normal distribution
set.seed(2971)
N <- 100
a <- rnorm(N, 178, 20)
b <- rlnorm(N, 0, 1)

plot(NULL, xlim = range(d2$weight), ylim = c(-100, 400), 
     xlab = "weight", ylab = "height")
abline(h = 0, lty = 2)
abline(h = 272, lty = 1, lwd = 0.5)
mtext("b ~ dnorm(0, 10)")
xbar <- mean(d2$weight)
for(i in 1:N) curve(a[i] + b[i] * (x - xbar), 
                    from = min(d2$weight), to = max(d2$weight),  
                    add = TRUE, col = col.alpha("black", 0.2))

# Find the posterior
# Define the average weight
xbar <- mean(d2$weight)
# Fit the model 
m4.3 <- quap(
  alist(
    height ~ dnorm(mu, sigma), 
    mu <- a + b * (weight - xbar), 
    a ~ dnorm(178, 20), 
    b ~ dlnorm(0, 1), 
    sigma ~ dunif(0, 50)
  ), data = d2)
# Fit another model 
m4.3b <- quap(
  alist(
    height ~ dnorm(mu, sigma), 
    mu <- a + exp(log_b) * (weight - xbar), 
    a ~ dnorm(178, 20), 
    log_b ~ dnorm(0, 1), 
    sigma ~ dunif(0, 50)
  ), data = d2)
# Summarize posterior (tables of marginal distributions)
precis(m4.3)
round(vcov(m4.3), 3)

# Plot posterior inference against the data
plot(height ~ weight, data = d2, col = rangi2)
post <- extract.samples(m4.3)
a_map <- mean(post$a)
b_map <- mean(post$b)
curve(a_map + b_map * (x - xbar), add = TRUE)

# Adding uncertainty around the mean
post <- extract.samples(m4.3)
post[1:5, ]
# Re-estimate the model
N <- 10
dN <- d2[1:N, ]
mN <- quap(
  alist(
    height ~ dnorm(mu, sigma), 
    mu <- a + b * (weight - mean(weight)), 
    a ~ dnorm(178, 20), 
    b ~ dlnorm(0, 1), 
    sigma ~ dunif(0, 50)
  ), data = dN)
# Plot 20 of these lines 
# Extract samples from the posterior
post <- extract.samples(mN, n = 20)
# Display raw data and sample size
plot(dN$weight, dN$height, 
     xlim = range(d2$weight), ylim = range(d2$height), 
     col = rangi2, xlab = "weight", ylab = "height")
mtext(concat("N = ", N))
# Plot the lines, with transparency
for(i in 1:20)
  curve(post$a[i] + post$b[i] * (x - mean(dN$weight)), 
    col = col.alpha("black", 0.3), add = TRUE)

# Compute and plot the density of the vector of predicted means
post <- extract.samples(m4.3)
mu_at_50 <- post$a + post$b * (50 - xbar)
dens(mu_at_50, col = rangi2, lwd = 2,  xlab = "mu|weight=50")
PI(mu_at_50, prob = 0.89)

# Repeat the above calculation using link function
mu <- link(m4.3)
str(mu)

# Compute a distribution of mu for each weight value
# Define sequence of weights to compute predictions for 
# These values will be on the horizontal axis
weight.seq <- seq(from = 25, to = 70, by = 1)
# Use link to compute mu 
# for each sample from posterior
# and for each weight in weight.seq
mu <- link(m4.3, data = data.frame(weight = weight.seq))
str(mu)

# Plot a distribution for each height
# Use type = "n" to hide raw data
plot(height ~ weight, d2, type = "n")
# Loop over samples and plot each mu value
for(i in 1:100)
  points(weight.seq, mu[i, ], pch = 16, col = col.alpha(rangi2, 0.1))

# Summarize the distribution of mu
mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI, prob = 0.89)

# Plot the summary tables
# Plot the raw data
# Fading out points to make line and interval more visible 
plot(height ~ weight, data = d2, col = col.alpha(rangi2, 0.5))
# Plot the MAP line, aka mean mu for each weight
lines(weight.seq, mu.mean)
# Plot a shaded region for 89% PI
shade(mu.PI, weight.seq)

# How link works
# Repeat the same without the link function
post <- extract.samples(m4.3)
mu.link <- function(weight) post$a + post$b * (weight - xbar)
weight.seq <- seq(from = 25, to = 70, by = 1)
mu <- sapply(weight.seq, mu.link)
mu.mean <- apply(mu, 2, mean)
mu.CI <- apply(mu, 2, PI, prob = 0.89)

# Simulate prediction intervals
sim.height <- sim(m4.3, data = list(weight = weight.seq))
str(sim.height)
# Summarize it
height.PI <- apply(sim.height, 2, PI, prob = 0.89)
# Create a summary plot
# Plot raw data
plot(height ~ weight, d2, col = col.alpha(rangi2, 0.5))
# Draw MAP line
lines(weight.seq, mu.mean)
# Draw HPDI version for line
shade(mu.CI, weight.seq)
# Draw PI region for simulated weights
shade(height.PI, weight.seq)

# Increase the number of samples to get rid of rough intervals
sim.height <- sim(m4.3, data = list(weight = weight.seq), n = 1e4)
height.PI <- apply(sim.height, 2, PI, prob = 0.89)

# How sim works
post <- extract.samples(m4.3)
weight.seq <- 25:70
sim.height <- sapply(weight.seq, function(weight)
  rnorm(
    n = nrow(post), 
    mean = post$a + post$b * (weight - xbar), 
    sd = post$sigma))
height.PI <- apply(sim.height, 2, PI, prob = 0.89)

### Curves from lines: 4.64-










































