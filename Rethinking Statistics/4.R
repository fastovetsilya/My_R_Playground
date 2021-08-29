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


### Building linear model

# Plot adult height and weight against one another
data(Howell1)
d <- Howell1
d2 <- d[d$age >= 18, ]
plot(d2$height, d2$weight)






























