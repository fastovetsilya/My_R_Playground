library(rethinking)
library(dagitty)

### Spurious association: 5.1-5.6 

# Load the data and standardize the variables 
data(WaffleDivorce)
d <- WaffleDivorce
d$D <- standardize(d$Divorce)
d$M <- standardize(d$Marriage)
d$A <- standardize(d$MedianAgeMarriage)

# Create a linear model 
sd(d$ MedianAgeMarriage)
m5.1 <- quap(
  alist(
    D ~ dnorm(mu, sigma), 
    mu <- a + bA * A, 
    a ~ dnorm(0, 0.2), 
    bA ~ dnorm(0, 0.5), 
    sigma ~ dexp(1)
  ), data = d)

# Simulate from the priors
set.seed(10)
prior <- extract.prior(m5.1)
mu <- link(m5.1, post = prior, data = list(A = c(-2, 2)))
plot(NULL, xlim = c(-2, 2), ylim = c(-2, 2), 
     xlab = "Median age marriage (std)", 
     ylab = "Divorce rate (std)")
for(i in 1:50) lines(c(-2, 2), mu[i, ], col = col.alpha("black", 0.4))

# Posterior predictions 
# Compute percentile interval of mean
A_seq <- seq(from = -3, to = 3.2, length.out = 30)
mu <- link(m5.1, data = list(A = A_seq))
mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI)
# Plot it all 
plot(D ~ A, data = d, col = rangi2)
lines(A_seq, mu.mean, lwd = 2)

# Fit similar model for the relationship with Marriage rate
m5.2 <- quap(
  alist(
    D ~ dnorm(mu, sigma), 
    mu <- a + bM * M, 
    a ~ dnorm(0, 0.2), 
    bM ~ dnorm(0, 0.5), 
    sigma ~ dexp(1)
  ), data = d)
# Compute percentile interval of mean
M_seq <- seq(from = -3, to = 3.2, length.out = 30)
mu <- link(m5.2, data = list(M = M_seq))
mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI)
# Plot it all 
plot(D ~ M, data = d, col = rangi2)
lines(M_seq, mu.mean, lwd = 2)

### Directed acyclic graphs
# Drawing a DAG
dag5.1 <- dagitty("dag{A -> D; A -> M; M -> D}")
coordinates(dag5.1) <- list(x = c(A = 0, D = 1, M = 2),
                            y = c(A = 0, D = 1, M = 0))
drawdag(dag5.1)

# Find the implications for DAG
DMA_dag2 <- dagitty("dag{D <- A -> M}")
impliedConditionalIndependencies(DMA_dag2)
DMA_dag1 <- dagitty("dag{D <- A -> M -> D}")
impliedConditionalIndependencies(DMA_dag1)


### Multiple regression: 5.10-5.27
# Approximate posterior with quap
m5.3 <- quap(
  alist(
    D ~ dnorm(mu, sigma), 
    mu <- a + bM * M + bA * A, 
    a ~ dnorm(0, 0.2), 
    bM ~ dnorm(0, 0.5), 
    bA ~ dnorm(0, 0.5), 
    sigma ~ dexp(1)
  ), data = d)
precis(m5.3)

# Plot
plot(coeftab(m5.1, m5.2, m5.3), par = c("bA", "bM"))

# Simulating the divorce example
N <- 50
age <- rnorm(N)
mar <- rnorm(N, -age)
div <- rnorm(N, age)

# Predictor residual plots
# Fit the model
m5.4 <- quap(
  alist(
    M ~ dnorm(mu, sigma), 
    mu <- a + bAM * A, 
    a ~ dnorm(0, 0.2), 
    bAM ~ dnorm(0, 0.5), 
    sigma ~ dexp(1)
  ), data = d)
# Compute residuals
mu <- link(m5.4)
mu_mean <- apply(mu, 2, mean)
mu_resid <- d$M - mu_mean

# Posterior predictive check
# Simulate predictions, averaging over posterior
# Call link without specifying new data, so it uses original data
mu <- link(m5.3)
# Summarize samples across cases 
mu_mean <- apply(mu, 2, mean)
mu_PI <- apply(mu, 2, PI)
# Simulate observations
# Again no new data, use original data
D_sim <- sim(m5.3, n=1e4)
D_PI <- apply(D_sim, 2, PI)
# Plot predictions against observed
plot(mu_mean ~ d$D, col = rangi2, ylim = range(mu_PI), 
     xlab = "Observed divorce", ylab = "Predicted divorce")
abline(a = 0, b = 1, lty = 2)
for(i in 1:nrow(d)) lines(rep(d$D[i], 2), mu_PI[, i], col = rangi2)
# Identify the outlier states
identify(x = d$D, y = mu_mean, labels = d$Loc)

# Simulate spurious association
N <- 100
x_real <- rnorm(N)
x_spur <- rnorm(N, x_real)
y <- rnorm(N, x_real)
d <- data.frame(y, x_real, x_spur)

# Counterfactual plots
data(WaffleDivorce)
d <- list()
d$A <- standardize(WaffleDivorce$MedianAgeMarriage)
d$D <- standardize(WaffleDivorce$Divorce)
d$M <- standardize(WaffleDivorce$Marriage)
# Run two regressions at the same time
m5.3_A <- quap(
  alist(
    ## A -> D <- M
    D ~ dnorm(mu, sigma), 
    mu <- a + bM * M + bA * A, 
    a ~ dnorm(0, 0.2), 
    bM ~ dnorm(0, 0.5), 
    bA ~ dnorm(0, 0.5), 
    sigma ~ dexp(1), 
    ## A -> M
    M ~ dnorm(mu_M, sigma_M), 
    mu_M <- aM + bAM * A, 
    aM <- dnorm(0, 0.2), 
    bAM ~ dnorm(0, 0.5), 
    sigma_M ~ dexp(1)
  ), data = d)
precis(m5.3_A)

# Manipulate A
A_seq <- seq(from = -2, to = 2, length.out = 30)
# Prep data
sim_dat <- data.frame(A = A_seq)
# Simulate M and then D, using A_seq
s <- sim(m5.3_A, data = sim_dat, vars = c("M", "D"))
# Plot the predictions
plot(sim_dat$A, colMeans(s$D), ylim = c(-2, 2), type = "l", 
     xlab = "Manipulated A", ylab = "Counterfactual D")
shade(apply(s$D, 2, PI), sim_dat$A)
mtext("Total counterfactual effect of A on D")

# Expected causal effect of increasing median age at marriage from 20 to 30
# New data frame, standardized to mean 26.1 and std dev 1.24
sim2_dat <- data.frame(A = (c(20, 30) - 26.1) / 1.24)
s2 <- sim(m5.3_A, data = sim2_dat, vars = c("M", "D"))
mean(s2$D[, 2] - s2$D[, 1])

# Simulate a counterfactual for an average state A = 0, and see what changing M does
sim_dat <- data.frame(M = seq(from = -2, to = 2, length.out = 30), A = 0)
s <- sim(m5.3_A, data = sim_dat, vars = "D")

plot(sim_dat$M, colMeans(s), ylim = c(-2, 2), type = "l", 
     xlab = "manipulated M", ylab = "counterfactual D")
shade(apply(s, 2, PI), sim_dat$M)
mtext("Total counterfactual effects of M on D")

# Simulating counterfactuals
# Define a range of values to assign to A
A_seq <- seq(from = -2, to = 2, length.out = 30)
# Extract the posterior samples
post <- extract.samples(m5.3_A)
M_sim <- with(post, sapply(1:30, 
              function(i) rnorm(1e3, aM + bAM * A_seq[i], sigma_M)))
# Simulate D
D_sim <- with(post, sapply(1:30, 
                           function(i) rnorm(1e3, a + bA * A_seq[i] + bM * M_sim[, i], sigma)))


### Masked relationships: 5.28 - 
# Load milk data
data(milk)
d <- milk
str(d)

# Standardize variables
d$K <- standardize(d$kcal.per.g)
d$N <- standardize(d$neocortex.perc)
d$M <- standardize(log(d$mass))

# Fit quap model 
m5.5_draft <- quap(
  alist(
    K ~ dnorm(mu, sigma), 
    mu <- a + bN * N, 
    a ~ dnorm(0, 1), 
    bN ~ dnorm(0, 1), 
    sigma ~ dexp(1)
  ), data = d)

# Remove missing rows
dcc <- d[complete.cases(d$K, d$N, d$M), ]

# Fit quap model on the new dataframe
m5.5_draft <- quap(
  alist(
    K ~ dnorm(mu, sigma), 
    mu <- a + bN * N, 
    a ~ dnorm(0, 1), 
    bN ~ dnorm(0, 1), 
    sigma ~ dexp(1)
  ), data = dcc)

# Simulate and plot 50 prior regression lines
prior <- extract.prior(m5.5_draft)
xseq <- c(-2, 2)
mu <- link(m5.5_draft, post = prior, data = list(N = xseq))
plot(NULL, xlim = xseq, ylim = xseq)
for(i in 1:50) lines(xseq, mu[i, ], col = col.alpha("black", 0.3))

# Build a new model with another prior
m5.5 <- quap(
  alist(
    K ~ dnorm(mu, sigma), 
    mu <- a + bN * N, 
    a ~ dnorm(0, 0.2), 
    bN ~ dnorm(0, 0.5), 
    sigma ~ dexp(1)
  ), data = dcc
)

# Look at the posterior
precis(m5.5)

























































































































































