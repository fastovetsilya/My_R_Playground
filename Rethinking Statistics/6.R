library(rethinking)
library(dagitty)

# Simulated science distortion
set.seed(1914)
N <- 200 # Num grant proposals
p <- 0.1 # Proportion to select
# uncorrelated newsworthiness and trustworthiness
nw <- rnorm(N)
tw <- rnorm(N)
# Select top 10% of combined scores
s <- nw + tw # total score
q <- quantile(s, 1 - p) # Top 10% threshold
selected <- ifelse(s >= q, TRUE, FALSE)
cor(tw[selected], nw[selected])

### Multicollinearity: 6.2 - 6.12
# Simulate height and leg length of 100 individuals
N <- 100 # Number of individuals
set.seed(909)
height <- rnorm(N, 10, 2) # Sim total height of each
leg_prop <- runif(N, 0.4, 0.5) # Leg as a proportion of height
leg_left <- leg_prop * height + # Sim left leg as a propotion + error
  rnorm(N, 0, 0.02) 
leg_right <- leg_prop * height + # Sim right leg as a proportion + error
  rnorm(N, 0, 0.02)
d <- data.frame(height, leg_left, leg_right)

# Fit the model 
m6.1 <- quap(
  alist(
    height ~ dnorm(mu, sigma), 
    mu <- a + bl * leg_left + br * leg_right, 
    a ~ dnorm(10, 100), 
    bl ~ dnorm(2, 10), 
    br ~ dnorm(2, 10), 
    sigma ~ dexp(1)
  ), data = d)
precis(m6.1)
plot(precis(m6.1))

# Look at the posterior
post <- extract.samples(m6.1)
plot(bl ~ br, post, col = col.alpha(rangi2, 0.1), pch = 16)

# Compute and plot the posterior of the sum of the predictors
sum_blbr <- post$bl + post$br
dens(sum_blbr, col = rangi2, lwd = 2, xlab = "sum of bl and br")

# Fit the regression with one predictor variable 
m6.2 <- quap(
  alist(
    height ~ dnorm(mu, sigma), 
    mu <- a + bl * leg_left, 
    a ~ dnorm(10, 100), 
    bl ~ dnorm(2, 10), 
    sigma ~ dexp(1)
  ), data = d)
precis(m6.2)

# Multicollinear milk 
data(milk)
d <- milk
d$K <- standardize(d$kcal.per.g)
d$F <- standardize(d$perc.fat)
d$L <- standardize(d$perc.lactose)

# Build bivariate regression models
# kcal.per.g regressed on perc.fat
m6.3 <- quap(
  alist(
    K ~ dnorm(mu, sigma), 
    mu <- a + bF * F, 
    a ~ dnorm(0, 0.2), 
    bF ~ dnorm(0, 0.5), 
    sigma ~ dexp(1)
  ), data = d)
# kcal.per.g regressed on perc.lactose
m6.4 <- quap(
  alist(
    K ~ dnorm(mu, sigma), 
    mu <- a + bL * L, 
    a ~ dnorm(0, 0.2), 
    bL ~ dnorm(0, 0.5), 
    sigma ~ dexp(1)
  ), data = d)

precis(m6.3)
precis(m6.4)

# Build the models with both predictors
m6.5 <- quap(
  alist(
    K ~ dnorm(mu, sigma), 
    mu <- a + bF * F + bL * L, 
    a ~ dnorm(0, 0.2), 
    bF ~ dnorm(0, 0.5), 
    bL ~ dnorm(0, 0.5), 
    sigma ~ dexp(1)
  ), data = d)
precis(m6.5)

# Create a pairs plot
pairs(~ kcal.per.g + perc.fat + perc.lactose, data = d, col = rangi2)

# Simulate collinearity
data(milk)
d <- milk
sim.coll <- function(r = 0.9) {
  d$x <- rnorm(nrow(d), mean = r * d$perc.fat, 
      sd = sqrt(1 - r ^ 2))
  m <- lm(kcal.per.g ~ perc.fat + x, data = d)
  sqrt(diag(vcov(m)))[2] # stddev of parameter
}
rep.sim.coll <- function(r = 0.9, n = 100) {
  stddev <- replicate(n, sim.coll(r))
  mean(stddev)
}
r.seq <- seq(from = 0, to = 0.99, by = 0.01)
stddev <- sapply(r.seq, function(z) rep.sim.coll(r = z, n = 100))
plot(stddev ~ r.seq, type = "l", col = rangi2, lwd = 2, xlab = "correlation")


### Post-treatment bias: 6.13 - 

# Simulate some data with post-treatment variable
set.seed(71)
# Number of plants
N <- 100
# Simulate initial heights
h0 <- rnorm(N, 10, 2)
# Assign treatments and simulate fungus and growth
treatment <- rep(0:1, each = N / 2)
fungus <- rbinom(N, size = 1, prob = 0.5 - treatment * 0.4)
h1 <- h0 + rnorm(N, 5 - 3 * fungus)
# Compose a clean data frame
d <- data.frame(h0 = h0, h1 = h1, treatment = treatment, fungus = fungus)
precis(d)

# Simulated prior distribution
sim_p <- rlnorm(1e4, 0, 0.25)
precis(data.frame(sim_p))
# Fit the model 
m6.6 <- quap(
  alist(
    h1 ~ dnorm(mu, sigma), 
    mu <- h0 * p, 
    p ~ dnorm(0, 0.25), 
    sigma ~ dexp(1)
  ), data = d)
precis(m6.6)

# Include treatment and fungus variables
m6.7 <- quap(
  alist(
    h1 ~ dnorm(mu, sigma), 
    mu <- h0 * p, 
    p <- a + bt * treatment + bf * fungus, 
    a ~ dlnorm(0, 0.2), 
    bt ~ dnorm(0, 0.5), 
    bf ~ dnorm(0, 0.5), 
    sigma ~ dexp(1)
  ), data = d)
precis(m6.7)

# Omit post-trearment variable fungus
m6.8 <- quap(
  alist(
    h1 ~ dnorm(mu, sigma),  
    mu <- h0 * p, 
    p <- a + bt * treatment, 
    a ~ dlnorm(0, 0.2), 
    bt ~ dnorm(0, 0.5), 
    sigma ~ dexp(1)
  ), data = d)
precis(m6.8)

# Fungus and d-separation
plant_dag <- dagitty("dag {
                     H_0 -> H_1
                     F -> H_1
                     T -> F
                     }")
coordinates(plant_dag) <- list(x = c(H_0 = 0, T = 2, F = 1.5, H_1 = 1), 
                               y = c(H_0 = 0, T = 0, F = 0, H_1 = 0))
drawdag(plant_dag)













































































































































































































































