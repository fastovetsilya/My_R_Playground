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


### Post-treatment bias: 6.13 - 6.20

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

impliedConditionalIndependencies(plant_dag)

#Include F in the model, then rerun the models above
set.seed(71)
N <- 1000
h0 <- rnorm(N, 10, 2)
treatment <- rep(0:1, each = N / 2)
M <- rbern(N)
fungus <- rbinom(N, size = 1, prob = 0.5 - treatment * 0.4 + 0.4 * M)
h1 <- h0 + rnorm(N, 5 + 3 * M)
d2 <- data.frame(h0 = h0, h1 = h1, treatment = treatment, fungus = fungus)

m6.7.d2 <- quap(
  alist(
    h1 ~ dnorm(mu, sigma), 
    mu <- h0 * p, 
    p <- a + bt * treatment + bf * fungus, 
    a ~ dlnorm(0, 0.2), 
    bt ~ dnorm(0, 0.5), 
    bf ~ dnorm(0, 0.5), 
    sigma ~ dexp(1)
  ), data = d2)
precis(m6.7.d2)

m6.8.d2 <- quap(
  alist(
    h1 ~ dnorm(mu, sigma),  
    mu <- h0 * p, 
    p <- a + bt * treatment, 
    a ~ dlnorm(0, 0.2), 
    bt ~ dnorm(0, 0.5), 
    sigma ~ dexp(1)
  ), data = d2)
precis(m6.8.d2)

### Collider bias: 6.21 - 6.28

d <- sim_happiness(seed = 1977, N_years = 1000)
precis(d)
# Rescale data
d2 <- d[d$age>17, ] # Only adults
d2$A <- (d2$age - 18) / (65 - 18)
# Approximate the posterior
d2$mid <- d2$married + 1
m6.9 <- quap(
  alist(
    happiness ~ dnorm(mu, sigma), 
    mu <- a[mid] + bA * A, 
    a[mid] ~ dnorm(0, 1), 
    bA ~ dnorm(0, 2), 
    sigma ~ dexp(1)
  ), data = d2
)
precis(m6.9, depth = 2)
# Model that omits marriage status
m6.10 <- quap(
  alist(
    happiness ~ dnorm(mu, sigma), 
    mu <- a + bA * A, 
    a ~ dnorm(0, 1), 
    bA ~ dnorm(0, 2), 
    sigma ~ dexp(1)
  ), data = d2
)
precis(m6.10)

# Simulate unobserved influence (grandparents)
N <- 200 # Number of grandparent-parent-child triads
b_GP <- 1 # Direct effect of G on P
b_GC <- 0 # Direct effect of G on C
b_PC <- 1 # Direct effect of P on C
b_U <- 2 # Direct effect of U on P and C
# Draw random observations
set.seed(1)
U <- 2 * rbern(N, 0.5) - 1
G <- rnorm(N)
P <- rnorm(N, b_GP * G + b_U * U)
C <- rnorm(N, b_PC * P + b_GC * G + b_U * U)
d <- data.frame(C = C, P = P, G = G, U = U)
# Build the model
m6.11 <- quap(
  alist(
    C ~ dnorm(mu, sigma), 
    mu <- a + b_PC * P + b_GC * G, 
    a ~ dnorm(0, 1),  
    c(b_PC, b_GC) ~ dnorm(0, 1), 
    sigma ~ dexp(1)
  ), data = d
)
precis(m6.11)
# Regression that conditions also on U
m6.12 <- quap(
  alist(
    C ~ dnorm(mu, sigma), 
    mu <- a + b_PC * P + b_GC * G + b_U * U, 
    a ~ dnorm(0, 1), 
    c(b_PC, b_GC, b_U) ~ dnorm(0, 1), 
    sigma ~ dexp(1)
  ), data = d
)
precis(m6.12)

### Confronting confounding: 6.29 - 6.31
# Block the backdoor
dag_6.1 <- dagitty("dag{
    U [unobserved]
    X -> Y
    X <- U <- A -> C -> Y
    U -> B <- C
}")
adjustmentSets(dag_6.1, exposure = "X", outcome = "Y")
# Trace backwards, starting at W and ending at D
dag_6.2 <- dagitty("dag {
    A -> D
    A -> M -> D
    A <- S -> M
    S -> W -> D
}")
adjustmentSets(dag_6.2, exposure = "W", outcome = "D")
impliedConditionalIndependencies(dag_6.2)
