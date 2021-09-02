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







































































