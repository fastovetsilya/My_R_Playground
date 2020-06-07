# Different combinations of nutrient elements 
x <- matrix(1:24, ncol=6)
cmb <- combn(ncol(x), 2)
r1 <- apply(cmb, 2, function(j) x[, j[1]]/x[, j[2]])
r2 <- apply(cmb, 2, function(j) x[, j[2]]/x[, j[1]])
cbind(r1, r2)

# Another way
pairwise.ratios <- function(x, prefix="probeset", char=":"){
  n <- ncol(x)
  cn <- colnames(x)
  if(length(cn) == 0){
    cn <- gsub(" ", "0", formatC(seq.int(n), width=nchar(n)))
    cn <- paste(prefix, cn, sep="")
  }
  cmb <- combn(n, 2)
  r1 <- apply(cmb, 2, function(j) x[, j[1]]/x[, j[2]])
  r2 <- apply(cmb, 2, function(j) x[, j[2]]/x[, j[1]])
  colnames(r1) <- apply(cmb, 2, function(j) paste(cn[j], collapse=char))
  colnames(r2) <- apply(cmb, 2, function(j) paste(cn[rev(j)], collapse=char))
  cbind(r1, r2)[, order(c(colnames(r1), colnames(r2)))]
} 
m1 <- matrix(1:24, ncol=6)
pairwise.ratios(m1)
