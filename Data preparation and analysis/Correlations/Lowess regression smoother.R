library('readxl')

Data <- read_excel('Correlations_La_Ce_rastv_pochv.xlsx', sheet = 1)
Datamat <- as.matrix(Data)
Datasample <- na.exclude(Datamat[, c('Ce_rastv_RootLength', 'Ce_rastv_FAC')])

# Cross-validation function (not used here)
# https://stats.stackexchange.com/questions/2002/how-do-i-decide-what-span-to-use-in-loess-regression-in-r
locv1 <- function(x1, y1, nd, span, ntrial)
{
  locvgcv <- function(sp, x1, y1)
  {
    nd <- length(x1)
    
    assign("data1", data.frame(xx1 = x1, yy1 = y1))
    fit.lo <- loess(yy1 ~ xx1, data = data1, span = sp, family = "gaussian", degree = 2, surface = "direct")
    res <- residuals(fit.lo)
    
    dhat2 <- function(x1, sp)
    {
      nd2 <- length(x1)
      diag1 <- diag(nd2)
      dhat <- rep(0, length = nd2)
      
      for(jj in 1:nd2){
        y2 <- diag1[, jj]
        assign("data1", data.frame(xx1 = x1, yy1 = y2))
        fit.lo <- loess(yy1 ~ xx1, data = data1, span = sp, family = "gaussian", degree = 2, surface = "direct")
        ey <- fitted.values(fit.lo)
        dhat[jj] <- ey[jj]
      }
      return(dhat)
    }
    
    dhat <- dhat2(x1, sp)
    trhat <- sum(dhat)
    sse <- sum(res^2)
    
    cv <- sum((res/(1 - dhat))^2)/nd
    gcv <- sse/(nd * (1 - (trhat/nd))^2)
    
    return(gcv)
  }
  
  gcv <- lapply(as.list(span1), locvgcv, x1 = x1, y1 = y1)
  #cvgcv <- unlist(cvgcv)
  #cv <- cvgcv[attr(cvgcv, "names") == "cv"]
  #gcv <- cvgcv[attr(cvgcv, "names") == "gcv"]
  
  return(gcv)
}


# Performing cross-validation (custom)

x <- Datasample[, 'Ce_rastv_RootLength']
y <- Datasample[, 'Ce_rastv_FAC']
df <- data.frame(x, y)
set.seed(1)

span.seq <- seq(from = 0.5, to = 0.95, by = 0.01) # explores range of spans
k <- 10 # number of folds
folds <- sample(x = 1:k, size = length(x), replace = TRUE)
cv.error.mtrx <- matrix(rep(x = NA, times = k * length(span.seq)), 
                        nrow = length(span.seq), ncol = k)

for(i in 1:length(span.seq)) {
  for(j in 1:k) {
    loess.fit <- loess(formula = y ~ x, data = df[folds != j, ], span = span.seq[i])
    preds <- predict(object = loess.fit, newdata = df[folds == j, ])
    cv.error.mtrx[i, j] <- mean((df$y[folds == j] - preds)^2, na.rm = TRUE)
    # some predictions result in `NA` because of the `x` ranges in each fold
  }
}

cv.errors <- rowMeans(cv.error.mtrx)

best.span.i <- which.min(cv.errors)
best.span.i
span.seq[best.span.i]

plot(x = span.seq, y = cv.errors, type = "l", main = "CV Plot")
points(x = span.seq, y = cv.errors, 
       pch = 20, cex = 0.75, col = "blue")
points(x = span.seq[best.span.i], y = cv.errors[best.span.i], 
       pch = 20, cex = 1, col = "red")

best.loess.fit <- loess(formula = y ~ x, data = df) #,span.seq[best.span.i])

x.seq <- seq(from = min(x), to = max(x), length = 100)

plot(x = df$x, y = df$y, main = "Ce_rastvor", 
     xlab = 'Root length', ylab = 'FAC')

best.loess.fit_pred <- predict(object = best.loess.fit, 
                               newdata = data.frame(x = x.seq), se=T)
lines(x = x.seq, y = best.loess.fit_pred$fit, 
      col = "red", lwd = 2)

lines(x = x.seq, y = best.loess.fit_pred$fit + qt(0.975,best.loess.fit_pred$df) * 
        best.loess.fit_pred$se, col = "blue", lty = 2)

lines(x = x.seq, y = best.loess.fit_pred$fit - qt(0.975,best.loess.fit_pred$df) *
        best.loess.fit_pred$se, col = "blue", lty = 2)


