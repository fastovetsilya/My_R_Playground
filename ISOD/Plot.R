deltafun <- function (x,c)
{
  if (x < c) y = 0
  if (x >= c) y = 1
  return(y)
}
x = seq(0,1,0.001)
y1 = numeric(length(x))
for (i in 1:length(y1)) y1[i] = deltafun(x[i], 0.3)

plot(x, y1, type ='l', xlab = 'Коэффициент корреляции', 
     ylab = 'Коэффициент знака обеспеченности', cex.lab = 1.5,
     cex.axis = 1.5, lwd = 2)

curve(x^2, xlab = 'Коэффициент корреляции', ylab = 'Коэффициент знака обеспеченности'
      ,cex.lab = 1.5, cex.axis = 1.5, add = T, lty=2, lwd = 2)
