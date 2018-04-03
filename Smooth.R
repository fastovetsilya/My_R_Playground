library('readxl')
Data<-read_excel('Smooth.xlsx')
m.deriv<-cbind(Data$ml[2:length(Data$ml)], diff(Data$pH))
spline<-smooth.spline(m.deriv[,1],m.deriv[,2],cv=TRUE)
kernel.box<-ksmooth(m.deriv[,1],m.deriv[,2],kernel='box',n.points=10^2)
kernel.norm<-ksmooth(m.deriv[,1],m.deriv[,2],kernel='normal',n.points=10^2)


sol.spline<-numeric(length(predict(spline, deriv=1)$x))
for (i in 2:length(predict(spline, deriv=1)$x))
{
fit<-lm(predict(spline,deriv=1)$y[(i-1):i] ~ predict(spline, deriv=1)$x[(i-1):i])
y.0<-as.numeric(-coef(fit)[1]/coef(fit)[2])
if(y.0>=predict(spline, deriv=1)$x[(i-1)] && y.0<=predict(spline, deriv=1)$x[i])
{
sol.spline[i]<-y.0
}
}
sol.spline[sol.spline!=0]

sol.kern<-numeric(length(diff(kernel.norm$y)))
for (i in 2:length(diff(kernel.norm$y)))
{
  fit.k<-lm(diff(kernel.norm$y)[(i-1):i] ~ kernel.norm$x[i:(i+1)])
  y.0k<-as.numeric(-coef(fit.k)[1]/coef(fit.k)[2])
  if(y.0k>=kernel.norm$x[i] && y.0k<=kernel.norm$x[(i+1)])
  {
    sol.kern[i]<-y.0k
  }
}
sol.kern[sol.kern!=0]

plot(Data$ml,Data$pH, pch='.',cex=3, main='integral plot',xlab='ml',ylab='pH')
plot(m.deriv[,1],m.deriv[,2],pch='.',cex=3, main='derivative plot', xlab='ml',ylab='^pH')
lines(spline, col='red')
lines(kernel.box, col='green')
lines(kernel.norm, col='blue')
abline(v=sol.spline[sol.spline!=0],col='red',lty=2)
abline(v=sol.kern[sol.kern!=0],col='blue',lty=2)
legend('topleft',c('Cubic spline','Nadaraya-Watson, Box Kernel','Nadaraya-Watson, Normal Kernel'),lty=c(1,1,1),lwd=c(1,1,1),col=c('red','green','blue'),cex=0.5)
plot(predict(spline, deriv=1)$x, predict(spline,deriv=1)$y, type='o', 
     main='second derivative plot:spline', xlab='ml', ylab='^^pH')
abline(v=sol.spline[sol.spline!=0],col='red',lty=2)
abline(h=0, col='green')
plot(kernel.norm$x[2:length(kernel.norm$x)],diff(kernel.norm$y), type='o',
     main='sec deriv plot: normal kernel (Nadaraya-Watson)', xlab='ml', ylab='^^pH')
abline(v=sol.kern[sol.kern!=0],col='blue',lty=2)
abline(h=0, col='green')


if(length(sol.spline[sol.spline!=0])>10)
{
  
rm(list=ls())
  
library('readxl')
Data<-read_excel('Smooth.xlsx')
m.deriv<-cbind(Data$ml[2:length(Data$ml)], diff(Data$pH))
sol.spline<-numeric()

for (k in seq(from=100, to=2, by=-1))
{
spline<-smooth.spline(m.deriv[,1],m.deriv[,2],cv=FALSE,df=k)
sol.spline<-numeric(length(predict(spline, deriv=1)$x))
for (i in 2:length(predict(spline, deriv=1)$x))
{
  fit<-lm(predict(spline,deriv=1)$y[(i-1):i] ~ predict(spline, deriv=1)$x[(i-1):i])
  y.0<-as.numeric(-coef(fit)[1]/coef(fit)[2])
  if(y.0>=predict(spline, deriv=1)$x[(i-1)] && y.0<=predict(spline, deriv=1)$x[i])
  {
    sol.spline[i]<-y.0
  }

}
if(length(sol.spline[sol.spline!=0])<10) break
}
kernel.norm<-ksmooth(m.deriv[,1],m.deriv[,2],kernel='normal',n.points=10^2)
kernel.box<-ksmooth(m.deriv[,1],m.deriv[,2],kernel='box',n.points=10^2)
sol.kern<-numeric(length(diff(kernel.norm$y)))
for (i in 2:length(diff(kernel.norm$y)))
{
  fit.k<-lm(diff(kernel.norm$y)[(i-1):i] ~ kernel.norm$x[i:(i+1)])
  y.0k<-as.numeric(-coef(fit.k)[1]/coef(fit.k)[2])
  if(y.0k>=kernel.norm$x[i] && y.0k<=kernel.norm$x[(i+1)])
  {
    sol.kern[i]<-y.0k
  }
}
sol.kern[sol.kern!=0]
plot(Data$ml,Data$pH, pch='.',cex=3, main='integral plot',xlab='ml',ylab='pH')
plot(m.deriv[,1],m.deriv[,2],pch='.',cex=3, main='corrected derivative plot', xlab='ml',ylab='^pH')
lines(spline, col='red')
lines(kernel.box, col='green')
lines(kernel.norm, col='blue')
abline(v=sol.spline[sol.spline!=0],col='red',lty=2)
abline(v=sol.kern[sol.kern!=0],col='blue',lty=2)
legend('topleft',c('Cubic spline','Nadaraya-Watson, Box Kernel','Nadaraya-Watson, Normal Kernel'),lty=c(1,1,1),lwd=c(1,1,1),col=c('red','green','blue'),cex=0.5)
plot(predict(spline, deriv=1)$x, predict(spline,deriv=1)$y, type='o', 
     main='corrected second derivative plot:spline', xlab='ml', ylab='^^pH')
abline(v=sol.spline[sol.spline!=0],col='red',lty=2)
abline(h=0, col='green')
plot(kernel.norm$x[2:length(kernel.norm$x)],diff(kernel.norm$y), type='o',
     main='corrected sec deriv plot: Nadaraya-Watson', xlab='ml', ylab='^^pH')
abline(v=sol.kern[sol.kern!=0],col='blue',lty=2)
abline(h=0, col='green')
}
sol.spline[sol.spline!=0]
sol.kern[sol.kern!=0]









