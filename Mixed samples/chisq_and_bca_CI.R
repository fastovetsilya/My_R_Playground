#library('boot')
mu<-0
sigma<-1
alpha=0.6
x<-c(9.55,10.31,11.55,8.81)
s<-sd(x)
n<-length(x)
sqrt((n-1)*s^2/qchisq(alpha/2,n-1,lower.tail = FALSE))#upper chi
sqrt((n-1)*s^2/qchisq(1-alpha/2,n-1,lower.tail = FALSE))#lower.chi

std<-function(data,i)
{
  return(mad(data[i]))
}
bootresult<-boot(x,std,R=10000)
boot.ci(bootresult, conf = 1-alpha,type='bca')