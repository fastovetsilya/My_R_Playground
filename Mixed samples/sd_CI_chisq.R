mu<-0
sigma<-1
alpha=0.1
N<-10^4
lower.conf<-numeric(N)
upper.conf<-numeric(N)
for(i in 1:N)
{
x<-rnorm(8,mu,sigma)
s<-sd(x)
n<-length(x)
lower.conf[i]<-sqrt((n-1)*s^2/qchisq(alpha/2,n-1,lower.tail = FALSE))
upper.conf[i]<-sqrt((n-1)*s^2/qchisq(1-alpha/2,n-1,lower.tail = FALSE))
}
length(upper.conf[upper.conf<sigma])/N
hist(lower.conf)
hist(upper.conf)
sd(upper.conf)

