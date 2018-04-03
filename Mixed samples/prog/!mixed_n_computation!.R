a<-c(0.21,0.10,0.13,0.56)#reference sample
no_of_comparisons<-6 #how many comparisons?
alpha<-0.05 #at what alpha to compare??
difference<-1.05 #difference in means
n<-length(a)
alpha<-1-(1-alpha)^(1/no_of_comparisons)#(Sidak correction of alpha)
if(n<8)
{
  print('Sample size less than 8! Calculating chi-squared approximated sigma...')
  s<-sd(a)
  lower.conf<-sqrt((n-1)*s^2/qchisq(0.05,n-1,lower.tail = FALSE))
  upper.conf<-sqrt((n-1)*s^2/qchisq(1-0.05,n-1,lower.tail = FALSE))
}
if (n>=8)
{
  print('Initial sample size sufficient. Estimating sigma nonparametrically with bias-corrected accelerated CI for mad...')
  library('boot')
  std<-function(data,i)
  {
    return(mad(data[i]))
  }
  bootresult<-boot(a,std,R=10^5)
  boot.ci(bootresult, conf = 0.9,type='bca')
  lower.conf<-boot.ci(bootresult, conf = 1-alpha,type='bca')$bca[4]
  upper.conf<-boot.ci(bootresult, conf = 1-alpha,type='bca')$bca[5]
}
sigma<-upper.conf
x1<-numeric()
x2<-numeric()
dx<-numeric()
print('Calculating sample size...Please wait')
for(i in 1:500)
{
for (k in 1:10^4)
{
x1<-rnorm(i,mean=0,sd=sigma)
x2<-rnorm(i,mean=0+difference,sd=sigma)
dx[k]<-mean(x2)-mean(x1)
}
if (quantile(dx,alpha/2)>0) break
}
hist(dx)
abline(v=quantile(dx,alpha/2),lty=2,col='blue')
if (quantile(dx,alpha/2)<=0)
{
  print('!!!Warning!!! There are >500 samples required!!! Please loosen the requirements and repeat calculations!')
}
print(c('Sidak-adjusted alpha is:',alpha)) #Sidak-adjusted alpha
print(c('Alpha-quantile is:',quantile(dx,alpha/2)))
print(c('Estimated sigma is:',sigma)) #estimated sigma
print(c('Estimated number of samples is:',i)) #estimated number of samples (two-sided test)

  
  
  
  
  