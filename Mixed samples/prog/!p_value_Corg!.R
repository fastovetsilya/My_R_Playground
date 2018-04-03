library('readxl')
setwd("C:/Users/Илья Фастовец/Desktop")

###INPUT FROM HERE
a<-c(3.1, 4.6, 4.3, 5.8, 4.6, 3.4, 4.7, 5.1)#reference sample
no_of_comparisons<-1 #how many comparisons?
alpha<-0.05 #at what alpha to compare??
mean1<-3 #sample mean 1(smaller)
mean2<-5 #sample mean 2(larger)
difference<-mean2-mean1
n0<-25 #sample size in groups to compare
###END OF INPUT

print('Fastovets Minimum Credence test of significance for mixed samples')
print('Adapted for carbon analysis in fractions')
sa1<-numeric()
sa2<-numeric()
if(mean1<=0.46)
{
  sa1<-0.1
}
if(mean1>0.46 && mean1<=0.7)
{
  sa1<-0.12
}
if(mean1>0.7 && mean1<=1.7)
{
  sa1<-0.19
}
if(mean1>1.7 && mean1<=5.2)
{
  sa1<-0.21
}
if(mean1>5.2)
{
  sa1<-0.41
}
if(mean2<=0.46)
{
  sa2<-0.1
}
if(mean2>0.46 && mean2<=0.7)
{
  sa2<-0.12
}
if(mean2>0.7 && mean2<=1.7)
{
  sa2<-0.19
}
if(mean2>1.7 && mean2<=5.2)
{
  sa2<-0.21
}
if(mean2>5.2)
{
  sa2<-0.41
}
n<-length(a)#length of reference sample
alpha<-1-(1-alpha)^(1/no_of_comparisons)#(Sidak correction of alpha)
if(n<8)
{
  print('Reference sample size less than 8! Calculating chi-squared approximated sigma...')
  s<-sd(a)
  lower.conf<-sqrt((n-1)*s^2/qchisq(0.05,n-1,lower.tail = FALSE))
  upper.conf<-sqrt((n-1)*s^2/qchisq(1-0.05,n-1,lower.tail = FALSE))
}
if (n>=8)
{
  print('Reference sample size sufficient. Estimating sigma nonparametrically with bias-corrected accelerated CI for mad...')
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
for (i in 1:10^5)
{
x1<-rnorm(n0,mean=rnorm(1,0,sa1),sd=sigma)
x2<-rnorm(n0,mean=rnorm(1,0,sa2)+difference,sd=sigma)
dx[i]<-mean(x2)-mean(x1)
}
hist(dx)
abline(v=0,col='red')
p.value<-mean(dx<=0)*2
print(c('Analytical sd around mean 1 is:',sa1))
print(c('Analytical sd around mean 2 is:',sa2))
print(c('Sigma is:',sigma))
print(c('p.value is:',p.value)) #two-sided p
print(c('Sidak-adjusted reference alpha is:',alpha)) #Sidak-adjusted reference alpha
print(c('Approximate Bonferroni recalculation of alpha for inference',p.value*no_of_comparisons)) # approximate(re-calculated after Bonferroni) p-value
if (p.value>alpha)
{
  print("NOT SIGNIFICANT")
}
if (p.value<=alpha && p.value*no_of_comparisons>=0.001)
{
  print("SIGNIFICANT")
}
if (p.value*no_of_comparisons<0.001)
{
  print('HIGHLY SIGNIFICANT')
}








