#g1<-rnorm(25,50,5)
#g2<-rnorm(25,60,5)
#t.test(g1,g2,var.equal=TRUE)
#wilcox.test(g1,g2)
library('readxl')
setwd("C:/Users/Илья Фастовец/Desktop")
Data<-read_excel('Data.xlsx')
g1<-Data$g1
g2<-Data$g2

a<-g1#reference sample
no_of_comparisons<-1 #how many comparisons?
alpha<-0.05 #at what alpha to compare??
mean1<-mean(g1) #sample mean 1(smaller)
mean2<-mean(g2) #sample mean 2(larger)
difference<-mean2-mean1
n0<-25 #sample size in groups to compare
n<-length(a)#length of reference sample
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
for (i in 1:10^5)
{
x1<-rnorm(n0,mean=0,sd=sigma)
x2<-rnorm(n0,mean=0+difference,sd=sigma)
dx[i]<-mean(x2)-mean(x1)
}
hist(dx)
abline(v=0,col='red')
p.value<-mean(dx<=0)*2
print(c('p.value is:',p.value)) #two-sided p
print(c('Sidak-adjusted reference alpha is:',alpha)) #Sidak-adjusted reference alpha
print(c('Approximate Bonferroni recalculation of alpha for inference',p.value*no_of_comparisons)) # approximate(re-calculated after Bonferroni) p-value
if (p.value>alpha)
{
  print("NOT SIGNIFICANT")
}
if (p.value<=alpha)
{
  print("SIGNIFICANT")
}

g1<-as.data.frame(g1)
g2<-as.data.frame(g2)







