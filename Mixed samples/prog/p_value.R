
###INPUT FROM HERE

a<- # Reference sample
no_of_comparisons<-1 # How many comparisons?
alpha<-0.05 # Overall confidence level
mean1<-3 # Sample mean 1(smaller)
mean2<-4.8 # Sample mean 2(larger)
n1<-7 # Sample size of group 1 
n2<-6 # Sample size of group 2 
method<-0 # Use '0' if the study is conducted in one lab, '1' and metrological data if in different labs
Nboot<-10^5 # Number of resamples in bootstrap 10^5 is recommended
Nsim<-10^5 #Number of resamples in simulation 10^5 is recommended
# Use the input below only if method=1
Prec.Data<-read_excel('DataGOST.xlsx')# Data on metrological characteristics 
element<-'K2O' # Chemical element

###END OF INPUT

print('Fastovets Minimum Credence test of significance for mixed samples')
# Extracting metrological data
if(method==1)
{
print('Adapted for RFA')
Prec.Data<-as.matrix(Prec.Data)
C<-mean1
sa1<-eval(parse(text=Prec.Data[3,element]))/100*C
C<-mean2
sa2<-eval(parse(text=Prec.Data[3,element]))/100*C
if (mean1<as.numeric(Prec.Data[1,element]) | mean1>as.numeric(Prec.Data[2,element]))
{
  print('!!!Mean 1 is out of precision range!!!')
}
if (mean2<as.numeric(Prec.Data[1,element]) | mean2>as.numeric(Prec.Data[2,element]))
{
  print('!!!Mean 2 is out of precision range!!!')
}
}
# Estimating sigma
difference<-mean2-mean1
n<-length(a)
alpha<-1-(1-alpha)^(1/no_of_comparisons)
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
  bootresult<-boot(a,std,R=Nboot)
  boot.ci(bootresult, conf = 0.9,type='bca')
  lower.conf<-boot.ci(bootresult, conf = 1-alpha,type='bca')$bca[4]
  upper.conf<-boot.ci(bootresult, conf = 1-alpha,type='bca')$bca[5]
}
sigma<-upper.conf
# Simulation
x1<-numeric()
x2<-numeric()
dx<-numeric()
if (method==1)
{
print('Corrected for metrological data method is used')
 for (i in 1:Nsim)
 {
 x1<-rnorm(n1,mean=rnorm(1,0,sa1),sd=sigma)
 x2<-rnorm(n2,mean=rnorm(1,0,sa2)+difference,sd=sigma)
 dx[i]<-mean(x2)-mean(x1)
 }
}
if (method==0)
{
print('Simple method is used')
 for (i in 1:Nsim)
 {
 x1<-rnorm(n1,mean=0,sd=sigma)
 x2<-rnorm(n2,mean=0+difference,sd=sigma)
 dx[i]<-mean(x2)-mean(x1)
 }
}
# Graph and inference
hist(dx)
abline(v=0,col='red')
p.value<-mean(dx<=0)*2
if(method==1)
{
print(paste0('Precision function is: ', Prec.Data[3,element]))
print(paste0('Analytical sd around mean 1 is: ',sa1))
print(paste0('Analytical sd around mean 2 is: ',sa2))
}
print(paste0('Sigma is: ',sigma))
print(paste0('p.value is: ',p.value)) 
print(paste0('Sidak-adjusted reference alpha is: ',alpha)) 
if (p.value*no_of_comparisons<1)
print(paste0('Approximate Bonferroni recalculation of alpha for inference is: ',p.value*no_of_comparisons)) 
if (p.value*no_of_comparisons>1)
  print(paste0('Approximate Bonferroni recalculation of alpha for inference is: ',1)) 
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








