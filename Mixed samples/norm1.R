mu1<-0
sigma1<-1
n1<-5
mu2<-0
sigma2<-1
n2<-5
sigma<-(sigma1+sigma2)/2
Nboot<-1000
Nreboot<-100
sdboot<-numeric(Nboot)
qsd<-numeric(Nreboot)
qsd.max<-numeric(Nreboot)
c.moresigma<-numeric(Nreboot)
ptm <- proc.time()
for(k in 1:Nreboot)
{
 x1<-rnorm(n1,mean=mu1,sd=sigma1)
 x2<-rnorm(n2,mean=mu2,sd=sigma2)
 for(i in 1:Nboot)
 {
  x1.boot<-sample(x1,replace=TRUE)
  x2.boot<-sample(x2,replace=TRUE)
  sdboot[i]<-(sd(x1.boot)+sd(x2.boot))/2
 }
 if (k<=10)
   hist(sdboot)
qsd[k]<-as.numeric(quantile(sdboot,0.99))
qsd.max[k]<-max(sdboot)
c.moresigma[k]<-length(sdboot[sdboot>sigma])/Nboot
}
proc.time() - ptm #time
hist(qsd)
reliability<-length(qsd[qsd<sigma])/Nreboot
reliability
reliability.max<-length(qsd.max[qsd.max<sigma])/Nreboot
reliability.max
1-quantile(c.moresigma,0.05)