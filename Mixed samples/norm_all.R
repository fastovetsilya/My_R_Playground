library ('boot')
#shape1<-20
#shape2<-1.5
#distr<-rbeta(10^7,shape1=shape1,shape2=shape2)
#mu<-shape1/(shape1+shape2)
#sigma<-sqrt(shape1*shape2/((shape1+shape2)^2*(shape1+shape2+1)))
#hist(distr)
mu<-0
sigma<-1
Nboot<-10000
Nreboot<-10000
qsd.perc<-numeric(Nreboot)
qsd.bca<-numeric(Nreboot)
qsd.max<-numeric(Nreboot)
c.moresigma<-numeric(Nreboot)
std<-function(data,i)
{
  return(mad(data[i]))
}
a<-8
b<-10
c<-15
d<-17
e<-25
f<-50
g<-150
h<-0
reliability.perc<-numeric(a+b+c+d+e+f+h)
reliability.bca<-numeric(a+b+c+d+e+f+h)
reliability.max<-numeric(a+b+c+d+e+f+h)
sd.qsd.max<-numeric(a+b+c+d+e+f+h)
msg<-numeric(a+b+c+d+e+f+h)
ptm <- proc.time()# timer start
for (c in c(a,b,c,d,e,f,g))
{
for (k in 1:Nreboot)
  {
  x<-rnorm(c,mu,sigma)
  #x<-rbeta(8,shape1=shape1,shape2=shape2)
  bootresult<-boot(x,std,R=Nboot)
  U.bca<-boot.ci(bootresult, conf = 0.90,type='bca')$bca[5]
  U.perc<-boot.ci(bootresult, conf = 0.90,type='perc')$perc[5]
  qsd.bca[k]<-U.bca
  qsd.perc[k]<-U.perc
  qsd.max[k]<-max(bootresult$t)
  c.moresigma[k]<-length(bootresult$t[bootresult$t>sigma])/Nboot
  if (k<=10)
  hist(bootresult$t)
}
reliability.perc[c]<-length(qsd.perc[qsd.perc<sigma])/Nreboot
reliability.bca[c]<-length(qsd.bca[qsd.bca<sigma])/Nreboot
reliability.max[c]<-length(qsd.max[qsd.max<sigma])/Nreboot
sd.qsd.max[c]<-sd(qsd.max)
msg[c]<-1-quantile(c.moresigma,0.05)
}
proc.time() - ptm #timer stop
msg[msg!=0]
reliability.perc[reliability.perc!=0]
reliability.bca[reliability.bca!=0]
reliability.max[reliability.max!=0]
sd.qsd.max[sd.qsd.max!=0]
hist(qsd.bca)
hist(qsd.max)