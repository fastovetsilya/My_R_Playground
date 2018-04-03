Nboot<-1000
x<-rnorm(8,0,1)
dx<-1.5
std<-function(x,i)
{
  return(mad(x[i]))
}
bootresult<-boot(x,std,R=Nboot)$t
sigma<-max(bootresult)
#sigma<-boot.ci(bootresult, conf = 0.90,type='bca')$bca[5]
n<-2*(qnorm(0.975)+qnorm(0.8))^2/(dx/sigma)^2
n
power.t.test(power = .8, delta = 1.5,sd=sd(x))$n;
