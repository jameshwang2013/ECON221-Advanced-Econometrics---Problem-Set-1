#Set work directory
setwd("/Users/jameswang/Documents/Dropbox/Notability/ECON221/P.Set_1")
getwd()

#import data
data <- read.table("PS0_2013_data.txt", header=TRUE)
attach(data)
hist(XLAG1)

###fit 8 different models
fit1 <- lm(Y~XLAG1)                             #linear
fit2 <- lm(Y~XLAG1+I(XLAG1^2))                  #second-order polynomial
fit3 <- lm(I(Y^2)~XLAG1)                        #Y^2
fit4 <- lm(log(Y)~XLAG1)                        #log(Y)

data.ts = ts(data$Y)
t = seq(1,length(data.ts))
fit5 = lm(data.ts~t+XLAG1)                      #time-series

log.data.ts = log(data.ts)
fit6 = lm(log.data.ts~t+XLAG1)                  #time-series w/log(Y)

q1.mat = matrix(data = rep(c(rep(1,1),rep(0,3)),20),nrow=80,ncol=1)
q2.mat = matrix(data = rep(c(rep(0,1),rep(1,1),rep(0,2)),20),nrow=80,ncol=1)
q3.mat = matrix(data = rep(c(rep(0,2),rep(1,1),rep(0,1)),20),nrow=80,ncol=1)
q4.mat = matrix(data = rep(c(rep(0,3),rep(1,1)),20),nrow=80,ncol=1)

q1.temp = rbind(q1.mat,c(1))
q2.temp = rbind(q2.mat,c(0))
q3.temp = rbind(q3.mat,c(0))
q4.temp = rbind(q4.mat,c(0))

q1 = rbind(q1.temp,c(0))
q2 = rbind(q2.temp,c(1))
q3 = rbind(q3.temp,c(0))
q4 = rbind(q4.temp,c(0))

n = length(log.data.ts)
fit7 = lm(log.data.ts~t+XLAG1+q1+q2+q3+q4)      #time series + seasonal dummies 

#perform model selection, use Bayesian Information Criterion (BIC/SIC)
BIC(fit1)     #234.7273
BIC(fit2)     #239.0946
BIC(fit3)     #907.4008
BIC(fit4)     #-323.8864 ***best
BIC(fit5)     #239.1339
BIC(fit6)     #-319.4803 ***second-best
BIC(fit7)     #-308.3462 ***third-best

#summary statistics for best models & fit1
summary(fit1)
summary(fit4)
summary(fit6)
summary(fit7)

#plot the best models
plot(XLAG1,Y,ylab="Healthcare Quarterly Earnings",xlab="Lagged-X Variable", main="Healthcare Quarterly Earnings on Lagged-X Variable")
abline(fit1,col="red")

plot(XLAG1,log(Y),ylab="Y",xlab="XLAG1")
abline(fit4,col="red")

plot.ts(t,log.data.ts,xlab="Time",ylab="Healthcare Quarterly Earnings",main="Healthcare Quarterly Earnings Time Series (FIT6)")
lines(fitted(fit6),col="red")

plot.ts(t,log.data.ts,,xlab="Time",ylab="Healthcare Quarterly Earnings",main="Healthcare Quarterly Earnings Time Series (FIT7)")
lines(fitted(fit7),col="red")

#plot residuals
fit.res1 <- resid(fit1)
plot(XLAG1, fit.res1, ylab="Residuals", xlab="XLAG1",main="Residual Plot of FIT1 Residuals")
abline(0,0)
hist(fit.res1,main="Residual Histogram of FIT1")

fit.res4 <- resid(fit4)
plot(XLAG1, fit.res4, ylab="Residuals", xlab="XLAG1")
abline(0,0)
hist(fit.res4,main="Residual Histogram of FIT4")

fit.res6 <- resid(fit6)
plot(XLAG1, fit.res6, ylab="Residuals", xlab="XLAG1")
abline(0,0)
hist(fit.res6)

fit.res7 <- resid(fit7)
plot(XLAG1, fit.res7, ylab="Residuals", xlab="XLAG1")
abline(0,0)
hist(fit.res7)

#plot normal qq plots
qqPlot(fit1,main="Normal QQ Plot for FIT1 Residuals")
qqPlot(fit4,main="Normal QQ Plot for FIT4 Residuals")
qqPlot(fit6,main="Normal QQ Plot for FIT6 Residuals")
qqPlot(fit7,main="Normal QQ Plot for FIT7 Residuals")

#plot leverage plot ... identify influential points
plot(fit1) #see 4th plot
plot(fit4) #see 4th plot
plot(fit6) #see 4th plot
plot(fit7) #see 4th plot

#conduct Bonferroni-adjusted outlier test ... test for outliers ... 
outlierTest(fit1) #p-value = .28974 ... no outliers
outlierTest(fit4) #p-value = .30527 ... no outliers
outlierTest(fit6) #p-value = .29547 ... no outliers
outlierTest(fit7) #p-value = .40510 ... no outliers

#conduct RESET test ... test for linearity
reset(fit1,power=2:3,type="regressor",data=data)  #.1222 not significant, functional form is linear
reset(fit4,power=2:3,type="regressor",data=data)  #.1284 not significant, functional form is linear
reset(fit6,power=2:3,type="regressor",data=data)  #.2568 not significant, functional form is linear
reset(fit7,power=2:3,type="regressor",data=data)  #.9396 not significant, functional form is linear

#conduct Durbin-Watson Test ... test for serial correlation
#install.packages("car")
#library(car)
dwt(fit1,simulate=TRUE)     #2.2982
dwt(fit4,simulate=TRUE)     #2.2912 ***best
dwt(fit6,simulate=TRUE)     #2.2913 ***second-best
dwt(fit7,simulate=TRUE)     #2.2978 ***third-best

#conduct Breusch-Pagan Test ... test for heteroskedasticity
#install.packages("lmtest")
#library(lmtest)
bptest(fit1)     #p-value = .2255 ... fail to reject null hypothesis of homoskedasticity
bptest(fit4)     #p-value = .1436 ... fail to reject null hypothesis of homoskedasticity
bptest(fit6)     #p-value = .3165 ... fail to reject null hypothesis of homoskedasticity
bptest(fit7)     #p-value = .4049 ... fail to reject null hypothesis of homoskedasticity

#conduct test ... test for structural change
#install.packages("strucchange")
#library(struccchange)
struc.change <- function(y,variables,data.set) {
  vec1 <- vector("numeric")
  vec2 <- logical()
   for (i in 1:82) {
      result <- sctest(y~variables,data=data.set,type="Chow",i)
      vec1 <- c(vec1,result$ p.value)
   }
  for (j in vec1)
  {
    if (j <= 0.05) {
      j = FALSE
      vec2 <- c(vec2,j)
    } else {
      j = TRUE
     vec2 <- c(vec2,j)
    }
  }
  return(vec2)
}
  
test1 <- struc.change(Y,c(XLAG1),data) #no rejections of the null...no structural change
test2 <- struc.change(log(Y),c(XLAG1),data) #no rejections of the null...no structural change
test3 <- struc.change(log(Y),c(XLAG1+t),data.ts) #no rejections of the null...no structural change
test4 <- struc.change(log(Y),c(XLAG1+t+q1+q2+q3+q4),data.ts) #no rejections of the null...no structural change

