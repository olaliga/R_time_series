library('TSA')
library('uroot')
library('urca')
library('MuMIn')
library('tseries')
library('PerformanceAnalytics')
library('xts')
library('rtsplot')
library('forecast')
library('lmtest')
library('lubridate')
library('compositions')


##### data section ######
data<-read.csv('躉售物價指數.csv',fileEncoding = 'big5')
data$value<-ts(data$value,start = c(1981,1),frequency = 12)
data$Month<-season(data$value)
train_data<-data[1:443,]

######data summary#########
plot(as.xts(data$value),col=1,main = 'Wholesale Index')
boxplot(data$value,main = 'Wholesale Index Box Plot')
hist(data$value,breaks = (40:63)*2,xlim = c(80,130),xlab = 'Value', ylab = 'Count', main = 'Wholesale Index Histogram (Breaks = 2)')

### 偵測 ####
######ACF#########
acf(data$value,lag.max = 100,main='Sample ACF') #要把ts去掉
######PACF#########
pacf(data$value,main='Sample PACF',lag.max = 100) #要把ts去掉
######EACF#########
eacf(data$value)
### BIC #######
res<-armasubsets(y=data$value,nar = 14, nma = 14,ar.method = 'ols')
plot(res)  

####### linear trend #####
modeldt<-data.frame(1:443,train_data$value)
colnames(modeldt)<-c('time','value')
model<-lm(modeldt$value~modeldt$time,data = modeldt)
summary(model)
acf(model$residuals,main='ACF of Trend line',lag.max = 100)
qqnorm((model$residuals-mean(model$residuals))/sd(model$residuals))
abline(a=0,b=1)
##########season#######
model_season = lm(train_data$value~train_data$Month)
summary(model_season)

####### model AR #########
modelar_1 = arima(train_data$value,order=c(2,0,0),method = 'ML')
coeftest(modelar_1)
summary(modelar_1)
plot(modelar_1)

######### model ARI(2,1) #######
modelar_2 = arima(train_data$value,order=c(2,1,0), method = 'ML')
modelar_2$aic
summary(modelar_2)
coeftest(modelar_2)
plot(modelar_2)
par(mfrow=c(1,1))
##### residual ######
plot((modelar_2$residuals-mean(modelar_2$residuals))/sd(modelar_2$residuals),main='Standardize Residuals',xlab='Time',ylab='')
qqnorm((modelar_2$residuals-mean(modelar_2$residuals))/sd(modelar_2$residuals))
abline(a=0,b=1)
acf(modelar_2$residuals,lag.max = 100,main=' ARI(2,1) Residual ACF')
pacf(modelar_2$residuals,lag.max = 100,main=' ARI(2,1) Residual PACF')
shapiro.test(modelar_2$residuals)


ljung_box_plot<-rep(0,100)
for(i in 1:100){
  ljung_box_plot[i]<-Box.test(modelar_2$residuals, lag = i, type = 'Ljung-Box')$p.value
}
plot(ljung_box_plot,main='p values for Ljung-Box statistic',ylab='p value',xlab='lag')
abline(h=0.05,col='blue')


ljung_plot<-function(x,t,name){
  ljung_box_plot<-rep(0,t)
  for(i in 1:t){
    ljung_box_plot[i]<-Box.test(x, lag = i, type = 'Ljung-Box')$p.value
  }
  plot(ljung_box_plot,main=paste('p values for Ljung-Box statistic-',name),ylab='p value',xlab='lag',ylim=c(0,1))
  abline(h=0.05,col='blue')
}
remove(ljung_box_plot)

######## model ARIma(2,1,q) #######
modelarma_211 = arima(train_data$value,order=c(2,1,1), method = 'ML')
modelarma_212 = arima(train_data$value,order=c(2,1,2), method = 'ML')
modelarma_213 = arima(train_data$value,order=c(2,1,3), method = 'ML')
coeftest(modelarma_211)
coeftest(modelarma_212)
coeftest(modelarma_213)
modelarma_211$aic
modelarma_212$aic
modelarma_213$aic
modelarmafix_212 = arima(train_data$value,order=c(2,1,2),fixed = c(NA,NA,0,NA), method = 'ML')
modelarmafix_213 = arima(train_data$value,order=c(2,1,3),fixed = c(0,NA,NA,NA,NA), method = 'ML')
modelarmafix_212$aic
modelarmafix_213$aic


coeftest(modelarma_211)
coeftest(modelarmafix_212)
coeftest(modelarmafix_213)

par(mfrow=c(1,1))
plot(modelarma_211)
plot(modelarmafix_212)
plot(modelarmafix_213)

par(mfrow=c(3,1))
acf(modelarma_211$residuals,lag.max = 100,main='ARIMA(2,1,1) ACF')
acf(modelarmafix_212$residuals,lag.max = 100,main='ARIMA(2,1,2) ACF')
acf(modelarmafix_213$residuals,lag.max = 100,main='ARIMA(2,1,3) ACF')

par(mfrow=c(3,1))
pacf(modelarma_211$residuals,lag.max = 100,main='ARIMA(2,1,1) PACF')
pacf(modelarmafix_212$residuals,lag.max = 100,main='ARIMA(2,1,2) PACF')
pacf(modelarmafix_213$residuals,lag.max = 100,main='ARIMA(2,1,3) PACF')

par(mfrow=c(3,1))
ljung_plot(modelarma_211$residuals,100,'ARIMA(2,1,1)')
ljung_plot(modelarmafix_212$residuals,100,'ARIMA(2,1,2)')
ljung_plot(modelarmafix_213$residuals,100,'ARIMA(2,1,3)')

par(mfrow=c(3,1))
qqnorm(modelarma_211$residuals)
abline(a=0,b=1)
qqnorm(modelarmafix_212$residuals)
abline(a=0,b=1)
qqnorm(modelarmafix_213$residuals)
abline(a=0,b=1)

stdize<-function(x){
  return((x-mean(x))/sd(x))
}

par(mfrow=c(3,1))
plot(stdize(modelarma_211$residuals),main='Standardize Residuals-ARIMA(2,1,1)',xlab='Time',ylab='')
plot(stdize(modelarmafix_212$residuals),main='Standardize Residuals-ARIMA(2,1,2)',xlab='Time',ylab='')
plot(ts(stdize(modelarmafix_213$residuals),start=c(1981,1),frequency = 12),main='Standardize Residuals-ARIMA(2,1,3)',xlab='Time',ylab='')
rs<-ts(stdize(modelarmafix_213$residuals),start=c(1981,1),frequency = 12)
absrs<-data.frame(rs,abs(rs),time(rs))
rs[which(abs(rs)>3)]


modelarmafix_213$loglik
summary(modelarmafix_213)

###### outlier #####

interest<-read.csv('interest.csv',fileEncoding = 'big5')[1:461,]
data4<-read.csv("01~18_出口.csv",fileEncoding = 'big5')
data5<-read.csv("Asia_US.csv",fileEncoding = 'big5')

exit_ratio<-xts(data4[,(2*(2:8)+1)], order.by=seq(as.Date("2001-12-31"), length = 18, by = "years"))
colnames(exit_ratio)<-c('Asia','Europe','Africa','North America','Central America','Latin America','Oceania')
plot(xts(exit_ratio),plot.type="s",auto.legend = TRUE)
addLegend(legend.loc="topright", legend.names=colnames(exit_ratio),lty = 1)

Asia_us_ratio<-xts(data5[,(2*(1:14)+1)], order.by=seq(as.Date("2001-12-31"), length = 18, by = "years"))
plot(xts(Asia_us_ratio),plot.type="s")
addLegend(legend.loc="topright", legend.names=colnames(Asia_us_ratio),lty = 1,family ='NotoSans-Medium')

plot(exit_ratio[,2])
plot(x=2001:2018,y=data4[,2],type='l',col=1)
plot(x=2001:2018,y=data4[,4],type='l',col=1)
plot(x=2001:2018,y=data4[,2],type='l',col=1)
plot(x=2001:2018,y=data4[,2],type='l',col=1)
plot(x=2001:2018,y=data4[,2],type='l',col=1)

plot(modelar_3)
acf(modelar_3$residuals,lag.max = 100)
pacf(modelar_3$residuals,lag.max = 100)
qqnorm(modelar_3$residuals)
abline(a=0,b=1)
shapiro.test(modelar_3$residuals)
plot(modelar_3$residuals)
tsdiag(modelar_3)


plot(modelar_2,n.ahead=18,type='b')
data$value<-ts(data$value)
lines(data$value,col='red')

tS3 <- 100 * cumulated(LPP2005REC[, 1:3])
str(xreg)
xreg$TWSE<-as.numeric(xreg$TWSE)
########## newmodel
####full##
xreg<-read.csv(file = "xreg_variable.csv")
modelarmafixxreg_213 = arima(train_data$value,order=c(2,1,3),xreg = xreg[1:443,2:10], method = 'ML')
summary(modelarmafixxreg_213)
modelarmafixxreg_213$residuals
coeftest(modelarmafixxreg_213)
rsxregfull<-data.frame(data$Time[1:443],modelarmafixxreg_213$residuals,abs(modelarmafixxreg_213$residuals))

plot(modelarmafixxreg_213,xreg=xreg[444:461,2:10],n.ahead=18)
plot(modelarmafixxreg_213$residuals)
acf(modelarmafixxreg_213$residuals,lag.max = 400)
acf(modelarmafix_213$residuals,lag.max = 400)
pacf(modelarmafixxreg_213$residuals,lag.max = 400)
pacf(modelarmafix_213$residuals,lag.max = 400)
modelarmafixxreg_213 = arima(train_data$value,order=c(2,1,3),fixed = c(0,NA,NA,NA,NA,NA,NA,NA,0,NA,NA,NA,NA,NA),xreg = xreg[1:443,2:10], method = 'ML')

### 窮舉全部模型 ###
all<-function(x){
    ls<-c(intToBits(x)[1],intToBits(x)[2],intToBits(x)[3],intToBits(x)[4],intToBits(x)[5],intToBits(x)[6],intToBits(x)[7],intToBits(x)[8],intToBits(x)[9])
    fix<-c(0,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
    fix[which(ls==0)+5]<-0
    arima(train_data$value,order=c(2,1,3),fixed = fix,xreg = xreg[1:443,2:10], method = 'ML')
}

all_selection<-lapply(c(0:511),all)
aic<-lapply(1:512,function(x){return(all_selection[[x]]$aic)})
aic<-unlist(aic) 
which(aic==min(aic))
head(aic[order(aic)])
aic<-data.frame(aic)


all_selection[[306]]$coef
acf(all_selection[[306]]$residuals,lag.max = 400)
pacf(all_selection[[306]]$residuals,lag.max = 400)
ljung_plot(all_selection[[306]]$residuals,100,'best')
qqnorm(all_selection[[306]]$residuals)
abline(a=0,b=1)
qqnorm(modelarma_211$residuals)

par(mfrow=c(3,1))
acf(all_selection[[306]]$residuals,lag.max = 400,main='Aic 1st Model ACF')
acf(all_selection[[305]]$residuals,lag.max = 400,main='Aic 2nd Model ACF')
acf(all_selection[[308]]$residuals,lag.max = 400,main='Aic 3rd Model ACF')

pacf(all_selection[[306]]$residuals,lag.max = 400,main='Aic 1st Model PACF')
pacf(all_selection[[305]]$residuals,lag.max = 400,main='Aic 2nd Model PACF')
pacf(all_selection[[308]]$residuals,lag.max = 400,main='Aic 1st Model PACF')

qqnorm(modelarmafix_213$residuals)
abline(a=0,b=1)
qqnorm(all_selection[[306]]$residuals,main='Aic 1st Model QQ-plot')
abline(a=0,b=1)
qqnorm(all_selection[[305]]$residuals,main='Aic 2nd Model QQ-plot')
abline(a=0,b=1)
qqnorm(all_selection[[308]]$residuals,main='Aic 3rd Model QQ-plot')
abline(a=0,b=1)


plot(stdize(all_selection[[306]]$residuals),main='Standardize Residuals-Aic 1st Model',xlab='Time',ylab='')
plot(stdize(all_selection[[305]]$residuals),main='Standardize Residuals-Aic 2nd Model',xlab='Time',ylab='')
plot(stdize(all_selection[[308]]$residuals),main='Standardize Residuals-Aic 3rd Model',xlab='Time',ylab='')


ljung_plot(all_selection[[306]]$residuals,100,'Aic 1st Model')
ljung_plot(all_selection[[305]]$residuals,100,'Aic 2nd Model')
ljung_plot(all_selection[[308]]$residuals,100,'Aic 3rd Model')


residual<-data.frame(data$Time[1:443],all_selection[[308]]$residuals,abs(all_selection[[308]]$residuals))

##predict
par(mfrow=c(2,1))
plot(modelarmafix_213,n.ahead=18,main='Forecast without Xreg')
lines(as.numeric(data$value),col=2)
plot(all_selection[[306]],newxreg = xreg[444:461,2:10],n.ahead=18,main='Forecast with Xreg')
lines(as.numeric(data$value),col=2)
