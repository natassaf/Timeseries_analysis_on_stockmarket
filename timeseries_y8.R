library("xlsx")
library(urca)
library(fGarch) 
library("rugarch")

temp <- data.frame(matrix(1,ncol = 1, nrow =177 ))
colnames(temp) <- "x16"

df <- read.xlsx("data_Assignment_2017_2018.xlsx",sheetIndex = 1,header = TRUE,stringsAsFactors=FALSE)
df<-sapply(df,as.numeric)
df<-as.data.frame(df)
ydescription<-df[1,]
df<-df[-1,]
df<-df[,-1]

#Creating a new dataframe with yt+1 and xt on the same row
ytrain<-df[c(2:178),c(1:9)]
xtrain<-df[c(1:177),c(10:24)]

mydf<-cbind(ytrain,xtrain,temp) #transformed dataframe excluding values of 2005 but including yof january 2005 as we have x for december 2004

xtest<-df[c(178:188),c(10:24)]
ytest<-df[c(178:189),8]
#Taking y1 -> convert in timeseries
y <- mydf$y8
j=ts(y, frequency=12, start = c(1990,4),end=c(2004,12))
plot(j,type="l", col='red', main="Time Series plot of HFRI ", ylab="Monthly returns")


#Testing Stationarity
m=ar(j)
res<-ur.df(j,type="none",lag=m$order-1)
summary(res)
print("Series is stationary as pvalue<0.05 and also series hdo not have a unit root test and they have drift")


#Autocorrelation problem 
Box.test(j, lag = 12, type = c("Ljung-Box"))#pvalue significantly low so there is autocorrelation problem

acf(j, 177, main="ACF of J&J")
acf(j, 48, main="ACF of J&J")#MA(6)

pacf(j,177,main="PACF")#Partial Autocorrelation problem!
pacf(j,48,main="PACF")#Indicator for AR(6)

#MODELING AUTOCORRELATION ISSUE
ar13<-arima(j,order=c(13,0,0),fixed=c(0,0,0,0,0,NA,0,0,0,0,0,0,NA,NA),transform.pars = FALSE)# lag 6 and lag 13 is significanit
ar13$aic
acf(ar13$residuals,50)
pacf(ar13$residuals,50)

ma6<-arima(j,order=c(0,0,6),fixed=c(0,0,0,0,0,NA,NA))#lag6 is significant with aic=-1168.78
ma6$aic
acf(ma6$residuals,50)
pacf(ma6$residuals,50)

arma66<-arima(j,order=c(6,0,6),fixed=c(0,0,0,0,0,NA,0,0,0,0,0,NA,NA),transform.pars = FALSE)#et-6 is insignificant so this is not a good model
acf(arma66$residuals,50)
pacf(arma66$residuals,50)
arma66$aic#-1170

#I choose ar13 as aic is much smaller and the ACF,PACF plots are better fitted
AIC(ar13)#-1174.478
BIC(ar13)#-1161.774
ar13fit<-ar13


#HETEROSCEDASTICITY ISSUE
acf(j^2,50)
pacf(j^2,50)
Box.test(j^2,lag=12,type="Ljung")#Ho rejected so  there is heteroscedasticity problem


#*****************************************************************************************************

#Create regression model.
y8res <- lm(y8~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12+x13+x14+x15,mydf)
y8res
AIC(y8res)
regression1<-y8res

#Check Autocorrelation of the residuals and heteroscedasticity
lmres<-residuals(y8res)
acf(lmres,50)#Small issue at lag 6
pacf(lmres,50)#Small issue at lag 6,13
acf(lmres^2,50) #Heteroscedasticity ok
pacf(lmres^2,50)#Heteroscedasticity problem on lag 2 
Box.test(lmres^2,lag=12,type="Ljung")#Ho not rejected, no heterosceasticity issue


#AUTOCORRELATION MODELLING
#AR(13)
ar13lm<-arima(lmres,order=c(13,0,0),fixed=c(0,0,0,0,0,NA,0,0,0,0,0,0,NA,NA),transform.pars = FALSE)
ar13lm$aic#-1191.583
acf(ar13lm$residuals,50)#No issue
Box.test(ar13lm$residuals,lag=12,type="Ljung")#H0 not rejected:No autocorrelation problem
pacf(ar13lm$residuals,50)#No issue
Box.test(ar13lm$residuals^2,lag=12,type="Ljung")#pvalue>0.05 no heteroscedasticity problem

#ARMA(6,6)
arma66lm<-arima(j,order=c(6,0,6),fixed=c(0,0,0,0,0,NA,NA,0,0,0,0,NA,NA),transform.pars = FALSE)#for AR significant lag=6 and for MA sign. lag=1, lag=6
arma66lm$aic#-886.2289
acf(arma66lm$residuals,50)#No issue
Box.test(arma66lm$residuals,lag=12,type="Ljung")#H0 not rejected:No autocorrelation problem
pacf(arma66lm$residuals,50)#No issue
Box.test(arma66lm$residuals^2,lag=12,type="Ljung")#pvalue>0.05 no heteroscedasticity problem

#As AR(13) has better aic i choose AR(13)

#HETEROSCEDASTICITY MODELING AFTER ARIMA
ar_residuals<-ar13lm$residuals
acf(ar_residuals^2,60)#problem at lag 9
pacf(ar_residuals^2,60)#problem at lag 11
Box.test(ar_residuals^2,lag=50,type="Ljung")#There is no heteroscedasticity proble as pvalue<0.05

#NORMALITY
shapiro.test(ar_residuals)#H0 is not rejected so it is normal


#FIX HETEROSCEDASTICITY INCLUDING ARMA MODEL
ext <- data.matrix(mydf[,c(11:25)]) 
m.spec <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(11,9)), mean.model = list(armaOrder = c(13, 0), include.mean = FALSE, external.regressors = ext), distribution.model = "norm", start.pars = NULL,fixed.pars = list() )
mfit <- ugarchfit(data = j, spec = m.spec, solver = "hybrid", fit.control = list(scale = 1)) 
mfit#I check here that the only significant parameter of garch is alpha10,alpha11,beta5 so i will now create GARCH(11,5) but keep only apha11,beta5 and from ARMA(13,0) will keep onle lags 6,13

plot(mfit)


ext <- data.matrix(mydf[,c(11:25)]) 
m.spec <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(11,5)), mean.model = list(armaOrder = c(13, 0), include.mean = FALSE, external.regressors = ext), distribution.model = "norm", start.pars = NULL,
                     fixed.pars = list(ar1=0,ar2=0,ar3=0,ar4=0,ar5=0,ar7=0,ar8=0,ar9=0,ar10=0,ar11=0,ar12=0,alpha1=0,alpha2=0,alpha3=0,alpha4=0,alpha5=0,alpha6=0,alpha7=0,alpha8=0,alpha9=0,beta1=0,beta2=0,beta3=0,beta4=0) )
mfit1 <- ugarchfit(data = j, spec = m.spec, solver = "hybrid", fit.control = list(scale = 1)) 

#RETRIVE STANDARDIZED RESIDUALS 
residuals1<-ts(residuals(mfit1,standardize=TRUE),freq=1)

#AUTOCORRELATION:
acf(residuals1, 100)#Autocorrelation ok
pacf(residuals1,100)#Autocorrelation ok
Box.test(residuals1,lag=50,type="Ljung")

#HETEROSCEDASTICITY
acf(residuals1^2, 100)#Heteroscedasticity fixed
pacf(residuals1^2,100)
Box.test(residuals1^2,lag=50,type="Ljung")

#NORMALITY
plot(mfit1)#type 9
shapiro.test((residuals1))#Ho not rejected residuals are normal!


#REGRESSION PARAMETERS REMOVE
mfit1#x12 is the least significant
aik1<- -6.65
bic1<- -6.28
#********************REMOVING X13*************************************************
#*********************************************************************************
y8res <- lm(y8~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12+x14+x15,mydf)

lmres<-residuals(y8res,standardized=TRUE)
acf(lmres,60)#issue at lag 6
pacf(lmres,60)#Small issue at lag 6,13
acf(lmres^2,60) #Heteroscedasticity problem 
pacf(lmres^2,60)#Heteroscedasticity problem
Box.test(j^2,lag=60,type="Ljung")#Ho not rejected, no heterosceasticity issue

#AUTOCORRELATION MODELLING
#AR(13)
ar13lm<-arima(lmres,order=c(13,0,0),fixed=c(0,0,0,0,0,NA,0,0,0,0,0,0,NA,NA),transform.pars = FALSE)
ar13lm$aic#-1190
acf(ar13lm$residuals,50)#No issue
Box.test(ar13lm$residuals,lag=12,type="Ljung")#H0 not rejected:No autocorrelation problem
pacf(ar13lm$residuals,50)#No issue
Box.test(ar13lm$residuals^2,lag=12,type="Ljung")#pvalue>0.05 no heteroscedasticity problem

#ARMA(6,6)
arma66lm<-arima(j,order=c(6,0,6),fixed=c(0,0,0,0,0,NA,NA,0,0,0,0,NA,NA),transform.pars = FALSE)#for AR significant lag=6 and for MA sign. lag=1, lag=6
arma66lm$aic#-1169
acf(arma66lm$residuals,50)#No issue
Box.test(arma66lm$residuals,lag=12,type="Ljung")#H0 not rejected:No autocorrelation problem
pacf(arma66lm$residuals,50)#No issue
Box.test(arma66lm$residuals^2,lag=12,type="Ljung")#pvalue>0.05 no heteroscedasticity problem

#As AR(13) has better aic i choose AR(13)

#HETEROSCEDASTICITY MODELING AFTER ARIMA
residuals2_prior<-residuals(ar13lm,standarized)
ar_residuals<-residuals2_prior
acf(ar_residuals^2,60)#problem at lag 10
pacf(ar_residuals^2,60)#problem at lag 10
Box.test(ar_residuals^2,lag=60,type="Ljung")#There is heteroscedasticity proble as pvalue<0.05

#NORMALITY
shapiro.test(ar_residuals)#H0 is not rejected so it is normal

#FORMULA AR-GRACH-LM:
# ARMA(13,0) AND GARCH(11,12) keeping in a2 a10,a11,b2,b10,b12 and use normal distribution as it fits my data better
ext <- data.matrix(mydf[,c(11:25)]) 
m.spec <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(11,12)), mean.model = list(armaOrder = c(13, 0), include.mean = FALSE, external.regressors = ext), distribution.model = "norm", start.pars = NULL,
                     fixed.pars = list(mxreg13=0,
                                       ar1=0,ar2=0,ar3=0,ar4=0,ar5=0,ar7=0,ar8=0,ar9=0,ar10=0,ar11=0,ar12=0,
                                       alpha1=0,alpha3=0,alpha4=0,alpha5=0,alpha6=0,alpha7=0,alpha8=0,alpha9=0,
                                       beta1=0,beta3=0,beta4=0,beta5=0,beta6=0,beta7=0,beta8=0,beta9=0,beta11=0) )
mfit2<- ugarchfit(data = j, spec = m.spec, solver = "hybrid", fit.control = list(scale = 1)) 

#RETRIVE STANDARDIZED RESIDUALS 
residuals2<-ts(residuals(mfit2,standardize=TRUE),freq=1)

#AUTOCORRELATION OK:
acf(residuals2, 100)#Autocorrelation ok
pacf(residuals2,100)#Autocorrelation ok
Box.test(residuals2,lag=100,type="Ljung")

#HETEROSCEDASTICITY OK
acf(residuals2^2, 100)#Heteroscedasticity fixed
pacf(residuals2^2,100)
Box.test(residuals2^2,lag=100,type="Ljung")#Heteroscedasticity fixed

#NORMALITY OK
plot(mfit2)#type 9
shapiro.test((residuals2))#Ho not rejected residuals are normal!

#REGRESSION PARAMETERS REMOVE
mfit2#x6 is the least significant
aik2<- -6.65
bic2<- -6.30
#********************REMOVING X6*************************************************
#*********************************************************************************
y8res <- lm(y8~x1+x2+x3+x4+x5+x7+x8+x9+x10+x11+x12+x14+x15,mydf)

lmres<-residuals(y8res,standardized=TRUE)
acf(lmres,60)#Small issue at lag 6
pacf(lmres,60)#Small issue at lag 6,13
acf(lmres^2,60) #Heteroscedasticity problem 
pacf(lmres^2,60)#Heteroscedasticity problem
Box.test(j^2,lag=60,type="Ljung")#Ho not rejected, no heterosceasticity issue

#AUTOCORRELATION MODELLING
#AR(13)
ar13lm<-arima(lmres,order=c(13,0,0),fixed=c(0,0,0,0,0,NA,0,0,0,0,0,0,NA,NA),transform.pars = FALSE)
ar13lm$aic#-1191
acf(ar13lm$residuals,50)#No issue
Box.test(ar13lm$residuals,lag=12,type="Ljung")#H0 not rejected:No autocorrelation problem
pacf(ar13lm$residuals,50)#No issue
Box.test(ar13lm$residuals^2,lag=12,type="Ljung")#pvalue>0.05 no heteroscedasticity problem

#ARMA(6,6)
arma66lm<-arima(j,order=c(6,0,6),fixed=c(0,0,0,0,0,NA,NA,0,0,0,0,NA,NA),transform.pars = FALSE)#for AR significant lag=6 and for MA sign. lag=1, lag=6
arma66lm$aic#-1169
acf(arma66lm$residuals,50)#No issue
Box.test(arma66lm$residuals,lag=12,type="Ljung")#H0 not rejected:No autocorrelation problem
pacf(arma66lm$residuals,50)#No issue
Box.test(arma66lm$residuals^2,lag=12,type="Ljung")#pvalue>0.05 no heteroscedasticity problem

#As AR(13) has better aic i choose AR(13)

#HETEROSCEDASTICITY MODELING AFTER ARIMA
residuals3_prior<-residuals(ar13lm,standarized)
ar_residuals<-residuals3_prior
acf(ar_residuals^2,60)#problem at lag 11,13
pacf(ar_residuals^2,60)#problem at lag 11
Box.test(ar_residuals^2,lag=60,type="Ljung")#There is heteroscedasticity proble as pvalue<0.05

#NORMALITY
shapiro.test(ar_residuals)#H0 is not rejected so it is normal

#FORMULA AR-GRACH-LM:
# ARMA(13,0) AND GARCH(11,12) keeping in a11,a13,b11 and use normal distribution as it fits my data better
ext <- data.matrix(mydf[,c(11:25)]) 
m.spec <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(13,11)), mean.model = list(armaOrder = c(13, 0), include.mean = FALSE, external.regressors = ext), distribution.model = "norm", start.pars = NULL,
                     fixed.pars = list(mxreg13=0,mxreg6=0,
                                       ar1=0,ar2=0,ar3=0,ar4=0,ar5=0,ar7=0,ar8=0,ar9=0,ar10=0,ar11=0,ar12=0,
                                       alpha1=0,alpha2=0,alpha3=0,alpha4=0,alpha5=0,alpha6=0,alpha7=0,alpha8=0,alpha9=0,alpha10=0,
                                       beta1=0,beta2=0,beta3=0,beta4=0,beta5=0,beta6=0,beta7=0,beta8=0,beta9=0,beta10=0) )
mfit3<- ugarchfit(data = j, spec = m.spec, solver = "hybrid", fit.control = list(scale = 1)) 

#RETRIVE STANDARDIZED RESIDUALS 
residuals3<-ts(residuals(mfit3,standardize=TRUE),freq=1)

#AUTOCORRELATION:OK
acf(residuals3, 100)#Autocorrelation ok
pacf(residuals3,100)#Autocorrelation ok

#HETEROSCEDASTICITY:OK
acf(residuals3^2, 100)#Heteroscedasticity fixed
pacf(residuals3^2,100)
Box.test(residuals3^2,lag=100,type="Ljung")#Heteroscedasticity fixed

#NORMALITY 0K
plot(mfit3)#type 9
shapiro.test((residuals3))#Ho not rejected residuals are normal!

#REGRESSION PARAMETERS REMOVE
mfit3#x9 is the least significant
aik3<- -6.6919
bic3<- -6.3868
#********************REMOVING X9*************************************************
#*********************************************************************************
y8res <- lm(y8~x1+x2+x3+x4+x5+x7+x8+x10+x11+x12+x14+x15,mydf)

lmres<-residuals(y8res,standardized=TRUE)
acf(lmres,60)#Small issue at lag 6
pacf(lmres,60)#Small issue at lag 6,13
acf(lmres^2,60) #Heteroscedasticity problem 
pacf(lmres^2,60)#Heteroscedasticity problem
Box.test(j^2,lag=60,type="Ljung")#Ho not rejected, no heterosceasticity issue

#AUTOCORRELATION MODELLING
#AR(13)
ar13lm<-arima(lmres,order=c(13,0,0),fixed=c(0,0,0,0,0,NA,0,0,0,0,0,0,NA,NA),transform.pars = FALSE)
ar13lm$aic#-1187
acf(ar13lm$residuals,50)#No issue
Box.test(ar13lm$residuals,lag=12,type="Ljung")#H0 not rejected:No autocorrelation problem
pacf(ar13lm$residuals,50)#No issue
Box.test(ar13lm$residuals^2,lag=12,type="Ljung")#pvalue>0.05 no heteroscedasticity problem

#HETEROSCEDASTICITY MODELING AFTER ARIMA
residuals4_prior<-residuals(ar13lm,standarized)
ar_residuals<-residuals4_prior
acf(ar_residuals^2,60)#problem at lag 10
pacf(ar_residuals^2,60)#problem at lag 10
Box.test(ar_residuals^2,lag=60,type="Ljung")#There is heteroscedasticity proble as pvalue<0.05

#NORMALITY
shapiro.test(ar_residuals)#H0 is not rejected so it is normal

#FORMULA AR-GRACH-LM:
# ARMA(13,0) AND GARCH(11,12) keeping in a11,b11 and use normal distribution as it fits my data better
ext <- data.matrix(mydf[,c(11:25)]) 
m.spec <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(11,11)), mean.model = list(armaOrder = c(13, 0), include.mean = FALSE, external.regressors = ext), distribution.model = "norm", start.pars = NULL,
                     fixed.pars = list(mxreg13=0,mxreg6=0,mxreg9=0,
                                       ar1=0,ar2=0,ar3=0,ar4=0,ar5=0,ar7=0,ar8=0,ar9=0,ar10=0,ar11=0,ar12=0,
                                       alpha1=0,alpha2=0,alpha3=0,alpha4=0,alpha5=0,alpha6=0,alpha7=0,alpha8=0,alpha9=0,alpha10=0,
                                       beta1=0,beta2=0,beta3=0,beta4=0,beta5=0,beta6=0,beta7=0,beta8=0,beta9=0,beta10=0) )
mfit4<- ugarchfit(data = j, spec = m.spec, solver = "hybrid", fit.control = list(scale = 1)) 

#RETRIVE STANDARDIZED RESIDUALS 
residuals4<-ts(residuals(mfit4,standardize=TRUE),freq=1)

#AUTOCORRELATION:OK
acf(residuals4, 100)#Autocorrelation ok
pacf(residuals4,100)#Autocorrelation ok
Box.test(residuals4,lag=100,type="Ljung")

#HETEROSCEDASTICITY:OK
acf(residuals4^2, 100)#Heteroscedasticity fixed
pacf(residuals4^2,100)
Box.test(residuals4^2,lag=100,type="Ljung")#Heteroscedasticity fixed

#NORMALITY 0K
plot(mfit4)#type 9
shapiro.test((residuals4))#Ho not rejected residuals are normal!

#REGRESSION PARAMETERS REMOVE
mfit4#x4 is the least significant
aic<- -6.4404
bic<- -6.1533
#********************REMOVING X4*************************************************
#*********************************************************************************
y8res <- lm(y8~x1+x2+x3+x5+x7+x8+x10+x11+x12+x14+x15,mydf)


lmres<-residuals(y8res,standardized=TRUE)
acf(lmres,60)#Small issue at lag 6
pacf(lmres,60)#Small issue at lag 6,13
acf(lmres^2,60) #Heteroscedasticity problem 
pacf(lmres^2,60)#Heteroscedasticity problem
Box.test(j^2,lag=60,type="Ljung")#Ho not rejected, no heterosceasticity issue

#AUTOCORRELATION MODELLING
#AR(13)
ar13lm<-arima(lmres,order=c(13,0,0),fixed=c(0,0,0,0,0,NA,0,0,0,0,0,0,NA,NA),transform.pars = FALSE)
ar13lm$aic#-1187
acf(ar13lm$residuals,50)#No issue
Box.test(ar13lm$residuals,lag=12,type="Ljung")#H0 not rejected:No autocorrelation problem
pacf(ar13lm$residuals,50)#No issue
Box.test(ar13lm$residuals^2,lag=12,type="Ljung")#pvalue>0.05 no heteroscedasticity problem

#HETEROSCEDASTICITY MODELING AFTER ARIMA
residuals5_prior<-residuals(ar13lm,standarized)
ar_residuals<-residuals5_prior
acf(ar_residuals^2,60)#problem at lag 10
pacf(ar_residuals^2,60)#problem at lag 10
Box.test(ar_residuals^2,lag=60,type="Ljung")#There is heteroscedasticity proble as pvalue<0.05

#NORMALITY
shapiro.test(ar_residuals)#H0 is not rejected so it is normal


#FORMULA AR-GRACH-LM:
# ARMA(13,0) AND GARCH(11,12) keeping in a10,b10 and use normal distribution as it fits my data better
ext <- data.matrix(mydf[,c(11:25)]) 
m.spec <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(10,10)), mean.model = list(armaOrder = c(13, 0), include.mean = FALSE, external.regressors = ext), distribution.model = "norm", start.pars = NULL,
                     fixed.pars = list(mxreg13=0,mxreg6=0,mxreg9=0,mxreg4=0,
                                       ar1=0,ar2=0,ar3=0,ar4=0,ar5=0,ar7=0,ar8=0,ar9=0,ar10=0,ar11=0,ar12=0,
                                       alpha1=0,alpha2=0,alpha3=0,alpha4=0,alpha5=0,alpha6=0,alpha7=0,alpha8=0,alpha9=0,
                                       beta1=0,beta2=0,beta3=0,beta4=0,beta5=0,beta6=0,beta7=0,beta8=0,beta9=0) )
mfit5<- ugarchfit(data = j, spec = m.spec, solver = "hybrid", fit.control = list(scale = 1)) 

#RETRIEVE STANDARDIZED RESIDUALS 
residuals5<-ts(residuals(mfit5,standardize=TRUE),freq=1)

#AUTOCORRELATION:OK
acf(residuals5, 50)#Autocorrelation ok
pacf(residuals5,50)#Autocorrelation ok
Box.test(residuals5,lag=50,type="Ljung")

#HETEROSCEDASTICITY:OK
acf(residuals5^2, 50)#Heteroscedasticity NOT FIXED
pacf(residuals5^2,50)
Box.test(residuals5^2,lag=50,type="Ljung")#Heteroscedasticity NOT FIXED

#NORMALITY 0K
plot(mfit5)#type 9
shapiro.test((residuals5))#Ho not rejected residuals are normal!


#REGRESSION PARAMETERS REMOVE
mfit5 #x3 is the least significant
aik<--6.4404
Bic<- -6.1533
#********************REMOVING X3*************************************************
#*********************************************************************************
y8res <- lm(y8~x1+x2+x5+x7+x8+x10+x11+x12+x14+x15,mydf)


lmres<-residuals(y8res,standardized=TRUE)
acf(lmres,60)#Small issue at lag 6
pacf(lmres,60)# issue at lag 6,13
acf(lmres^2,60) #Heteroscedasticity problem 
pacf(lmres^2,60)#Heteroscedasticity problem
Box.test(j^2,lag=60,type="Ljung")#Ho not rejected, no heterosceasticity issue

#AUTOCORRELATION MODELLING
#AR(13)
ar13lm<-arima(lmres,order=c(13,0,0),fixed=c(0,0,0,0,0,NA,0,0,0,0,0,0,NA,NA),transform.pars = FALSE)
ar13lm$aic#-1187
acf(ar13lm$residuals,50)#No issue
Box.test(ar13lm$residuals,lag=12,type="Ljung")#H0 not rejected:No autocorrelation problem
pacf(ar13lm$residuals,50)#No issue
Box.test(ar13lm$residuals^2,lag=12,type="Ljung")#pvalue>0.05 no heteroscedasticity problem

#HETEROSCEDASTICITY MODELING AFTER ARIMA
residuals6_prior<-residuals(ar13lm,standarized)
ar_residuals<-residuals6_prior
acf(ar_residuals^2,60)#problem at lag 10
pacf(ar_residuals^2,60)#problem at lag 11
Box.test(ar_residuals^2,lag=60,type="Ljung")#There is heteroscedasticity proble as pvalue<0.05

#NORMALITY
shapiro.test(ar_residuals)#H0 is not rejected so it is normal

#FORMULA AR-GRACH-LM:
# ARMA(13,0) AND GARCH(11,12) keeping in a2,a11 ,b2,b10 and use normal distribution as it fits my data better
ext <- data.matrix(mydf[,c(11:25)]) 
m.spec <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(11,10)), mean.model = list(armaOrder = c(13, 0), include.mean = FALSE, external.regressors = ext), distribution.model = "norm", start.pars = NULL,
                     fixed.pars = list(mxreg13=0,mxreg6=0,mxreg9=0,mxreg4=0,mxreg3=0,
                                       ar1=0,ar2=0,ar3=0,ar4=0,ar5=0,ar7=0,ar8=0,ar9=0,ar10=0,ar11=0,ar12=0,
                                       alpha1=0,alpha3=0,alpha4=0,alpha5=0,alpha6=0,alpha7=0,alpha8=0,alpha9=0,alpha10=0,
                                       beta1=0,beta3=0,beta4=0,beta5=0,beta6=0,beta7=0,beta8=0,beta9=0) )
mfit6<- ugarchfit(data = j, spec = m.spec, solver = "hybrid", fit.control = list(scale = 1)) 

#RETRIEVE STANDARDIZED RESIDUALS 
residuals6<-ts(residuals(mfit6,standardize=TRUE),freq=1)

#AUTOCORRELATION:OK
acf(residuals6, 50)#Autocorrelation problem at lag 25(insignificant-small)
pacf(residuals6,50)#Autocorrelation ok
Box.test(residuals6,lag=100,type="Ljung")#Autocorellation ok

#HETEROSCEDASTICITY:OK
acf(residuals6^2, 50)#Heteroscedasticity fixed
pacf(residuals6^2,50)
Box.test(residuals6^2,lag=50,type="Ljung")#Heteroscedasticity fixed

#NORMALITY 0K
plot(mfit6)#type 9
shapiro.test((residuals6))#Ho not rejected residuals are normal!

#REGRESSION PARAMETERS REMOVE
mfit6 #x11 is the least significant
aik6<--6.6851
bic6<--6.3800
#********************REMOVING X11*************************************************
#*********************************************************************************
y8res <- lm(y8~x1+x2+x5+x7+x8+x10+x12+x14+x15,mydf)

lmres<-residuals(y8res,standardized=TRUE)
acf(lmres,60)#Issue at lag 6
pacf(lmres,60)#Issue at lag 6
acf(lmres^2,60) #Heteroscedasticity problem 
pacf(lmres^2,60)#Heteroscedasticity problem
Box.test(j^2,lag=60,type="Ljung")#Ho not rejected, no heterosceasticity issue

#AUTOCORRELATION MODELLING
#AR(13)
ar13lm<-arima(lmres,order=c(13,0,0),fixed=c(0,0,0,0,0,NA,0,0,0,0,0,0,NA,NA),transform.pars = FALSE)
ar13lm$aic#-1186
acf(ar13lm$residuals,50)#No issue
Box.test(ar13lm$residuals,lag=12,type="Ljung")#H0 not rejected:No autocorrelation problem
pacf(ar13lm$residuals,50)#No issue
Box.test(ar13lm$residuals^2,lag=12,type="Ljung")#pvalue>0.05 no heteroscedasticity problem

#HETEROSCEDASTICITY MODELING AFTER ARIMA
residuals7_prior<-residuals(ar13lm,standarized)
ar_residuals<-residuals7_prior
acf(ar_residuals^2,60)#problem at lag 11
pacf(ar_residuals^2,60)#problem at lag 11
Box.test(ar_residuals^2,lag=60,type="Ljung")#There is heteroscedasticity proble as pvalue<0.05

#NORMALITY
shapiro.test(ar_residuals)#H0 is not rejected so it is normal

#FORMULA AR-GRACH-LM:
# ARMA(13,0) AND GARCH(11,12) keeping in a11 ,b11 and use normal distribution as it fits my data better
ext <- data.matrix(mydf[,c(11:25)]) 
m.spec <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(11,11)), mean.model = list(armaOrder = c(13, 0), include.mean = FALSE, external.regressors = ext), distribution.model = "norm", start.pars = NULL,
                     fixed.pars = list(mxreg13=0,mxreg6=0,mxreg9=0,mxreg4=0,mxreg3=0,mxreg11=0,
                                       ar1=0,ar2=0,ar3=0,ar4=0,ar5=0,ar7=0,ar8=0,ar9=0,ar10=0,ar11=0,ar12=0,
                                       alpha1=0,alpha2=0,alpha3=0,alpha4=0,alpha5=0,alpha6=0,alpha7=0,alpha8=0,alpha9=0,alpha10=0,
                                       beta1=0,beta2=0,beta3=0,beta4=0,beta5=0,beta6=0,beta7=0,beta8=0,beta9=0,beta10=0) )
mfit7<- ugarchfit(data = j, spec = m.spec, solver = "hybrid", fit.control = list(scale = 1)) 

#RETRIEVE STANDARDIZED RESIDUALS 
residuals7<-ts(residuals(mfit7,standardize=TRUE),freq=1)

#AUTOCORRELATION:OK
acf(residuals7, 50)#Autocorrelation problem at lag 25(insignificant-small)
pacf(residuals7,50)#Autocorrelation ok
Box.test(residuals7,lag=100,type="Ljung")#Autocorellation ok

#HETEROSCEDASTICITY:OK
acf(residuals7^2, 100)#Heteroscedasticity fixed
pacf(residuals7^2,100)
Box.test(residuals6^2,lag=100,type="Ljung")#Heteroscedasticity fixed

#NORMALITY 0K
plot(mfit7)#type 9
shapiro.test((residuals7))#Ho not rejected residuals are normal!

#REGRESSION PARAMETERS REMOVE
mfit7 #x12 is the least significant
aik7= -6.7189
bic7= -6.4677
#********************REMOVING X12*************************************************
#*********************************************************************************
y8res <- lm(y8~x1+x2+x5+x7+x8+x10+x14+x15,mydf)

lmres<-residuals(y8res,standardized=TRUE)
acf(lmres,60)#Small issue at lag 5
pacf(lmres,60)#issue at lag 6
acf(lmres^2,60) 
pacf(lmres^2,60)
Box.test(j^2,lag=60,type="Ljung")#Ho not rejected, no heterosceasticity issue

#AUTOCORRELATION MODELLING
#AR(13) This time i will keep also ar5 ,ar13
ar13lm<-arima(lmres,order=c(13,0,0),fixed=c(0,0,0,0,NA,NA,0,0,0,0,0,0,NA,NA),transform.pars = FALSE)
ar13lm$aic
acf(ar13lm$residuals,50)#No issue
Box.test(ar13lm$residuals,lag=12,type="Ljung")#H0 not rejected:No autocorrelation problem
pacf(ar13lm$residuals,50)#No issue
Box.test(ar13lm$residuals^2,lag=50,type="Ljung")#pvalue>0.05 no heteroscedasticity problem

#HETEROSCEDASTICITY MODELING AFTER ARIMA
residuals8_prior<-residuals(ar13lm,standarized)
ar_residuals<-residuals8_prior
acf(ar_residuals^2,60)#problem at lag 11
pacf(ar_residuals^2,60)#problem at lag 9,10
Box.test(ar_residuals^2,lag=60,type="Ljung")#There is heteroscedasticity proble as pvalue<0.05

#NORMALITY
shapiro.test(ar_residuals)#H0 is not rejected so it is normal


#FORMULA AR-GRACH-LM:
# ARMA(13,0) keeping also ar5 ar6,ar13 AND GARCH(11,13) keeping in a9,a11,a13,b9,b10,b11,b13 and use normal distribution as it fits my data better
ext <- data.matrix(mydf[,c(11:25)]) 
m.spec <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(13,13)), mean.model = list(armaOrder = c(13, 0), include.mean = FALSE, external.regressors = ext), distribution.model = "norm", start.pars = NULL,
                     fixed.pars = list(mxreg13=0,mxreg6=0,mxreg9=0,mxreg4=0,mxreg3=0,mxreg11=0,mxreg12=0,
                                       ar1=0,ar2=0,ar3=0,ar4=0,ar7=0,ar8=0,ar9=0,ar10=0,ar11=0,ar12=0,
                                       alpha1=0,alpha2=0,alpha3=0,alpha4=0,alpha5=0,alpha6=0,alpha7=0,alpha8=0,alpha10=0,alpha12=0,
                                       beta1=0,beta2=0,beta3=0,beta4=0,beta5=0,beta6=0,beta7=0,beta8=0,beta9=0,beta12=0) )
mfit8<- ugarchfit(data = j, spec = m.spec, solver = "hybrid", fit.control = list(scale = 1)) 

#RETRIEVE STANDARDIZED RESIDUALS 
residuals8<-ts(residuals(mfit8,standardize=TRUE),freq=1)

#AUTOCORRELATION:OK
acf(residuals8, 50)#Autocorrelation problem at lag 25(insignificant-small)
pacf(residuals8,50)#Autocorrelation ok
Box.test(residuals8,lag=100,type="Ljung")#Autocorellation ok

#HETEROSCEDASTICITY: NOT FIXED EVEN THOUGH I INCLUDED GARCH (13,13)
acf(residuals8^2, 50)
pacf(residuals8^2,50)
Box.test(residuals8^2,lag=50,type="Ljung")

#NORMALITY 0K
plot(mfit8)#type 9 It is not a best fit but it is ok for now considering shapiro test was fine
shapiro.test((residuals8))#Ho not rejected residuals are normal!

#REGRESSION PARAMETERS REMOVE
mfit8 #x2 is the least significant
aik8= -6.5952
Bic8= -6.2722

#********************REMOVING X2*************************************************
#*********************************************************************************
y8res <- lm(y8~x1+x5+x7+x8+x10+x14+x15,mydf)

lmres<-residuals(y8res,standardized=TRUE)
acf(lmres,60)#Small issue at lag 5
pacf(lmres,60)#issue at lag 6
acf(lmres^2,60) 
pacf(lmres^2,60)
Box.test(j^2,lag=60,type="Ljung")#Ho not rejected, no heterosceasticity issue

#AUTOCORRELATION MODELLING
#AR(13) This time i will keep also ar12 
ar13lm<-arima(lmres,order=c(13,0,0),fixed=c(0,0,0,0,0,NA,0,0,0,0,0,NA,NA,NA),transform.pars = FALSE)
ar13lm$aic#-1177
acf(ar13lm$residuals,50)
Box.test(ar13lm$residuals,lag=12,type="Ljung")#H0 not rejected:No autocorrelation problem
pacf(ar13lm$residuals,50)
Box.test(ar13lm$residuals^2,lag=50,type="Ljung")#pvalue>0.05 no heteroscedasticity problem

#HETEROSCEDASTICITY MODELING AFTER ARIMA
residuals9_prior<-residuals(ar13lm,standarized)
ar_residuals<-residuals9_prior
acf(ar_residuals^2,60)#problem at lag 11
pacf(ar_residuals^2,60)#problem at lag 11
Box.test(ar_residuals^2,lag=60,type="Ljung")#There is heteroscedasticity proble as pvalue<0.05

#NORMALITY
shapiro.test(ar_residuals)#H0 is not rejected so it is normal

#FORMULA AR-GRACH-LM:
# ARMA(13,0) keeping also ar6,ar13 AND GARCH(11,11) keeping in a2,a11,b2,b11 and use normal distribution as it fits my data better
ext <- data.matrix(mydf[,c(11:25)]) 
m.spec <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(11,11)), mean.model = list(armaOrder = c(13, 0), include.mean = FALSE, external.regressors = ext), distribution.model = "norm", start.pars = NULL,
                     fixed.pars = list(mxreg13=0,mxreg6=0,mxreg9=0,mxreg4=0,mxreg3=0,mxreg11=0,mxreg12=0,mxreg2=0,
                                       ar1=0,ar2=0,ar3=0,ar4=0,ar5=0,ar7=0,ar8=0,ar9=0,ar10=0,ar11=0,ar12=0,
                                       alpha1=0,alpha2=0,alpha3=0,alpha4=0,alpha5=0,alpha6=0,alpha7=0,alpha8=0,alpha10=0,
                                       beta1=0,beta2=0,beta3=0,beta4=0,beta5=0,beta6=0,beta7=0,beta8=0,beta9=0,beta10=0) )
mfit9<- ugarchfit(data = j, spec = m.spec, solver = "hybrid", fit.control = list(scale = 1)) 

#RETRIEVE STANDARDIZED RESIDUALS 
residuals9<-ts(residuals(mfit9,standardize=TRUE),freq=1)

#AUTOCORRELATION:OK
acf(residuals9, 50)#Autocorrelation problem at lag 25(insignificant-small)
pacf(residuals9,50)#Autocorrelation ok
Box.test(residuals9,lag=100,type="Ljung")#Autocorellation ok

#HETEROSCEDASTICITY:OK
acf(residuals9^2,50)#Heteroscedasticity fixed
pacf(residuals9^2,50)
Box.test(residuals9^2,lag=100,type="Ljung")#Heteroscedasticity fixed

#NORMALITY 0K
plot(mfit9)#type 9
shapiro.test((residuals9))#Ho not rejected residuals are normal!

#REGRESSION PARAMETERS REMOVE
mfit9 #x5 is the least significant
aik9=-6.72
bic9= -6.48

#********************REMOVING X5*************************************************
#*********************************************************************************
y8res <- lm(y8~x1+x7+x8+x10+x14+x15,mydf)


lmres<-residuals(y8res,standardized=TRUE)
acf(lmres,60)#issue at lag 5
pacf(lmres,60)#issue at lag 6,13
acf(lmres^2,60) #small issue at lag 2,10
pacf(lmres^2,60)#small issue at lag 1,10
Box.test(j^2,lag=60,type="Ljung")#Ho not rejected, no heterosceasticity issue

#AUTOCORRELATION MODELLING
#AR(13) This time i will keep also ar12 
ar13lm<-arima(lmres,order=c(13,0,0),fixed=c(0,0,0,0,0,NA,0,0,0,0,0,NA,NA,NA),transform.pars = FALSE)
ar13lm$aic#-1178
acf(ar13lm$residuals,50)
Box.test(ar13lm$residuals,lag=12,type="Ljung")#H0 not rejected:No autocorrelation problem
pacf(ar13lm$residuals,50)
Box.test(ar13lm$residuals^2,lag=50,type="Ljung")#pvalue>0.05 no heteroscedasticity problem

#HETEROSCEDASTICITY MODELING AFTER ARIMA
residuals10_prior<-residuals(ar13lm,standarized)
ar_residuals<-residuals10_prior
acf(ar_residuals^2,60)#problem at lag 2,11
pacf(ar_residuals^2,60)#problem at lag 2,11
Box.test(ar_residuals^2,lag=60,type="Ljung")#There is heteroscedasticity proble as pvalue<0.05

#NORMALITY
shapiro.test(ar_residuals)#H0 is not rejected so it is normal

#FORMULA AR-GRACH-LM:
# ARMA(13,0) keeping also ar6,ar13 AND GARCH(11,11) keeping in a2,a11,b2,b11 and use normal distribution as it fits my data better
ext <- data.matrix(mydf[,c(11:25)]) 
m.spec <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(11,11)), mean.model = list(armaOrder = c(13, 0), include.mean = FALSE, external.regressors = ext), distribution.model = "norm", start.pars = NULL,
                     fixed.pars = list(mxreg13=0,mxreg6=0,mxreg9=0,mxreg4=0,mxreg3=0,mxreg11=0,mxreg12=0,mxreg2=0,mxreg5=0,
                                       ar1=0,ar2=0,ar3=0,ar4=0,ar5=0,ar7=0,ar8=0,ar9=0,ar10=0,ar11=0,ar12=0,
                                       alpha1=0,alpha2=0,alpha3=0,alpha4=0,alpha5=0,alpha6=0,alpha7=0,alpha8=0,alpha10=0,
                                       beta1=0,beta2=0,beta3=0,beta4=0,beta5=0,beta6=0,beta7=0,beta8=0,beta9=0,beta10=0) )
mfit10<- ugarchfit(data = j, spec = m.spec, solver = "hybrid", fit.control = list(scale = 1)) 

#RETRIEVE STANDARDIZED RESIDUALS 
residuals10<-ts(residuals(mfit10,standardize=TRUE),freq=1)

#AUTOCORRELATION:OK
acf(residuals10, 50)#Autocorrelation problem at lag 25(insignificant-small)
pacf(residuals10,50)#Autocorrelation ok
Box.test(residuals10,lag=100,type="Ljung")#Autocorellation ok

#HETEROSCEDASTICITY:OK
acf(residuals10^2,50)#Heteroscedasticity fixed
pacf(residuals10^2,50)
Box.test(residuals10^2,lag=100,type="Ljung")#Heteroscedasticity fixed

#NORMALITY 0K
plot(mfit10)#type 9
shapiro.test((residuals10))#Ho not rejected residuals are normal!

#REGRESSION PARAMETERS REMOVE
mfit10 #x7 is the least significant

aik= -6.7201
bic= -6.4868


#********************REMOVING X7*************************************************
#*********************************************************************************
y8res <- lm(y8~x1+x8+x10+x14+x15,mydf)


lmres<-residuals(y8res,standardized=TRUE)
acf(lmres,60)#issue at lag 5
pacf(lmres,60)#issue at lag 5
acf(lmres^2,60) #small issue at lag 2,101
Box.test(j^2,lag=60,type="Ljung")#Ho not rejected, no heterosceasticity issue

#AUTOCORRELATION MODELLING
#AR(13) This time i will keep also ar12 
ar13lm<-arima(lmres,order=c(13,0,0),fixed=c(0,0,0,0,0,NA,0,0,0,0,0,NA,NA,NA),transform.pars = FALSE)
ar13lm$aic#-1178
acf(ar13lm$residuals,50)
Box.test(ar13lm$residuals,lag=12,type="Ljung")#H0 not rejected:No autocorrelation problem
pacf(ar13lm$residuals,50)
Box.test(ar13lm$residuals^2,lag=50,type="Ljung")#pvalue>0.05 no heteroscedasticity problem

#HETEROSCEDASTICITY MODELING AFTER ARIMA
residuals11_prior<-residuals(ar13lm,standarized)
ar_residuals<-residuals11_prior
acf(ar_residuals^2,60)#problem at lag 11
pacf(ar_residuals^2,60)#problem at lag 11
Box.test(ar_residuals^2,lag=60,type="Ljung")#There is heteroscedasticity problem as pvalue<0.05

#NORMALITY
shapiro.test(ar_residuals)#H0 is not rejected so it is normal

#FORMULA AR-GRACH-LM:
# ARMA(13,0) keeping also ar6,ar13 AND GARCH(11,11) keeping in a2,a11,b2,b11 and use normal distribution as it fits my data better
ext <- data.matrix(mydf[,c(11:25)]) 
m.spec <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(11,11)), mean.model = list(armaOrder = c(13, 0), include.mean = FALSE, external.regressors = ext), distribution.model = "norm", start.pars = NULL,
                     fixed.pars = list(mxreg13=0,mxreg6=0,mxreg9=0,mxreg4=0,mxreg3=0,mxreg11=0,mxreg12=0,mxreg2=0,mxreg5=0,mxreg7=0,
                                       ar1=0,ar2=0,ar3=0,ar4=0,ar5=0,ar7=0,ar8=0,ar9=0,ar10=0,ar11=0,ar12=0,
                                       alpha1=0,alpha2=0,alpha3=0,alpha4=0,alpha5=0,alpha6=0,alpha7=0,alpha8=0,alpha10=0,
                                       beta1=0,beta2=0,beta3=0,beta4=0,beta5=0,beta6=0,beta7=0,beta8=0,beta9=0,beta10=0) )
mfit11<- ugarchfit(data = j, spec = m.spec, solver = "hybrid", fit.control = list(scale = 1)) 

#RETRIEVE STANDARDIZED RESIDUALS 
residuals11<-ts(residuals(mfit11,standardize=TRUE),freq=1)

#AUTOCORRELATION:OK
acf(residuals11, 50)#Autocorrelation problem at lag 25(insignificant-small)
pacf(residuals11,50)#Autocorrelation ok
Box.test(residuals11,lag=100,type="Ljung")#Autocorellation ok

#HETEROSCEDASTICITY:OK
acf(residuals11^2,50)#Heteroscedasticity fixed
pacf(residuals11^2,50)
Box.test(residuals11^2,lag=100,type="Ljung")#Heteroscedasticity fixed

#NORMALITY 0K
plot(mfit11)#type 9
shapiro.test((residuals11))#Ho not rejected residuals are normal!

#REGRESSION PARAMETERS REMOVE
mfit11 #x15 is the least significant

aik11=-6.7409
bic11=-6.5435

#********************REMOVING X15*************************************************
#*********************************************************************************
y8res <- lm(y8~x1+x8+x10+x14,mydf)


lmres<-residuals(y8res,standardized=TRUE)
acf(lmres,60)#issue at lag 5
pacf(lmres,60)#issue at lag 5
acf(lmres^2,60) #small issue at lag 2,101
Box.test(j^2,lag=60,type="Ljung")#Ho not rejected, no heterosceasticity issue

#AUTOCORRELATION MODELLING
#AR(13) This time i will keep also ar12 
ar13lm<-arima(lmres,order=c(13,0,0),fixed=c(0,0,0,0,0,NA,0,0,0,0,0,NA,NA,NA),transform.pars = FALSE)
ar13lm$aic#-1178
acf(ar13lm$residuals,50)
Box.test(ar13lm$residuals,lag=12,type="Ljung")#H0 not rejected:No autocorrelation problem
pacf(ar13lm$residuals,50)
Box.test(ar13lm$residuals^2,lag=50,type="Ljung")#pvalue<0.05 heteroscedasticity problem

#HETEROSCEDASTICITY MODELING AFTER ARIMA
residuals12_prior<-residuals(ar13lm,standarized)
ar_residuals<-residuals12_prior
acf(ar_residuals^2,60)#problem at lag 11
pacf(ar_residuals^2,60)#problem at lag 10
Box.test(ar_residuals^2,lag=60,type="Ljung")#There is heteroscedasticity problem as pvalue<0.05

#NORMALITY
shapiro.test(ar_residuals)#H0 is not rejected so it is normal

#FORMULA AR-GRACH-LM:
# ARMA(13,0) keeping also ar6,ar13 AND GARCH(12,13) keeping in a10,a12,b11,b13 and use normal distribution as it fits my data better
ext <- data.matrix(mydf[,c(11:25)]) 
m.spec <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(12,13)), mean.model = list(armaOrder = c(13, 0), include.mean = FALSE, external.regressors = ext), distribution.model = "norm", start.pars = NULL,
                     fixed.pars = list(mxreg13=0,mxreg6=0,mxreg9=0,mxreg4=0,mxreg3=0,mxreg11=0,mxreg12=0,mxreg2=0,mxreg5=0,mxreg7=0,mxreg15=0,
                                       ar1=0,ar2=0,ar3=0,ar4=0,ar5=0,ar7=0,ar8=0,ar9=0,ar10=0,ar11=0,ar12=0,
                                       alpha1=0,alpha2=0,alpha3=0,alpha4=0,alpha5=0,alpha6=0,alpha7=0,alpha8=0,alpha9=0,
                                       beta1=0,beta2=0,beta3=0,beta4=0,beta5=0,beta6=0,beta7=0,beta8=0,beta9=0,beta10=0,beta12=0) )
mfit12<- ugarchfit(data = j, spec = m.spec, solver = "hybrid", fit.control = list(scale = 1)) 

#RETRIEVE STANDARDIZED RESIDUALS 
residuals12<-ts(residuals(mfit12,standardize=TRUE),freq=1)

#AUTOCORRELATION:OK
acf(residuals12, 50)#Autocorrelation problem at lag 25(insignificant-small)
pacf(residuals12,50)#Autocorrelation ok
Box.test(residuals12,lag=100,type="Ljung")#Autocorellation ok

#HETEROSCEDASTICITY:OK
acf(residuals12^2,50)#Heteroscedasticity fixed
pacf(residuals12^2,50)
Box.test(residuals12^2,lag=100,type="Ljung")#Heteroscedasticity fixed

#NORMALITY 0K
plot(mfit12)#type 9
shapiro.test((residuals12))#Ho not rejected residuals are normal!

#REGRESSION PARAMETERS REMOVE
mfit12 #x10 is the least significant

aik12=-6.73
bic12=-6.52

#********************REMOVING X15*************************************************
#*********************************************************************************
y8res <- lm(y8~x1+x8+x14,mydf)


lmres<-residuals(y8res,standardized=TRUE)
acf(lmres,60)
pacf(lmres,60)
acf(lmres^2,60) 
Box.test(j^2,lag=60,type="Ljung")

#AUTOCORRELATION MODELLING
#AR(13) This time i will keep also ar12 
ar13lm<-arima(lmres,order=c(13,0,0),fixed=c(0,0,0,0,0,NA,0,0,0,0,0,NA,NA,NA),transform.pars = FALSE)
ar13lm$aic#-1180
acf(ar13lm$residuals,50)
Box.test(ar13lm$residuals,lag=12,type="Ljung")#H0 not rejected:No autocorrelation problem
pacf(ar13lm$residuals,50)
Box.test(ar13lm$residuals^2,lag=50,type="Ljung")#pvalue<0.05 heteroscedasticity problem

#HETEROSCEDASTICITY MODELING AFTER ARIMA
residuals13_prior<-residuals(ar13lm,standarized)
ar_residuals<-residuals13_prior
acf(ar_residuals^2,60)#problem at lag 10
pacf(ar_residuals^2,60)#problem at lag 10
Box.test(ar_residuals^2,lag=60,type="Ljung")#There is heteroscedasticity problem as pvalue<0.05

#NORMALITY
shapiro.test(ar_residuals)#H0 is not rejected so it is normal

#FORMULA AR-GRACH-LM:
# ARMA(13,0) keeping also ar6,ar13 AND GARCH(12,13) keeping in a10,a12,b11,b13 and use normal distribution as it fits my data better
ext <- data.matrix(mydf[,c(11:25)]) 
m.spec <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(11,11)), mean.model = list(armaOrder = c(13, 0), include.mean = FALSE, external.regressors = ext), distribution.model = "norm", start.pars = NULL,
                     fixed.pars = list(mxreg13=0,mxreg6=0,mxreg9=0,mxreg4=0,mxreg3=0,mxreg11=0,mxreg12=0,mxreg2=0,mxreg5=0,mxreg7=0,mxreg15=0,mxreg10=0,
                                       ar1=0,ar2=0,ar3=0,ar4=0,ar5=0,ar7=0,ar8=0,ar9=0,ar10=0,ar11=0,ar12=0,
                                       alpha1=0,alpha2=0,alpha3=0,alpha4=0,alpha5=0,alpha6=0,alpha7=0,alpha8=0,alpha9=0,
                                       beta1=0,beta2=0,beta3=0,beta4=0,beta5=0,beta6=0,beta7=0,beta8=0,beta9=0) )
mfit13<- ugarchfit(data = j, spec = m.spec, solver = "hybrid", fit.control = list(scale = 1)) 

#RETRIEVE STANDARDIZED RESIDUALS 
residuals13<-ts(residuals(mfit13,standardize=TRUE),freq=1)

#AUTOCORRELATION:OK
acf(residuals13, 50)#Autocorrelation problem at lag 25(insignificant-small)
pacf(residuals13,50)#Autocorrelation ok
Box.test(residuals13,lag=100,type="Ljung")#Autocorellation ok

#HETEROSCEDASTICITY:OK
acf(residuals13^2,50)#Heteroscedasticity fixed
pacf(residuals13^2,50)
Box.test(residuals13^2,lag=100,type="Ljung")#Heteroscedasticity fixed

#NORMALITY 0K
plot(mfit13)#type 9
shapiro.test((residuals13))#Ho not rejected residuals are normal!

#REGRESSION PARAMETERS REMOVE
mfit13 #x10 is the least significant

aik13=-6.75
bic13=-6.56

#x1 is on the limit of being signficant so I will keep it in my model

#*******QUESTION 5*********************************************************************************

#FORECASTS WITH AR13 of Question 1:
forc1<-predict(ar13fit,12)

#FORECASTS WITH REGRESSION of Question 2:
forecast2 <- data.frame(x1 = df[178:189, 10], x2 = df[178:189, 11], x3 = df[178:189, 12], 
                          x4 = df[178:189, 13], x5 = df[178:189, 14], x6 = df[178:189, 15], 
                          x7 = df[178:189, 16], x8 = df[178:189, 17], x9 = df[178:189, 18], 
                          x10 = df[178:189, 19], x11 = df[178:189, 20], x12 = df[178:189, 11], 
                          x13 = df[178:189, 22], x14 = df[178:189, 23], x15 = df[178:189, 24])

forc2 <- predict(regression1,forecast2,12)
forc2

#FORECAST WITH LM-ARMA-GARCH MODEL
forc3 <- ugarchforecast(mfit13,  external.forecasts = list(mregfor = ext),n.ahead = 12)
forc3<-fitted(forc3)
forc3<-ts(forc3, frequency=12, start = c(2005,1), end = c(2005,12))
print(forc3)

#****************************************************************************
#MSE
#MSE OF LM-GARCH-ARMA IN Q3:
mse3=0
for(i in c(1:length(forc3))){
  mse3=mse3+(ytest[i]-forc3[i])**2
  
}
cat("MSE=",mse3)

#HIT RATIO of LM_ARMA-GARCH MODEL
count=0
all=0
for(i in c(1:length(forc3))){
  if ((forc3[i]>0 & ytest[i]>0)|(forc3[i]<0 & ytest[i]<0)){
    count=count+1
  }
  all=all+1
}
hit_ratio3<-count/all
cat("hit_ratio is:",hit_ratio3)



