# Silvia Goñi Mendia
# Advanced Time Series Analysis Project 2021

rm(list = ls())

###############################
# Upload and pre-process data #
###############################

# upload unemployment data
# Source: https://www.bde.es/webbde/es/estadis/infoest/series/be2402.csv
data1 = read.csv("unemployment.csv")
data1 = data1[6:84, 1:2]
colnames(data1) = c("time","unemployment")
data1$unemployment = as.numeric(data1$unemployment)

# upload tourism data
# Source: https://www.bde.es/webbde/es/estadis/infoest/series/ie0301.csv
data3 = read.csv("tourism.csv")
data3 = data3[447:683, c(1,15)]
colnames(data3) = c("time","tourism")
data3$tourism = as.numeric(data3$tourism)/1000000

# upload gdp data, in millions of Euros
# Source: https://datosmacro.expansion.com/pib/espana?anio=2021
data4 = read.csv("pib2.csv",sep = ";")
data4 = data4[2:80,]
colnames(data4) = c("year","quarter","gdp")
data4$gdp = as.numeric(data4$gdp)/10000

attach(data1)

###################################
# Univariate time series analysis #
###################################

# Declare variable "unemployment" as time series
unemployment_ts = ts(unemployment, frequency = 4, start = c(2002,1))

# Make a plot
par(mfrow=c(1,1))
ts.plot(unemployment_ts)
# Trend: No clear trend, positive, then negative
# Seasonality: some seasonal trend can be seen
# Stationarity: it does not look stationary

# Basic descriptive statistics
summary(unemployment_ts)

# Correlogram
acf(unemployment_ts)
# Clear autocorrelations until lag 3.5

# Take the log
logunemployment_ts = log(unemployment_ts)
ts.plot(logunemployment_ts)
# This series still does not look stationary
acf(logunemployment_ts)

# Check whether the trend in the log-transformed series is deterministic or stochastic
library(CADFtest)
dim(data1)
max.lag = round(sqrt(79))
CADFtest(logunemployment_ts, type = "trend", criterion= "BIC", max.lag.y=max.lag)
# The p-value = 0.94 > 5%, so we do not reject H0 and conclude that there is a stochastic trend

# Construct the series in log-differences 
dlogunemployment_ts = diff(log(unemployment_ts))
ts.plot(dlogunemployment_ts)
acf(dlogunemployment_ts)
# Some autocorrelations can be appreciated

# Seasonality
monthplot(dlogunemployment_ts)
# Seasonality can be appreciated, Q2 < Q3 < Q4 < Q1

# Unit root test to formally check stationarity
library(CADFtest)
max.lag = round(sqrt(length(dlogunemployment_ts))) 
CADFtest(dlogunemployment_ts, type= "drift", criterion= "BIC", max.lag.y=max.lag)
# The p-value = 0.06 > 5%, so we cannot reject H0 and conclude that the ts is not stationary

# Take seasonal differences 
dslogunemployment_ts = diff(diff(log(unemployment_ts)), lag = 4)
ts.plot(dslogunemployment_ts)
acf(dslogunemployment_ts)
# Some autocorrelations can be appreciated, corresponding to Q1 and Q3 in each lag

# Unit root test to formally check stationarity
library(CADFtest)
max.lag = round(sqrt(length(dlogunemployment_ts))) 
CADFtest(dslogunemployment_ts, type= "drift", criterion= "BIC", max.lag.y=max.lag)
# The p-value = 0.00 < 5%, so we reject H0 and conclude that the ts is stationary

# Correlogram and partial correlogram 
par(mfrow=c(2,1))
acf(dslogunemployment_ts)  #MA(?)
pacf(dslogunemployment_ts) #AR(1)
# There are some significant correlations and partial correlations 
# The time series is not white noise

# Ljung-Box test to check whether the ts is white noise
Box.test(dslogunemployment_ts, lag = max.lag, type = "Ljung-Box")
# The p-value = 0.00 < 5%, we reject H0 and conclude that the ts is not white noise

#####################
# LINEAR REGRESSION # 
#####################

# log(population_ts) = β0 + β1TREND + β2Q1 + β3Q2 + β4Q3 + ε
dim(data1) 
# we have 79 rows, however, it is not a divisible number so we eliminate the first data point
datab = data1[2:79,]
unemployment2_ts = ts(datab$unemployment, frequency = 4, start = c(2002,1))

TREND = 1:78
Q1 <- rep(c(1,0,0),26)
Q2 <- rep(c(0,1,0),26)
Q3 <- rep(c(0,0,1),26)
fit <- lm(log(unemployment2_ts) ~ TREND + Q1 + Q2 + Q3)
summary(fit)

# Q3 is omitted to avoid perfect multi-collinearity. 
# Assuming that the model is valid, we can interpret the output:
# R2 = 0.285, that is 28.5% of the variance of log(unemployment) is explained by the regressors.
# F-statistic is associated to p-value = 0.00 < 5%, thus we reject H0 and conclude that the regressors are jointly significant.
# The p-value = 0.00 < 5%, thus we reject H0 and conclude that TREND is significant.
# Every year "unemployment" decreases with 0.81% on average (due to TREND).

# Linear regression residual plot
ts.plot(fit$residuals)
Box.test(fit$residuals, lag = max.lag, type = "Ljung-Box")

# The residuals are clearly not white noise.
# The p-value = 0.00 < 5%. We reject H0 and conclude that the residuals are not white noise. The model is not valid.

############################
# Compare different models # 
############################

# Seasonal ARIMA(1,1,0)(0,1,0)
fit_sar <- arima(logunemployment_ts, order = c(1,1,0), seasonal = list(order = c(0,1,0))) 
fit_sar

# Check significance of terms
abs(fit_sar$coef/sqrt(diag(fit_sar$var.coef)))

# Check the validity of the model (plot, correlogram and Q-test of residuals)
ts.plot(fit_sar$residuals)
acf(fit_sar$residuals)
Box.test(fit_sar$residuals, lag = max.lag, type = "Ljung-Box")
# The p-value = 0.06 > 5%. We cannot reject H0 and conclude that the residuals are white noise. The model is valid.

# Seasonal ARIMA(0,1,4)(0,1,0) 
fit_sma4 <- arima(logunemployment_ts, order = c(0,1,4), seasonal = list(order = c(0,1,0))) 
fit_sma4

# Check significance of terms
abs(fit_sma4$coef/sqrt(diag(fit_sma4$var.coef)))
# The 4th term is not significant

# Check the validity of the model (plot, correlogram and Q-test of residuals)
ts.plot(fit_sma4$residuals)
acf(fit_sma4$residuals)
Box.test(fit_sma4$residuals, lag = max.lag, type = "Ljung-Box")
# The p-value = 0.66 > 5%. We cannot reject H0 and conclude that the residuals are white noise. 
# The model is valid, but as the last term is not significant, we will exclude it in the next model

# Seasonal ARIMA(0,1,3)(0,1,0) 
fit_sma3 <- arima(logunemployment_ts, order = c(0,1,3), seasonal = list(order = c(0,1,0))) 
fit_sma3

# Check significance of terms
abs(fit_sma3$coef/sqrt(diag(fit_sma3$var.coef)))

# Check the validity of the model (plot, correlogram and Q-test of residuals)
ts.plot(fit_sma3$residuals)
acf(fit_sma3$residuals)
Box.test(fit_sma3$residuals, lag = max.lag, type = "Ljung-Box")
# The p-value = 0.66 > 5%. We cannot reject H0 and conclude that the residuals are white noise. The model is valid.

# Seasonal ARIMA(0,1,2)(0,1,0) 
fit_sma2 <- arima(logunemployment_ts, order = c(0,1,2), seasonal = list(order = c(0,1,0))) 
fit_sma2

# Check significance of terms
abs(fit_sma2$coef/sqrt(diag(fit_sma2$var.coef)))

# Check the validity of the model (plot, correlogram and Q-test of residuals)
ts.plot(fit_sma2$residuals)
acf(fit_sma2$residuals)
Box.test(fit_sma2$residuals, lag = max.lag, type = "Ljung-Box")
# The p-value = 0.09 > 5%. We cannot reject H0 and conclude that the residuals are white noise. The model is valid.


# Seasonal ARIMA(0,1,1)(0,1,0) 
fit_sma1 <- arima(logunemployment_ts, order = c(0,1,1), seasonal = list(order = c(0,1,0))) 
fit_sma1

# Check significance of terms
abs(fit_sma1$coef/sqrt(diag(fit_sma1$var.coef)))
# All ARMA terms are significant.

# Check the validity of the model (plot, correlogram and Q-test of residuals)
ts.plot(fit_sma1$residuals)
acf(fit_sma1$residuals)
Box.test(fit_sma1$residuals, lag = max.lag, type = "Ljung-Box")
# The p-value = 0.03 < 5%. We reject H0 and conclude that the residuals are not white noise. The model is not valid.

# Valid models: SARIMA(1,1,0), SARIMA(0,1,2), SARIMA(0,1,3)

# Akaike Information Criterion
AIC(fit_sma3) # -275.84
AIC(fit_sma2) # -259.87
AIC(fit_sar) # -262.54

# Schwarz Information Criterion
AIC(fit_sma3, k=log(79)) # -266.36
AIC(fit_sma2, k=log(79)) # -252.76
AIC(fit_sar, k=log(79)) # -257.80

# Order of preference: SARIMA(0,1,3)(0,1,0), SARIMA(1,1,0)(0,1,0), SARIMA(0,1,2)(0,1,0)

############
# FORECAST #
############

# Forecast for the next 12 quarters using a SARIMA(0,1,3)(0,1,0) model
myforecastSMA3 = predict(fit_sma3,n.ahead=12)
expected = myforecastSMA3$pred
lower = myforecastSMA3$pred-qnorm(0.975)*myforecastSMA3$se;
upper = myforecastSMA3$pred+qnorm(0.975)*myforecastSMA3$se;
cbind(lower,expected,upper)

# Plot forecasted values and prediction interval
par(mfrow=c(1,1))
plot.ts(logunemployment_ts,xlim=c(2015,2026),ylim=c(0,4))
lines(expected,col="red")
lines(lower,col="blue")
lines(upper,col="blue")

# Forecast for the next 12 quarters using a SARIMA(1,1,0)(0,1,0) model
myforecastSAR = predict(fit_sar,n.ahead=12)
expected= myforecastSAR$pred
lower = myforecastSAR$pred-qnorm(0.975)*myforecastSAR$se;
upper = myforecastSAR$pred+qnorm(0.975)*myforecastSAR$se;
cbind(lower,expected,upper)

# Plot forecasted values and prediction interval
plot.ts(logunemployment_ts,xlim=c(2015,2026),ylim=c(0,4))
lines(expected,col="red")
lines(lower,col="blue")
lines(upper,col="blue")

# Absolute value loss
y = logunemployment_ts
S = round(0.75*length(y));h=1;
error1.h = c()
for (i in S:(length(y)-h))
{
  mymodel.sub = arima(y[1:i], order = c(1,1,0),seasonal=c(0,1,0))
  predict.h = predict(mymodel.sub,n.ahead=h)$pred[h]
  error1.h = c(error1.h,y[i+h]-predict.h)
}
error2.h = c()
for (i in S:(length(y)-h))
{
  mymodel.sub = arima(y[1:i], order = c(0,1,3),seasonal=c(0,1,0))
  predict.h = predict(mymodel.sub,n.ahead=h)$pred[h]
  error2.h = c(error2.h,y[i+h]-predict.h)
}
# Mean Absolute Forecast Error
mean(abs(error1.h));mean(abs(error2.h))
# The SARIMA(0,1,3)(0,1,0) has the lowest MAE

# Diebold-Mariano test
dm.test(error1.h,error2.h,h=h,power=1)
# We obtain a p-value = 0.26 > 5%, so we reject H0 and conclude that the forecast 
# performance of the two models, using the absolute value loss, is not significantly different.

# Squared value loss
error1.h = c()
for (i in S:(length(y)-h))
{
  mymodel.sub = arima(y[1:i], order = c(1,1,0),seasonal=c(0,1,0))
  predict.h = predict(mymodel.sub,n.ahead=h)$pred[h]
  error1.h = c(error1.h,y[i+h]-predict.h)
}
error2.h = c()
for (i in S:(length(y)-h))
{
  mymodel.sub = arima(y[1:i], order = c(0,1,3),seasonal=c(0,1,0))
  predict.h = predict(mymodel.sub,n.ahead=h)$pred[h]
  error2.h = c(error2.h,y[i+h]-predict.h)
}

# Compute the Mean Squared Forecast Error
mean(error1.h^2);mean(error2.h^2)
# The SARIMA(0,1,3)(0,1,0) has the lowest MAE

# Diebold-Mariano test
dm.test(error1.h,error2.h,h=h,power=2)
# We obtain a p-value = 0.11 > 5%, so we reject H0 and conclude that the forecast
# performance of the two models, using the squared value loss, is not significantly different.

#####################################
# Multivariate time series analysis #
#####################################

# Introduce other variables from the data set
gdp_ts <- ts(data4$gdp, frequency = 4,start = c(2002,1))
ts.plot(gdp_ts)

tourism_ts <- ts(data3$tourism, frequency = 12,start = c(2002,1))
tourism_ts <- aggregate(tourism_ts, nfrequency = 4) # change frequency of time series to quarterly
ts.plot(tourism_ts)

# Multiple graph
par(mfrow=c(2,1))
ts.plot(unemployment_ts,gdp_ts,tourism_ts,col=c("black","red","green"))
# It is very hard to see any type of relationship between the variables, except in 2020

# log(population_ts) = β0 + β1log(gdp_ts) + β2log(tourism_ts) + ε
fit2 <- lm(log(unemployment_ts) ~ log(gdp_ts) + log(tourism_ts))
summary(fit2)
# β1 = 0.89: if gdp increases by 1%, then unemployment increases by 0.89% on average
# β2 = 0.03: if tourism increases by 1%, then unemployment increases by 0.03% on average

ts.plot(fit2$residuals)
acf(fit2$residuals)
Box.test(fit2$residuals, lag = max.lag, type = "Ljung-Box")
# The residuals still show significant autocorrelations
# The p-value = 0.00 < 5%. We reject H0 and conclude that the residuals are not white noise. The model is not valid.
# The residuals are not white noise and the model cannot be validated