library('forecast')

Data <- read.csv('HistoricalQuotes.csv')
y <- rev(Data$close)[850:length(rev(Data$close))]


# #fit_arima <- Arima(y, order = c(0, 1, 1), seasonal = c(1, 1, 1), xreg = NULL,
#                    include.mean = T, include.drift = T)
                   
fit_arima <- auto.arima(y, stepwise = FALSE, approximation = FALSE, ic='aic', 
                        trace = T)

plot(forecast(fit_arima, h = 100))

fit_nnetar <- nnetar(y, size = 10, repeats = 20)
plot(forecast(fit_nnetar, h=100))
     