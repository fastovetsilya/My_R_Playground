# This tool is currently under development


library('forecast')
Data <- read.csv('LSTM_Predicted.csv')
fit_arima <- auto.arima(Data, d = NA, D = NA, max.p = 15, max.q = 15, max.P = 15,
                        max.Q = 15, max.order = 50, max.d = 15, max.D = 15, start.p = 1,
                        start.q = 1, start.P = 1, start.Q = 1, stationary = FALSE,
                        seasonal = TRUE, ic = 'aic', stepwise = FALSE,
                        trace = TRUE, approximation = FALSE,
                        truncate = NULL, xreg = NULL, test = 'kpss',
                        seasonal.test = "ocsb", allowdrift = TRUE, allowmean = TRUE,
                        lambda = NULL, biasadj = FALSE, parallel = TRUE, num.cores = 2)
plot(forecast(fit_arima, h = 100), xlim = c(150, nrow(Data)+20), 
     ylim = c(200,250))
fit_nnetar <- nnetar(Data[,1], size=20, repeats = 20)
plot(forecast(fit_nnetar, h = 100))
