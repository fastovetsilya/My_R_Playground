# This tool is currently under development


library('forecast')
Data <- read.csv('LSTM_Predicted.csv')
Max_params <- 10
fit_arima <- auto.arima(Data, d = NA, D = NA, max.p = Max_params, max.q = Max_params, 
                        max.P = Max_params, max.Q = Max_params, max.order = 50, 
                        max.d = Max_params, max.D = Max_params, 
                        start.p = 1,
                        start.q = 1, start.P = 1, start.Q = 1, stationary = FALSE,
                        seasonal = FALSE, ic = 'aic', stepwise = FALSE,
                        trace = TRUE, approximation = FALSE,
                        truncate = NULL, xreg = NULL, test = 'kpss',
                        seasonal.test = "ocsb", allowdrift = TRUE, allowmean = TRUE,
                        lambda = NULL, biasadj = FALSE, parallel = TRUE, num.cores = 2)
plot(forecast(fit_arima, h = 100), xlim = c(210, nrow(Data)+20), 
     ylim = c(200, 260))
fit_nnetar <- nnetar(Data[,1], size=20, repeats = 20)
plot(forecast(fit_nnetar, h = 30))
  