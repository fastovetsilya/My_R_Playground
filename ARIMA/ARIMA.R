# This tool is currently under development


library('forecast')

Data <- read.csv('HistoricalQuotes.csv')
Data <- Data[-1,]
y <- rev(Data$close)[950:length(rev(Data$close))]
x <- as.numeric(rev(Data$volume)[1:length(rev(Data$volume))])

# #fit_arima <- Arima(y, order = c(0, 1, 1), seasonal = c(1, 1, 1), xreg = NULL,
#                    include.mean = T, include.drift = T)
                   
fit_arima <- auto.arima(y, d = NA, D = NA, max.p = 5, max.q = 5, max.P = 5,
                        max.Q = 5, max.order = 10, max.d = 5, max.D = 5, start.p = 1,
                        start.q = 1, start.P = 1, start.Q = 1, stationary = FALSE,
                        seasonal = TRUE, ic = 'aic', stepwise = FALSE,
                        trace = TRUE, approximation = FALSE,
                        truncate = NULL, xreg = NULL, test = 'kpss',
                        seasonal.test = c("ocsb", "ch"), allowdrift = TRUE, allowmean = TRUE,
                        lambda = NULL, biasadj = FALSE, parallel = TRUE, num.cores = 2)

plot(forecast(fit_arima, h = 100))


fit_nnetar <- nnetar(y, size = 10, repeats = 20)

     