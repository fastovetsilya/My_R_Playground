# Load the library
library('forecast')

# Initialize parameters
Nsteps <- 15
Max_params <- 10
D <- NA # 1 or NA
Crit <- 'aicc' # 'aic', 'bic', or 'aicc'
root_test <- 'kpss' # 'kpss, 'adf', or 'pp'
seas_test <- 'ocsb' # 'ocsb', or 'ch'
Crit_point_1 <- 218.39
Crit_point_2 <- 37.7

# Load the data
Data <- read.csv('LSTM_Predicted.csv', header = FALSE)
D_r <- as.data.frame(Data[1:(nrow(Data)-15),])

# Compute processor time
ptm <- proc.time()

# Define custom ARIMA function
Custom_arima <- function(Data, Max_params, D)
{
fit_arima <- auto.arima(Data, d = NA, D = D, max.p = Max_params, max.q = Max_params, 
                        max.P = Max_params, max.Q = Max_params, max.order = 50, 
                        max.d = Max_params, max.D = Max_params, 
                        start.p = 1,
                        start.q = 1, start.P = 1, start.Q = 1, stationary = FALSE,
                        seasonal = TRUE, ic = Crit, stepwise = FALSE,
                        trace = TRUE, approximation = FALSE,
                        truncate = NULL, xreg = NULL, test = root_test,
                        seasonal.test = seas_test, allowdrift = TRUE, allowmean = TRUE,
                        lambda = NULL, biasadj = FALSE, parallel = TRUE, num.cores = 2)
return(fit_arima)
}

# Run ARIMA
arima_test <- Custom_arima(D_r, Max_params, D)
arima_predict <- Custom_arima(Data, Max_params, D)

# Run simple neural network
nnetar_test <- nnetar(D_r[,1], p = arima_test$arma[1], P = arima_test$arma[4], 
                      size = 10)
nnetar_predict <- nnetar(Data[,1], p = arima_predict$arma[1], P = arima_predict$arma[4], 
                         size = 10)

# Forecast Nsteps steps ahead
forecast_arimatest <- forecast(arima_test, h = Nsteps)
forecast_arimapredict <- forecast(arima_predict, h = Nsteps)
forecast_nnetartest <- forecast(nnetar_test, h = Nsteps)
forecast_nnetarpredict <- forecast(nnetar_predict, h =Nsteps)

# Plot neural network predictions
plot(forecast_nnetartest, 
     xlim = c(100, nrow(Data) + Nsteps)) 
lines(Data)
plot(forecast_nnetarpredict, 
     xlim = c(100, nrow(Data) + Nsteps)) 

# Plot ARIMA predictions
plot(forecast_arimatest, 
     xlim = c(100, nrow(Data) + Nsteps)) 
lines(Data)
plot(forecast_arimapredict, 
     xlim = c(100, nrow(Data)+Nsteps)) 
abline(h = Crit_point_1, col = 'green')
abline(h = Crit_point_2, col = 'red')

# Compute lower and upper bounds for the forecast
Lower <- as.numeric(forecast_arimapredict$lower[,2])
Upper <- as.numeric(forecast_arimapredict$upper[,2])

# Compute the risks as the proportion of the confidence interval tha covers critical points
Crit_1_risk <- signif((sum(Crit_point_1 - Lower[Lower < Crit_point_1])/
               sum(Upper - Lower)) * 100, 4)
Crit_2_risk <- signif((sum(Crit_point_2 - Lower[Lower < Crit_point_2])/
                sum(Upper-Lower)) * 100, 4)

# Print the results
print(proc.time() - ptm)
print(paste0('Risk of Crit_point_1 is: ', Crit_1_risk, ' %'))
print(paste0('Risk of Crit_point_2 is: ', Crit_2_risk, ' %'))