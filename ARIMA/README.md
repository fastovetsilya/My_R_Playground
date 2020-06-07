## LSTM-ARIMA for time series forcasting

This algorithm uses autimatic ARIMA to predict time series data. The hyperparameters
 of ARIMA are computed by exhaustive search and thus take some time.
The algorithm may be employed to predict stock market data (prices etc.).
One of my ideas is to apply this automatic ARIMA model to the representation
 of the time series obtained from LSTM (see my Python Playground). By doing
 so, we attempt to reduce the influence of stochastic processes in the market
 on the final prediction. The advantage of using ARIMA here is that it
 1) simple and 2) it calculates confidence bounds of the prediction to
 assess the risk of investments.
The algorithm also computes a simple neural network prediction based on the optimal parameters from ARIMA