# implied_volatility-yahoo_finance_data.R

################################################################################
# Heuristically trying to replicate characteristics of the volatility surface
# by calculating implied volatility for SPX (European) call options as of for
# 1) Various expiration dates and fixed strike price
# 2) Various strike prices and fixed expiration date (for options "in" and "out" the money)
################################################################################
# We manually extracted option prices from yahoo-finance on 24.05.2024. Since
# the examined options were not that liquid and the last traded price was not 
# even between current Bid and Ask, we proxied fair option prices as arithmetic 
# mean of Bid and Ask
################################################################################
# We use SPX closing price S=5267 on 23.05.2024
# Risk free interest rate in the US: r = 0.0525
# SPX dividend yield as of March 2024 is 0.0135
################################################################################


library(RQuantLib)

################################################################################
# 1) Various expiration dates, fixed strike price K = 5200
################################################################################

# Time to maturity
tau = c(7,14,21,28,42,68,147,175,210,238,273) 

# Option prices
Val = c((82.70+83.80),(96.60+98.30),(111.6+112.50),(120.8+121.7),
        (139.1+143.5),(176.1+176.9),(268.5+270.7), (304.9+307.6), 
        (340.5+342.3), (369.10+372.00),(401.50+405.1)) /2

# Expiration dates
dates = c("2024-05-31", "2024-06-07", "2024-06-14", "2024-06-21", 
          "2024-07-05", "2024-07-31", "2024-10-18", "2024-11-15",
          "2024-12-20", "2025-01-17", "2025-02-21")

ImpVol = numeric(11)

# Calculating implied volatilty
for (i in 1:length(tau)){
  ImpVol[i] = EuropeanOptionImpliedVolatility(type="call", value=Val[i], underlying = 5267, 
                                              strike = 5200, maturity = tau[i]/365, 
                                              riskFreeRate = 0.0525, dividendYield = 0.0135, volatility = 0.10)
}

# Saving as table
as.table(cbind(dates,tau,ImpVol))


################################################################################
# 2) Various strike prices, fixed expiration date 31.05.2024
################################################################################

# Option "in the money" ########################################################

# Strike prices
K = c(5200, 5150, 5100, 5050, 5000, 4500)

# Option prices
Val = c((89.9+90.5), (137.9 +136.5),(185.6+187),(232.2+240.6), (281.7+290),(780.5+788.9) )/2

ImpVol = numeric(length(K))
# vol = 15.9, 20.13, 24.51, 30.9, 35.16, 72

# Calculating implied volatilty
for (i in 1:length(K)){
  ImpVol[i] = EuropeanOptionImpliedVolatility(type="call", value=Val[i], underlying = 5267, 
                                              strike = K[i], maturity = 7/365, 
                                              riskFreeRate = 0.0525, dividendYield = 0.0135, volatility = 0.10)
}

ImpVol
# Saving as table
as.table(cbind(K,Val,ImpVol))

# Options "out of the money" ###################################################

# Strike prices
K = c(5300,5350, 5400, 5450, 5500, 6000)

# Option prices
Val = c( (19+19.4), (4.5+4.7), (0.7+0.8), (0.25+0.15), (0.05+0.20),(0+0.05))/2
# vol = 10, 9.25, 8, 9.8, 11.74,27

# Calculating implied volatilty
for (i in 1:length(K)){
  ImpVol[i] = EuropeanOptionImpliedVolatility(type="call", value=Val[i], underlying = 5267, 
                                              strike = K[i], maturity = 7/365, 
                                              riskFreeRate = 0.0525, dividendYield = 0.0135, volatility = 0.10)
}

ImpVol
# Saving as table
as.table(cbind(K,Val,ImpVol))
