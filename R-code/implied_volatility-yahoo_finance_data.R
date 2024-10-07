library(RQuantLib)
# Heuristically trying to replicate the Volatility Surface by calculating 
# Implied Volatility for SPX (European) Call Options as of 24.05.2024 for
# 1) Various Expiration Dates
# 2) Various Strike Prices

# Throughout we use the following data: 
# From Yahoo Finance: https://finance.yahoo.com/quote/%5ESPX/options
# SPX Closing Price on May 23th 2024: S = 5267
# Option Prices are provided there

# Risk free interest in US is currently about r = 0.0525
# c.f. https://www.tagesschau.de/wirtschaft/finanzen/fed-leitzins-146.html#:~:text=Hohe%20US%2DInflation%20Fed%20l%C3%A4sst%20Leitzins%20unver%C3%A4ndert&text=Angesichts%20der%20hartn%C3%A4ckig%20hohen%20Inflation,k%C3%B6nnen%20sich%20Gesch%C3%A4ftsbanken%20Zentralbankgeld%20leihen.

# SPX Dividend Yield is approximately (as of March 2024)
# 0.0135 
# c.f. https://ycharts.com/indicators/sp_500_dividend_yield#:~:text=Basic%20Info,month%20and%201.66%25%20last%20year.

# Remark: Since the specific listed Options are not that liquid and
# the last traded price is often not even in between Bid and Ask, we 
# approximate the Option Value as the arithmetic mean of Bid and Ask.

################################################################################
# 1) Various Expiration Dates
################################################################################
# Calculating Implied Volatility for various Expiration Dates
# for SPX European Call Options with Strike K = 5200

tau = c(7,14,21,28,42,68,147,175,210,238,273) # Time to maturity
Val = c((82.70+83.80),(96.60+98.30),(111.6+112.50),(120.8+121.7),
        (139.1+143.5),(176.1+176.9),(268.5+270.7), (304.9+307.6), 
        (340.5+342.3), (369.10+372.00),(401.50+405.1)) /2
dates = c("2024-05-31", "2024-06-07", "2024-06-14", "2024-06-21", 
          "2024-07-05", "2024-07-31", "2024-10-18", "2024-11-15",
          "2024-12-20", "2025-01-17", "2025-02-21")
ImpVol = numeric(11)

for (i in 1:length(tau)){
  ImpVol[i] = EuropeanOptionImpliedVolatility(type="call", value=Val[i], underlying = 5267, 
                                              strike = 5200, maturity = tau[i]/365, 
                                              riskFreeRate = 0.0525, dividendYield = 0.0135, volatility = 0.10)
}

as.table(cbind(dates,tau,ImpVol))
################################################################################
################################################################################

################################################################################
# 2) Various Strike Prices
################################################################################
# Calculating Implied Volatility for various Strike Prices
# for SPX European Call Options with expiration date "31.05.2024"

# Option "in the money"
K = c(5200, 5150, 5100, 5050, 5000, 4500)
Val = c((89.9+90.5), (137.9 +136.5),(185.6+187),(232.2+240.6), (281.7+290),(780.5+788.9) )/2
ImpVol = numeric(length(K))
# vol = 15.9, 20.13, 24.51, 30.9, 35.16, 72

for (i in 1:length(K)){
  ImpVol[i] = EuropeanOptionImpliedVolatility(type="call", value=Val[i], underlying = 5267, 
                                              strike = K[i], maturity = 7/365, 
                                              riskFreeRate = 0.0525, dividendYield = 0.0135, volatility = 0.10)
}
ImpVol
as.table(cbind(K,Val,ImpVol))

# Options "out of the money"
K = c(5300,5350, 5400, 5450, 5500, 6000)
Val = c( (19+19.4), (4.5+4.7), (0.7+0.8), (0.25+0.15), (0.05+0.20),(0+0.05))/2
# vol = 10, 9.25, 8, 9.8, 11.74,27

for (i in 1:length(K)){
  ImpVol[i] = EuropeanOptionImpliedVolatility(type="call", value=Val[i], underlying = 5267, 
                                              strike = K[i], maturity = 7/365, 
                                              riskFreeRate = 0.0525, dividendYield = 0.0135, volatility = 0.10)
}
ImpVol
as.table(cbind(K,Val,ImpVol))
