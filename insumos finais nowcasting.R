# Insumos Finais - Relatório - Nowcasting
# Biblioteca
library(meedr)
library(zoo)
library(lubridate)
library(sidrar)
library(zoo)
library(tidyverse)
library(forecast)
library(dplyr)
library(HDeconometrics)
library(glmnet)
library(Matrix)
library(neuralnet)
library(randomForest)
library(e1071)
library(ggplot2)
library(rbcb)
library(cli)
library(GetBCBData)
library(RColorBrewer)


# Preparações

## get QoQ from YoY Number
get_qoq <- function(yoy){
  library(sidrar)
  library(lubridate)
  library(stats)
  library(dplyr)
  library(seasonal)
  CNT = sidrar::get_sidra(api = "/t/1620/n1/all/v/all/p/all/c11255/90707/d/v583%202")
  pib = CNT  %>% 
    dplyr::select(`Trimestre (Código)`, Valor)
  colnames(pib) <- c("trimestre", "valor")
  pib.ts = ts(pib$valor,
              start = c(as.numeric(substr(pib$trimestre,1,4))[1],
                        as.numeric(substr(pib$trimestre,5,6))[1]), 
              end = c(as.numeric(substr(pib$trimestre,1,4))[nrow(pib)],
                      as.numeric(substr(pib$trimestre,5,6))[nrow(pib)]), frequency = 4)
  cresce = (1+yoy/100)*pib.ts[length(pib.ts)-3]
  Fpib = c(pib.ts,round(cresce, digits = 2))
  Fpib.ts = ts(Fpib, start = c(1996,01), frequency = 4)
  Fpib.meu.dessaz = seas(x = Fpib.ts, 
                         outlier.types = "all",
                         pickmdl.method = "best",
                         pickmdl.identify = "all",
                         forecast.maxlead = 6,
                         forecast.maxback = 0,
                         forecast.print = "none",
                         estimate = "",
                         x11.savelog = "q",
                         automdl = NULL)
  round(100*(Fpib.meu.dessaz$data[,1][nrow(Fpib.meu.dessaz$data)]/Fpib.meu.dessaz$data[,1][nrow(Fpib.meu.dessaz$data)-1]-1), digits = 1)
}

# transforms non detrended series into YoY rate of Change
yoy = function(series){
  index = 1:length(series)
  for(i in index){
    YoY <- series[12+index]/series[index]-1
  }
  return(ts(na.omit(YoY),start = c(start(series)[1]+1,start(series)[2]),frequency = 12))
}

# transforms quarterly into monthly data, with interpolation
qtr2month <- function(x, reference_month = 3, interpolation = FALSE){
  
  if(!reference_month %in% c(1,2,3)){
    stop("The reference_month should be 1,2 or 3")
  }
  
  if(!is.ts(x)){
    stop("x should be a ts object")
  }
  
  if(!is.null(dim(x))){
    stop("x should be a single ts object")
  }
  
  data_q <- zoo::as.Date(x)
  data_m <- seq(data_q[1], data_q[length(data_q)], by = 'months')
  out_x <- ts(rep(NA,length(data_m)),
              start =  as.numeric(c(substr(data_q[1],1,4), substr(data_q[1],6,7))),
              frequency = 12)
  out_x[data_m %in% data_q] <- x
  
  if(reference_month %in% c(2,3)){
    out_x <- stats::lag(out_x, -(reference_month-1))
    data_q <- zoo::as.Date(out_x)
    data_m <- seq(data_q[1], data_q[length(data_q)], by = 'months')
  }
  
  if(interpolation){
    xout <- zoo::as.yearmon(data_m)
    out_x <- stats::approx(out_x, xout = xout, method = "linear")
    out_x <- ts(out_x$y, start = out_x$x[1], end = out_x$x[length(out_x$x)], frequency = 12)
  }
  
  # output
  return(out_x)
}

# balances panel function, balances unbalanced matrix using statistical tecniques
balanced <- function(df) {
  df %>% 
    mutate_all(function(col) {
      if(all(is.na(col))) {
        return(col)
      } else if(any(is.na(col))) {
        col <- na.approx(col, na.rm = FALSE) # interpola as séries
        fit <- auto.arima(na.omit(col)) # omite NA's antes de fitar
        col[is.na(col)] <- forecast(fit, h = sum(is.na(col)))$mean
      }
      return(col)
    })
}

# convert quarterly dates into monthly ones
convert_date_pib <- function(quarterly_dates){
  # Initialize an empty vector to store the monthly dates
  monthly_dates <- c()
  
  # Loop over each quarterly date
  for (quarterly_date in quarterly_dates) {
    # Extract the year and the quarter
    year <- substr(quarterly_date, 1, 4)
    quarter <- as.numeric(substr(quarterly_date, 5, 6))
    
    # Calculate the first month of the quarter
    first_month <- (quarter - 1) * 3 + 1
    
    # Add each month of the quarter to the monthly_dates vector
    for (i in 0:2) {
      month <- sprintf("%02d", first_month + i)  # Make sure the month is two digits
      date_str <- paste0(year, "-", month, "-01")
      monthly_dates <- c(monthly_dates, date_str)
    }
  }
  
  # Convert the date strings to Date objects
  monthly_dates <- as.Date(monthly_dates)
  
  return(monthly_dates)
}

# convert sidra string into date format
convert_date_sidra <- function(date_str){
  as.Date(paste0(substr(date_str, 1, 4), "-", substr(date_str, 5, 6), "-01"))
}

# gets expectative given reference
exp <- function(ref){
  meedr::get_quarterly(indicator = "PIB Total",
                       first_date = Sys.Date()-30*365,
                       reference_date = ref)
}

# gets today expectative
GDP.now.expec <- function(trihoje){
  if (trihoje == "1") {
    exp <- function(ref){
      meedr::get_quarterly(indicator = "PIB Total",
                           first_date = Sys.Date()-30*365,
                           reference_date = ref)
    }
    PIB.Expect <- c()
    for (j in 2002:as.numeric(format(Sys.Date(), "%Y"))) {
      for (i in 1:4) {
        gdp_val <- exp(paste0(i,"/",j))$median[1]
        PIB.Expect[paste0(i,"/",j)] <- gdp_val
      }
    }
    EXPEC.PIB.HOJE <- PIB.Expect[1:(length(PIB.Expect)-3)]
    
  }
  if (trihoje == "2") {
    exp <- function(ref){
      meedr::get_quarterly(indicator = "PIB Total",
                           first_date = Sys.Date()-30*365,
                           reference_date = ref)
    }
    PIB.Expect <- c()
    for (j in 2002:as.numeric(format(Sys.Date(), "%Y"))) {
      for (i in 1:4) {
        gdp_val <- exp(paste0(i,"/",j))$median[1]
        PIB.Expect[paste0(i,"/",j)] <- gdp_val
      }
    }
    EXPEC.PIB.HOJE <- PIB.Expect[1:(length(PIB.Expect)-2)]
  }
  if (trihoje == "3") {
    exp <- function(ref){
      meedr::get_quarterly(indicator = "PIB Total",
                           first_date = Sys.Date()-30*365,
                           reference_date = ref)
    }
    PIB.Expect <- c()
    for (j in 2002:as.numeric(format(Sys.Date(), "%Y"))) {
      for (i in 1:4) {
        gdp_val <- exp(paste0(i,"/",j))$median[1]
        PIB.Expect[paste0(i,"/",j)] <- gdp_val
      }
    }
    EXPEC.PIB.HOJE <- PIB.Expect[1:(length(PIB.Expect)-1)]
  }
  if (trihoje == "4") {
    exp <- function(ref){
      meedr::get_quarterly(indicator = "PIB Total",
                           first_date = Sys.Date()-30*365,
                           reference_date = ref)
    }
    PIB.Expect <- c()
    for (j in 2002:as.numeric(format(Sys.Date(), "%Y"))) {
      for (i in 1:4) {
        gdp_val <- exp(paste0(i,"/",j))$median[1]
        PIB.Expect[paste0(i,"/",j)] <- gdp_val
      }
    }
    EXPEC.PIB.HOJE <- PIB.Expect
  }
  return(EXPEC.PIB.HOJE)
}

# date of recent quarter 
trihoje = quarter(Sys.Date())

# Algorithm to assemble time series of expectations
PIB.Expect <- c()
for (j in 2002:as.numeric(format(Sys.Date(), "%Y"))) {
  for (i in 1:4) {
    gdp_val <- exp(paste0(i,"/",j))$median[1]
    PIB.Expect[paste0(i,"/",j)] <- gdp_val
  }
}
print(PIB.Expect)
gdp.expec = GDP.now.expec(trihoje)
gdp.expec.ts = round(ts(gdp.expec,start = c(2002,1), frequency = 4), digits = 2)
monthly.gdp.expec.ts = qtr2month(gdp.expec.ts, reference_month = 3, interpolation = TRUE)
time.M.expec = as.Date(time(monthly.gdp.expec.ts))
plot(monthly.gdp.expec.ts)


# Main Code
PIB = get_sidra(api = "/t/5932/n1/all/v/6561/p/all/c11255/90707/d/v6561%201")
pib = ts(PIB$Valor, start = c(1996,01), frequency = 4)
datas.pib = PIB$`Trimestre (Código)`
pib.m = qtr2month(pib, reference_month = 3, interpolation = TRUE)
PIB.M = ts(c(pib.m, rep(tail(pib.m)[6],2)), start = c(1996,01), frequency = 12)
PIM = get_sidra(api = "/t/8888/n1/all/v/11602/p/all/c544/129314/d/v11602%201")
pim = ts(PIM$Valor, start = c(2002,01), frequency = 12)
PMCA = get_sidra(api = "/t/8881/n1/all/v/11709/p/all/c11046/56736/d/v11709%201")
pmca = ts(PMCA$Valor, start = c(2003,01), frequency = 12)
PMS = get_sidra(api = "/t/5906/n1/all/v/11624/p/all/c11046/56726/d/v11624%201")
pms = ts(PMS$Valor[13:nrow(PMS)], start = c(2012,01), frequency = 12)
IBC_BR = rbcb::get_series(c(IBC_BR = 24363), start_date = "01-01-2003")
datas.ibc = IBC_BR$date
novo_caged = ipeadatar::ipeadata(code = "CAGED12_SALDON12")
antigo_caged = ipeadatar::ipeadata(code = "CAGED12_SALDO12")
CAGED <- bind_rows(antigo_caged,novo_caged)
caged.a12 = CAGED %>%
  mutate(cumulative_12_months = rollapply(value, width = 12, FUN = sum, align = "right", fill = NA))
caged = ts(caged.a12$cumulative_12_months, start = c(1999,05), frequency = 12)
datas.caged = caged.a12$date
nuci <- rbcb::get_series(c(nuci = 24352), start_date = "01-01-2003")

# dataframes
df_pib <- data.frame(date = as.Date(time(PIB.M)), pib = PIB.M)
df_pim <- data.frame(date = as.Date(time(pim)), pim = PIM$Valor)
df_pmca <- data.frame(date = as.Date(time(pmca)), pmca = PMCA$Valor)
df_ibc <- data.frame(date = datas.ibc[13:length(datas.ibc)], ibc_br = round(100*yoy(IBC_BR$IBC_BR),digits = 2))
df_expec <- data.frame(date = time.M.expec, expec = monthly.gdp.expec.ts)
df_caged <- data.frame(date = datas.caged, caged_a12 = caged/10000000)
df_nuci <- data.frame(date = nuci$date[13:length(nuci$date)], nuci = yoy(nuci$nuci))
df_pms <- data.frame(date = as.Date(time(pms)), pms = pms)


# Balancear Painel
df_all <- Reduce(function(df1, df2) merge(df1, df2, by = "date", all = TRUE),
                 list(df_pib, 
                      df_pim, 
                      df_pmca, 
                      df_ibc,
                      df_expec,
                      df_caged,
                      df_nuci,
                      df_pms))
tail(df_all)
DF <- df_all[97:nrow(df_all),]
df <- DF[,-c(2)]
tail(df)
# balancear painel
balanced.DF <- balanced(df)
tail(balanced.DF)



# set de treino e set de teste




