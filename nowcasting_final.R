# BIBLIOTECAS (DEPOIS VAI TER QUE COLOCAR O COMANDO PARA INSTALAR AUTOMÁTICO SE NÃO TIVER O PACOTE (VITOR))
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
library(neuralnet)
library(randomForest)
library(e1071)
library(ggplot2)
library(RColorBrewer)

# FUNÇÕES
get_qoq <- function(yoy){
  library(sidrar)
  library(lubridate)
  library(stats)
  library(seasonal)
  CNT = sidrar::get_sidra(api = "/t/1620/n1/all/v/all/p/all/c11255/90707/d/v583%202")
  pib = CNT  %>% 
    select(`Trimestre (Código)`, Valor)
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
                         transform.function = "auto",
                         regression.aictest = c("td","easter"),
                         regression.variables = c("ls2008.04","TC2020.2","Easter[1]","TD"),
                         arima.model = "(4 0 0)(0 1 1)",
                         pickmdl.method = "best",
                         pickmdl.identify = "all",
                         outlier.types = "all",
                         forecast.maxlead = 6,
                         forecast.maxback = 0,
                         x11.savelog = "q")
  round(100*(Fpib.meu.dessaz$data[,1][nrow(Fpib.meu.dessaz$data)]/Fpib.meu.dessaz$data[,1][nrow(Fpib.meu.dessaz$data)-1]-1), digits = 1)
}
yoy = function(series){
  index = 1:length(series)
  for(i in index){
    YoY <- series[12+index]/series[index]-1
  }
  return(ts(na.omit(YoY),start = c(start(series)[1]+1,start(series)[2]),frequency = 12))
}
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
convert_date_sidra <- function(date_str){
  as.Date(paste0(substr(date_str, 1, 4), "-", substr(date_str, 5, 6), "-01"))
}
exp <- function(ref){
  meedr::get_quarterly(indicator = "PIB Total",
                       first_date = Sys.Date()-30*365,
                       reference_date = ref, be_quiet = TRUE)
}
GDP.now.expec <- function(trihoje){
  if (trihoje == "1") {
    exp <- function(ref){
      meedr::get_quarterly(indicator = "PIB Total",
                           first_date = Sys.Date()-30*365,
                           reference_date = ref,
                           be_quiet = TRUE)
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
                           reference_date = ref,
                           be_quiet = TRUE)
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
                           reference_date = ref,
                           be_quiet = TRUE)
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
                           reference_date = ref,
                           be_quiet = TRUE)
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
trihoje = quarter(Sys.Date())
# PREPARAR EXPECTATIVAS
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
# EXPECTATIVAS MENSAIS
plot(monthly.gdp.expec.ts)
# CÓDIGO PRINCIPAL
# BASE DE DADOS - enxuta
PIB = get_sidra(api = "/t/5932/n1/all/v/6561/p/all/c11255/90707/d/v6561%201")
pib = ts(PIB$Valor, start = c(1996,01), frequency = 4)
datas.pib = PIB$`Trimestre (Código)`
pib.m = qtr2month(pib, reference_month = 3, interpolation = TRUE)
PIB.M = ts(c(pib.m, rep(tail(pib.m)[6],2)), start = c(1996,01), frequency = 12)
PIM = get_sidra(api = "/t/8888/n1/all/v/11602/p/all/c544/129314/d/v11602%201")
pim = ts(PIM$Valor, start = c(2002,01), frequency = 12)
PMCA = get_sidra(api = "/t/8881/n1/all/v/11709/p/all/c11046/56736/d/v11709%201")
pmca = ts(PMCA$Valor, start = c(2003,01), frequency = 12)
PMS = get_sidra(api = "/t/8163/n1/all/v/11624/p/all/c11046/56726/c1274/56703/d/v11624%201")
datas.pms = PMS$`Mês (Código)`
IBC_BR = rbcb::get_series(c(IBC_BR = 24363), start_date = "01-01-2003")
datas.ibc = IBC_BR$date
# dataframes
df_pib <- data.frame(date = as.Date(time(PIB.M)), pib = PIB.M)
df_pim <- data.frame(date = as.Date(time(pim)), pim = PIM$Valor)
df_pmca <- data.frame(date = as.Date(time(pmca)), pmca = PMCA$Valor)
df_ibc <- data.frame(date = datas.ibc[13:length(datas.ibc)], ibc_br = round(100*yoy(IBC_BR$IBC_BR),digits = 2))
df_expec <- data.frame(date = time.M.expec, expec = monthly.gdp.expec.ts)
# merge que funfa
df_all <- Reduce(function(df1, df2) merge(df1, df2, by = "date", all = TRUE),
                 list(df_pib, 
                      df_pim, 
                      df_pmca, 
                      df_ibc, 
                      df_expec))

# confeccionar base final 
DF <- df_all[97:nrow(df_all),]
df <- DF[,-c(2)]
balanced.DF <- balanced(na.omit(df))
# treino e teste
treino = balanced.DF %>%
  filter(date < "2021-01-01")
teste = balanced.DF %>%
  filter(date > "2021-01-01")
# treino e teste pib
treino.PIB = DF[,1:2] %>%
  filter(date < "2021-01-01")
teste.PIB = DF[,1:2] %>%
  filter(date > "2021-01-01")
# modelos
# LASSO Projeção é o spot 27 do predict.lasso
LASSO = ic.glmnet(treino[,-1], treino.PIB[,-1], crit = "bic")
predict.lasso = predict(LASSO,newdata =teste[,-1])
plot(teste.PIB[,2], type = "l")
lines(predict.lasso, col = 2)
# adaLASSO
tau=1
first.step.coef=coef(LASSO)[-1]
penalty.factor=abs(first.step.coef+1/sqrt(nrow(teste)))^(-tau)
adalasso=ic.glmnet(treino[,2:ncol(treino)],treino.PIB$pib,crit="bic",penalty.factor=penalty.factor)
pred.adalasso=predict(adalasso,newdata=teste[,2:ncol(teste)])
# Ridge Regression
ridge.fit = cv.glmnet(as.matrix(treino[,-1]), treino.PIB[,-1], type = "mse", alpha = 0, family = "gaussian")
ridge.predict = predict(ridge.fit, s = ridge.fit$lambda.1se, newx = as.matrix(teste[,-1]))
# elastic net
en.fit = cv.glmnet(as.matrix(treino[,-1]), treino.PIB[,-1], type = "deviance", alpha = 0.5, family = "gaussian")
en.predicted = predict(en.fit, newx = as.matrix(teste[,-1]))
# redes neurais
# nn1
NN1.fit = neuralnet::neuralnet(treino.PIB[,-1] ~ ., data = as.matrix(treino[,-1]), hidden = 2, err.fct = "sse", threshold = 0.01, linear.output = TRUE)
predict.NN1 = predict(NN1.fit, newdata = teste[,-1])
# nn2
NN2.fit = neuralnet::neuralnet(treino.PIB[,-1] ~ ., data = as.matrix(treino[,-1]), hidden = 3, err.fct = "sse", threshold = 0.05, linear.output = TRUE)
predict.NN2 = predict(NN2.fit, newdata = teste[,-1])
# nn3
NN3.fit = neuralnet::neuralnet(treino.PIB[,-1] ~ ., data = as.matrix(treino[,-1]), hidden = c(4,2), err.fct = "sse", threshold = 0.1, linear.output = TRUE)
predict.NN3 = predict(NN3.fit, newdata = teste[,-1])
# tem que transformar no long o dataframe
df_pib.teste <- data.frame(date = teste.PIB[,1], pib = teste.PIB[,2])
df_predict_lasso <- data.frame(date = teste.PIB[,1], lasso = c(predict.lasso,NA,NA))
df_predict_adalasso <- data.frame(date = teste.PIB[,1], adaLASSO = c(pred.adalasso,NA,NA))
df_predict_ridge <- data.frame(date = teste.PIB[,1], RIDGE = c(ridge.predict,NA,NA))
df_predict_en <- data.frame(date = teste.PIB[,1], EN = c(en.predicted,NA,NA))
df_predict_nn1 <- data.frame(date = teste.PIB[,1], NN1 = c(predict.NN1,NA,NA))
df_predict_nn2 <- data.frame(date = teste.PIB[,1], NN2 = c(predict.NN2,NA,NA))
df_predict_nn3 <- data.frame(date = teste.PIB[,1], NN3 = c(predict.NN3,NA,NA))
plota.1 = bind_cols(df_pib.teste, df_predict_lasso[,2], df_predict_adalasso[,2], df_predict_ridge[,2], df_predict_en[,2], df_predict_nn1[,2], df_predict_nn2[,2], df_predict_nn3[,2])
colnames(plota.1) <- c("date","pib","lasso", "adaLASSO", "ridge","en", "NN1", "NN2","NN3")
# tentar plot style
plota.1 %>%
  ggplot2::ggplot() +
  ggplot2::aes(x = date) +
  ggplot2::geom_line(aes(y = pib, color = "(%)YQoQ"), linetype = "solid") +
  ggplot2::geom_line(aes(y = lasso, color = "LASSO"), linetype = "dashed") +
  ggplot2::geom_line(aes(y = adaLASSO, color = "adaLASSO"), linetype = "dashed") +
  ggplot2::geom_line(aes(y = ridge, color = "ridge"), linetype = "dashed") +
  ggplot2::geom_line(aes(y = en, color = "en"), linetype = "dashed") +
  ggplot2::geom_line(aes(y = NN1, color = "NN1"), linetype = "dashed") +
  ggplot2::geom_line(aes(y = NN2, color = "NN2"), linetype = "dashed") +
  ggplot2::geom_line(aes(y = NN3, color = "NN3"), linetype = "dashed") +
  ggplot2::scale_color_manual(values = c("(%)YQoQ" = "black", "LASSO" = "red", "adaLASSO" = "blue", "ridge" = "orange", "en" = "purple", "NN1" = "green", "NN2" = "pink", "NN3" = "grey")) +
  ggplot2::theme_grey() +
  ggplot2::labs(
    title    = paste("BRGDP IDB Nowcasting"),
    subtitle = "CSC/CBR intelectual property, Selected Machine Learning Models",
    x        = "Date",
    y        = "%QoQSA",
    caption  = "IDB Nowcasting",
    color    = "Line"
  )
predictions <- data.frame(
  Model = c("LASSO", "adaLASSO", "RIDGE", "Elastic Net", "NN1", "NN2", "NN3", "average"),
  GDP_YQoQ = c(predict.lasso[length(predict.lasso)],
                 pred.adalasso[length(pred.adalasso)],
                 ridge.predict[length(ridge.predict)],
                 en.predicted[length(en.predicted)],
                 predict.NN1[length(predict.NN1)],
                 predict.NN2[length(predict.NN2)],
                 predict.NN3[length(predict.NN3)],
                 mean(predict.lasso[length(predict.lasso)],
                      pred.adalasso[length(pred.adalasso)],
                      ridge.predict[length(ridge.predict)],
                      en.predicted[length(en.predicted)],
                      predict.NN1[length(predict.NN1)],
                      predict.NN2[length(predict.NN2)],
                      predict.NN3[length(predict.NN3)]))
)
ggplot(predictions, aes(x = Model, y = GDP_YQoQ, fill = Model)) +
  geom_bar(stat = "identity", color = "black") +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(x = "Model", y = "GDP_YQoQ", title = "BRGDP Nowcasting 2Q2023 (YQoQ)") +
  geom_text(aes(label = round(GDP_YQoQ, 2)), vjust = -0.3)
# dessaz e novo gráfico
predictions <- data.frame(
  Model = c("LASSO", "adaLASSO", "RIDGE", "Elastic Net", "NN1", "NN2", "NN3", "average"),
  GDP_YQoQ = c(get_qoq(predict.lasso[length(predict.lasso)]),
               get_qoq(pred.adalasso[length(pred.adalasso)]),
               get_qoq(ridge.predict[length(ridge.predict)]),
               get_qoq(en.predicted[length(en.predicted)]),
               get_qoq(predict.NN1[length(predict.NN1)]),
               get_qoq(predict.NN2[length(predict.NN2)]),
               get_qoq(predict.NN3[length(predict.NN3)]),
               get_qoq(mean(predict.lasso[length(predict.lasso)],
                    pred.adalasso[length(pred.adalasso)],
                    ridge.predict[length(ridge.predict)],
                    en.predicted[length(en.predicted)],
                    predict.NN1[length(predict.NN1)],
                    predict.NN2[length(predict.NN2)],
                    predict.NN3[length(predict.NN3)])))
)
ggplot(predictions, aes(x = Model, y = GDP_YQoQ, fill = Model)) +
  geom_bar(stat = "identity", color = "black") +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(x = "Model", y = "GDP_YQoQ", title = "BRGDP Nowcasting 2Q2023 (YQoQ)") +
  geom_text(aes(label = round(GDP_YQoQ, 2)), vjust = -0.3)

# plot de médio prazo - manual esta parte ainda tem que fazer função para automatizar

Values <- c(PIB$Valor[length(PIB$Valor)],PIB.Expect[length(PIB.Expect)-2], PIB.Expect[length(PIB.Expect)-1],PIB.Expect[length(PIB.Expect)])
medio.prazo <- data.frame(
  Data = c("1Qrt","2Qrt","3Qrt","4Qrt", "2023"),
  Value = c(PIB$Valor[length(PIB$Valor)],PIB.Expect[length(PIB.Expect)-2], PIB.Expect[length(PIB.Expect)-1],PIB.Expect[length(PIB.Expect)],mean(Values)))  
# Order the plots
medio.prazo$Data <- factor(medio.prazo$Data, levels = c("1Qrt","2Qrt","3Qrt","4Qrt", "2023"))
# Create plot
medio.prazo %>%
  ggplot(aes(x=Data, y=Value)) +
  geom_bar(aes(fill=Data), stat="identity", color="black", size=0.25) +
  scale_fill_brewer(palette="Set2") +
  geom_text(aes(label=round(Value, digits = 2)), vjust=-0.3, size=3.5) +
  theme_minimal() +
  labs(x="Quarters", y="YQoQ (%)", title="Real GDP YQoQ (%) Growth in 2023 Quarters") +
  theme(
    plot.title = element_text(hjust = 0.5, face="bold", size=14),
    axis.text = element_text(size=12),
    axis.title = element_text(size=14, face="bold"),
    legend.position = "none"
  )


