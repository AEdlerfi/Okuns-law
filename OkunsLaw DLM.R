#--------------------------------------------------------------------------------------------------------------------------
# Okuns law and potential output - recreation of Lancaster and Tulips 2015 RDP
#--------------------------------------------------------------------------------------------------------------------------

setwd("C:/Users/aelde/OneDrive/Documents/GitHub/Okuns law")

library(dlm)
library(dplyr)
library(tidyverse)
library(readxl)


#--------------------------------------------------------------------------------------------------------------------------
# import data 
#--------------------------------------------------------------------------------------------------------------------------

# Data series used:
#   Unemployment rate
#   real unit labour costs (authors estimates - COE+TLS/(GDP*Deflator)) - data is spliced with rulc from NA growth
#   Real gdp


tldata <- read_excel("TLdata2015.xlsx", sheet = "Sheet1", range = "A1:D224")


#--------------------------------------------------------------------------------------------------------------------------
# data preperation
#--------------------------------------------------------------------------------------------------------------------------

dur <- tldata$ur-lag(tldata$ur)                                          # first difference of the unemployment rate
d2lgdp <- 400*(log(tldata$gdp)-log(lag(tldata$gdp,2)))/2                 # Two-quarter GDP growth, annualised log changes
d2lrulc <- 400*(log(tldata$rulc)-log(lag(tldata$rulc,2)))/2              # Two-quarter rulc growth, annualised log changes

mod_data <- data.frame(
  dur = dur,
  
  d2lgdp = d2lgdp,
  
  d2lrulc = d2lrulc
)

mod_data <- mod_data[-c(1:2),]

mod_data <- mod_data %>% 
  mutate(dur_lag1 = lag(dur),
         d2lrulc_lag2 =lag(d2lrulc,2)) %>% 
  filter(!is.na(d2lrulc_lag2))


#--------------------------------------------------------------------------------------------------------------------------
# constant coefficients model
#-------------------------------------------------------------------------------------------------------------------------

  
const_coef <- nls(formula = dur~ b1*dur_lag1 + b2*(d2lgdp-b0) + b3*d2lrulc_lag2 ,
    start = list(b0 =0.1, b1=0.1, b2=0.1, b3=0.1),
    data = mod_data) 



#--------------------------------------------------------------------------------------------------------------------------
# tvp model
#--------------------------------------------------------------------------------------------------------------------------

# Construct DLM

OkunsDLM <- dlm(
  
  
  FF = matrix(c(1,1,1),ncol = 3, byrow = TRUE),
  
  V = matrix(1),
  
  GG = matrix(c(1,0,0,
                0,1,0,
                0,0,1), ncol = 3, byrow = TRUE),
  
  W =  matrix(c(1,0,0,
                    0,1,0,
                    0,0,1), ncol = 3, byrow = TRUE),
  
  JFF = matrix(c(1,2,0),ncol = 3, byrow = TRUE),
  
  X = cbind(mod_data$dur_lag1,mod_data$d2lgdp),
  
  m0 = c(0,0,0),
  
  C0 = matrix(c(1e+07,0,0,
                    0,1e+07,0,
                    0,0,1e+07), ncol = 3, byrow = TRUE)
  
)

lag1 <- mod_data$dur[1]
beta1 <- 0.04


buildOkuns <- function(p){
  
  
  V(OkunsDLM)  <- exp(p[1])
  
  GG(OkunsDLM)[1,1]  <- 1
  
  GG(OkunsDLM)[2,2]  <- p[2] 
  
  GG(OkunsDLM)[3,3]  <- -1*p[2] 
  
  W(OkunsDLM)[1,1] <- exp(p[3])
  
  W(OkunsDLM)[2,2] <- 0
  
  W(OkunsDLM)[3,3] <- exp(p[4])
  
  m0(OkunsDLM) <- c(0,0,4)
  
  C0(OkunsDLM)[1,1] <- 1
  
  C0(OkunsDLM)[3,3] <- 5
  
  
  return(OkunsDLM)

  }



okuns.est <-  dlmMLE(y = mod_data$dur, parm = c(-1.4,-0.049,-4,-3), build = buildOkuns)


OkunsDLM1 <- buildOkuns(okuns.est$par)

smoothed <- dlmSmooth(y = mod_data$dur, mod = OkunsDLM1)