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
# tvp model full model - dur = alpha*dur(-1)+ beta(dgdp-potential) + gamma*wages
#--------------------------------------------------------------------------------------------------------------------------

# Construct DLM

OkunsDLMfm <- dlm(
  
  
  FF = matrix(c(1,1,1,1),ncol = 4, byrow = TRUE),
  
  V = matrix(1),
  
  GG = matrix(c(1,0,0,0,
                0,1,0,0,
                0,0,1,0,
                0,0,0,1), ncol = 4, byrow = TRUE),
  
  W =  matrix(c(1,0,0,0,
                0,1,0,0,
                0,0,1,0,
                0,0,0,1), ncol = 4, byrow = TRUE),
  
  JFF = matrix(c(1,2,3,0),ncol = 4, byrow = TRUE),
  
  X = cbind(mod_data$dur_lag1,mod_data$d2lgdp, mod_data$d2lrulc_lag2),
  
  m0 = c(0,0,0,0),
  
  C0 = matrix(c(1e+07,0,0,0,
                0,1e+07,0,0,
                0,0,1e+07,0,
                0,0,0,1e+07), ncol = 4, byrow = TRUE)
  
)


buildOkunsFM <- function(p){
  
  #FF(OkunsDLMfm)[4] <-  

  V(OkunsDLMfm)  <- exp(p[2])
  
  GG(OkunsDLMfm)[1,1]  <- 1
  
  GG(OkunsDLMfm)[2,2]  <- 1
  
  GG(OkunsDLMfm)[3,3]  <- 1 
  
  GG(OkunsDLMfm)[4,4]  <- 1
  
  W(OkunsDLMfm)[1,1] <- exp(p[3])
  
  W(OkunsDLMfm)[2,2] <- 0
  
  W(OkunsDLMfm)[3,3] <- 0
  
  W(OkunsDLMfm)[4,4] <- exp(p[4])
  
  m0(OkunsDLMfm) <- c(0,0,0,p[1]*-4)
  
  C0(OkunsDLMfm)[1,1] <- 1
  
  C0(OkunsDLMfm)[4,4] <- 5
  
  
  return(OkunsDLMfm)
  
}



okuns.estfm <-  dlmMLE(y = mod_data$dur, parm = c(-0.049,-1.4,-6,-5), build = buildOkunsFM)


OkunsDLM1fm <- buildOkunsFM(okuns.estfm$par)


#--------------------------------------------------------------------------------------------------------------------------
# filtered and smoothed estimates
#--------------------------------------------------------------------------------------------------------------------------

filtered.fm <- dlmFilter(y = mod_data$dur, mod = OkunsDLM1fm)

smoothed <- dlmSmooth(y = mod_data$dur, mod = OkunsDLM1fm)

residuals(filtered)


data.frame(GDPgrowth = mod_data$d2lgdp[-1],
           Potential = dropFirst(filtered.fm$a[,4])/(-1*dropFirst(filtered.fm$a[,2])),
           Date = seq(as.Date("1960-12-01"),as.Date("2015-03-01"), by = "quarter")) %>% 
  gather(Var, Val, -Date) %>% 
  ggplot()+
  geom_line(aes(Date,Val, colour = Var))


States <- data.frame(Date = seq(as.Date("1960-12-01"),as.Date("2015-03-01"), by = "quarter"), 
                     GDPgrowth = mod_data$d2lgdp[-1],
                     Potential = dropFirst(filtered.fm$a[,4])/(-1*dropFirst(filtered.fm$a[,2])),
                     alpha = dropFirst(filtered.fm$a[,1]),
                     beta = dropFirst(filtered.fm$a[,2]),
                     gamma = dropFirst(filtered.fm$a[,3])
                     
                     
                     )


#--------------------------------------------------------------------------------------------------------------------------
# recreating estiamtes in table 1
#--------------------------------------------------------------------------------------------------------------------------

data.frame(alpha = last(States$alpha),
           beta = last(States$beta),
           `Potential output growth` = last(States$Potential),
           gamma = last(States$gamma),
           `Okuns coef` = 4*last(States$beta)/(1-last(States$alpha))
)
