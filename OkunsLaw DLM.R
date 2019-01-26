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


tldata <- read_excel("TLdata2018.xlsx", sheet = "Sheet1", range = "A1:D238")


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

beta.start <- as.vector(coef(const_coef)[3])

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
  
  V(OkunsDLMfm)  <- exp(p[2])
  
  GG(OkunsDLMfm)[1,1]  <- 1
  
  GG(OkunsDLMfm)[2,2]  <- 1
  
  GG(OkunsDLMfm)[3,3]  <- 1 
  
  GG(OkunsDLMfm)[4,4]  <- 1
  
  W(OkunsDLMfm)[1,1] <- exp(p[3])
  
  W(OkunsDLMfm)[2,2] <- 0
  
  W(OkunsDLMfm)[3,3] <- 0
  
  W(OkunsDLMfm)[4,4] <- exp(p[4])
  
  m0(OkunsDLMfm) <- c(0,0,0,p[1]*4)
  
  C0(OkunsDLMfm)[1,1] <- 1
  
  C0(OkunsDLMfm)[4,4] <- 5
  
  
  return(OkunsDLMfm)
  
}



okuns.estfm <-  dlmMLE(y = mod_data$dur, parm = c(beta.start,-1.4,-6,-5), build = buildOkunsFM)


OkunsDLM1fm <- buildOkunsFM(okuns.estfm$par)


#--------------------------------------------------------------------------------------------------------------------------
# filtered and smoothed estimates
#--------------------------------------------------------------------------------------------------------------------------

filtered.fm <- dlmFilter(y = mod_data$dur, mod = OkunsDLM1fm)

smoothed <- dlmSmooth(y = mod_data$dur, mod = OkunsDLM1fm)

variances <-  dlmSvd2var(filtered.fm$U.C,filtered.fm$D.C) %>%
  lapply(sqrt) %>%
  sapply(diag) %>%
  t() %>% 
  data.frame()

variances$Date <- seq(as.Date("1960-06-01"),as.Date("2018-09-01"), by = "quarter")

Upperb <- filtered.fm$m[-1,2]+ variances[-1,2]#*qnorm(0.05,lower = FALSE)
Lowerb <- filtered.fm$m[-1,2]- variances[-1,2]#*qnorm(0.05,lower = FALSE)


Upperbp <- filtered.fm$m[-1,4]+variances[-1,4]#*qnorm(0.05,lower = FALSE)
Lowerbp <- filtered.fm$m[-1,4]-variances[-1,4]#*qnorm(0.05,lower = FALSE)



data.frame(GDPgrowth = mod_data$d2lgdp,
           Potential = dropFirst(filtered.fm$m[,4])/(-1*dropFirst(filtered.fm$m[,2])),
           upper = Upperbp/(abs(Upperb)),
           lower = Lowerbp/(abs(Lowerb)),
           Date = seq(as.Date("1960-09-01"),as.Date("2018-09-01"), by = "quarter")) %>% 
  gather(Var, Val, -Date) %>% 
  filter(Var != "GDPgrowth" & Date >= "1979-12-01") %>%
  spread(Var, Val) %>% 
  ggplot()+
  geom_line(aes(Date,Potential), colour = "red", size = 1)+
  geom_ribbon(aes(x = Date, ymin = lower, ymax = upper) , alpha = 0.2)+
  theme_bw()+
  theme(legend.position = "none")+
  ylim(0.5,8.5)+
  ylab("")+
  xlab("")+
  ggtitle("Potential output growth",subtitle = "annualised % change")+
  annotate("segment", x=ymd("2010-06-01"), xend =ymd("2012-06-01") , y= 6, yend= 4.5, arrow = arrow(), colour = "dark grey" )+
  annotate("text", x=ymd("2010-06-01") , y= 6.5, label = "+/- One S.D", colour = "dark grey")+
  annotate("segment", x=ymd("2012-06-01"), xend =ymd("2017-06-01") , y= 1, yend= 2.5, arrow = arrow(), colour = "red" )+
  annotate("text", x=ymd("2012-06-01") , y=0.75 , label = "Potenatial output growth", colour = "red")


data.frame(GDPgrowth = mod_data$d2lgdp,
           Potential = dropFirst(filtered.fm$m[,4])/(-1*dropFirst(filtered.fm$m[,2])),
           upper = Upperbp/(abs(Upperb)),
           lower = Lowerbp/(abs(Lowerb)),
           Date = seq(as.Date("1960-09-01"),as.Date("2018-09-01"), by = "quarter")) %>% 
  gather(Var, Val, -Date) %>% 
  filter(!Var %in% c("upper","lower") & Date >= "1979-12-01") %>%
  spread(Var, Val) %>% 
  ggplot()+
  geom_line(aes(Date,GDPgrowth), colour = "black", size = 0.5)+
  geom_line(aes(Date,Potential), colour = "red", size = 1)+
  theme_bw()+
  theme(legend.position = "none")+
  ylab("")+
  xlab("")+
  ggtitle("Potential and actual output growth",subtitle = "annualised % change")+
  annotate("text", x=ymd("2011-06-01") , y= 5.5, label = "Real GDP growth", colour = "black")+
  annotate("text", x=ymd("2012-06-01") , y=0 , label = "Potential output growth", colour = "red")



States <- data.frame(Date = seq(as.Date("1960-09-01"),as.Date("2018-09-01"), by = "quarter"), 
                     GDPgrowth = mod_data$d2lgdp,
                     Potential = dropFirst(filtered.fm$m[,4])/(-1*dropFirst(filtered.fm$m[,2])),
                     alpha = dropFirst(filtered.fm$m[,1]),
                     beta = dropFirst(filtered.fm$m[,2]),
                     gamma = dropFirst(filtered.fm$m[,3])
                     
                     
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

