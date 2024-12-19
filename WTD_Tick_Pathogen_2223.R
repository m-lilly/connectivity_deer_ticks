##Analysis of nymphal tick abundance and pathogen infection

setwd("/users/marielilly/Desktop/R/UWIN_MH/Ch1_Manuscript_Code")
tick <- read.csv("WTD_Ticks_Pathogen_2022_2023.csv", header = T)

######Correlation tests for tick hazard predictor variables######
library(corrplot)
library(dplyr)
library(tidyr)
cor_tick <- tick %>%
  dplyr::select(Connect300v1_ParkBuffer1000, 
         Prop_BarrenLand_ParkBuffer1000, 
         PCT_Impervious_Park1000, 
         PCT_TCC_PARK_BUFFER1000, 
         HU_Density_Km_Park1000, 
         PCT_TCC_PARK,
         Occu_Predicted_Connectivity)

corrplot(cor(cor_tick), method = 'number')

cor_tick500 <- tick %>%
  dplyr::select(Connect300v1_ParkBuffer1000, 
                Prop_BarrenLand_ParkBuffer500, 
                PCT_Impervious_Park500, 
                PCT_TCC_PARK_BUFFER500, 
                HU_Density_Km_Park500, 
                PCT_TCC_PARK,
                Occu_Predicted_Connectivity)

corrplot(cor(cor_tick500), method = 'number')


cor_tick100 <- tick %>%
  dplyr::select(Connect300v1_ParkBuffer1000, 
                Prop_BarrenLand_ParkBuffer100, 
                PCT_Impervious_Park100, 
                PCT_TCC_PARK_BUFFER100, 
                HU_Density_Km_Park100, 
                PCT_TCC_PARK,
                Occu_Predicted_Connectivity)

corrplot(cor(cor_tick100), method = 'number')

#####Model construction for predictors of tick abundance######
#Remove pathogen data from dataframe to prevent NA issue at sites where no ticks were collected and is therefore NA for pathogen status
tick <- subset(tick, select = -c(Tested, DIN_Bburg, DIN_Babesia, DIN_Ap))

library(glmmTMB)
library(MuMIn)


####Buffer distances evaluation####
##Compare different buffers of each variable to first select buffer distance

##Barren land 500
M_Barren1000 <- glmmTMB(IS_N ~ Prop_BarrenLand_ParkBuffer1000 + offset(log(M_Dragged)) + (1|Site) + (1|Year), family = nbinom2(link="log"), data = tick) 
summary(M_Barren1000)

M_Barren500 <- glmmTMB(IS_N ~ Prop_BarrenLand_ParkBuffer500 + offset(log(M_Dragged)) + (1|Site) + (1|Year), family = nbinom2(link="log"), data = tick) 
summary(M_Barren500)

M_Barren100 <- glmmTMB(IS_N ~ Prop_BarrenLand_ParkBuffer100 + offset(log(M_Dragged)) + (1|Site) + (1|Year), family = nbinom2(link="log"), data = tick) 
summary(M_Barren100)

cand.mod <- list(M_Barren1000, M_Barren500, M_Barren100)
cand.names <- c("1000", "500", "100")
aictab(cand.set = cand.mod, cand.names)
##RESULTS! Basically all the same, buffer 500 has slightly lower AIC but none significant on their own

##TCC 1000
M_TCC1000 <- glmmTMB(IS_N ~ PCT_TCC_PARK_BUFFER1000 + offset(log(M_Dragged)) + (1|Site) + (1|Year), family = nbinom2(link="log"), data = tick) 
summary(M_TCC1000)

M_TCC500 <- glmmTMB(IS_N ~ PCT_TCC_PARK_BUFFER500 + offset(log(M_Dragged)) + (1|Site) + (1|Year), family = nbinom2(link="log"), data = tick) 
summary(M_TCC500)

M_TCC100 <- glmmTMB(IS_N ~ PCT_TCC_PARK_BUFFER100 + offset(log(M_Dragged)) + (1|Site) + (1|Year), family = nbinom2(link="log"), data = tick) 
summary(M_TCC100)

cand.mod <- list(M_TCC1000, M_TCC500, M_TCC100)
cand.names <- c("1000", "500", "100")
aictab(cand.set = cand.mod, cand.names)
##RESULTS: 1000 lowest AIC, all significant

##Impervious 500
M_I1000 <- glmmTMB(IS_N ~ PCT_Impervious_Park1000 + offset(log(M_Dragged)) + (1|Site) + (1|Year), family = nbinom2(link="log"), data = tick) 
summary(M_I1000)

M_I500 <- glmmTMB(IS_N ~ PCT_Impervious_Park500 + offset(log(M_Dragged)) + (1|Site) + (1|Year), family = nbinom2(link="log"), data = tick) 
summary(M_I500)
##Lowest AIC

M_I100 <- glmmTMB(IS_N ~ PCT_Impervious_Park100 + offset(log(M_Dragged)) + (1|Site) + (1|Year), family = nbinom2(link="log"), data = tick) 
summary(M_I100)

cand.mod <- list(M_I1000, M_I500, M_I100)
cand.names <- c("1000", "500", "100")
aictab(cand.set = cand.mod, cand.names)
##RESULTS: 500 lowest AIC but only 1 point lower than 1000, all significant


##Housing Unit Density 100
M_HU1000 <- glmmTMB(IS_N ~ HU_Density_Km_Park1000 + offset(log(M_Dragged)) + (1|Site) + (1|Year), family = nbinom2(link="log"), data = tick) 
summary(M_HU1000)

M_HU500 <- glmmTMB(IS_N ~ HU_Density_Km_Park500 + offset(log(M_Dragged)) + (1|Site) + (1|Year), family = nbinom2(link="log"), data = tick) 
summary(M_HU500)

M_HU100 <- glmmTMB(IS_N ~ HU_Density_Km_Park100 + offset(log(M_Dragged)) + (1|Site) + (1|Year), family = nbinom2(link="log"), data = tick) 
summary(M_HU100)
##Lowest AIC, not quite significant. None signfiicant

#####Occupancy model outputs  as predictor of ticks####
##Connectivity occupancy model
M_tick_deer <- glmmTMB(IS_N ~ Occu_Predicted_Connectivity + offset(log(M_Dragged)) + (1|Site) + (1|Year), family = nbinom2(link="log"), data = tick) 
summary(M_tick_deer)
##LOWEST AIC

##Model averaged occupancy model
M_tick_deer3 <- glmmTMB(IS_N ~ Occu_Predicted_ModelAverage + offset(log(M_Dragged)) + (1|Site) + (1|Year), family = nbinom2(link="log"), ziformula = ~1, data = tick) 
summary(M_tick_deer3)

##Tree canopy cover occupancy model
M_tick_deer6 <- glmmTMB(IS_N ~ Occu_Predicted_TCC + offset(log(M_Dragged)) + (1|Site) + (1|Year), family = nbinom2(link="log"), data = tick)
summary(M_tick_deer6)

##Impervious surface occupancy model
M_tick_deer7 <- glmmTMB(IS_N ~ Occu_Predicted_Impervious + offset(log(M_Dragged)) + (1|Site) + (1|Year), family = nbinom2(link="log"), data = tick)
summary(M_tick_deer7)


#####Building full models to dredge with best buffer distances#####

###Variables of interest with best buffers and model dredging

#Connect
M_connect <- glmmTMB(IS_N ~ Connect300v1_ParkBuffer1000 + Prop_BarrenLand_ParkBuffer500  + HU_Density_Km_Park100 + PCT_TCC_PARK + Park_Size_Ha + offset(log(M_Dragged)) + (1|Site) + (1|Year), family = nbinom2(link="log"), data = tick, na.action="na.fail")
summary(M_connect)
dredge(M_connect)
##Best model includes connectivity alone (AICc 424.4), next best connectivity and TCC within park, next best connectivity and HU

#Impervious
M_I <- glmmTMB(IS_N ~ PCT_Impervious_Park500 + Prop_BarrenLand_ParkBuffer500  + HU_Density_Km_Park100 + PCT_TCC_PARK + Park_Size_Ha + offset(log(M_Dragged)) + (1|Site) + (1|Year), family = nbinom2(link="log"), data = tick, na.action="na.fail") 
summary(M_I)
dredge(M_I)
##Best model includes Impervious alone (AICc 406.5) followed by Imp, TCC within park, Prop barren land (AICc 406.9), followed by IMP + TCC Park (AICc 407.1) followed by HU + IMP + TCC PARK (AIC 407.2)

#Occupancy
M_occu <- glmmTMB(IS_N ~ Occu_Predicted_Connectivity + Prop_BarrenLand_ParkBuffer500  + HU_Density_Km_Park100 + PCT_TCC_PARK + Park_Size_Ha + offset(log(M_Dragged)) + (1|Site) + (1|Year), family = nbinom2(link="log"), data = tick, na.action="na.fail") 
summary(M_occu)
dredge(M_occu)
##Best model includes Occupancy and TCC within park (404.6), followed by Occupancy +  TCC within park + Park size (404.9), followed by Occupancy +  TCC within park + barren land (406.1), followed by full model not far

#TCC
M_TCC <- glmmTMB(IS_N ~ PCT_TCC_PARK_BUFFER1000 + Prop_BarrenLand_ParkBuffer500  + HU_Density_Km_Park100 + PCT_TCC_PARK + Park_Size_Ha + offset(log(M_Dragged)) + (1|Site) + (1|Year), family = nbinom2(link="log"), data = tick, na.action="na.fail") 
summary(M_TCC)
dredge(M_TCC)
##Best model includes TCC and barren land (409.4), followed by TCC barren land and TCC within park


#####Build best models based on model dredging and compare#####
#Connect
M_connect1 <- glmmTMB(IS_N ~ Connect300v1_ParkBuffer1000 + offset(log(M_Dragged)) + (1|Site) + (1|Year), family = nbinom2(link="log"), data = tick) 
summary(M_connect1)

M_connect2 <- glmmTMB(IS_N ~ Connect300v1_ParkBuffer1000 + PCT_TCC_PARK + offset(log(M_Dragged)) + (1|Site) + (1|Year), family = nbinom2(link="log"), data = tick) 
summary(M_connect2)

#Impervious
M_I1 <- glmmTMB(IS_N ~ PCT_Impervious_Park500 + offset(log(M_Dragged)) + (1|Site) + (1|Year), family = nbinom2(link="log"), data = tick) 
summary(M_I)

M_I2 <- glmmTMB(IS_N ~ PCT_Impervious_Park500 + Prop_BarrenLand_ParkBuffer500  + PCT_TCC_PARK + offset(log(M_Dragged)) + (1|Site) + (1|Year), family = nbinom2(link="log"), data = tick) 
#summary(M_I2)

M_I3 <- glmmTMB(IS_N ~ PCT_Impervious_Park500 + PCT_TCC_PARK + offset(log(M_Dragged)) + (1|Site) + (1|Year), family = nbinom2(link="log"), data = tick) 
#summary(M_I3)

#Occupancy
M_occu1 <- glmmTMB(IS_N ~ Occu_Predicted_Connectivity + PCT_TCC_PARK + offset(log(M_Dragged)) + (1|Site) + (1|Year), family = nbinom2(link="log"), data = tick) 
summary(M_occu1)
r2(M_occu1)

M_occu2 <- glmmTMB(IS_N ~ Occu_Predicted_Connectivity + PCT_TCC_PARK + Park_Size_Ha + offset(log(M_Dragged)) + (1|Site) + (1|Year), family = nbinom2(link="log"), data = tick) 
summary(M_occu2)

M_occu3 <- glmmTMB(IS_N ~ Occu_Predicted_Connectivity + PCT_TCC_PARK + Prop_BarrenLand_ParkBuffer500 + offset(log(M_Dragged)) + (1|Site) + (1|Year), family = nbinom2(link="log"), data = tick) 
summary(M_occu3)

#TCC
M_TCC1 <- glmmTMB(IS_N ~ PCT_TCC_PARK_BUFFER1000 + Prop_BarrenLand_ParkBuffer500 + offset(log(M_Dragged)) + (1|Site) + (1|Year), family = nbinom2(link="log"), data = tick) 
summary(M_TCC1)

M_TCC2 <- glmmTMB(IS_N ~ PCT_TCC_PARK_BUFFER1000 + Prop_BarrenLand_ParkBuffer500 + PCT_TCC_PARK + offset(log(M_Dragged)) + (1|Site) + (1|Year), family = nbinom2(link="log"), data = tick) 
summary(M_TCC2)


M_occu_only <- glmmTMB(IS_N ~ Occu_Predicted_Connectivity + offset(log(M_Dragged)) + (1|Site) + (1|Year), family = nbinom2(link="log"), data = tick) 
summary(M_occu_only)
##Higher AIC than occupancy + other variables

####Model comparison####
library(AICcmodavg) # For GOF tests
cand.mod <- list(M_connect1, M_connect2, M_I1, M_I2, M_I3, M_occu1, M_occu2, M_occu3, M_TCC1, M_TCC2)
cand.names <- c("Connect", "Connect + TCC Park", "Impervious", "Impervious + Barren + TCC Park", "Impervious + TCC Park", "Occu + TCC Park", "Occu + TCC Park + Park size", "Occu + TCC Park + Barren", "TCC + Barren", "TCC + Barren + TCC Park")
aictab(cand.set = cand.mod, cand.names)

summary(M_occu1)
summary(M_occu3)


##########RESULTs:
##BEST model is Occu + TCC PARK, followed closely by Occu + TCC Park + Park size



######Spatial autocorrelation test#######
##Spatial autocorrelation? 
###MORANS I TEST for ticks
#calculate nearest neighbors to account for spatial correlation
library(tidyverse)  # Modern data science workflow
library(spdep)
library(spatialreg)
library(sf)
tick.coords <- cbind(tick$Longitude, tick$Latitude)
tick.5nn <- knearneigh(tick.coords, k=5, longlat = TRUE)
tick.5nn.nb <- knn2nb(tick.5nn)
plot(tick.5nn.nb, tick.coords)
listw <- nb2listw(tick.5nn.nb)


moran.mc(residuals(M_occu1), listw, 1000, zero.policy = T)

###NOT SIGNIFICANT, so do not need to do spatial model



#####Tick abundance figures######
library(ggiraph)
library(ggiraphExtra)
library(plyr)
library(MASS)
library(ggeffects)


#simplified model for visualization
M_tick_deer_p <- glm.nb(IS_N ~ Occu_Predicted_Connectivity, data = tick) 
summary(M_tick_deer_p)
pred_td2 <- ggpredict(M_tick_deer_p)
##Cut off at Y=300
plot(pred_td2, colors = "purple4",
     show_data = TRUE,
     alpha = 0.55,
     dot_alpha = 0.45,
     dot_size=4,
     line_size = 1)+
  ylim(0,300)+
  theme_classic(base_size = 16)+
  labs(title = NULL, x = "White-tailed deer occupancy", y = "Nymphal blacklegged tick density")



#Figure with tick and percent canopy cover within park
M_tick_TCC_p <- glm.nb(IS_N ~ PCT_TCC_PARK, data = tick) 
summary(M_tick_TCC_p)

#Cut off at Y=300
plot(pred_td2, colors = "darkgreen",
     show_data = TRUE,
     alpha = 0.55,
     dot_alpha = 0.45,
     dot_size=4,
     line_size = 1)+
  ylim(0,500)+
  theme_classic(base_size = 16)+
  labs(title = NULL, x = "Percent tree canopy cover within greenspace", y = "Nymphal blacklegged tick density")


######## 
###Model construction for predictors of pathogen infected tick abundance######
#re-import data to include pathogen variables
tick <- read.csv("WTD_Ticks_Pathogen_2022_2023.csv", header = T)
library(performance)

#Pathogen summaries
sum(tick2$DIN_Bburg)/sum(tick2$IS_N)
#0.2866712
sum(tick2$DIN_Ap)/sum(tick2$IS_N)
#0.03353756
sum(tick2$DIN_Babesia)/sum(tick2$IS_N)
#0.144808


###Borrelia burgdorferi models
M_occu_din <- glmmTMB(DIN_Bburg ~ Occu_Predicted_Connectivity + PCT_TCC_PARK + offset(log(M_Dragged))+ (1|Site) + (1|Year), family = nbinom2(link="log"), ziformula = ~1, data = tick)
summary(M_occu_din)
r2(M_occu_din)

M_occu_din3 <- glmmTMB(DIN_Bburg ~ Occu_Predicted_Connectivity + Park_Size_Ha + offset(log(M_Dragged))+ (1|Site) + (1|Year), family = nbinom2(link="log"), ziformula = ~1, data = tick)
summary(M_occu_din3)

hist(tick$DIN_Bburg)
M_occu_din2 <- glmmTMB(DIN_Bburg ~ Occu_Predicted_Connectivity + offset(log(M_Dragged))+ (1|Site) + (1|Year), family = nbinom2(link="log"), ziformula = ~1, data = tick)
summary(M_occu_din2)
r2(M_occu_din2)
#best model by 1.9 AIC score

M_occu_din2 <- glmmTMB(DIN_Bburg ~ Occu_Predicted_Connectivity + offset(log(M_Dragged))+ (1|Site) + (1|Year), family = nbinom2(link="log"), ziformula = ~1, data = tick)
summary(M_occu_din2)
r2(M_occu_din2)

M_connect_din <- glmmTMB(DIN_Bburg ~ Connect300v1_ParkBuffer1000 + offset(log(M_Dragged))+ (1|Site) + (1|Year), family = nbinom2(link="log"), ziformula = ~1, data = tick)
summary(M_connect_din)
r2(M_connect_din)

M_tcc_din <- glmmTMB(DIN_Bburg ~ PCT_TCC_PARK_BUFFER1000 + offset(log(M_Dragged))+ (1|Site) + (1|Year), family = nbinom2(link="log"), ziformula = ~1, data = tick)
summary(M_tcc_din)
r2(M_occu_din2)

M_imp_din <- glmmTMB(DIN_Bburg ~ PCT_Impervious_Park1000 + offset(log(M_Dragged))+ (1|Site) + (1|Year), family = nbinom2(link="log"), ziformula = ~1, data = tick)
summary(M_imp_din)

M_imp_din <- glmmTMB(DIN_Bburg ~ PCT_Impervious_Park1000 + offset(log(M_Dragged))+ (1|Site) + (1|Year), family = nbinom2(link="log"), ziformula = ~1, data = tick)
summary(M_imp_din)

M_huoccu_din <- glmmTMB(DIN_Bburg ~ Occu_Predicted_Connectivity + HU_Density_Km_Park1000 + offset(log(M_Dragged))+ (1|Site) + (1|Year), family = nbinom2(link="log"), ziformula = ~1, data = tick)
summary(M_huoccu_din)

M_hu_din <- glmmTMB(DIN_Bburg ~  HU_Density_Km_Park1000 + offset(log(M_Dragged))+ (1|Site) + (1|Year), family = nbinom2(link="log"), ziformula = ~1, data = tick)
summary(M_hu_din) # not significant


cand.mod <- list(M_occu_din3, M_imp_din, M_occu_din2, M_huoccu_din, M_connect_din, M_tcc_din)
#FROM ABOVE
cand.names <- c("Occupancy + Area", "Impervious", "Occupancy", "Housing Density + occupancy", "Connectivitiy", "TCC")
aictab(cand.set = cand.mod, cand.names)

###Babesia models
M_occu_Bab <- glmmTMB(DIN_Babesia ~ Occu_Predicted_Connectivity + offset(log(M_Dragged))+ (1|Site) + (1|Year), family = nbinom2(link="log"), ziformula = ~1, data = tick)
summary(M_occu_Bab)
r2(M_occu_Bab)

M_imp_Bab <- glmmTMB(DIN_Babesia ~ PCT_Impervious_Park1000 + offset(log(M_Dragged))+ (1|Site) + (1|Year), family = nbinom2(link="log"), ziformula = ~1, data = tick)
summary(M_imp_Bab)
r2(M_imp_Bab)

M_tcc_Bab <- glmmTMB(DIN_Babesia ~ PCT_TCC_PARK_BUFFER1000 + offset(log(M_Dragged))+ (1|Site) + (1|Year), family = nbinom2(link="log"), ziformula = ~1, data = tick)
summary(M_tcc_Bab)
r2(M_tcc_Bab)

M_hu_Bab <- glmmTMB(DIN_Babesia ~ HU_Density_Km_Park1000 + offset(log(M_Dragged))+ (1|Site) + (1|Year), family = nbinom2(link="log"), ziformula = ~1, data = tick)
summary(M_hu_Bab)
r2(M_hu_Bab)

M_connect_Bab <- glmmTMB(DIN_Babesia ~ Connect300v1_ParkBuffer1000 + offset(log(M_Dragged))+ (1|Site) + (1|Year), family = nbinom1(link="log"), ziformula = ~1, data = tick)
summary(M_connect_Bab)
r2(M_connect_Bab)

cand.mod <- list(M_occu_Bab, M_imp_Bab, M_tcc_Bab , M_hu_Bab, M_connect_Bab)
#FROM ABOVE
cand.names <- c("Occupancy", "Impervious", "TCC", "Housing Density", "Connectivitiy")
aictab(cand.set = cand.mod, cand.names)
##Impervious best model!


###Anaplasma models
M_occu_Ana <- glmmTMB(DIN_Ap ~ Occu_Predicted_Connectivity + offset(log(M_Dragged))+ (1|Site) + (1|Year), family = nbinom2(link="log"), ziformula = ~1, data = tick)
summary(M_occu_Ana)
r2(M_occu_Ana)

M_imp_Ana <- glmmTMB(DIN_Ap ~ PCT_Impervious_Park1000 + offset(log(M_Dragged))+ (1|Site) + (1|Year), family = nbinom2(link="log"), ziformula = ~1, data = tick)
summary(M_imp_Ana)
r2(M_imp_Ana)

M_tcc_Ana <- glmmTMB(DIN_Ap ~ PCT_TCC_PARK_BUFFER1000 + offset(log(M_Dragged))+ (1|Site) + (1|Year), family = nbinom2(link="log"), ziformula = ~1, data = tick)
summary(M_tcc_Ana)
r2(M_tcc_Ana)

M_hu_Ana <- glmmTMB(DIN_Ap ~ HU_Density_Km_Park1000 + offset(log(M_Dragged))+ (1|Site) + (1|Year), family = nbinom2(link="log"), ziformula = ~1, data = tick)
summary(M_hu_Ana)
r2(M_hu_Ana)

M_connect_Ana <- glmmTMB(DIN_Ap ~  Connect300v1_ParkBuffer1000 + offset(log(M_Dragged))+ (1|Site) + (1|Year), family = nbinom2(link="log"), ziformula = ~1, data = tick)
summary(M_connect_Ana)
r2(M_connect_Ana)

cand.mod <- list(M_occu_Ana, M_imp_Ana, M_tcc_Ana , M_hu_Ana, M_connect_Ana)
#FROM ABOVE
cand.names <- c("Occupancy", "Impervious", "TCC", "Housing Density", "Connectivitiy")
aictab(cand.set = cand.mod, cand.names)
##Impervious best model

######Pathogen infected tick abundance figures######

##Deer occupancy and Borrelia burgdorferi infected nymphal tick abundance
M_din_deer_p <- glm.nb(DIN_Bburg ~ Occu_Predicted_Connectivity, data = tick) 
summary(M_din_deer_p)

pred_dd <- ggpredict(M_din_deer_p)

plot(pred_dd, colors = "green",
     show_data = TRUE,
     alpha = 0.55,
     dot_alpha = 0.45,
     dot_size=4,
     line_size = 1)+
  theme_classic(base_size = 16)+
  ylim(0, 100)+
  labs(title = NULL, x = "White-tailed deer occupancy", y = "Boreelia burgdorferi infected nymph abundance")


##Impervious surface and Borrelia burgdorferi infected nymphal tick abundance
M_dinb_imp_p <- glm.nb(DIN_Bburg ~ PCT_Impervious_Park1000, data = tick) 
summary(M_dinb_imp_p)

pred_dd <- ggpredict(M_dinb_imp_p)

plot(pred_dd, colors = "green",
     show_data = TRUE,
     alpha = 0.55,
     dot_alpha = 0.45,
     dot_size=4,
     line_size = 1)+
  theme_classic(base_size = 16)+
  ylim(0, 100)+
  labs(title = NULL, x = "Percent impervious surface", y = "Density of infected nymphs")

##Impervious surface and Babesia infected nymphal tick abundance
M_din_bab_p <- glm.nb(DIN_Babesia ~ PCT_Impervious_Park1000, data = tick) 
summary(M_din_bab_p)

pred_dd <- ggpredict(M_din_bab_p)

plot(pred_dd, colors = "cyan2",
     show_data = TRUE,
     alpha = 0.55,
     dot_alpha = 0.45,
     dot_size=4,
     line_size = 1)+
  theme_classic(base_size = 16)+
  ylim(0, 100)+
  labs(title = NULL, x = "Percent impervious surface", y = "Density of Babesia infected nymphs")

M_din_bab_p <- glm.nb(DIN_Babesia ~ PCT_Impervious_Park1000, data = tick) 
summary(M_din_bab_p)

##Deer occupancy and Babesia infected nymphal tick abundance
M_din_bab_p <- glm.nb(DIN_Babesia ~ Occu_Predicted_Connectivity, data = tick) 
pred_dd <- ggpredict(M_din_bab_p)

plot(pred_dd, colors = "cyan2",
     show_data = TRUE,
     alpha = 0.55,
     dot_alpha = 0.45,
     dot_size=4,
     line_size = 1)+
  theme_classic(base_size = 16)+
  ylim(0, 100)+
  labs(title = NULL, x = "White-tailed deer occupancy", y = "Density of Babesia infected nymphs")

##Impervious surface and Babesia infected nymphal tick abundance
M_din_ap_p <- glm.nb(DIN_Ap ~ PCT_Impervious_Park1000, data = tick) 
summary(M_din_ap_p)

pred_dd <- ggpredict(M_din_ap_p)

plot(pred_dd, colors = "magenta",
     show_data = TRUE,
     alpha = 0.55,
     dot_alpha = 0.45,
     dot_size=4,
     line_size = 1)+
  theme_classic(base_size = 16)+
  ylim(0, 100)+
  labs(title = NULL, x = "Percent impervious surface", y = "Density of Anaplasma infected nymphs")

##Deer occupancy and Anaplasma infected nymphal tick abundance
M_din_ap_p <- glm.nb(DIN_Ap ~ Occu_Predicted_Connectivity, data = tick) 
summary(M_din_ap_p)

pred_dd <- ggpredict(M_din_ap_p)

plot(pred_dd, colors = "magenta",
     show_data = TRUE,
     alpha = 0.55,
     dot_alpha = 0.45,
     dot_size=4,
     line_size = 1)+
  theme_classic(base_size = 16)+
  ylim(0, 100)+
  labs(title = NULL, x = "White-tailed deer occupancy", y = "Density of Anaplasma infected nymphs")
