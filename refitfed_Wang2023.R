#################################
##########################the data from Wang han is 3 clim.trait. Rdata
##the data set name climate.trait)
load("~/E/eechange project/analysis/work2024/LMA_wang2023/Data/3 clim.trait.RData")
climate.trait$fD <- ifelse(climate.trait$DecEv == "D", climate.trait$f, NA)
climate.trait$fE <- ifelse(climate.trait$DecEv == "E", climate.trait$f, NA)

####fit the linear model based on Wang 2023 and control the other three are not changed

fixed_effectE <-0.5*log(climate.trait$IPAR,base=exp(1))+0.25*log(climate.trait$fE,base=exp(1))-0.01*climate.trait$T_g_m
modelE <- lm(log(LMA,base=exp(1)) ~ fixed_effectE+log(alpha_m,base=exp(1)),data=climate.trait)
# View summary of the model
summary(modelE)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                  3.14125    0.22111  14.207  < 2e-16 ***
#   fixed_effectE                1.48329    0.21124   7.022 3.17e-12 ***
#   log(alpha_m, base = exp(1)) -0.85784    0.07371 -11.639  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.5652 on 1685 degrees of freedom
# (3640 observations deleted due to missingness)
# Multiple R-squared:  0.3076,	Adjusted R-squared:  0.3067 
# F-statistic: 374.2 on 2 and 1685 DF,  p-value: < 2.2e-16
fixed_effectD <-1*log(climate.trait$IPAR,base=exp(1))+1*log(climate.trait$fE,base=exp(1))-0.05*climate.trait$T_g_m
modelD<-lm(log(LMA,base=exp(1)) ~ fixed_effectD+log(alpha_m,base=exp(1)),data=climate.trait)
summary(modelD)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                  3.99303    0.09494  42.058  < 2e-16 ***
#   fixed_effectD                0.45928    0.06182   7.429 1.72e-13 ***
#   log(alpha_m, base = exp(1)) -0.98992    0.05921 -16.720  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.5643 on 1685 degrees of freedom
# (3640 observations deleted due to missingness)
# Multiple R-squared:  0.3099,	Adjusted R-squared:  0.3091 
# F-statistic: 378.3 on 2 and 1685 DF,  p-value: < 2.2e-16
#########################################rework all the thing 
##########################################use the climate data included in in CPTD
#########################################only the LMA use the CPTD
#################################################################################
library(dplyr)
library(ggplot2)
library(tidyverse)
Traitdata<-read.csv("/Users/echo/E/eechange project/analysis/work2024/ecy2091-sup-0002-datas1/Hard Traits.csv")
sitedata<-read.csv("/Users/echo/E/eechange project/analysis/work2024/ecy2091-sup-0002-datas1/Sites.csv")
####I read the High res climate.csv, because the original High-res climate has format issue
climate <- read_csv("/Users/echo/E/eechange project/analysis/work2024/ecy2091-sup-0002-datas1/High res climate.csv")
siteandsample<-read.csv("/Users/echo/E/eechange project/analysis/work2024/ecy2091-sup-0002-datas1/Species translations.csv")

####check the Northeast China Transect (NECT) in sitedata is from 1-33 in Site.ID
siteNECT <- sitedata %>% filter(grepl("NECT", Site.Name))
write.csv(siteNECT,"/Users/echo/E/eechange project/analysis/work2024/ReworkALL/Sites_NECT.csv")
climateNECT<-subset(climate, `Site ID` >= 1 & `Site ID` <= 33)
ClimatecombineNECT<-data.frame(cbind(siteNECT,climateNECT))
write.csv(ClimatecombineNECT,"/Users/echo/E/eechange project/analysis/work2024/ReworkALL/NECT_Climate.csv")

siteandsampleNECT <- subset(siteandsample, Site.ID >= 1 & Site.ID <= 33) ##794
###combine the data of siteandsampleNECT and Traitdata
########combine the ClimatecombineNECT and siteandsampleNECT, where the ClimatecombineNECT is 33 obs, and siteandsampleNECT is 794, this is because  each site has have many sample data, so the Site.ID in siteandsampleNECT is repetitive 
########
NECT_trait_climate <- siteandsampleNECT %>%
  left_join(Traitdata, by = "SAMPLE.ID")

NECT_sample_climate <- NECT_trait_climate %>%
  left_join(ClimatecombineNECT, by = "Site.ID")
write.csv(NECT_sample_climate,"/Users/echo/E/eechange project/analysis/work2024/ReworkALL/NECT_Climate_trait.csv")
##################################################
##the database contains information on total annual photosynthetically active radiation during the growing season when mean daily temperatures are >0°C (PAR0), 
##the daily mean photosynthetically active radiation during the growing season (mPAR0), 
##growing degree days above a baseline of 0°C (GDD0), 
##the daily mean temperature during the growing season (mGDD0), 
##the ratio of actual to equilibrium evapotranspiration (α),
##a moisture index (MI)
##Defining the (thermal) growing season as the period when mean quasi-daily temperature (interpolated from monthly data) is above 0°C, 
##we calculated the ratio of growing-season length to the number of days in the year
wangCPT<-climate.trait %>% filter(grepl("CPT", source))

# Filter rows where the 'Site' column contains values from 1 to 33
wangCPT_filtered <- wangCPT %>%
  filter(site >= 1 & site <= 33)

wangCPT_filtered$fD <- ifelse(wangCPT_filtered$DecEv == "D", wangCPT_filtered$f, NA)
wangCPT_filtered$fE <- ifelse(wangCPT_filtered$DecEv == "E", wangCPT_filtered$f, NA)

lnlma_evergreen<-0.5*log(wangCPT_filtered$IPAR,base=exp(1))+0.25*log(wangCPT_filtered$fE,base=exp(1))-0.01*wangCPT_filtered$T_g_m-0.86*log(wangCPT_filtered$alpha_m,base=exp(1))+3.14
lma_evergreen<-exp(lnlma_evergreen)

lnlmanew1<-log(wangCPT_filtered$IPAR,base=exp(1))+log(wangCPT_filtered$fD,base=exp(1))-0.05*wangCPT_filtered$T_g_m-0.99*log(wangCPT_filtered$alpha_m,base=exp(1))+3.99
lma_D<-exp(lnlmanew1)
###Evergreen



newdata<-data.frame(cbind(wangCPT_filtered,lma_D,lma_evergreen))

newdata$lmamix <- ifelse(is.na(newdata$lma_evergreen), newdata$lma_D, newdata$lma_evergreen)

observedlma1<-newdata$LMA
prelma<-newdata$lmamix
lnobservedlma1<-log(observedlma1,base=exp(1))
lnprelma<-log(prelma,base=exp(1))
Tg_df<-data.frame(cbind(newdata$site,lnobservedlma1,lnprelma))
Tg_df <- na.omit(Tg_df)
lmaobsmean<-aggregate(Tg_df$lnobservedlma1, by=list(Tg_df$V1),mean)
lma.sd<-aggregate(Tg_df$lnobservedlma1, by=list(Tg_df$V1),FUN=sd)
lnpre<-aggregate(Tg_df$lnprelma, by=list(Tg_df$V1),mean)
map<-unique(Tg_df$V1)
ClimatecombineNECT<-read.csv("/Users/echo/E/eechange project/analysis/work2024/ReworkALL/NECT_Climate.csv")

select<-data.frame(ClimatecombineNECT$Site.ID,ClimatecombineNECT$MAP..mm.)
select_filtered <- select %>%
  filter(!ClimatecombineNECT.Site.ID %in% c(11, 12, 13, 14, 19))

lmadata<-data.frame(cbind(select_filtered$ClimatecombineNECT.MAP..mm.,lnpre$x,lma.sd$x,lmaobsmean$x))
names(lmadata) <- c("map", "lnprelma", "lma.sd","lmaobsmean")

plot1<- ggplot(data=lmadata,aes(x=map,y=lmaobsmean))+geom_point() + 
  geom_smooth(aes(y=lmaobsmean,x=map),method = 'lm') +
  geom_errorbar(aes(ymax=lmaobsmean+lma.sd,ymin=lmaobsmean-lma.sd),position=position_dodge(0.9),width=1.0)+
  xlab(expression('MAP'*~'('*'mm'*')')) + ylab(expression('lnLMA'*~'('*g*'/'*m^2*')')) +
  geom_point(aes(x=map,y=lnprelma),color="red") + 
  geom_smooth(aes(y=lnprelma,x=map, col='red'),method = 'lm') +
  theme_bw()+
  theme(legend.position="none",axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),
        axis.text=element_text(size=rel(1.5)),panel.grid.major=element_line(color="white"),panel.grid.minor=element_line(color="white"))+
  annotate("text", x=693, y=7, label="(c)",size=6.5)
plot1
############################################################re-fitted the LMA model
NECT_sample_climate<-read.csv("/Users/echo/E/eechange project/analysis/work2024/ReworkALL/NECT_Climate_trait.csv")
test<-data.frame(cbind(NECT_sample_climate$LMA..kg.m2.,NECT_sample_climate$Site.ID))
test<-na.omit(test)##498

modelevergreen<-lm(log(wangCPT_filtered$LMA) ~ log(wangCPT_filtered$IPAR,base=exp(1))+log(wangCPT_filtered$fE,base=exp(1))+wangCPT_filtered$T_g_m+log(wangCPT_filtered$alpha_m,base=exp(1)))

summary(modelevergreen) ##keep to use the Wang 2023


modeldeciduous<-lm(log(wangCPT_filtered$LMA)~log(wangCPT_filtered$IPAR,base=exp(1))+log(wangCPT_filtered$fD,base=exp(1))+wangCPT_filtered$T_g_m+log(wangCPT_filtered$alpha_m,base=exp(1)))
summary(modeldeciduous)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                   2.42461    1.34002   1.809 0.072461 .  
# log(wangCPT_filtered$IPAR, base = exp(1))     1.40981    0.67744   2.081 0.039184 *  
#   log(wangCPT_filtered$fD, base = exp(1))       3.30219    0.96246   3.431 0.000784 ***
#   wangCPT_filtered$T_g_m                       -0.05027    0.03470  -1.449 0.149553    
# log(wangCPT_filtered$alpha_m, base = exp(1)) -0.47774    0.28923  -1.652 0.100743    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.3295 on 145 degrees of freedom
# (348 observations deleted due to missingness)
# Multiple R-squared:  0.5338,	Adjusted R-squared:  0.521 
# F-statistic: 41.51 on 4 and 145 DF,  p-value: < 2.2e-16

lnlma_evergreen<-0.5*log(wangCPT_filtered$IPAR,base=exp(1))+0.25*log(wangCPT_filtered$fE,base=exp(1))-0.01*wangCPT_filtered$T_g_m-0.86*log(wangCPT_filtered$alpha_m,base=exp(1))+3.14
lma_E<-exp(lnlma_evergreen)

lnlma_deciduous<-log(1.41*wangCPT_filtered$IPAR,base=exp(1))+3.30*log(wangCPT_filtered$fD,base=exp(1))-0.05*wangCPT_filtered$T_g_m-0.48*log(wangCPT_filtered$alpha_m,base=exp(1))+2.42
lma_D<-exp(lnlma_deciduous)

lnlma_deciduous <- predict(modeldeciduous, data = wangCPT_filtered)

# Convert predictions back to the original scale
lma_D <- exp(lnlma_deciduous)
lma_D<-data.frame(lma_D)
lma_D <- cbind(Data.ID = rownames(lma_D), lma_D)
lma_D$Data.ID<-as.numeric(lma_D$Data.ID)

newdata<-data.frame(cbind(wangCPT_filtered,lma_E))
newdata$Data.ID <- seq_len(nrow(newdata))


newdata <- newdata %>%
  left_join(lma_D, by = "Data.ID")

newdata$lmamix <- ifelse(is.na(newdata$lma_E), newdata$lma_D, newdata$lma_E)

observedlma1<-newdata$LMA
prelma<-newdata$lmamix
lnobservedlma1<-log(observedlma1,base=exp(1))
lnprelma<-log(prelma,base=exp(1))
Tg_df<-data.frame(cbind(newdata$site,lnobservedlma1,lnprelma))
Tg_df <- na.omit(Tg_df)
lmaobsmean<-aggregate(Tg_df$lnobservedlma1, by=list(Tg_df$V1),mean)
lma.sd<-aggregate(Tg_df$lnobservedlma1, by=list(Tg_df$V1),FUN=sd)
lnpre<-aggregate(Tg_df$lnprelma, by=list(Tg_df$V1),mean)
map<-unique(Tg_df$V1)
##1  2  3  4  5  6  7  8  9 10 15 16 17 18 20 21 22 23 24 25 26 27 28 29 30 31 32 33
select<-data.frame(ClimatecombineNECT$Site.ID,ClimatecombineNECT$MAP..mm.)
select_filtered <- select %>%
  filter(!ClimatecombineNECT.Site.ID %in% c(11, 12, 13, 14, 19))

lmadata<-data.frame(cbind(select_filtered$ClimatecombineNECT.MAP..mm.,lnpre$x,lma.sd$x,lmaobsmean$x))
names(lmadata) <- c("map", "lnprelma", "lma.sd","lmaobsmean")

plot2<- ggplot(data=lmadata,aes(x=map,y=lmaobsmean))+geom_point() + 
  geom_smooth(aes(y=lmaobsmean,x=map),method = 'lm') +
  geom_errorbar(aes(ymax=lmaobsmean+lma.sd,ymin=lmaobsmean-lma.sd),position=position_dodge(0.9),width=1.0)+
  xlab(expression('MAP'*~'('*'mm'*')')) + ylab(expression('lnLMA'*~'('*g*'/'*m^2*')')) +
  geom_point(aes(x=map,y=lnprelma),color="red") + 
  geom_smooth(aes(y=lnprelma,x=map, col='red'),method = 'lm') +
  theme_bw()+
  theme(legend.position="none",axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),
        axis.text=element_text(size=rel(1.5)),panel.grid.major=element_line(color="white"),panel.grid.minor=element_line(color="white"))+
  annotate("text", x=693, y=6, label="(c)",size=6.5)
plot2
write.csv(lmadata,"/Users/echo/E/eechange project/analysis/work2024/LMA_wang2023/predicted_LMA.csv")
############################################
############################################################re-fitted the LMA model evergreen use Wang2023, decedeous use self data fitted
NECT_sample_climate<-read.csv("/Users/echo/E/eechange project/analysis/work2024/ReworkALL/NECT_Climate_trait.csv")
test<-data.frame(cbind(NECT_sample_climate$LMA..kg.m2.,NECT_sample_climate$Site.ID))
test<-na.omit(test)##498

modelevergreen<-lm(log(wangCPT_filtered$LMA) ~ log(wangCPT_filtered$IPAR,base=exp(1))+log(wangCPT_filtered$fE,base=exp(1))+wangCPT_filtered$T_g_m+log(wangCPT_filtered$alpha_m,base=exp(1)))

summary(modelevergreen) ##keep to use the Wang 2023 there is NA


modeldeciduous<-lm(log(wangCPT_filtered$LMA)~log(wangCPT_filtered$IPAR,base=exp(1))+log(wangCPT_filtered$fD,base=exp(1))+wangCPT_filtered$T_g_m+log(wangCPT_filtered$alpha_m,base=exp(1)))
summary(modeldeciduous)
lnlma_evergreen<-0.5*log(wangCPT_filtered$IPAR,base=exp(1))+0.25*log(wangCPT_filtered$fE,base=exp(1))-0.01*wangCPT_filtered$T_g_m-0.27*log(wangCPT_filtered$alpha_m,base=exp(1))+3.78
lma_E<-exp(lnlma_evergreen)

lnlma_deciduous<-log(1.41*wangCPT_filtered$IPAR,base=exp(1))+3.30*log(wangCPT_filtered$fD,base=exp(1))-0.05*wangCPT_filtered$T_g_m-0.48*log(wangCPT_filtered$alpha_m,base=exp(1))+2.42
lma_D<-exp(lnlma_deciduous)

lnlma_deciduous <- predict(modeldeciduous, data = wangCPT_filtered)

# Convert predictions back to the original scale
lma_D <- exp(lnlma_deciduous)
lma_D<-data.frame(lma_D)
lma_D <- cbind(Data.ID = rownames(lma_D), lma_D)
lma_D$Data.ID<-as.numeric(lma_D$Data.ID)

newdata<-data.frame(cbind(wangCPT_filtered,lma_E))
newdata$Data.ID <- seq_len(nrow(newdata))


newdata <- newdata %>%
  left_join(lma_D, by = "Data.ID")

newdata$lmamix <- ifelse(is.na(newdata$lma_E), newdata$lma_D, newdata$lma_E)

observedlma1<-newdata$LMA
prelma<-newdata$lmamix
lnobservedlma1<-log(observedlma1,base=exp(1))
lnprelma<-log(prelma,base=exp(1))
Tg_df<-data.frame(cbind(newdata$site,lnobservedlma1,lnprelma))
Tg_df <- na.omit(Tg_df)
lmaobsmean<-aggregate(Tg_df$lnobservedlma1, by=list(Tg_df$V1),mean)
lma.sd<-aggregate(Tg_df$lnobservedlma1, by=list(Tg_df$V1),FUN=sd)
lnpre<-aggregate(Tg_df$lnprelma, by=list(Tg_df$V1),mean)
map<-unique(Tg_df$V1)
##1  2  3  4  5  6  7  8  9 10 15 16 17 18 20 21 22 23 24 25 26 27 28 29 30 31 32 33
select<-data.frame(ClimatecombineNECT$Site.ID,ClimatecombineNECT$MAP..mm.)
select_filtered <- select %>%
  filter(!ClimatecombineNECT.Site.ID %in% c(11, 12, 13, 14, 19))

lmadata<-data.frame(cbind(select_filtered$ClimatecombineNECT.MAP..mm.,lnpre$x,lma.sd$x,lmaobsmean$x))
names(lmadata) <- c("map", "lnprelma", "lma.sd","lmaobsmean")

plot2<- ggplot(data=lmadata,aes(x=map,y=lmaobsmean))+geom_point() + 
  geom_smooth(aes(y=lmaobsmean,x=map),method = 'lm') +
  geom_errorbar(aes(ymax=lmaobsmean+lma.sd,ymin=lmaobsmean-lma.sd),position=position_dodge(0.9),width=1.0)+
  xlab(expression('MAP'*~'('*'mm'*')')) + ylab(expression('lnLMA'*~'('*g*'/'*m^2*')')) +
  geom_point(aes(x=map,y=lnprelma),color="red") + 
  geom_smooth(aes(y=lnprelma,x=map, col='red'),method = 'lm') +
  theme_bw()+
  theme(legend.position="none",axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),
        axis.text=element_text(size=rel(1.5)),panel.grid.major=element_line(color="white"),panel.grid.minor=element_line(color="white"))+
  annotate("text", x=693, y=6, label="(c)",size=6.5)
plot2
write.csv(lmadata,"/Users/echo/E/eechange project/analysis/work2024/LMA_wang2023/predicted_LMA_2025.csv")##

##################################Narea
vamaxdata<-read.csv("/Users/echo/E/eechange project/analysis/work2024/ReworkALL/chelsa_2006/vcmax_threemodel_obsvcmax.csv")
##fiited the model narea
# Create a new data frame combining relevant columns
data_combined <- data.frame(
  Narea = NECT_sample_climate$Narea..g.m2.,
  LMA = NECT_sample_climate$LMA..kg.m2. * 1000,  # Multiply by 1000 for correct units
  vcmax25_Own = vamaxdata$vcmax25_Own
)
data_combined <- na.omit(data_combined)
lnneramodel <- lm(Narea ~ LMA + vcmax25_Own-1, data = data_combined)
summary(lnneramodel)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# LMA         0.0212564  0.0008031  26.469  < 2e-16 ***
#   vcmax25_Own 0.0016956  0.0005230   3.242  0.00129 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.4243 on 376 degrees of freedom
# Multiple R-squared:  0.9282,	Adjusted R-squared:  0.9278 
# F-statistic:  2431 on 2 and 376 DF,  p-value: < 2.2e-16
vcmax25data<-data.frame(cbind(vamaxdata$Site_ID,vamaxdata$MAP_mm,vamaxdata$vcmax25_Own,NECT_sample_climate$Narea..g.m2.,log(NECT_sample_climate$Narea..g.m2.)))
names(vcmax25data) <- c("X", "MAP_mm", "vcmax25_Own","Narea..g.m2.","lnNarea..g.m2.")
vcmax25data<-na.omit(vcmax25data)

vcmax25<-aggregate(vcmax25data$vcmax25_Own, by=list(vcmax25data$X),mean)
lnNarea.sd<-aggregate(vcmax25data$lnNarea..g.m2., by=list(vcmax25data$X),FUN=sd)
lnNarea.mean<-aggregate(vcmax25data$lnNarea..g.m2., by=list(vcmax25data$X),mean)
MAP<-aggregate(vcmax25data$MAP_mm, by=list(vcmax25data$X),mean)

Madata<-data.frame(cbind(newdata$site,newdata$Rain,newdata$lmamix))
names(Madata) <- c("X","MAP", "lmamix")
Madata<-na.omit(Madata)
Ma<-aggregate(Madata$lmamix, by=list(Madata$X),mean)

datapre<-data.frame(cbind(MAP,vcmax25$x,lnNarea.mean$x,lnNarea.sd$x))
names(datapre) <- c("Site.ID", "MAP_mm", "vcmax25_Own","lnNarea_gm2_mean.","lnNarea_sd")

datapre <- datapre %>%
  filter(!Site.ID %in% c(11, 12, 13, 14, 19))

datapre<-data.frame(cbind(datapre,Ma))

datapre$predictednarea<-0.0212*datapre$x+0.0017*datapre$vcmax25_Own 
datapre$ln_prenaera<-log(datapre$predictednarea)

n4<- ggplot(data=datapre,aes(x=MAP_mm,y=lnNarea_gm2_mean.))+geom_point() + 
  geom_smooth(aes(y=lnNarea_gm2_mean.,x=MAP_mm),method = 'lm') +
  geom_errorbar(aes(ymax=lnNarea_gm2_mean.+lnNarea_sd,ymin=lnNarea_gm2_mean.-lnNarea_sd),position=position_dodge(0.9),width=1.0)+
  xlab(expression('MAP'*~'('*'mm'*')')) + ylab(expression('lnNarea'*~'('*g*'/'*m^2*')')) +
  geom_point(aes(x=MAP_mm,y=ln_prenaera),color="red") + 
  geom_smooth(aes(y=ln_prenaera,x=MAP_mm, col='red'),method = 'lm') +
  theme_bw()+
  theme(legend.position="none",axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),
        axis.text=element_text(size=rel(1.5)),panel.grid.major=element_line(color="white"),panel.grid.minor=element_line(color="white"))+
  annotate("text", x=693, y=1.5, label="(c)",size=6.5)
n4
write.csv(datapre,"/Users/echo/E/eechange project/analysis/work2024/LMA_wang2023/predicted_Narea.csv")
write.csv(datapre,"/Users/echo/E/eechange project/analysis/work2024/LMA_wang2023/predicted_Narea_2025.csv")##this is the lma evergreen use all the equation from Wang 2023



