#########################################################################Fig 1
NECT_sample_climate<-read.csv("/Users/echo/E/eechange project/analysis/work2024/ReworkALL/chelsa_2006/Climate_trait_beta.csv")
pa1 <- ggplot(NECT_sample_climate,aes(x=NECT_sample_climate$MAP..mm., y=tg_chelsa_2006)) + geom_point() + 
  geom_smooth(method = 'lm', formula = y ~ x) + 
  ylab(expression('Tg'*~'('*degree*'C'*')')) + 
  xlab(expression('MAP'*~'('*'mm'*')'))+
  theme_bw()+
  theme( axis.title.x=element_text(size=15),axis.title.y=element_text(size=13),
         axis.text=element_text(size=rel(1.5)),panel.grid.major=element_line(color="white"),
         panel.grid.minor=element_line(color="white")) +
  annotate("text", x=693, y=11.2, label="(a)",size=6.5)
pa1
pa2 <- ggplot(NECT_sample_climate,aes(x=MAP..mm., y=NECT_sample_climate$vpd_chelsa_2006)) + geom_point() + 
  geom_smooth(method = 'lm', formula = y ~ x) + 
  ylab(expression('VPD'*~'('*'Pa'*')')) + 
  xlab(expression('MAP'*~'('*'mm'*')'))+
  annotate("text", x=693, y=800, label="(b)",size=6.5)+
  theme_bw()+
  theme( axis.title.x=element_text(size=15),axis.title.y=element_text(size=15),
         axis.text=element_text(size=rel(1.5)),panel.grid.major=element_line(color="white"),
         panel.grid.minor=element_line(color="white")) 
pa2
pa3 <- ggplot(NECT_sample_climate,aes(x=MAP..mm., y=ppfd_chelsa_2006)) + geom_point() + 
  geom_smooth(method = 'lm', formula = y ~ x) + 
  ylab(expression('PPFD'*~'('*'mol'*~m^-2*day^-1*')')) + 
  xlab(expression('MAP'*~'('*'mm'*')'))+
  theme_bw()+
  theme( axis.title.x=element_text(size=15),axis.title.y=element_text(size=15),
         axis.text=element_text(size=rel(1.5)),panel.grid.major=element_line(color="white"),
         panel.grid.minor=element_line(color="white")) +
  annotate("text", x=693, y=34, label="(c)",size=6.5)
pa3
############################################This figure is not update from plot_all_figure.R
plot_grid(pa1,pa2,pa3,label_x = 0.2,nrow = 2)
ggsave(filename = "Fig1_environmental.jpg", path = "/Users/echo/E/eechange project/analysis/work2024/ReworkALL/chelsa_2006/plot")
###################################################
############################################################
lmtemp<-lm(tg_chelsa_2006~MAP..mm.,data=NECT_sample_climate)
summary(lmtemp)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  9.2850383  0.1325922  70.027  < 2e-16 ***
#   MAP..mm.    -0.0009442  0.0002862  -3.299  0.00101 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 1.29 on 792 degrees of freedom
# Multiple R-squared:  0.01356,	Adjusted R-squared:  0.01231 
# F-statistic: 10.88 on 1 and 792 DF,  p-value: 0.001013
lmvpd<-lm(vpd_chelsa_2006~MAP..mm.,data=NECT_sample_climate)
summary(lmvpd)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 766.46555    6.59642  116.19   <2e-16 ***
#   MAP..mm.     -0.42976    0.01424  -30.18   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 64.18 on 792 degrees of freedom
# Multiple R-squared:  0.535,	Adjusted R-squared:  0.5344 
# F-statistic: 911.1 on 1 and 792 DF,  p-value: < 2.2e-16
lmppfd<-lm(ppfd_chelsa_2006~MAP..mm.,data=NECT_sample_climate)
summary(lmppfd)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 35.8113356  0.0693944  516.05   <2e-16 ***
#   MAP..mm.    -0.0100699  0.0001498  -67.23   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.6752 on 792 degrees of freedom
# Multiple R-squared:  0.8509,	Adjusted R-squared:  0.8507 
# F-statistic:  4520 on 1 and 792 DF,  p-value: < 2.2e-16

chicombine<-read.csv("/Users/echo/E/eechange project/analysis/work2024/ReworkALL/chelsa_2006/Chi_threemodel.csv")
site_summary <- chicombine %>%
  group_by(Site_ID) %>%
  summarize(
    MAP_mm = mean(MAP_mm, na.rm = TRUE),  # Mean MAP for each site
    vpd_pa= mean(VPD_CHELSA_2006, na.rm = TRUE),
    Chi_Obs_Mean = mean(Chi_Obs, na.rm = TRUE),
    Chi_Obs_SD = sd(Chi_Obs, na.rm = TRUE),
    Chi_Constant_Mean = mean(Chi_Constant, na.rm = TRUE),
    Chi_Alienor_Mean = mean(Chi_Alienor, na.rm = TRUE),
    Chi_Own_Mean = mean(Chi_Own, na.rm = TRUE)
  )
# Load ggplot2
library(ggplot2)
##############################################Fig.S2, Same as the plot_all_figure.R, not update
chi1<-ggplot(site_summary, aes(x = MAP_mm)) +
  # Add points and error bars for Chi_Obs
  geom_point(aes(y = Chi_Obs_Mean, color = "Chi_Obs")) +
  geom_point(aes(y = Chi_Constant_Mean, color = "Chi_Constant")) +
  geom_point(aes(y = Chi_Alienor_Mean, color = "Chi_Alienor")) +
  geom_point(aes(y = Chi_Own_Mean, color = "Chi_Own")) +
  geom_errorbar(
    aes(
      ymin = Chi_Obs_Mean - Chi_Obs_SD,
      ymax = Chi_Obs_Mean + Chi_Obs_SD,
      color = "Chi_Obs"
    ),
    width = 10
  ) +
  # Add smooth curves for other chi types
  geom_smooth(aes(y = Chi_Constant_Mean, color = "Chi_Constant"), method = "lm", se = TRUE, size = 1) +
  geom_smooth(aes(y = Chi_Alienor_Mean, color = "Chi_Alienor"), method = "lm", se = TRUE, size = 1) +
  geom_smooth(aes(y = Chi_Own_Mean, color = "Chi_Own"), method = "lm", se = TRUE, size = 1) +
  geom_smooth(aes(y = Chi_Obs_Mean, color = "Chi_Obs"), method = "lm", se = TRUE, size = 1) +
  # Customize labels and axis titles
  ylab(expression(~chi)) +
  xlab(expression('MAP'*~'('*'mm'*')')) +
  # Use a minimal theme and customize appearance
  theme_bw()+
  scale_color_manual(
    values = c(
      "Chi_Obs" = "blue",
      "Chi_Constant" = "purple",
      "Chi_Alienor" = "#65bab7",
      "Chi_Own" = "red"
    )
  ) +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text = element_text(size = rel(1.5)),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_line(color = "white")
  ) #+
#annotate("text", x=693, y=0.9, label="(a)",size=6.5)
# Annotate the plot
chi1## save work2024--plot
############################################Fig S3
####################################################################
################### Iabs is taken to be 80% of IL, the estimated canopy-average PPFD, given by inversion of Beer's law as IL ≈ −k I0 fv/ln (1 – fv)
##########https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/1365-2745.13967
########Dong ning journal of ecology
NECT_sample_climate<-read.csv("/Users/echo/E/eechange project/analysis/work2024/ReworkALL/chelsa_2006/Climate_trait_beta.csv")
model_constant_df<-read.csv("/Users/echo/E/eechange project/analysis/work2024/ReworkALL/chelsa_2006/Pmodel_constant_output.csv")
sw<-read.csv("/Users/echo/E/eechange project/analysis/work2024/ReworkALL/chelsa_2006/rsds/sw_w_m2_NECT.csv")
Tg_K<-NECT_sample_climate$tg_chelsa_2006+273.15
fpapr<-NECT_sample_climate$fapar_MODIS_2006

Gstar <- pressure*4.275*exp((37830/8.31447)*(1/298.15 - 1/(Tg_K)))/101325    #Gstar in Pa [Bernacchi et al. (2001)], KC, KO and ??* equal to 79430, 36380 and 37830 J/mol# standard values of KC, KO and ??* at 25?? equals to 39.37, 27480 and 4.22 Pa. O
Ca_constant<-model_constant_df$ca
Ci_constant<-model_constant_df$ci
m_constant <- (Ci_constant-Gstar)/(Ci_constant+2*Gstar)

Kc <- 40.49*exp((79430/8.31447)*(1/298.15 - 1/(Tg_K)))                      #Kc in Pa -> 1
Oc <-0.21*pressure #Oc in Pa -> month/daily
Ko <-27840*exp((36380/8.31447)*(1/298.15 - 1/(Tg_K)))                       #Ko in Pa -> month/daily
K  <- Kc*(1+(Oc/Ko)) 
k=0.5
value<-1/(1-fpapr)

PPFD<-rowMeans(sw[, 2:13], na.rm = TRUE)*0.45/0.22 ##here is I0
I<-PPFD*k*fpapr/log(value)

Vcmax_constant <- 0.093*0.8*I*((Ci_constant+K)/(Ci_constant+2*Gstar))
Vcmax25_constant <- Vcmax_constant*exp((65330/8.31447)*((1/Tg_K)-(1/298.15)))
####*sqrt(1-(0.41/m_constant)^(2/3))
model_Alienor_df<-read.csv("/Users/echo/E/eechange project/analysis/work2024/ReworkALL/chelsa_2006/Pmodel_Alienor_output.csv")
Gstar <- pressure*4.275*exp((37830/8.31447)*(1/298.15 - 1/(Tg_K)))/101325    #Gstar in Pa [Bernacchi et al. (2001)], KC, KO and ??* equal to 79430, 36380 and 37830 J/mol# standard values of KC, KO and ??* at 25?? equals to 39.37, 27480 and 4.22 Pa. O
Ca_Alienor<-model_Alienor_df$ca
Ci_Alienor<-model_Alienor_df$ci
m_Alienor <- (Ci_Alienor-Gstar)/(Ci_Alienor+2*Gstar)

Kc <- 40.49*exp((79430/8.31447)*(1/298.15 - 1/(Tg_K)))                      #Kc in Pa -> 1
Oc <-0.21*pressure #Oc in Pa -> month/daily
Ko <-27840*exp((36380/8.31447)*(1/298.15 - 1/(Tg_K)))                       #Ko in Pa -> month/daily
K  <- Kc*(1+(Oc/Ko)) 


Vcmax_Alienor <- 0.093*I*((Ci_Alienor+K)/(Ci_Alienor+2*Gstar))
Vcmax25_Alienor <- Vcmax_Alienor*exp((65330/8.31447)*((1/Tg_K)-(1/298.15)))
#####*sqrt(1-(0.41/m_Alienor)^(2/3))
model_own_df<-read.csv("/Users/echo/E/eechange project/analysis/work2024/ReworkALL/chelsa_2006/Pmodel_ownbeta_output.csv")
Gstar <- pressure*4.275*exp((37830/8.31447)*(1/298.15 - 1/(Tg_K)))/101325    #Gstar in Pa [Bernacchi et al. (2001)], KC, KO and ??* equal to 79430, 36380 and 37830 J/mol# standard values of KC, KO and ??* at 25?? equals to 39.37, 27480 and 4.22 Pa. O
Ca_own<-model_own_df$ca
Ci_own<-model_own_df$ci
m_own <- (Ci_own-Gstar)/(Ci_own+2*Gstar)

Kc <- 40.49*exp((79430/8.31447)*(1/298.15 - 1/(Tg_K)))                      #Kc in Pa -> 1
Oc <-0.21*pressure #Oc in Pa -> month/daily
Ko <-27840*exp((36380/8.31447)*(1/298.15 - 1/(Tg_K)))                       #Ko in Pa -> month/daily
K  <- Kc*(1+(Oc/Ko)) 

Vcmax_own <-  0.093*0.8*I*((Ci_own+K)/(Ci_own+2*Gstar))
Vcmax25_own <- Vcmax_own*exp((65330/8.31447)*((1/Tg_K)-(1/298.15)))
##*sqrt(1-(0.41/m_own)^(2/3))
vcmaxcombine<-data.frame(cbind(NECT_sample_climate$Site.ID.1,NECT_sample_climate$MAP..mm.,Vcmax_constant,Vcmax_Alienor,Vcmax_own,Vcmax25_constant,Vcmax25_Alienor,Vcmax25_own))
colnames(vcmaxcombine) <- c(
  "Site_ID", 
  "MAP_mm", 
  "vcmax_Constant", 
  "vcmax_Alienor", 
  "vcmax_Own",
  "vcmax25_Constant", 
  "vcmax25_Alienor", 
  "vcmax25_Own"
)

vcmaxG<-stack("/Users/echo/E/eechange project/analysis/work2024/vcmax/GOME2_VcmaxTg_05deg.tif")
plot(vcmaxG)##unit:µmol m−2 s−1 Vcmax from GOME-2 SIF

site<-cbind(NECT_sample_climate$Longitude,NECT_sample_climate$Latitude)
obsvcmax<-raster::extract(vcmaxG,site,method = 'bilinear')
vcmaxcombine1<-data.frame(cbind(NECT_sample_climate$Site.ID.1,NECT_sample_climate$MAP..mm.,Vcmax_constant,Vcmax_Alienor,Vcmax_own,obsvcmax,Vcmax25_constant,Vcmax25_Alienor,Vcmax25_own))
colnames(vcmaxcombine1) <- c(
  "Site_ID", 
  "MAP_mm", 
  "vcmax_Constant", 
  "vcmax_Alienor", 
  "vcmax_Own",
  "Obs_vcmax",
  "vcmax25_Constant", 
  "vcmax25_Alienor", 
  "vcmax25_Own"
)

vcmaxcombine1$Vcmax25_obs <- vcmaxcombine1$Obs_vcmax*exp((65330/8.31447)*((1/Tg_K)-(1/298.15)))
# Summarize only necessary columns
site_summary2 <- vcmaxcombine1 %>%
  group_by(Site_ID) %>%
  summarize(
    MAP_mm = mean(MAP_mm, na.rm = TRUE),  # Mean MAP for each site
    vcmax_Constant_Mean = mean(vcmax_Constant, na.rm = TRUE),
    vcmax_Alienor_Mean = mean(vcmax_Alienor, na.rm = TRUE),
    vcmax_Own_Mean = mean(vcmax_Own, na.rm = TRUE),
    Obs_vcmax = mean(Obs_vcmax, na.rm = TRUE),
    vcmax25_Constant_Mean = mean(vcmax25_Constant, na.rm = TRUE),
    vcmax25_Alienor_Mean = mean(vcmax25_Alienor, na.rm = TRUE),
    vcmax25_Own_Mean = mean(vcmax25_Own, na.rm = TRUE),
    vcmax25_obs_Mean = mean(Vcmax25_obs, na.rm = TRUE)
  )


# Plot for Vcmax_Own
vcmax12 <- ggplot(site_summary2, aes(x = MAP_mm)) +
  geom_point(aes(y = vcmax_Own_Mean, color = "vcmax_Own"), size = 1) +
  geom_smooth(aes(y = vcmax_Own_Mean, color = "vcmax_Own"), method = "lm", se = TRUE, size = 1) +
  geom_point(aes(y = Obs_vcmax, color = "Obs_vcmax")) +
  # Add smooth curves for other vcmax types
  geom_smooth(aes(y = Obs_vcmax, color = "Obs_vcmax"), method = "lm", se = TRUE, size = 1) +
  ylab(expression('Vcmax'*~'('*mu*'mol'*m^-2*s^-1*')')) +
  xlab(expression('MAP'*~'('*'mm'*')')) +
  theme_bw() +
  scale_color_manual(
    values = c(
      "Obs_vcmax" = "blue",
      "vcmax_Own" = "red"
    )
  ) +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text = element_text(size = rel(1.5)),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_line(color = "white")
  ) +
  annotate("text", x = 693, y = 40, label = "(a)", size = 6.5)
vcmax12
vcmax22 <- ggplot(site_summary2, aes(x = MAP_mm)) +
  geom_point(aes(y = vcmax25_Own_Mean, color = "vcmax25_Own_Mean")) +
  geom_smooth(aes(y = vcmax25_Own_Mean, color = "vcmax25_Own_Mean"), method = "lm", se = TRUE, size = 1) +
  geom_point(aes(y = vcmax25_obs_Mean, color = "vcmax25_obs_Mean")) +
  # Add smooth curves for other vcmax types
  geom_smooth(aes(y = vcmax25_obs_Mean, color = "vcmax25_obs_Mean"), method = "lm", se = TRUE, size = 1) +
  ylab(expression('Vcmax25'*~'('*mu*'mol'*m^-2*s^-1*')')) +
  xlab(expression('MAP'*~'('*'mm'*')')) +
  theme_bw() +
  scale_color_manual(
    values = c(
      "vcmax25_obs_Mean" = "blue",
      "vcmax25_Own_Mean" = "red"
    )
  ) +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text = element_text(size = rel(1.5)),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_line(color = "white")
  ) +
  annotate("text", x = 730, y = 200, label = "(b)", size = 6.5)
vcmax22
# Combine plots
plot_grid(vcmax12, vcmax22, label_x = 0.2, nrow = 1)
write.csv(site_summary2,"/Users/echo/E/eechange project/analysis/work2024/ReworkALL/chelsa_2006/vcmax_vcmax25_site_obs_Ning2017_0.8IL.csv")
#################################################################Fig2
# Annotate the plot
f1 <- ggplot(data=site_summary,aes(x=MAP_mm,y=Chi_Obs_Mean))+geom_point() + 
  geom_smooth(aes(y=Chi_Obs_Mean,x=MAP_mm),method = 'lm') +
  #geom_bar(stat="identity",position="dodge")+
  geom_errorbar(aes(ymax=Chi_Obs_Mean+Chi_Obs_SD,ymin=Chi_Obs_Mean - Chi_Obs_SD),position=position_dodge(0.9),width=1.0)+
  #scale_fill_brewer(palette = "Set1")+
  xlab(expression('MAP'*~'('*'mm'*')')) + ylab(expression(~chi)) +
  geom_point(aes(x=MAP_mm,y=Chi_Own_Mean),color="red") + 
  geom_smooth(aes(y=Chi_Own_Mean,x=MAP_mm, col='red'),method = 'lm') +
  theme_bw()+
  theme(legend.position="none",axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),
        axis.text=element_text(size=rel(1.5)),panel.grid.major=element_line(color="white"),panel.grid.minor=element_line(color="white")) +
  annotate("text", x=693, y=0.9, label="(a)",size=6.5)
f1

lmadata<-read.csv("/Users/echo/E/eechange project/analysis/work2024/LMA_wang2023/predicted_LMA.csv")

f2<- ggplot(data=lmadata,aes(x=map,y=lmaobsmean))+geom_point() + 
  geom_smooth(aes(y=lmaobsmean,x=map),method = 'lm') +
  geom_errorbar(aes(ymax=lmaobsmean+lma.sd,ymin=lmaobsmean-lma.sd),position=position_dodge(0.9),width=1.0)+
  xlab(expression('MAP'*~'('*'mm'*')')) + ylab(expression('lnLMA'*~'('*g*'/'*m^2*')')) +
  geom_point(aes(x=map,y=lnprelma),color="red") + 
  geom_smooth(aes(y=lnprelma,x=map, col='red'),method = 'lm') +
  theme_bw()+
  theme(legend.position="none",axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),
        axis.text=element_text(size=rel(1.5)),panel.grid.major=element_line(color="white"),panel.grid.minor=element_line(color="white"))+
  annotate("text", x=693, y=5, label="(c)",size=6.5)
f2

##################################Narea
vamaxdata<-vcmaxcombine1
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
# LMA         0.0191092  0.0009714  19.672  < 2e-16 ***
#   vcmax25_Own 0.0040104  0.0007641   5.249 2.83e-07 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.4282 on 313 degrees of freedom
# Multiple R-squared:  0.9313,	Adjusted R-squared:  0.9308 
# F-statistic:  2121 on 2 and 313 DF,  p-value: < 2.2e-16
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

datapre$predictednarea<-0.0191*datapre$x+0.0040*datapre$vcmax25_Own 
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
  annotate("text", x=693, y=1.5, label="(d)",size=6.5)
n4
write.csv(datapre,"/Users/echo/E/eechange project/analysis/work2024/LMA_wang2023/predicted_Narea_NewVcmax25.csv")

plot_grid(f1,f2,n4,label_x = 0.2,nrow = 2)

vcmaxcombine2<-read.csv("/Users/echo/E/eechange project/analysis/work2024/ReworkALL/chelsa_2006/vcmax_threemodel_Ning2017_obsvcmax_newoption.csv")
vcmaxcombine2$Vcmax25_obs <- vcmaxcombine2$Obs_vcmax*exp((65330/8.31447)*((1/Tg_K)-(1/298.15)))
site_summary2 <- vcmaxcombine2 %>%
  group_by(Site_ID) %>%
  summarize(
    MAP_mm = mean(MAP_mm, na.rm = TRUE),  # Mean MAP for each site
    vcmax_Constant_Mean = mean(vcmax_Constant, na.rm = TRUE),
    vcmax_Alienor_Mean = mean(vcmax_Alienor, na.rm = TRUE),
    vcmax_Own_Mean = mean(vcmax_Own, na.rm = TRUE),
    Obs_vcmax = mean(Obs_vcmax, na.rm = TRUE),
    vcmax25_Constant_Mean = mean(vcmax25_Constant, na.rm = TRUE),
    vcmax25_Alienor_Mean = mean(vcmax25_Alienor, na.rm = TRUE),
    vcmax25_Own_Mean = mean(vcmax25_Own, na.rm = TRUE),
    vcmax25_obs_Mean = mean(Vcmax25_obs, na.rm = TRUE)
  )


# Plot for Vcmax_Own
vcmax12 <- ggplot(site_summary2, aes(x = MAP_mm)) +
  geom_point(aes(y = vcmax_Own_Mean, color = "red"), size = 1) +
  geom_smooth(aes(y = vcmax_Own_Mean, color = "red"), method = "lm", se = TRUE, size = 1) +
  geom_point(aes(y = Obs_vcmax)) +
  # Add smooth curves for other vcmax types
  geom_smooth(aes(y = Obs_vcmax), method = "lm", se = TRUE, size = 1) +
  ylab(expression('Vcmax'*~'('*mu*'mol'*m^-2*s^-1*')')) +
  xlab(expression('MAP'*~'('*'mm'*')')) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text = element_text(size = rel(1.5)),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_line(color = "white")
  ) +
  annotate("text", x = 693, y = 40, label = "(b)", size = 6.5)
vcmax12

plot_grid(f1,vcmax12,f2,n4,label_x = 0.2,nrow = 2)

lmobvcmax<-lm(Obs_vcmax~MAP_mm,data=site_summary2)
summary(lmobvcmax)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 41.084895   2.598058  15.814  < 2e-16 ***
#   MAP_mm      -0.031172   0.006312  -4.939 2.56e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 5.435 on 31 degrees of freedom
# Multiple R-squared:  0.4403,	Adjusted R-squared:  0.4223 
# F-statistic: 24.39 on 1 and 31 DF,  p-value: 2.559e-05
lmorevcmax<-lm(vcmax_Own_Mean~MAP_mm,data=site_summary2)
summary(lmorevcmax)
# Call:
#   lm(formula = vcmax_Own_Mean ~ MAP_mm, data = site_summary2)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -1.9654 -0.6573  0.1295  0.6272  2.1976 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 25.866483   0.495761   52.17  < 2e-16 ***
#   MAP_mm      -0.018185   0.001204  -15.10 7.75e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 1.037 on 31 degrees of freedom
# Multiple R-squared:  0.8803,	Adjusted R-squared:  0.8764 
# F-statistic:   228 on 1 and 31 DF,  p-value: 7.755e-16
site_summary2<-read.csv("/Users/echo/E/eechange project/analysis/work2024/ReworkALL/chelsa_2006/vcmax_vcmax25_site_obs_Ning2017_0.8IL.csv")

site_summary2$MAP_m <- site_summary2$MAP_mm / 1000
lmobvcmax1<-lm(Obs_vcmax~MAP_m,data=site_summary2)
summary(lmobvcmax1)
lmorevcmax1<-lm(vcmax_Own_Mean~MAP_m,data=site_summary2)
summary(lmorevcmax1)
# summary(lmobvcmax1)
# 
# Call:
#   lm(formula = Obs_vcmax ~ MAP_m, data = site_summary2)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -10.571  -3.405  -1.983   2.962  10.623 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   41.085      2.598  15.814  < 2e-16 ***
#   MAP_m        -31.172      6.312  -4.939 2.56e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 5.435 on 31 degrees of freedom
# Multiple R-squared:  0.4403,	Adjusted R-squared:  0.4223 
# F-statistic: 24.39 on 1 and 31 DF,  p-value: 2.559e-05
# 
# > summary(lmorevcmax1)
# 
# Call:
#   lm(formula = vcmax_Own_Mean ~ MAP_m, data = site_summary2)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -1.9654 -0.6573  0.1295  0.6272  2.1976 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  25.8665     0.4958   52.17  < 2e-16 ***
#   MAP_m       -18.1847     1.2044  -15.10 7.75e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 1.037 on 31 degrees of freedom
# Multiple R-squared:  0.8803,	Adjusted R-squared:  0.8764 
# F-statistic:   228 on 1 and 31 DF,  p-value: 7.755e-16


lmobservedx<-lm(Chi_Obs_Mean~MAP_mm,data=site_summary)
summary(lmobservedx)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 4.143e-01  1.746e-02   23.72  < 2e-16 ***
#   MAP_mm      5.306e-04  4.242e-05   12.51 1.19e-13 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.03653 on 31 degrees of freedom
# Multiple R-squared:  0.8346,	Adjusted R-squared:  0.8293 
# F-statistic: 156.4 on 1 and 31 DF,  p-value: 1.194e-13
lmpredicedx<-lm(Chi_Own_Mean~MAP_mm,data=site_summary)
summary(lmpredicedx)

# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 4.225e-01  1.718e-02   24.59  < 2e-16 ***
#   MAP_mm      4.974e-04  4.174e-05   11.92 4.15e-13 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.03594 on 31 degrees of freedom
# Multiple R-squared:  0.8208,	Adjusted R-squared:  0.8151 
# F-statistic:   142 on 1 and 31 DF,  p-value: 4.15e-13
lmpredicedxconstant<-lm(Chi_Constant_Mean~MAP_mm,data=site_summary)
summary(lmpredicedxconstant)
coefs <- summary(lmpredicedxconstant)$coefficients
round(coefs, 2)
# Call:
#   lm(formula = Chi_Constant_Mean ~ MAP_mm, data = site_summary)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.023508 -0.011168  0.003346  0.011776  0.020095 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 5.641e-01  6.549e-03  86.135  < 2e-16 ***
#   MAP_mm      9.135e-05  1.591e-05   5.742 2.56e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.0137 on 31 degrees of freedom
# Multiple R-squared:  0.5154,	Adjusted R-squared:  0.4997 
# F-statistic: 32.97 on 1 and 31 DF,  p-value: 2.564e-06
lmpredicedxAlienor<-lm(Chi_Alienor_Mean~MAP_mm,data=site_summary)
summary(lmpredicedxAlienor)
# Call:
#   lm(formula = Chi_Alienor_Mean ~ MAP_mm, data = site_summary)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.036488 -0.006147 -0.000255  0.008035  0.021894 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 0.5384246  0.0057614   93.45  < 2e-16 ***
#   MAP_mm      0.0001686  0.0000140   12.04 3.16e-13 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.01205 on 31 degrees of freedom
# Multiple R-squared:  0.8239,	Adjusted R-squared:  0.8183 
# F-statistic: 145.1 on 1 and 31 DF,  p-value: 3.163e-13

site_summary$MAP_m <- site_summary$MAP_mm / 1000
model1<-lm(Chi_Own_Mean ~ MAP_m, data = site_summary)
summary(model1)
model2<-lm(Chi_Constant_Mean ~ MAP_m, data = site_summary)
summary(model2)
model3<-lm(Chi_Alienor_Mean ~ MAP_m, data = site_summary)
summary(model3)
summary(model1)
model4<-lm(Chi_Obs_Mean~MAP_m,data=site_summary)
summary(model4)
# summary(model1)
# Call:
#   lm(formula = Chi_Own_Mean ~ MAP_m, data = site_summary)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.061729 -0.023342 -0.003368  0.023490  0.093072 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.42246    0.01718   24.59  < 2e-16 ***
#   MAP_m        0.49740    0.04174   11.92 4.15e-13 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.03594 on 31 degrees of freedom
# Multiple R-squared:  0.8208,	Adjusted R-squared:  0.8151 
# F-statistic:   142 on 1 and 31 DF,  p-value: 4.15e-13
# 
# > summary(model2)
# 
# Call:
#   lm(formula = Chi_Constant_Mean ~ MAP_m, data = site_summary)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.023508 -0.011168  0.003346  0.011776  0.020095 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 0.564118   0.006549  86.135  < 2e-16 ***
#   MAP_m       0.091353   0.015911   5.742 2.56e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.0137 on 31 degrees of freedom
# Multiple R-squared:  0.5154,	Adjusted R-squared:  0.4997 
# F-statistic: 32.97 on 1 and 31 DF,  p-value: 2.564e-06
# 
# > summary(model3)
# 
# Call:
#   lm(formula = Chi_Alienor_Mean ~ MAP_m, data = site_summary)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.036488 -0.006147 -0.000255  0.008035  0.021894 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 0.538425   0.005761   93.45  < 2e-16 ***
#   MAP_m       0.168583   0.013997   12.04 3.16e-13 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.01205 on 31 degrees of freedom
# Multiple R-squared:  0.8239,	Adjusted R-squared:  0.8183 
# F-statistic: 145.1 on 1 and 31 DF,  p-value: 3.163e-13
# summary(model4)
# 
# Call:
#   lm(formula = Chi_Obs_Mean ~ MAP_m, data = site_summary)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.07276 -0.02245 -0.00003  0.02580  0.09026 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.41427    0.01746   23.72  < 2e-16 ***
#   MAP_m        0.53059    0.04242   12.51 1.19e-13 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.03653 on 31 degrees of freedom
# Multiple R-squared:  0.8346,	Adjusted R-squared:  0.8293 
# F-statistic: 156.4 on 1 and 31 DF,  p-value: 1.194e-13

lmlma<-lm(lmaobsmean~map,data=lmadata)
summary(lmlma)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  4.7916920  0.1158783  41.351  < 2e-16 ***
#   map         -0.0018425  0.0002788  -6.608 5.23e-07 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.2395 on 26 degrees of freedom
# Multiple R-squared:  0.6268,	Adjusted R-squared:  0.6124 
# F-statistic: 43.66 on 1 and 26 DF,  p-value: 5.23e-07
lmlmapre<-lm(lnprelma~map,data=lmadata)
summary(lmlmapre)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  4.8219993  0.0603147   79.95  < 2e-16 ***
#   map         -0.0019833  0.0001451  -13.66 2.24e-13 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1247 on 26 degrees of freedom
# Multiple R-squared:  0.8778,	Adjusted R-squared:  0.8731 
# F-statistic: 186.7 on 1 and 26 DF,  p-value: 2.235e-13
lmadata$MAP_m <- lmadata$map / 1000
lma1<-lm(lmaobsmean~MAP_m,data=lmadata)
lma2<-lm(lnprelma~MAP_m,data=lmadata)
summary(lma1)
summary(lma2)
# summary(lma1)

# Call:
#   lm(formula = lmaobsmean ~ MAP_m, data = lmadata)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.35964 -0.15678 -0.04025  0.09782  0.73542 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   4.7917     0.1159  41.351  < 2e-16 ***
#   MAP_m        -1.8425     0.2788  -6.608 5.23e-07 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.2395 on 26 degrees of freedom
# Multiple R-squared:  0.6268,	Adjusted R-squared:  0.6124 
# F-statistic: 43.66 on 1 and 26 DF,  p-value: 5.23e-07
# 
# > summary(lma2)
# 
# Call:
#   lm(formula = lnprelma ~ MAP_m, data = lmadata)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.182859 -0.114540 -0.007858  0.090828  0.222942 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  4.82200    0.06031   79.95  < 2e-16 ***
#   MAP_m       -1.98328    0.14514  -13.66 2.24e-13 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1247 on 26 degrees of freedom
# Multiple R-squared:  0.8778,	Adjusted R-squared:  0.8731 
# F-statistic: 186.7 on 1 and 26 DF,  p-value: 2.235e-13

lmobservednarea<-lm(lnNarea_gm2_mean.~MAP_mm,data=datapre)
summary(lmobservednarea)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  1.1522381  0.1029181  11.196 1.92e-11 ***
#   MAP_mm      -0.0021951  0.0002477  -8.863 2.46e-09 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.2127 on 26 degrees of freedom
# Multiple R-squared:  0.7513,	Adjusted R-squared:  0.7418 
# F-statistic: 78.56 on 1 and 26 DF,  p-value: 2.456e-09
lmpredictednarea<-lm(ln_prenaera~MAP_mm,data=datapre)
summary(lmpredictednarea)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  1.046463   0.043629   23.98  < 2e-16 ***
#   MAP_mm      -0.001748   0.000105  -16.65 2.19e-15 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.09019 on 26 degrees of freedom
# Multiple R-squared:  0.9143,	Adjusted R-squared:  0.911 
# F-statistic: 277.2 on 1 and 26 DF,  p-value: 2.186e-15
datapre<-read.csv("/Users/echo/E/eechange project/analysis/work2024/LMA_wang2023/predicted_Narea_NewVcmax25.csv")

datapre$MAP_m <- datapre$MAP_mm / 1000
naera1<-lm(lnNarea_gm2_mean.~MAP_m,data=datapre)
nara2<-lm(ln_prenaera~MAP_m,data=datapre)
summary(naera1)
summary(nara2)

# summary(naera1)
# 
# Call:
#   lm(formula = lnNarea_gm2_mean. ~ MAP_m, data = datapre)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -0.3991 -0.1565 -0.0509  0.1053  0.4211 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   1.1522     0.1029  11.196 1.92e-11 ***
#   MAP_m        -2.1951     0.2477  -8.863 2.46e-09 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.2127 on 26 degrees of freedom
# Multiple R-squared:  0.7513,	Adjusted R-squared:  0.7418 
# F-statistic: 78.56 on 1 and 26 DF,  p-value: 2.456e-09
# 
# > summary(nara2)
# 
# Call:
#   lm(formula = ln_prenaera ~ MAP_m, data = datapre)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.166943 -0.070182 -0.002851  0.046217  0.186559 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  1.04646    0.04363   23.98  < 2e-16 ***
#   MAP_m       -1.74808    0.10498  -16.65 2.19e-15 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.09019 on 26 degrees of freedom
# Multiple R-squared:  0.9143,	Adjusted R-squared:  0.911 
# F-statistic: 277.2 on 1 and 26 DF,  p-value: 2.186e-15
###########################################################Fig 3

##lai https://land.copernicus.eu/global/products/lai
##lai for 2006 1km
##the datasets is the different datasets are period data. there are three datasets for different time
joined_data<-read.csv("/Users/echo/E/eechange project/analysis/work2024/faparmax_LAI/Pre_Obs_fapar_LAI")
k = 0.5
joined_data$LAI_obs<-(-1/k)*log(1-joined_data$fapar_max)
site_mean <- joined_data %>%
  group_by(X1) %>%
  summarize(across(where(is.numeric), mean, na.rm = TRUE))

pb1 <- ggplot(data=site_mean,aes(x=X3,y=fapar_max))+geom_point() + 
  geom_smooth(aes(y=fapar_max,x=X3),method = 'lm') +
  xlab(expression('MAP'*~'('*'mm'*')')) + ylab(expression('fAPAR')) +
  geom_point(aes(x=X3,y=fapar_max_sim),color="red") + 
  geom_smooth(aes(y=fapar_max_sim,x=X3, col='red'),method = 'lm') +
  theme_bw()+
  theme(legend.position="none",axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),
        axis.text=element_text(size=rel(1.5)),panel.grid.major=element_line(color="white"),panel.grid.minor=element_line(color="white"))+
  annotate("text", x=693, y=1, label="(a)",size=6.5)
pb1


#######################################LAI

pb2 <- ggplot(data=site_mean,aes(x=X3,y=LAI_obs))+geom_point() + 
  geom_smooth(aes(y=LAI_obs,x=X3),method = 'lm') +
  xlab(expression('MAP'*~'('*'mm'*')')) + ylab(expression('LAI')) +
  geom_point(aes(x=X3,y=lai),color="red") + 
  geom_smooth(aes(y=lai,x=X3, col='red'),method = 'lm') +
  theme_bw()+
  theme(legend.position="none",axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),
        axis.text=element_text(size=rel(1.5)),panel.grid.major=element_line(color="white"),panel.grid.minor=element_line(color="white"))+
  annotate("text", x=693, y=4, label="(b)",size=6.5)

pb2

lmfaparobs<-lm(fapar_max~X3,data=site_mean)
summary(lmfaparobs)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -1.330e-02  2.761e-02  -0.482    0.633    
# X3           1.292e-03  6.707e-05  19.265   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.05775 on 31 degrees of freedom
# Multiple R-squared:  0.9229,	Adjusted R-squared:  0.9204 
# F-statistic: 371.1 on 1 and 31 DF,  p-value: < 2.2e-16
lmfaparpre<-lm(fapar_max_sim~X3,data=site_mean)
summary(lmfaparpre)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 1.258e-01  3.546e-02   3.548  0.00126 ** 
#   X3          9.908e-04  8.614e-05  11.502 1.03e-12 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.07417 on 31 degrees of freedom
# Multiple R-squared:  0.8102,	Adjusted R-squared:  0.804 
# F-statistic: 132.3 on 1 and 31 DF,  p-value: 1.025e-12
lmlaiobs<-lm(LAI_obs~X3,data=site_mean)
summary(lmlaiobs)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -0.6968231  0.1507839  -4.621 6.34e-05 ***
#   X3           0.0057302  0.0003663  15.643 2.93e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.3154 on 31 degrees of freedom
# Multiple R-squared:  0.8876,	Adjusted R-squared:  0.8839 
# F-statistic: 244.7 on 1 and 31 DF,  p-value: 2.926e-16

lmlai<-lm(lai~X3,data=site_mean)
summary(lmlai)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -0.1338273  0.1505552  -0.889    0.381    
# X3           0.0044013  0.0003658  12.033 3.24e-13 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.3149 on 31 degrees of freedom
# Multiple R-squared:  0.8237,	Adjusted R-squared:  0.818 
# F-statistic: 144.8 on 1 and 31 DF,  p-value: 3.238e-13
###################################seasonal fapar and LAI--Boya version
site_mean<-read.csv("/Users/echo/E/eechange project/analysis/work2024/faparmax_LAI/Pre_Obs_fapar_LAI_Boyaseasonal.csv")

pb1 <- ggplot(data=site_mean,aes(x=MAP,y=fapar_obs))+geom_point() + 
  geom_smooth(aes(y=fapar_obs,x=MAP),method = 'lm') +
  xlab(expression('MAP'*~'('*'mm'*')')) + ylab(expression('fAPAR')) +
  geom_point(aes(x=MAP,y=fapar_pre),color="red") + 
  geom_smooth(aes(y=fapar_pre,x=MAP, col='red'),method = 'lm') +
  theme_bw()+
  theme(legend.position="none",axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),
        axis.text=element_text(size=rel(1.5)),panel.grid.major=element_line(color="white"),panel.grid.minor=element_line(color="white"))+
  annotate("text", x=693, y=0.6, label="(a)",size=6.5)
pb1

pb1 <- ggplot(data=site_mean,aes(x=MAP,y=fapar_obs))+geom_point(size = 2.8) + 
  geom_smooth(aes(y=fapar_obs,x=MAP),method = 'lm') +
  xlab(expression('MAP'*~'('*'mm'*')')) + ylab(expression('fAPAR')) +
  geom_point(aes(x=MAP,y=fapar_pre),color="red",size = 2.8) + 
  geom_smooth(aes(y=fapar_pre,x=MAP, col='red'),method = 'lm') +
  theme_bw()+
  theme(legend.position="none",axis.title.x=element_blank(),axis.title.y=element_text(size=29),
        axis.text=element_text(size=rel(2.5)),panel.grid.major=element_line(color="white"),panel.grid.minor=element_line(color="white"))
pb1
#######################################LAI

pb2 <- ggplot(data=site_mean,aes(x=MAP,y=LAI_obs))+geom_point() + 
  geom_smooth(aes(y=LAI_obs,x=MAP),method = 'lm') +
  xlab(expression('MAP'*~'('*'mm'*')')) + ylab(expression('LAI')) +
  geom_point(aes(x=MAP,y=LAI_pre),color="red") + 
  geom_smooth(aes(y=LAI_pre,x=MAP, col='red'),method = 'lm') +
  theme_bw()+
  theme(legend.position="none",axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),
        axis.text=element_text(size=rel(1.5)),panel.grid.major=element_line(color="white"),panel.grid.minor=element_line(color="white"))+
  annotate("text", x=693, y=1.5, label="(b)",size=6.5)

pb2

pb2 <- ggplot(data=site_mean,aes(x=MAP,y=LAI_obs))+geom_point(size = 2.8) + 
  geom_smooth(aes(y=LAI_obs,x=MAP),method = 'lm') +
  xlab(expression('MAP'*~'('*'mm'*')')) + ylab(expression('LAI')) +
  geom_point(aes(x=MAP,y=LAI_pre),color="red",,size = 2.8) + 
  geom_smooth(aes(y=LAI_pre,x=MAP, col='red'),method = 'lm') +scale_y_continuous(breaks = c(0.3, 0.6, 0.9, 1.2))+
  theme_bw()+
  theme(legend.position="none",axis.title.x=element_blank(),axis.title.y=element_text(size=29),
        axis.text=element_text(size=rel(2.5)),panel.grid.major=element_line(color="white"),panel.grid.minor=element_line(color="white"))

pb2
lmfaparobs<-lm(fapar_obs~MAP,data=site_mean)
summary(lmfaparobs)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 1.168e-02  1.867e-02   0.625    0.536    
# MAP         6.284e-04  4.537e-05  13.852 8.02e-15 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.03906 on 31 degrees of freedom
# Multiple R-squared:  0.8609,	Adjusted R-squared:  0.8564 
# F-statistic: 191.9 on 1 and 31 DF,  p-value: 8.023e-15
lmfaparpre<-lm(fapar_pre~MAP,data=site_mean)
summary(lmfaparpre)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -3.400e-03  1.671e-02  -0.203     0.84    
# MAP          6.808e-04  4.059e-05  16.773   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.03495 on 31 degrees of freedom
# Multiple R-squared:  0.9007,	Adjusted R-squared:  0.8975 
# F-statistic: 281.3 on 1 and 31 DF,  p-value: < 2.2e-16
lmlaiobs<-lm(LAI_obs~MAP,data=site_mean)
summary(lmlaiobs)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -0.072948   0.058849   -1.24    0.224    
# MAP          0.001762   0.000143   12.32 1.76e-13 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1231 on 31 degrees of freedom
# Multiple R-squared:  0.8304,	Adjusted R-squared:  0.825 
# F-statistic: 151.8 on 1 and 31 DF,  p-value: 1.761e-13
lmlai<-lm(LAI_pre~MAP,data=site_mean)
summary(lmlai)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -0.1026431  0.0528744  -1.941   0.0614 .  
# MAP          0.0018792  0.0001285  14.629 1.84e-15 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1106 on 31 degrees of freedom
# Multiple R-squared:  0.8735,	Adjusted R-squared:  0.8694 
# F-statistic:   214 on 1 and 31 DF,  p-value: 1.835e-15
site_mean$MAP_m <- site_mean$MAP / 1000
lmfaparobs1<-lm(fapar_obs~MAP_m,data=site_mean)
summary(lmfaparobs1)
lmfaparpre1<-lm(fapar_pre~MAP_m,data=site_mean)
summary(lmfaparpre1)
# summary(lmfaparobs1)
# 
# Call:
#   lm(formula = fapar_obs ~ MAP_m, data = site_mean)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.108768 -0.022559  0.003045  0.016700  0.075142 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.01168    0.01867   0.625    0.536    
# MAP_m        0.62840    0.04537  13.852 8.02e-15 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.03906 on 31 degrees of freedom
# Multiple R-squared:  0.8609,	Adjusted R-squared:  0.8564 
# F-statistic: 191.9 on 1 and 31 DF,  p-value: 8.023e-15
# 
# > summary(lmfaparpre1)
# 
# Call:
#   lm(formula = fapar_pre ~ MAP_m, data = site_mean)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.086425 -0.007835 -0.003051  0.008960  0.093524 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -0.00340    0.01671  -0.203     0.84    
# MAP_m        0.68079    0.04059  16.773   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.03495 on 31 degrees of freedom
# Multiple R-squared:  0.9007,	Adjusted R-squared:  0.8975 
# F-statistic: 281.3 on 1 and 31 DF,  p-value: < 2.2e-16
lmlaiobs1<-lm(LAI_obs~MAP_m,data=site_mean)
summary(lmlaiobs1)
lmlai2<-lm(LAI_pre~MAP_m,data=site_mean)
summary(lmlai2)
# summary(lmlaiobs1)
# 
# Call:
#   lm(formula = LAI_obs ~ MAP_m, data = site_mean)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.32707 -0.07866  0.00298  0.05460  0.24597 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -0.07295    0.05885   -1.24    0.224    
# MAP_m        1.76149    0.14297   12.32 1.76e-13 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1231 on 31 degrees of freedom
# Multiple R-squared:  0.8304,	Adjusted R-squared:  0.825 
# F-statistic: 151.8 on 1 and 31 DF,  p-value: 1.761e-13
# 
# > summary(lmlai2)
# 
# Call:
#   lm(formula = LAI_pre ~ MAP_m, data = site_mean)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.25442 -0.03649  0.00115  0.01686  0.33158 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -0.10264    0.05287  -1.941   0.0614 .  
# MAP_m        1.87921    0.12845  14.629 1.84e-15 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1106 on 31 degrees of freedom
# Multiple R-squared:  0.8735,	Adjusted R-squared:  0.8694 
# F-statistic:   214 on 1 and 31 DF,  p-value: 1.835e-15
###################plot GPP, The new GPP dataset, FLUXCOM-X (X-BASE), shows significant improvement over MTE.
model_season_df<-read.csv("/Users/echo/E/eechange project/analysis/work2024/GPP/gppobs_gCm2y_NECT_pmodel_FLUX_X.csv")
joined_data<-read.csv("/Users/echo/E/eechange project/analysis/work2024/faparmax_LAI/Pre_Obs_fapar_LAI")
model_season_df$site<-joined_data$X1
model_season_df$Precipitation<-joined_data$X3
model_season_df$Tg<-climatedata$tg_chelsa_2006
tmean<-read.csv("/Users/echo/E/eechange project/analysis/work2024/GPP/tmean_C_NECT_1987_2006.csv")
model_season_df$Tmean<-rowMeans(tmean[, 2:13], na.rm = TRUE)
gpp_mean <- model_season_df %>%
  group_by(site) %>%
  summarize(across(where(is.numeric), mean, na.rm = TRUE))


pb3 <- ggplot(data=gpp_mean,aes(x=Precipitation,y=gpp))+geom_point(color="red") + 
  geom_smooth(method = 'lm', formula = y ~ x, aes(color='red')) +
  xlab(expression('MAP'*~'('*'mm'*')')) + ylab(expression('GPP'*~'('*'gC'*~m^-1*d^-1*')')) +
  geom_point(aes(x=Precipitation,y=gppobs_FLUX_X)) + 
  geom_smooth(aes(x=Precipitation,y=gppobs_FLUX_X),method = 'lm') +
  theme_bw()+
  theme(legend.position="none",axis.title.x=element_text(size=15),axis.title.y=element_text(size=15),
        axis.text=element_text(size=rel(1.5)),panel.grid.major=element_line(color="white"),panel.grid.minor=element_line(color="white"))+
  annotate("text", x=655, y=5.5, label="(c)",size=6.5)
pb3

p4 <- ggplot(data=gpp_mean,aes(x=Precipitation,y=gpp*0.46))+geom_point(color="red") + 
  geom_smooth(method = 'lm', formula = y ~ x, aes(color='red')) +
  xlab(expression('MAP'*~'('*'mm'*')')) +  ylab(expression('NPP'*~'('*'gC'*~m^-1*d^-1*')'))  +
  geom_point(aes(x=Precipitation,y=gppobs_FLUX_X*0.46)) + 
  geom_smooth(aes(x=Precipitation,y=gppobs_FLUX_X*0.46),method = 'lm') +
  theme_bw()+
  theme(legend.position="none",axis.title.x=element_text(size=15),axis.title.y=element_text(size=15),
        axis.text=element_text(size=rel(1.5)),panel.grid.major=element_line(color="white"),panel.grid.minor=element_line(color="white"))+
  annotate("text", x=693, y=3, label="(d)",size=6.5)
p4
lmgpp<-lm(gpp~Precipitation,data=gpp_mean)
summary(lmgpp) 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   0.0270677  0.1462310   0.185    0.854    
# Precipitation 0.0055104  0.0003553  15.511  3.7e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.3059 on 31 degrees of freedom
# Multiple R-squared:  0.8859,	Adjusted R-squared:  0.8822 
# F-statistic: 240.6 on 1 and 31 DF,  p-value: 3.696e-16
lmgppobx<-lm(gppobs_FLUX_X~Precipitation,data=gpp_mean)
summary(lmgppobx) 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   -1.5513810  0.3638345  -4.264 0.000175 ***
#   Precipitation  0.0082369  0.0008839   9.319 1.68e-10 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.7611 on 31 degrees of freedom
# Multiple R-squared:  0.7369,	Adjusted R-squared:  0.7284 
# F-statistic: 86.84 on 1 and 31 DF,  p-value: 1.681e-10
lmnpp<-lm((gpp*0.46)~Precipitation,data=gpp_mean)
summary(lmnpp)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   0.0124512  0.0672663   0.185    0.854    
# Precipitation 0.0025348  0.0001634  15.511  3.7e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1407 on 31 degrees of freedom
# Multiple R-squared:  0.8859,	Adjusted R-squared:  0.8822 
# F-statistic: 240.6 on 1 and 31 DF,  p-value: 3.696e-16
lmnppobs<-lm((gppobs_FLUX_X*0.46)~Precipitation,data=gpp_mean)
summary(lmnppobs)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   -0.7136352  0.1673639  -4.264 0.000175 ***
#   Precipitation  0.0037890  0.0004066   9.319 1.68e-10 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.3501 on 31 degrees of freedom
# Multiple R-squared:  0.7369,	Adjusted R-squared:  0.7284 
# F-statistic: 86.84 on 1 and 31 DF,  p-value: 1.681e-10
wholedata<-read.csv("/Users/echo/E/eechange project/analysis/work2024/ReworkALL/chelsa_1987_2006/Climate_R_S_NECT.csv")
pb4<-ggplot(wholedata,aes(x=MAP..mm., y=log(R_S))) + geom_point() + 
  geom_smooth(method = 'lm', formula = y ~ x)+
  xlab(expression('MAP'*~'('*'mm'*')')) + ylab(expression('ln R:S'))+
  geom_point(aes(x=MAP..mm.,y=log(preR_S)),color="red") + 
  geom_smooth(aes(y=log(preR_S),x=MAP..mm., col='red'),method = 'lm') +
  theme_bw()+
  theme(legend.position="none",axis.title.x=element_text(size=15),axis.title.y=element_text(size=15),
        axis.text=element_text(size=rel(1.5)),panel.grid.major=element_line(color="white"),panel.grid.minor=element_line(color="white"))+
  annotate("text", x=655, y=0.6, label="(e)",size=6.5)
pb4

wholedata$sitegroup <- ifelse(wholedata$BGB.tree_gC.m.2 > 0, "treesite", "grasssite")

pb4group<-ggplot(wholedata,aes(x=MAP..mm., y=log(R_S))) + geom_point() + 
  #geom_smooth(method = 'lm', formula = y ~ x)+
  xlab(expression('MAP'*~'('*'mm'*')')) + ylab(expression('ln R:S'))+
  geom_point(aes(x=MAP..mm.,y=log(preR_S)),color="red") + 
  geom_smooth(aes(y = log(R_S), group = sitegroup), method = 'lm',color = "blue", formula = y ~ x, se = TRUE) +
  
  # Prediction regression line (red)
  geom_smooth(aes(y = log(preR_S), group = sitegroup), method = 'lm',color = "red", formula = y ~ x, se = TRUE) +
  
  #geom_smooth(aes(y=log(preR_S),x=MAP..mm., col='red'),method = 'lm') +
  theme_bw()+
  theme(legend.position="none",axis.title.x=element_text(size=15),axis.title.y=element_text(size=15),
        axis.text=element_text(size=rel(1.5)),panel.grid.major=element_line(color="white"),panel.grid.minor=element_line(color="white"))+
  annotate("text", x=655, y=0.6, label="(e)",size=6.5)
pb4group

wholedata<-read.csv("/Users/echo/E/eechange project/analysis/work2024/ReworkALL/chelsa_1987_2006/Climate_R_S_NECT.csv")
wholedata$preRMF<-wholedata$preR_S/(wholedata$preR_S+1)
wholedata$prebelow<-wholedata$preRMF*wholedata$total.biomass_gC.m.2

pbBB<-ggplot(wholedata,aes(x=MAP..mm., y=log(BGB_gC.m.2))) + geom_point() + 
  geom_smooth(aes(x=MAP..mm., y=log(BGB_gC.m.2)),method = 'lm', formula = y ~ x)+
  xlab(expression('MAP'*~'('*'mm'*')')) + ylab(expression('ln BGB'*~'('*'gC'*~m^-2*')'))+
  geom_point(aes(x=MAP..mm.,y=log(prebelow)),color="red") + 
  geom_smooth(aes(y=log(prebelow),x=MAP..mm., col='red'),method = 'lm') +
  theme_bw()+
  theme(legend.position="none",axis.title.x=element_text(size=15),axis.title.y=element_text(size=15),
        axis.text=element_text(size=rel(1.5)),panel.grid.major=element_line(color="white"),panel.grid.minor=element_line(color="white"))+
  annotate("text", x=655, y=12, label="(f)",size=6.5)
pbBB

RSimprove2<-read.csv("/Users/echo/E/eechange project/analysis/work2024/ReworkALL/chelsa_1987_2006/Climate_R_S_NECT_predictRS_biomass.csv")

pbRSgroup <- ggplot(RSimprove2, aes(x = MAP..mm.)) +
  # Points: Black for observations, red for predictions
  geom_point(aes(y = log(R_S)), color = "black") +
  geom_point(aes(y = pred_lnRSnew), color = "red") +
  
  # Overall regression line (black)
  #geom_smooth(aes(y = log(R_S)), method = 'lm', formula = y ~ x, color = "blue", se = TRUE) +
  
  # Regression lines by sitegroup
  geom_smooth(aes(y = log(R_S), group = sitegroup), method = 'lm',color = "blue", formula = y ~ x, se = TRUE) +
  
  # Prediction regression line (red)
  #geom_smooth(aes(y = pred_lnRSnew), method = 'lm', formula = y ~ x, color = "red", se = TRUE) +
  geom_smooth(aes(y = pred_lnRSnew, group = sitegroup), method = 'lm',color = "red", formula = y ~ x, se = TRUE) +
  
  # Labels and theme
  xlab(expression('MAP'*~'('*'mm'*')')) +
  ylab(expression('lnR:S')) +
  theme_bw() +
  theme(
    legend.position = "right",  # Keep legend for sitegroup
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text = element_text(size = rel(1.5)),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_line(color = "white")
  ) +
  annotate("text", x = 655, y = 1, label = "(e)", size = 6.5)

pbRSgroup

pbRSgroup <- ggplot(RSimprove2, aes(x = MAP..mm.)) +
  geom_point(aes(y = log(R_S)), color = "black",size = 2.8) +
  geom_point(aes(y = pred_lnRSnew), color = "red",size = 2.8) +
  geom_smooth(aes(y = log(R_S), group = sitegroup), method = 'lm',color = "blue", formula = y ~ x, se = TRUE) +
  geom_smooth(aes(y = pred_lnRSnew, group = sitegroup), method = 'lm',color = "red", formula = y ~ x, se = TRUE) +
  xlab(expression('MAP'*~'('*'mm'*')')) +
  ylab(expression('lnR:S')) +scale_y_continuous(breaks = c(-7, -5, -3, -1))+
  theme_bw() +
  theme(
    legend.position = "right",  # Keep legend for sitegroup
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 29),
    axis.text = element_text(size = rel(2.5)),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_line(color = "white")
  ) 
pbRSgroup

pbAGBgroup <- ggplot(RSimprove2, aes(x = MAP..mm.)) +
  # Points: Black for observations, red for predictions
  geom_point(aes(y = log(AGB_gC.m.2)), color = "black") +
  geom_point(aes(y = log(preAGB_new)), color = "red") +
  
  # Overall regression line (black)
  #geom_smooth(aes(y = log(BGB_gC.m.2)), method = 'lm', formula = y ~ x, color = "blue", se = TRUE) +
  
  # Regression lines by sitegroup
  geom_smooth(aes(y = log(AGB_gC.m.2), group = sitegroup), method = 'lm',color = "blue", formula = y ~ x, se = TRUE) +
  
  # Prediction regression line (red)
  #geom_smooth(aes(y = log(prebelow)), method = 'lm', formula = y ~ x, color = "red", se = TRUE) +
  geom_smooth(aes(y = log(preAGB_new), group = sitegroup), method = 'lm',color = "red", formula = y ~ x, se = TRUE) +
  
  # Labels and theme
  xlab(expression('MAP'*~'('*'mm'*')')) +
  ylab(expression('lnAGB'*~'('*'gC'*~m^-2*')')) +
  theme_bw() +
  theme(
    legend.position = "right",  # Keep legend for sitegroup
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text = element_text(size = rel(1.5)),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_line(color = "white")
  ) +
  annotate("text", x = 655, y = 11, label = "(f)", size = 6.5)

pbAGBgroup

pbAGBgroup <- ggplot(RSimprove2, aes(x = MAP..mm.)) +
  geom_point(aes(y = log(AGB_gC.m.2)), color = "black",size = 2.8) +
  geom_point(aes(y = log(preAGB_new)), color = "red",size = 2.8) +
  geom_smooth(aes(y = log(AGB_gC.m.2), group = sitegroup), method = 'lm',color = "blue", formula = y ~ x, se = TRUE) +
  geom_smooth(aes(y = log(preAGB_new), group = sitegroup), method = 'lm',color = "red", formula = y ~ x, se = TRUE) +
  
  # Labels and theme
  xlab(expression('MAP'*~'('*'mm'*')')) +
  ylab(expression('lnAGB'*~'('*'gC'*~m^-2*')')) +
  theme_bw() +
  theme(
    legend.position = "right",  # Keep legend for sitegroup
    axis.title.x =element_blank(),
    axis.title.y = element_text(size = 29),
    axis.text = element_text(size = rel(2.5)),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_line(color = "white")
  )
pbAGBgroup

pbBBgroup <- ggplot(RSimprove2, aes(x = MAP..mm.)) +
  # Points: Black for observations, red for predictions
  geom_point(aes(y = log(BGB_gC.m.2)), color = "black") +
  geom_point(aes(y = log(prebelow_new)), color = "red") +
  
  # Overall regression line (black)
  #geom_smooth(aes(y = log(BGB_gC.m.2)), method = 'lm', formula = y ~ x, color = "blue", se = TRUE) +
  
  # Regression lines by sitegroup
  geom_smooth(aes(y = log(BGB_gC.m.2), group = sitegroup), method = 'lm',color = "blue", formula = y ~ x, se = TRUE) +
  
  # Prediction regression line (red)
  #geom_smooth(aes(y = log(prebelow)), method = 'lm', formula = y ~ x, color = "red", se = TRUE) +
  geom_smooth(aes(y = log(prebelow_new), group = sitegroup), method = 'lm',color = "red", formula = y ~ x, se = TRUE) +
  
  # Labels and theme
  xlab(expression('MAP'*~'('*'mm'*')')) +
  ylab(expression('lnBGB'*~'('*'gC'*~m^-2*')')) +
  theme_bw() +
  theme(
    legend.position = "right",  # Keep legend for sitegroup
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text = element_text(size = rel(1.5)),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_line(color = "white")
  ) +
  annotate("text", x = 655, y = 10, label = "(g)", size = 6.5)

pbBBgroup

pbBBgroup <- ggplot(RSimprove2, aes(x = MAP..mm.)) +
  # Points: Black for observations, red for predictions
  geom_point(aes(y = log(BGB_gC.m.2)), color = "black",size = 2.8) +
  geom_point(aes(y = log(prebelow_new)), color = "red",size = 2.8) +
  geom_smooth(aes(y = log(BGB_gC.m.2), group = sitegroup), method = 'lm',color = "blue", formula = y ~ x, se = TRUE) +
  geom_smooth(aes(y = log(prebelow_new), group = sitegroup), method = 'lm',color = "red", formula = y ~ x, se = TRUE) +
  xlab(expression('MAP'*~'('*'mm'*')')) +
  ylab(expression('lnBGB'*~'('*'gC'*~m^-2*')')) +scale_y_continuous(breaks = c(0, 2, 4, 6, 8))+
  theme_bw() +
  theme(
    legend.position = "right",  # Keep legend for sitegroup
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 29),
    axis.text = element_text(size = rel(2.5)),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_line(color = "white")
  ) 
pbBBgroup

multifigure<-plot_grid(pb1,pb2,pb3,p4,label_x = 0.2,ncol=4,nrow=1)
multifigure1<-plot_grid(pbRSgroup,pbAGBgroup,pbBBgroup,label_x = 0.2,ncol=4,nrow=1)
plot_grid(multifigure,multifigure1,label_x = 0.2,nrow=2)

plot_grid(pb1,pb2,pb3,p4,pbRSgroup,pbAGBgroup,pbBBgroup,label_x = 0.2,ncol=4,nrow=2)
# woodbiomass<- RSimprove2 %>% 
#   filter(!is.na(mean_pretreebiomass_gCm2))  # Removes NA and NaN
# pbtotolbiomass<-ggplot(woodbiomass,aes(x=MAP..mm., y=log(total.biomass_gC.m.2))) + geom_point() + 
#   geom_smooth(aes(x=MAP..mm., y=log(total.biomass_gC.m.2)),method = 'lm', formula = y ~ x)+
#   xlab(expression('MAP'*~'('*'mm'*')')) + ylab(expression('lnTotal biomass'*~'('*'gC'*~m^-2*')'))+
#   geom_point(aes(x=MAP..mm.,y=log(pre_totalbiomass)),color="red") + 
#   geom_smooth(aes(y=log(pre_totalbiomass),x=MAP..mm., col='red'),method = 'lm') +
#   theme_bw()+
#   theme(legend.position="none",axis.title.x=element_text(size=15),axis.title.y=element_text(size=15),
#         axis.text=element_text(size=rel(1.5)),panel.grid.major=element_line(color="white"),panel.grid.minor=element_line(color="white"))
# pbtotolbiomass
lowBGB <- RSimprove2 %>% 
  filter(MAP..mm. < 420, !is.na(MAP..mm.))
BGBmodelobs_low<-lm(log(BGB_gC.m.2)~MAP..mm.,data=lowBGB)
summary(BGBmodelobs_low)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -0.3863034  0.1358089  -2.844  0.00464 ** 
#   MAP..mm.     0.0017523  0.0004062   4.314 1.95e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.8071 on 477 degrees of freedom
# Multiple R-squared:  0.03755,	Adjusted R-squared:  0.03553 
# F-statistic: 18.61 on 1 and 477 DF,  p-value: 1.95e-05
BGBmodelpre_low<-lm(log(prebelow_new)~MAP..mm.,data=lowBGB)
summary(BGBmodelpre_low)
# Call:
#   lm(formula = log(prebelow_new) ~ MAP..mm., data = lowBGB)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -2.2494 -0.2096  0.1895  0.3164  1.2150 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -0.5443943  0.1229199  -4.429 1.18e-05 ***
#   MAP..mm.     0.0023516  0.0003676   6.397 3.79e-10 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.7305 on 477 degrees of freedom
# Multiple R-squared:  0.079,	Adjusted R-squared:  0.07707 
# F-statistic: 40.92 on 1 and 477 DF,  p-value: 3.795e-10
AGBmodelobs_low<-lm(log(AGB_gC.m.2)~MAP..mm.,data=lowBGB)
summary(AGBmodelobs_low)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 3.4356138  0.0729776   47.08   <2e-16 ***
#   MAP..mm.    0.0057769  0.0002183   26.47   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.4337 on 477 degrees of freedom
# Multiple R-squared:  0.5949,	Adjusted R-squared:  0.5941 
# F-statistic: 700.5 on 1 and 477 DF,  p-value: < 2.2e-16
AGBmodelpre_low<-lm(log(preAGB_new)~MAP..mm.,data=lowBGB)
summary(AGBmodelpre_low)
# Call:
#   lm(formula = log(preAGB_new) ~ MAP..mm., data = lowBGB)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -1.01418 -0.26282  0.01165  0.25270  0.73771 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 3.4333873  0.0726543   47.26   <2e-16 ***
#   MAP..mm.    0.0057863  0.0002173   26.63   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.4318 on 477 degrees of freedom
# Multiple R-squared:  0.5978,	Adjusted R-squared:  0.597 
# F-statistic: 709.1 on 1 and 477 DF,  p-value: < 2.2e-16
RSmodelobs_low<-lm(log(R_S)~MAP..mm.,data=lowBGB)
summary(RSmodelobs_low)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -3.8219172  0.1546174 -24.719   <2e-16 ***
#   MAP..mm.    -0.0040246  0.0004624  -8.703   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.9189 on 477 degrees of freedom
# Multiple R-squared:  0.137,	Adjusted R-squared:  0.1352 
# F-statistic: 75.74 on 1 and 477 DF,  p-value: < 2.2e-16
RSmodelobs_lowpre<-lm(pred_lnRSnew~MAP..mm.,data=lowBGB)
summary(RSmodelobs_lowpre)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -3.9777816  0.1155047 -34.438   <2e-16 ***
#   MAP..mm.    -0.0034346  0.0003455  -9.942   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.6864 on 477 degrees of freedom
# Multiple R-squared:  0.1717,	Adjusted R-squared:  0.1699 
# F-statistic: 98.85 on 1 and 477 DF,  p-value: < 2.2e-16
lowBGB$MAP_m <- lowBGB$MAP..mm. / 1000
BGBmodelobs_low1<-lm(log(BGB_gC.m.2)~MAP_m,data=lowBGB)
summary(BGBmodelobs_low1)
BGBmodelpre_low1<-lm(log(prebelow_new)~MAP_m,data=lowBGB)
summary(BGBmodelpre_low1)
AGBmodelobs_low1<-lm(log(AGB_gC.m.2)~MAP_m,data=lowBGB)
summary(AGBmodelobs_low1)
AGBmodelpre_low1<-lm(log(preAGB_new)~MAP_m,data=lowBGB)
summary(AGBmodelpre_low1)
RSmodelobs_low1<-lm(log(R_S)~MAP_m,data=lowBGB)
summary(RSmodelobs_low1)
RSmodelobs_lowpre1<-lm(pred_lnRSnew~MAP_m,data=lowBGB)
summary(RSmodelobs_lowpre1)

# summary(BGBmodelobs_low1)
# 
# Call:
#   lm(formula = log(BGB_gC.m.2) ~ MAP_m, data = lowBGB)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -1.6775 -0.4440  0.1467  0.6779  1.2008 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  -0.3863     0.1358  -2.844  0.00464 ** 
#   MAP_m         1.7523     0.4062   4.314 1.95e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.8071 on 477 degrees of freedom
# Multiple R-squared:  0.03755,	Adjusted R-squared:  0.03553 
# F-statistic: 18.61 on 1 and 477 DF,  p-value: 1.95e-05
# 
# > summary(BGBmodelpre_low1)
# 
# Call:
#   lm(formula = log(prebelow_new) ~ MAP_m, data = lowBGB)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -2.2494 -0.2096  0.1895  0.3164  1.2150 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  -0.5444     0.1229  -4.429 1.18e-05 ***
#   MAP_m         2.3516     0.3676   6.397 3.79e-10 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.7305 on 477 degrees of freedom
# Multiple R-squared:  0.079,	Adjusted R-squared:  0.07707 
# F-statistic: 40.92 on 1 and 477 DF,  p-value: 3.795e-10
# 
# > summary(AGBmodelobs_low1)
# 
# Call:
#   lm(formula = log(AGB_gC.m.2) ~ MAP_m, data = lowBGB)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -1.01115 -0.27101  0.01313  0.25473  0.74019 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  3.43561    0.07298   47.08   <2e-16 ***
#   MAP_m        5.77692    0.21827   26.47   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.4337 on 477 degrees of freedom
# Multiple R-squared:  0.5949,	Adjusted R-squared:  0.5941 
# F-statistic: 700.5 on 1 and 477 DF,  p-value: < 2.2e-16
# 
# > summary(AGBmodelpre_low1)
# # 
# Call:
#   lm(formula = log(preAGB_new) ~ MAP_m, data = lowBGB)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -1.01418 -0.26282  0.01165  0.25270  0.73771 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  3.43339    0.07265   47.26   <2e-16 ***
#   MAP_m        5.78626    0.21730   26.63   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.4318 on 477 degrees of freedom
# Multiple R-squared:  0.5978,	Adjusted R-squared:  0.597 
# F-statistic: 709.1 on 1 and 477 DF,  p-value: < 2.2e-16
# > summary(RSmodelobs_low1)
# 
# Call:
#   lm(formula = log(R_S) ~ MAP_m, data = lowBGB)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -2.19717 -0.45039 -0.04239  0.68447  1.66000 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  -3.8219     0.1546 -24.719   <2e-16 ***
#   MAP_m        -4.0246     0.4624  -8.703   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.9189 on 477 degrees of freedom
# Multiple R-squared:  0.137,	Adjusted R-squared:  0.1352 
# F-statistic: 75.74 on 1 and 477 DF,  p-value: < 2.2e-16
# 
# > summary(RSmodelobs_lowpre1)
# 
# Call:
#   lm(formula = pred_lnRSnew ~ MAP_m, data = lowBGB)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -2.10032 -0.38662  0.02791  0.52668  1.03518 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  -3.9778     0.1155 -34.438   <2e-16 ***
#   MAP_m        -3.4346     0.3455  -9.942   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.6864 on 477 degrees of freedom
# Multiple R-squared:  0.1717,	Adjusted R-squared:  0.1699 
# F-statistic: 98.85 on 1 and 477 DF,  p-value: < 2.2e-16


highBGB <- RSimprove2 %>% 
  filter(MAP..mm. > 420, !is.na(MAP..mm.))
BGBmodelobs_high<-lm(log(BGB_gC.m.2)~MAP..mm.,data=highBGB)
summary(BGBmodelobs_high)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 7.4805613  0.4254959  17.581   <2e-16 ***
#   MAP..mm.    0.0006547  0.0006985   0.937    0.349    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.6796 on 313 degrees of freedom
# Multiple R-squared:  0.002799,	Adjusted R-squared:  -0.0003871 
# F-statistic: 0.8785 on 1 and 313 DF,  p-value: 0.3493
BGBmodelpre_high<-lm(log(prebelow_new)~MAP..mm.,data=highBGB)
summary(BGBmodelpre_high)
# Call:
#   lm(formula = log(prebelow_new) ~ MAP..mm., data = highBGB)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -1.97630 -0.79569  0.09071  0.77554  1.51591 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 6.4243999  0.5886901  10.913   <2e-16 ***
#   MAP..mm.    0.0023368  0.0009664   2.418   0.0162 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.9403 on 313 degrees of freedom
# Multiple R-squared:  0.01834,	Adjusted R-squared:  0.0152 
# F-statistic: 5.847 on 1 and 313 DF,  p-value: 0.01617
AGBmodelobs_high<-lm(log(AGB_gC.m.2)~MAP..mm.,data=highBGB)
summary(AGBmodelobs_high)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 7.112418   0.620746  11.458  < 2e-16 ***
#   MAP..mm.    0.003420   0.001019   3.356 0.000889 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.9915 on 313 degrees of freedom
# Multiple R-squared:  0.03473,	Adjusted R-squared:  0.03165 
# F-statistic: 11.26 on 1 and 313 DF,  p-value: 0.0008888

AGBmodelpre_high<-lm(log(preAGB_new)~MAP..mm.,data=highBGB)
summary(AGBmodelpre_high)
# Call:
#   lm(formula = log(preAGB_new) ~ MAP..mm., data = highBGB)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -1.8420 -0.6988  0.2392  1.0083  1.1533 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 7.4430982  0.5749900  12.945   <2e-16 ***
#   MAP..mm.    0.0029011  0.0009439   3.073   0.0023 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.9184 on 313 degrees of freedom
# Multiple R-squared:  0.02929,	Adjusted R-squared:  0.02619 
# F-statistic: 9.446 on 1 and 313 DF,  p-value: 0.002302
RSmodelobs_high<-lm(log(R_S)~MAP..mm.,data=highBGB)
summary(RSmodelobs_high)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.3681433  0.2189088   1.682   0.0936 .  
# MAP..mm.    -0.0027650  0.0003594  -7.694 1.87e-13 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.3496 on 313 degrees of freedom
# Multiple R-squared:  0.1591,	Adjusted R-squared:  0.1564 
# F-statistic:  59.2 on 1 and 313 DF,  p-value: 1.873e-13
RSmodelobs_highpre<-lm(pred_lnRSnew~MAP..mm.,data=highBGB)
summary(RSmodelobs_highpre)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -1.0186983  0.1830203  -5.566 5.61e-08 ***
#   MAP..mm.    -0.0005642  0.0003005  -1.878   0.0613 .  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.2923 on 313 degrees of freedom
# Multiple R-squared:  0.01114,	Adjusted R-squared:  0.007983 
# F-statistic: 3.527 on 1 and 313 DF,  p-value: 0.06131

highBGB$MAP_m <- highBGB$MAP..mm. / 1000
BGBmodelobs_high1<-lm(log(BGB_gC.m.2)~MAP_m,data=highBGB)
summary(BGBmodelobs_high1)
BGBmodelpre_high1<-lm(log(prebelow_new)~MAP_m,data=highBGB)
summary(BGBmodelpre_high1)
AGBmodelobs_high1<-lm(log(AGB_gC.m.2)~MAP_m,data=highBGB)
summary(AGBmodelobs_high1)
AGBmodelpre_high1<-lm(log(preAGB_new)~MAP_m,data=highBGB)
summary(AGBmodelpre_high1)
RSmodelobs_high1<-lm(log(R_S)~MAP_m,data=highBGB)
summary(RSmodelobs_high1)
RSmodelobs_highpre1<-lm(pred_lnRSnew~MAP_m,data=highBGB)
summary(RSmodelobs_highpre1)

# summary(BGBmodelobs_high1)
# 
# Call:
#   lm(formula = log(BGB_gC.m.2) ~ MAP_m, data = highBGB)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -1.4344 -0.5048  0.2313  0.5361  1.0263 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   7.4806     0.4255  17.581   <2e-16 ***
#   MAP_m         0.6547     0.6985   0.937    0.349    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.6796 on 313 degrees of freedom
# Multiple R-squared:  0.002799,	Adjusted R-squared:  -0.0003871 
# F-statistic: 0.8785 on 1 and 313 DF,  p-value: 0.3493
# 
# > summary(BGBmodelpre_high1)
# 
# Call:
#   lm(formula = log(prebelow_new) ~ MAP_m, data = highBGB)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -1.97630 -0.79569  0.09071  0.77554  1.51591 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   6.4244     0.5887  10.913   <2e-16 ***
#   MAP_m         2.3368     0.9664   2.418   0.0162 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.9403 on 313 degrees of freedom
# Multiple R-squared:  0.01834,	Adjusted R-squared:  0.0152 
# F-statistic: 5.847 on 1 and 313 DF,  p-value: 0.01617
# 
# > summary(AGBmodelobs_high1)
# 
# Call:
#   lm(formula = log(AGB_gC.m.2) ~ MAP_m, data = highBGB)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -1.9885 -0.8416  0.2861  1.0919  1.1747 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   7.1124     0.6207  11.458  < 2e-16 ***
#   MAP_m         3.4197     1.0190   3.356 0.000889 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.9915 on 313 degrees of freedom
# Multiple R-squared:  0.03473,	Adjusted R-squared:  0.03165 
# F-statistic: 11.26 on 1 and 313 DF,  p-value: 0.0008888
# 
# > summary(AGBmodelpre_high1)
# 
# Call:
#   lm(formula = log(preAGB_new) ~ MAP_m, data = highBGB)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -1.8420 -0.6988  0.2392  1.0083  1.1533 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   7.4431     0.5750  12.945   <2e-16 ***
#   MAP_m         2.9011     0.9439   3.073   0.0023 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.9184 on 313 degrees of freedom
# Multiple R-squared:  0.02929,	Adjusted R-squared:  0.02619 
# F-statistic: 9.446 on 1 and 313 DF,  p-value: 0.002302
# 
# > summary(RSmodelobs_high1)
# 
# Call:
#   lm(formula = log(R_S) ~ MAP_m, data = highBGB)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.55582 -0.14840 -0.05478  0.33686  0.61183 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   0.3681     0.2189   1.682   0.0936 .  
# MAP_m        -2.7650     0.3594  -7.694 1.87e-13 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.3496 on 313 degrees of freedom
# Multiple R-squared:  0.1591,	Adjusted R-squared:  0.1564 
# F-statistic:  59.2 on 1 and 313 DF,  p-value: 1.873e-13
# 
# > summary(RSmodelobs_highpre1)
# 
# Call:
#   lm(formula = pred_lnRSnew ~ MAP_m, data = highBGB)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.39331 -0.21707 -0.04753  0.12053  0.73312 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  -1.0187     0.1830  -5.566 5.61e-08 ***
#   MAP_m        -0.5642     0.3005  -1.878   0.0613 .  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.2923 on 313 degrees of freedom
# Multiple R-squared:  0.01114,	Adjusted R-squared:  0.007983 
# F-statistic: 3.527 on 1 and 313 DF,  p-value: 0.06131
######################################################
BGBmodelobs<-lm(log(BGB_gC.m.2)~MAP..mm.,data=wholedata)
summary(BGBmodelobs)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -5.8499182  0.1956649  -29.90   <2e-16 ***
#   MAP..mm.     0.0208879  0.0004223   49.46   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 1.904 on 792 degrees of freedom
# Multiple R-squared:  0.7554,	Adjusted R-squared:  0.7551 
# F-statistic:  2446 on 1 and 792 DF,  p-value: < 2.2e-16
BGBmodelpre<-lm(log(prebelow)~MAP..mm.,data=wholedata)
summary(BGBmodelpre)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -5.8361718  0.1920199  -30.39   <2e-16 ***
#   MAP..mm.     0.0208713  0.0004145   50.36   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 1.868 on 792 degrees of freedom
# Multiple R-squared:  0.762,	Adjusted R-squared:  0.7617 
# F-statistic:  2536 on 1 and 792 DF,  p-value: < 2.2e-16
##########################################

pbBBgroup <- ggplot(wholedata, aes(x = MAP..mm.)) +
  # Points: Black for observations, red for predictions
  geom_point(aes(y = log(BGB_gC.m.2)), color = "black") +
  geom_point(aes(y = log(prebelow)), color = "red") +
  
  # Overall regression line (black)
  #geom_smooth(aes(y = log(BGB_gC.m.2)), method = 'lm', formula = y ~ x, color = "blue", se = TRUE) +
  
  # Regression lines by sitegroup
  geom_smooth(aes(y = log(BGB_gC.m.2), group = sitegroup), method = 'lm',color = "blue", formula = y ~ x, se = TRUE) +
  
  # Prediction regression line (red)
  #geom_smooth(aes(y = log(prebelow)), method = 'lm', formula = y ~ x, color = "red", se = TRUE) +
  geom_smooth(aes(y = log(prebelow), group = sitegroup), method = 'lm',color = "red", formula = y ~ x, se = TRUE) +
  
  # Labels and theme
  xlab(expression('MAP'*~'('*'mm'*')')) +
  ylab(expression('lnBGB'*~'('*'gC'*~m^-2*')')) +
  theme_bw() +
  theme(
    legend.position = "right",  # Keep legend for sitegroup
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text = element_text(size = rel(1.5)),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_line(color = "white")
  ) +
  annotate("text", x = 655, y = 10, label = "(f)", size = 6.5)

pbBBgroup
########################root: shoot and BGB with temperature
tmean<-read.csv("/Users/echo/E/eechange project/analysis/work2024/GPP/tmean_C_NECT_1987_2006.csv")
wholedata$Tmean<-rowMeans(tmean[, 2:13], na.rm = TRUE)

pbBBTg<-ggplot(wholedata,aes(x=tg_1987_2006, y=log(BGB_gC.m.2))) + geom_point() + 
  geom_smooth(aes(x=tg_1987_2006, y=log(BGB_gC.m.2)),method = 'lm', formula = y ~ x)+
  xlab(expression('Tg'*" ("*paste(degree,C)*")")) + ylab(expression('ln BGB'*~'('*'gC'*~m^-2*')'))+
  geom_point(aes(x=tg_1987_2006,y=log(prebelow)),color="red") + 
  geom_smooth(aes(y=log(prebelow),x=tg_1987_2006, col='red'),method = 'lm') +
  theme_bw()+
  theme(legend.position="none",axis.title.x=element_text(size=15),axis.title.y=element_text(size=15),
        axis.text=element_text(size=rel(1.5)),panel.grid.major=element_line(color="white"),panel.grid.minor=element_line(color="white"))+
  annotate("text", x=10, y=12, label="(b)",size=6.5)
pbBBTg

pbBBTmean<-ggplot(wholedata,aes(x=Tmean, y=log(BGB_gC.m.2))) + geom_point() + 
  geom_smooth(aes(x=Tmean, y=log(BGB_gC.m.2)),method = 'lm', formula = y ~ x)+
  xlab(expression('Tmean'*" ("*paste(degree,C)*")")) + ylab(expression('ln BGB'*~'('*'gC'*~m^-2*')'))+
  geom_point(aes(x=Tmean,y=log(prebelow)),color="red") + 
  geom_smooth(aes(y=log(prebelow),x=Tmean, col='red'),method = 'lm') +
  theme_bw()+
  theme(legend.position="none",axis.title.x=element_text(size=15),axis.title.y=element_text(size=15),
        axis.text=element_text(size=rel(1.5)),panel.grid.major=element_line(color="white"),panel.grid.minor=element_line(color="white"))+
  annotate("text", x=8, y=12, label="(b)",size=6.5)
pbBBTmean

pbRSTg<-ggplot(wholedata,aes(x=tg_1987_2006, y=log(R_S))) + geom_point() + 
  geom_smooth(aes(x=tg_1987_2006, y=log(R_S)),method = 'lm', formula = y ~ x)+
  xlab(expression('Tg'*" ("*paste(degree,C)*")")) + ylab(expression('ln R:S'))+
  geom_point(aes(x=tg_1987_2006,y=log(preR_S)),color="red") + 
  geom_smooth(aes(y=log(preR_S),x=tg_1987_2006, col='red'),method = 'lm') +
  theme_bw()+
  theme(legend.position="none",axis.title.x=element_text(size=15),axis.title.y=element_text(size=15),
        axis.text=element_text(size=rel(1.5)),panel.grid.major=element_line(color="white"),panel.grid.minor=element_line(color="white"))+
  annotate("text", x=10, y=0.6, label="(a)",size=6.5)
pbRSTg

pbRSTmean<-ggplot(wholedata,aes(x=Tmean, y=log(R_S))) + geom_point() + 
  geom_smooth(aes(x=Tmean, y=log(R_S)),method = 'lm', formula = y ~ x)+
  xlab(expression('Tmean'*" ("*paste(degree,C)*")")) + ylab(expression('ln R:S'))+
  geom_point(aes(x=Tmean,y=log(preR_S)),color="red") + 
  geom_smooth(aes(y=log(preR_S),x=Tmean, col='red'),method = 'lm') +
  theme_bw()+
  theme(legend.position="none",axis.title.x=element_text(size=15),axis.title.y=element_text(size=15),
        axis.text=element_text(size=rel(1.5)),panel.grid.major=element_line(color="white"),panel.grid.minor=element_line(color="white"))+
  annotate("text", x=8, y=0.6, label="(a)",size=6.5)
pbRSTmean

Tg_RS_BGB<-plot_grid(pbRSTg,pbBBTg,label_x = 0.2,nrow=1)
Tg_RS_BGB
Tmean_RS_BGB<-plot_grid(pbRSTmean,pbBBTmean,label_x = 0.2,nrow=1)
Tmean_RS_BGB

#################################
library(dplyr)
library(stringr)

# Standardize lifeform categories
RSimprove3 <- RSimprove2 %>%
  mutate(
    Lifeform_Standardized = case_when(
      # Trees (any variation containing "tree")
      str_detect(life.form, regex("tree", ignore_case = TRUE)) ~ "Tree",
      
      # Shrubs (includes "shrub", "woody", "low shrub", "tall shrub")
      str_detect(life.form, regex("shrub|woody", ignore_case = TRUE)) ~ "Shrub",
      
      # Herbs (graminoids, forbs, grasses, non-woody)
      str_detect(life.form, regex("graminoid|forb|herb|grass|non.?woody", ignore_case = TRUE)) ~ "Herb",
      
      # Default to Herb if none of the above match
      TRUE ~ "Herb"
    )
  )

lifeform_summary <- RSimprove3 %>%
  group_by(Site, Lifeform_Standardized) %>%
  summarise(Count = n_distinct(Lifeform_Standardized)) %>%
  mutate(Proportion = Count / sum(Count))
# Plot
ggplot(lifeform_summary, aes(x = Site, y = Proportion, fill = Lifeform_Standardized)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("Tree" = "#1b9e77", "Shrub" = "#d95f02", "Herb" = "#7570b3")) +
  labs(title = "Community Structure by Site", y = "Proportion of Species", x = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Step 1: Calculate MAP mean per site (if not already done)
site_map <- RSimprove3 %>%
  group_by(Site) %>%
  summarise(MAP_mean = mean(MAP..mm., na.rm = TRUE))

# 1. Ensure data is summarized correctly
lifeform_summary <- RSimprove3 %>%
  group_by(Site, Lifeform_Standardized) %>%
  summarise(
    Count = n(),
    MAP_mean = mean(MAP..mm., na.rm = TRUE)  # Calculate MAP mean per site
  ) %>%
  mutate(Proportion = Count / sum(Count)) %>%
  ungroup() %>%
  arrange(MAP_mean) %>%  # Sort by MAP
  mutate(Site = factor(Site, levels = unique(Site)))  # Fix site order by MAP

# 2. Plot with MAP_mean as discrete categories
ggplot(lifeform_summary, 
       aes(x = factor(MAP_mean),  # Treat MAP_mean as categorical
           y = Proportion, 
           fill = Lifeform_Standardized)) +
  geom_col(position = "stack", width = 0.8) +  # Adjust width as needed
  scale_fill_manual(values = c("Tree" = "#1b9e77", "Shrub" = "#d95f02", "Herb" = "#7570b3")) +
  labs(
    x = "Mean Annual Precipitation (mm)",
    y = "Proportion of Species",
    fill = "Lifeform"
  ) +
  theme_minimal() +
  theme(legend.position = "top",
    # Increase text sizes
    legend.text = element_text(size = 20),        # Larger legend text
    legend.title = element_text(size = 14),       # Larger legend title
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # X-axis text
    axis.text.y = element_text(size = 10),        # Y-axis text
    axis.title.x = element_text(size = 18),       # X-axis label
    axis.title.y = element_text(size = 18),       # Y-axis label
    plot.title = element_text(size = 14)          # Plot title (if you add one)
  ) +
  scale_x_discrete(
    labels = ~ round(as.numeric(.x), 0)  # Show rounded MAP values on x-axis
  )

test <- filter(RSimprove3, MAP..mm. %in% c(384.15, 400.52))
#####################################################
NECT_sample_climate<-read.csv("/Users/echo/E/eechange project/analysis/work2024/ReworkALL/chelsa_2006/Climate_trait_beta.csv")
modelbeta<-lm(log(NECT_sample_climate$beta_own)~NECT_sample_climate$theta_m3m3,data=NECT_sample_climate)
summary(modelbeta)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                      3.1861     0.1515   21.03   <2e-16 ***
#   NECT_sample_climate$theta_m3m3   7.2918     0.5107   14.28   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.5884 on 315 degrees of freedom
# (477 observations deleted due to missingness)
# Multiple R-squared:  0.3929,	Adjusted R-squared:  0.391 
# F-statistic: 203.9 on 1 and 315 DF,  p-value: < 2.2e-16
pbbetatheta<-ggplot(NECT_sample_climate,aes(x=theta_m3m3, y=log(beta_own))) + geom_point() + 
  geom_smooth(aes(x=theta_m3m3, y=log(beta_own)),method = 'lm', formula = y ~ x)+
  xlab(expression(theta*" (m"^3*"/m"^3*")")) + ylab(expression('ln'*beta))+
  theme_bw()+
  theme(legend.position="none",axis.title.x=element_text(size=15),axis.title.y=element_text(size=15),
        axis.text=element_text(size=rel(1.5)),panel.grid.major=element_line(color="white"),panel.grid.minor=element_line(color="white"))
pbbetatheta
#################################
library(dplyr)

# Step 1: Calculate mean and SE/SD per site
NECT_means <- NECT_sample_climate %>%
  group_by(Site.ID) %>%  # Replace "site" with your actual grouping column
  summarise(
    mean_theta = mean(theta_m3m3, na.rm = TRUE),
    mean_logbeta = mean(log(beta_own), na.rm = TRUE),
    se_logbeta = sd(log(beta_own), na.rm = TRUE) / sqrt(n())  # Standard error
  )
modelbetasite<-lm(log(NECT_means$mean_logbeta)~NECT_means$mean_theta,data=NECT_means)
summary(modelbetasite)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)            1.27818    0.05772   22.14  < 2e-16 ***
#   NECT_means$mean_theta  1.29843    0.20641    6.29 5.38e-07 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.08218 on 31 degrees of freedom
# Multiple R-squared:  0.5607,	Adjusted R-squared:  0.5465 
# F-statistic: 39.57 on 1 and 31 DF,  p-value: 5.385e-07
# Step 2: Plot with error bars
pbbetatheta <- ggplot(NECT_means, aes(x = mean_theta, y = mean_logbeta)) +
  geom_point(size = 1.5) +  # Plot site means
  geom_errorbar(aes(ymin = mean_logbeta - se_logbeta, ymax = mean_logbeta + se_logbeta), 
                width = 0.005) +  # Vertical error bars
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE) +  # Keep regression line
  xlab(expression(theta*" (m"^3*"/m"^3*")")) + 
  ylab(expression('ln'*beta)) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text = element_text(size = rel(1.5)),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_line(color = "white")
  )

print(pbbetatheta)
############################################
chicombine<-read.csv("/Users/echo/E/eechange project/analysis/work2024/ReworkALL/chelsa_2006/Chi_threemodel.csv")
site_summary <- chicombine %>%
  group_by(Site_ID) %>%
  summarize(
    MAP_mm = mean(MAP_mm, na.rm = TRUE),  # Mean MAP for each site
    vpd_pa= mean(VPD_CHELSA_2006, na.rm = TRUE),
    Chi_Obs_Mean = mean(Chi_Obs, na.rm = TRUE),
    Chi_Obs_SD = sd(Chi_Obs, na.rm = TRUE),
    Chi_Constant_Mean = mean(Chi_Constant, na.rm = TRUE),
    Chi_Alienor_Mean = mean(Chi_Alienor, na.rm = TRUE),
    Chi_Own_Mean = mean(Chi_Own, na.rm = TRUE)
  )
# Load ggplot2
library(ggplot2)


chi1<-ggplot(site_summary, aes(x = MAP_mm)) +
  # Add points and error bars for Chi_Obs
  geom_point(aes(y = Chi_Constant_Mean, color = "Chi_Constant")) +
  geom_point(aes(y = Chi_Own_Mean, color = "Chi_Own")) +
  # Add smooth curves for other chi types
  geom_smooth(aes(y = Chi_Constant_Mean, color = "Chi_Constant"), method = "lm", se = TRUE, size = 1) +
  geom_smooth(aes(y = Chi_Own_Mean, color = "Chi_Own"), method = "lm", se = TRUE, size = 1) +
  # Customize labels and axis titles
  ylab(expression(~chi)) +
  xlab(expression('MAP'*~'('*'mm'*')')) +
  # Use a minimal theme and customize appearance
  theme_bw()+
  scale_color_manual(
    values = c(
      "Chi_Constant" = "purple",
      "Chi_Own" = "red"
    )
  ) +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text = element_text(size = rel(1.5)),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_line(color = "white")
  ) #+
#annotate("text", x=693, y=0.9, label="(a)",size=6.5)
# Annotate the plot
chi1## save work2024--plot

###########################################BPE
###################plot GPP, The new GPP dataset, FLUXCOM-X (X-BASE), shows significant improvement over MTE. and also change 0.46 change as BPE
BPE_final<-stack("/Users/echo/PhD_third/CUE_Collalti/result_data/BPE_pre_Ev_De.tif")
model_season_df<-read.csv("/Users/echo/E/eechange project/analysis/work2024/GPP/gppobs_gCm2y_NECT_pmodel_FLUX_X.csv")
joined_data<-read.csv("/Users/echo/E/eechange project/analysis/work2024/faparmax_LAI/Pre_Obs_fapar_LAI")
model_season_df$site<-joined_data$X1
model_season_df$Precipitation<-joined_data$X3
model_season_df$Tg<-NECT_sample_climate$tg_chelsa_2006
tmean<-read.csv("/Users/echo/E/eechange project/analysis/work2024/GPP/tmean_C_NECT_1987_2006.csv")
model_season_df$Tmean<-rowMeans(tmean[, 2:13], na.rm = TRUE)
csite<-cbind(NECT_sample_climate$Longitude,NECT_sample_climate$Latitude)
model_season_df$BPE<-raster::extract(BPE_final,csite,method = 'bilinear')
gpp_mean <- model_season_df %>%
  group_by(site) %>%
  summarize(across(where(is.numeric), mean, na.rm = TRUE))


pb3 <- ggplot(data=gpp_mean,aes(x=Precipitation,y=gpp))+geom_point(color="red") + 
  geom_smooth(method = 'lm', formula = y ~ x, aes(color='red')) +
  xlab(expression('MAP'*~'('*'mm'*')')) + ylab(expression('GPP'*~'('*'gC'*~m^-1*d^-1*')')) +
  geom_point(aes(x=Precipitation,y=gppobs_FLUX_X)) + 
  geom_smooth(aes(x=Precipitation,y=gppobs_FLUX_X),method = 'lm') +
  theme_bw()+
  theme(legend.position="none",axis.title.x=element_text(size=15),axis.title.y=element_text(size=15),
        axis.text=element_text(size=rel(1.5)),panel.grid.major=element_line(color="white"),panel.grid.minor=element_line(color="white"))+
  annotate("text", x=655, y=5.5, label="(c)",size=6.5)
pb3
pb3 <- ggplot(data=gpp_mean,aes(x=Precipitation,y=gpp))+geom_point(color="red",size = 2.8) + 
  geom_smooth(method = 'lm', formula = y ~ x, aes(color='red')) +
  xlab(expression('MAP'*~'('*'mm'*')')) + ylab(expression('GPP'*~'('*'gC'*~m^-1*d^-1*')')) +
  geom_point(aes(x=Precipitation,y=gppobs_FLUX_X),size = 2.8) + 
  geom_smooth(aes(x=Precipitation,y=gppobs_FLUX_X),method = 'lm') +
  theme_bw()+
  theme(legend.position="none",axis.title.x=element_blank(),axis.title.y=element_text(size=29),
        axis.text=element_text(size=rel(2.5)),panel.grid.major=element_line(color="white"),panel.grid.minor=element_line(color="white"))
pb3


p4 <- ggplot(data=gpp_mean,aes(x=Precipitation,y=gpp*BPE))+geom_point(color="red") + 
  geom_smooth(method = 'lm', formula = y ~ x, aes(color='red')) +
  xlab(expression('MAP'*~'('*'mm'*')')) +  ylab(expression('NPP'*~'('*'gC'*~m^-1*d^-1*')'))  +
  geom_point(aes(x=Precipitation,y=gppobs_FLUX_X*BPE)) + 
  geom_smooth(aes(x=Precipitation,y=gppobs_FLUX_X*BPE),method = 'lm') +
  theme_bw()+
  theme(legend.position="none",axis.title.x=element_text(size=15),axis.title.y=element_text(size=15),
        axis.text=element_text(size=rel(1.5)),panel.grid.major=element_line(color="white"),panel.grid.minor=element_line(color="white"))+
  annotate("text", x=693, y=3, label="(d)",size=6.5)
p4

p4 <- ggplot(data=gpp_mean,aes(x=Precipitation,y=gpp*BPE))+geom_point(color="red",size = 2.8) + 
  geom_smooth(method = 'lm', formula = y ~ x, aes(color='red')) +
  xlab(expression('MAP'*~'('*'mm'*')')) +  ylab(expression('NPP'*~'('*'gC'*~m^-1*d^-1*')'))  +
  geom_point(aes(x=Precipitation,y=gppobs_FLUX_X*BPE),size = 2.8) + 
  geom_smooth(aes(x=Precipitation,y=gppobs_FLUX_X*BPE),method = 'lm') +
  theme_bw()+
  theme(legend.position="none",axis.title.x=element_blank(),axis.title.y=element_text(size=29),
        axis.text=element_text(size=rel(2.5)),panel.grid.major=element_line(color="white"),panel.grid.minor=element_line(color="white"))
p4
lmgpp<-lm(gpp~Precipitation,data=gpp_mean)
summary(lmgpp) 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   0.0270677  0.1462310   0.185    0.854    
# Precipitation 0.0055104  0.0003553  15.511  3.7e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.3059 on 31 degrees of freedom
# Multiple R-squared:  0.8859,	Adjusted R-squared:  0.8822 
# F-statistic: 240.6 on 1 and 31 DF,  p-value: 3.696e-16
lmgppobx<-lm(gppobs_FLUX_X~Precipitation,data=gpp_mean)
summary(lmgppobx) 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   -1.5513810  0.3638345  -4.264 0.000175 ***
#   Precipitation  0.0082369  0.0008839   9.319 1.68e-10 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.7611 on 31 degrees of freedom
# Multiple R-squared:  0.7369,	Adjusted R-squared:  0.7284 
# F-statistic: 86.84 on 1 and 31 DF,  p-value: 1.681e-10
lmnpp<-lm((gpp*BPE)~Precipitation,data=gpp_mean)
summary(lmnpp)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   -0.245909   0.088936  -2.765   0.0101 *  
#   Precipitation  0.002607   0.000217  12.012 2.43e-12 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1721 on 27 degrees of freedom
# (4 observations deleted due to missingness)
# Multiple R-squared:  0.8424,	Adjusted R-squared:  0.8365 
# F-statistic: 144.3 on 1 and 27 DF,  p-value: 2.426e-12
lmnppobs<-lm((gppobs_FLUX_X*BPE)~Precipitation,data=gpp_mean)
summary(lmnppobs)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   -0.8634087  0.1846289  -4.676 7.27e-05 ***
#   Precipitation  0.0038556  0.0004506   8.558 3.59e-09 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.3573 on 27 degrees of freedom
# (4 observations deleted due to missingness)
# Multiple R-squared:  0.7306,	Adjusted R-squared:  0.7206 
# F-statistic: 73.23 on 1 and 27 DF,  p-value: 3.588e-09
gpp_mean$MAP_m <- gpp_mean$Precipitation / 1000

lmgppobx1<-lm(gppobs_FLUX_X~MAP_m,data=gpp_mean)
summary(lmgppobx1) 
lmgpp1<-lm(gpp~MAP_m,data=gpp_mean)
summary(lmgpp1) 
lmnppobs1<-lm((gppobs_FLUX_X*BPE)~MAP_m,data=gpp_mean)
summary(lmnppobs1)
lmnpp1<-lm((gpp*BPE)~MAP_m,data=gpp_mean)
summary(lmnpp1)

# summary(lmgppobx1) 

# Call:
#   lm(formula = gppobs_FLUX_X ~ MAP_m, data = gpp_mean)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -0.9840 -0.5936 -0.1786  0.7754  1.8296 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  -1.5514     0.3638  -4.264 0.000175 ***
#   MAP_m         8.2369     0.8839   9.319 1.68e-10 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.7611 on 31 degrees of freedom
# Multiple R-squared:  0.7369,	Adjusted R-squared:  0.7284 
# F-statistic: 86.84 on 1 and 31 DF,  p-value: 1.681e-10
# 
# > summary(lmgpp1) 
# 
# Call:
#   lm(formula = gpp ~ MAP_m, data = gpp_mean)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.80547 -0.12584 -0.04371  0.15449  0.75712 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.02707    0.14623   0.185    0.854    
# MAP_m        5.51038    0.35525  15.511  3.7e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.3059 on 31 degrees of freedom
# Multiple R-squared:  0.8859,	Adjusted R-squared:  0.8822 
# F-statistic: 240.6 on 1 and 31 DF,  p-value: 3.696e-16
# 
# > summary(lmnppobs1)
# 
# Call:
#   lm(formula = (gppobs_FLUX_X * BPE) ~ MAP_m, data = gpp_mean)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -0.4152 -0.2950 -0.1412  0.3349  0.7247 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  -0.8634     0.1846  -4.676 7.27e-05 ***
#   MAP_m         3.8556     0.4506   8.558 3.59e-09 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.3573 on 27 degrees of freedom
# (4 observations deleted due to missingness)
# Multiple R-squared:  0.7306,	Adjusted R-squared:  0.7206 
# F-statistic: 73.23 on 1 and 27 DF,  p-value: 3.588e-09
# 
# > summary(lmnpp1)
# 
# Call:
#   lm(formula = (gpp * BPE) ~ MAP_m, data = gpp_mean)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.31869 -0.10315 -0.05288  0.10615  0.38493 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -0.24591    0.08894  -2.765   0.0101 *  
#   MAP_m        2.60697    0.21703  12.012 2.43e-12 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1721 on 27 degrees of freedom
# (4 observations deleted due to missingness)
# Multiple R-squared:  0.8424,	Adjusted R-squared:  0.8365 
# F-statistic: 144.3 on 1 and 27 DF,  p-value: 2.426e-12
####################################################################replot the figure with the China map and trsect and the raster of preciptation
# ---- Packages ----
library(sf)
library(ggplot2)
library(cowplot)
library(rnaturalearth)
library(rnaturalearthdata)

# ---- Paths (edit to your actual files) ----
china_province_shp <- "~/E/eechange project/data/Archive/China_Province.shp"
china_boundary_shp <- "~/E/eechange project/data/Archive/China_Boundary_Nineline.shp"
sites_csv          <- "~/E/eechange project/data/MAP/extract_example.csv"

# ---- Data ----
# If your shapefiles are available, use them:
china_prov <- st_read(china_province_shp, quiet = TRUE)
china_line <- st_read(china_boundary_shp, quiet = TRUE)

# Fallback (if shapefiles aren’t available): outline from Natural Earth
# china_outline <- rnaturalearth::ne_countries(country = "China", returnclass = "sf")

# Sites (assumes columns named longitude, latitude)
sta <- read.csv(sites_csv, stringsAsFactors = FALSE)
sta_33 <- sta[1:33, ]
sites_sf <- st_as_sf(sta_33, coords = c("longitude","latitude"), crs = 4326)

# ---- Whole-China map with transect points ----
p_china <- ggplot() +
  geom_sf(data = china_prov, fill = "white", color = "grey20", linewidth = 0.2) +
  geom_sf(data = china_line, color = "grey10", linewidth = 0.3) +
  geom_sf(data = sites_sf, color = "black", size = 2) +
  coord_sf(xlim = c(70, 140), ylim = c(0, 55), expand = FALSE) +
  theme_void() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4)
  )
p_china



# install.packages(c("terra","sf","rnaturalearth","rnaturalearthdata","tidyterra","ggplot2","viridis"))
library(terra)
library(sf)
library(rnaturalearth)
library(tidyterra)
library(ggplot2)
library(viridis)
chelsa_dir<-'/Users/echo/E/eechange project/analysis/work2024/ReworkALL/chelsa_2006/pr'
files <- list.files(chelsa_dir, pattern = "^CHELSA_pr_\\d+_2006_V\\.2\\.1\\.tif$", full.names = TRUE)
mons  <- as.integer(gsub(".*pr_(\\d+)_.*", "\\1", basename(files)))
files <- files[order(mons)]        # ensure 1..12

r_stack <- rast(files)             # SpatRaster with 12 layers
pr_mean <- app(r_stack, sum, na.rm = TRUE)      # mean monthly precipitation (your requested "mean")

# (If you actually want ANNUAL TOTAL instead of mean, use this one line:)
# pr_mean <- app(r_stack, sum, na.rm = TRUE)
pr_meanmm<-pr_mean/100 ##kg m-2 month1/100 kg/m²/month is numerically equivalent to mm/month
mask <- terra::rast("/Users/echo/PhD_second_year/Net carbon profit/0.05ddegree/Net_carbon_profit_mean0.05_LAImax_growing.nc")
pr_meanmm_aligned <- resample(pr_meanmm, mask, method = "bilinear")
pr_meanmm1 <- pr_meanmm_aligned[mask > 0]
pr_land <- mask(pr_meanmm_aligned, mask)
plot(pr_land)

# 1) Bring in China’s boundary
china_sf <- ne_countries(country = "China", scale = "large", returnclass = "sf")

# 2) Match CRS to your raster and convert to terra vect
china_sf  <- st_transform(china_sf, crs(pr_land))
china_vect <- vect(china_sf)

# 3) Clip/mask the raster to China → everything outside becomes NA (white)
pr_china <- pr_land |> crop(china_vect) |> mask(china_vect)

# (Optional) 4) Crop to the same NE-China window as your example map
#    Adjust these bounds if you want it tighter/looser.
#    xmin, xmax, ymin, ymax  (lon/lat in degrees)
bb_ne <- ext(111, 140, 30, 55)
pr_ne <- crop(pr_china, bb_ne)

# 5) Plot (choose one object to plot: pr_china or pr_ne)

# Convert cropped raster to data frame
pr_df <- as.data.frame(pr_ne, xy = TRUE, na.rm = FALSE)
colnames(pr_df) <- c("x", "y", "MAP")
# Plot
ggplot(pr_df) +
  geom_raster(aes(x = x, y = y, fill = MAP)) +
  scale_fill_viridis_c(name = "MAP (mm)", na.value = "white")+
  geom_sf(data = sites_sf, shape = 21, size = 2.6, color = "black", fill = "black", stroke = 0.5) +
  xlab(expression('Longitude')) +  ylab(expression('Latitude'))  +
  coord_sf(expand = FALSE) +
  theme_bw()+
  theme(legend.position="right",axis.title.x=element_text(size=15),axis.title.y=element_text(size=15),
        axis.text=element_text(size=rel(1.5)),panel.grid.major=element_line(color="white"),panel.grid.minor=element_line(color="white"))

######################################################change scale color
# Crop to 40–50N latitude (and keep longitudes as before)
# Crop to 40–50N latitude
bb_ne <- ext(111, 140, 40, 50)
pr_ne <- crop(pr_china, bb_ne)

# Convert cropped raster to data frame
pr_df <- as.data.frame(pr_ne, xy = TRUE, na.rm = FALSE)
colnames(pr_df) <- c("x", "y", "MAP")


# Re-bin MAP with thresholds starting at 200 (everything below is one bin)
pr_df$MAP_class <- cut(
  pr_df$MAP,
  breaks = c(0,100, 200, 300, 400, 500, 600, 700, 800, 900, Inf),
  labels = c("0-100","200","300","400","500","600","700","800","900",">900"),
  right = TRUE, include.lowest = TRUE
)

# Smooth yellow → blue palette (10 colors)
library(scales)

my_colors <- colorRampPalette(c("yellow", "orchid", "blue"))(11)

# Ensure legend order is exactly as labels above
legend_levels <- c("0-100","200","300","400","500","600","700","800","900",">900")
pr_df$MAP_class <- factor(pr_df$MAP_class, levels = legend_levels)

ggplot(pr_df) +
  geom_raster(aes(x = x, y = y, fill = MAP_class)) +
  scale_fill_manual(
    values = my_colors,
    name   = "MAP (mm)",
    limits = legend_levels,         # order in legend
    breaks = legend_levels,         # show all entries
    na.translate = FALSE            # <-- removes the "NA" box from legend
  ) +
  geom_sf(data = sites_sf, shape = 21, size = 2.5,
          color = "black", fill = "black", stroke = 0.5) +
  xlab("Longitude") + ylab("Latitude") +
  coord_sf(expand = FALSE) +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme_bw() +
  theme(
    plot.background  = element_rect(fill = "white", colour = NA),  # white outside
    panel.background = element_rect(fill = "white", colour = NA),  # white panel
    legend.background = element_rect(fill = "white", colour = NA),
    legend.key = element_rect(fill = "white", colour = NA),
    legend.position = "right",
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text    = element_text(size = rel(1.5)),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_line(color = "white")
  )

##################################################
ndata<-read.csv("~/E/eechange project/data/predict 15n/Hard Traits15N(with_predictingN)_ftmean.csv")
nitrogendata<-read.csv("~/E/eechange project/data/nitrogen_rework/Hard Traits15N(repeatftemp_frunoff)_ftmean.csv")
fTemp<-nitrogendata$ftemp_new
fR<-nitrogendata$frunoff_new
nitrogenisotope_model<--7.6265 +(10.4682/(1+1.1405*(fR/fTemp)))
datamodel<-data.frame(cbind(ndata$MAP,ndata$d15N.14N,nitrogenisotope_model))
datamodel<-na.omit(datamodel)
nobsmean<-aggregate(datamodel$V2, by=list(datamodel$V1),mean)
n.sd<-aggregate(datamodel$V2, by=list(datamodel$V1),FUN=sd)
nlnpre<-aggregate(datamodel$nitrogenisotope_model, by=list(datamodel$V1),mean)

ndata1<-data.frame(cbind(nlnpre$Group.1,nlnpre$x,n.sd$x,nobsmean$x))
names(ndata1) <- c("map", "lnpren", "n.sd","nobsmean")
ndata1$MAP_m <- ndata1$map / 1000
lm15N1<-lm(nobsmean~MAP_m,data=ndata1)
summary(lm15N1)

lm15Npre1<-lm(lnpren~MAP_m,data=ndata1)
summary(lm15Npre1)

# summary(lm15N1)

# Call:
#   lm(formula = nobsmean ~ MAP_m, data = ndata1)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -2.8082 -1.1376 -0.4313  0.4770  5.6153 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   5.0995     0.8483   6.011 1.19e-06 ***
#   MAP_m       -12.3217     2.0609  -5.979 1.30e-06 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 1.774 on 31 degrees of freedom
# Multiple R-squared:  0.5356,	Adjusted R-squared:  0.5206 
# F-statistic: 35.75 on 1 and 31 DF,  p-value: 1.303e-06
# 
# > summary(lm15Npre1)
# 
# Call:
#   lm(formula = lnpren ~ MAP_m, data = ndata1)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -0.6244 -0.2654  0.1369  0.2341  0.3984 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   4.9605     0.1443   34.38   <2e-16 ***
#   MAP_m       -12.1863     0.3505  -34.77   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.3018 on 31 degrees of freedom
# Multiple R-squared:  0.975,	Adjusted R-squared:  0.9742 
# F-statistic:  1209 on 1 and 31 DF,  p-value: < 2.2e-16

