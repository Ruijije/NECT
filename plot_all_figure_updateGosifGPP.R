##################################
library(png)
library(raster)
library(ggplot2)
library(dplyr)
img <- readPNG("/Users/echo/E/eechange project/analysis/work2024/plot/GPP_Fig3.png")
dim(img)
h_px <- dim(img)[1]
w_px <- dim(img)[2]

#############################the csv is produced from reading laptop, i copy the code in folder
model_season_df<-read.csv("/Users/echo/E/eechange project/analysis/work2024/GPP/gppobs_gCm2y_NECT_pmodel_FLUX_X_GOSIF_PML15.csv")

gpp_mean <- model_season_df %>%
  group_by(site) %>%
  summarize(across(where(is.numeric), mean, na.rm = TRUE))


pb3 <- ggplot(data=gpp_mean,aes(x=Precipitation,y=gpp))+geom_point(color="red",size = 2.8) + 
  geom_smooth(method = 'lm', formula = y ~ x, aes(color='red')) +
  xlab(expression('MAP'*~'('*'mm'*')')) + ylab(expression('GPP'*~'('*'gC'*~m^-1*d^-1*')')) +
  geom_point(aes(x=Precipitation,y=gpp_gosif),size = 2.8) + 
  geom_smooth(aes(x=Precipitation,y=gpp_gosif),method = 'lm') +
  theme_bw()+
  theme(legend.position="none",axis.title.x=element_blank(),axis.title.y=element_text(size=29),
        axis.text=element_text(size=rel(2.5)),panel.grid.major=element_line(color="white"),panel.grid.minor=element_line(color="white"))
pb3


p4 <- ggplot(data=gpp_mean,aes(x=Precipitation,y=gpp*BPE))+geom_point(color="red",size = 2.8) + 
  geom_smooth(method = 'lm', formula = y ~ x, aes(color='red')) +
  xlab(expression('MAP'*~'('*'mm'*')')) +  ylab(expression('NPP'*~'('*'gC'*~m^-1*d^-1*')'))  +
  geom_point(aes(x=Precipitation,y=gpp_gosif*BPE),size = 2.8) + 
  geom_smooth(aes(x=Precipitation,y=gpp_gosif*BPE),method = 'lm') +
  ylim(-1,2)+
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
lmgppobx<-lm(gpp_gosif~Precipitation,data=gpp_mean)
summary(lmgppobx) 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   -1.2508983  0.2448769  -5.108 1.57e-05 ***
#   Precipitation  0.0070838  0.0005949  11.908 4.24e-13 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.5122 on 31 degrees of freedom
# Multiple R-squared:  0.8206,	Adjusted R-squared:  0.8148 
# F-statistic: 141.8 on 1 and 31 DF,  p-value: 4.241e-13

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
lmnppobs<-lm((gpp_gosif*BPE)~Precipitation,data=gpp_mean)
summary(lmnppobs)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   -0.6909977  0.1375427  -5.024 2.86e-05 ***
#   Precipitation  0.0032386  0.0003357   9.649 3.04e-10 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.2662 on 27 degrees of freedom
# (4 observations deleted due to missingness)
# Multiple R-squared:  0.7752,	Adjusted R-squared:  0.7669 
# F-statistic:  93.1 on 1 and 27 DF,  p-value: 3.041e-10
gpp_mean$MAP_m <- gpp_mean$Precipitation / 1000

lmgppobx1<-lm(gpp_gosif~MAP_m,data=gpp_mean)
summary(lmgppobx1) 
lmgpp1<-lm(gpp~MAP_m,data=gpp_mean)
summary(lmgpp1) 
lmnppobs1<-lm((gpp_gosif*BPE)~MAP_m,data=gpp_mean)
summary(lmnppobs1)
lmnpp1<-lm((gpp*BPE)~MAP_m,data=gpp_mean)
summary(lmnpp1)
# lmgppobx1<-lm(gpp_gosif~MAP_m,data=gpp_mean)
# > summary(lmgppobx1) 
# 
# Call:
#   lm(formula = gpp_gosif ~ MAP_m, data = gpp_mean)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.73079 -0.39811 -0.09983  0.42850  1.21114 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  -1.2509     0.2449  -5.108 1.57e-05 ***
#   MAP_m         7.0838     0.5949  11.908 4.24e-13 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.5122 on 31 degrees of freedom
# Multiple R-squared:  0.8206,	Adjusted R-squared:  0.8148 
# F-statistic: 141.8 on 1 and 31 DF,  p-value: 4.241e-13
# 
# > lmgpp1<-lm(gpp~MAP_m,data=gpp_mean)
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
# > lmnppobs1<-lm((gpp_gosif*BPE)~MAP_m,data=gpp_mean)
# > summary(lmnppobs1)
# 
# Call:
#   lm(formula = (gpp_gosif * BPE) ~ MAP_m, data = gpp_mean)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.32247 -0.21706 -0.09906  0.23988  0.51187 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  -0.6910     0.1375  -5.024 2.86e-05 ***
#   MAP_m         3.2386     0.3357   9.649 3.04e-10 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.2662 on 27 degrees of freedom
# (4 observations deleted due to missingness)
# Multiple R-squared:  0.7752,	Adjusted R-squared:  0.7669 
# F-statistic:  93.1 on 1 and 27 DF,  p-value: 3.041e-10
# 
# > lmnpp1<-lm((gpp*BPE)~MAP_m,data=gpp_mean)
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
