#########################################################################Fig 1
NECT_sample_climate<-read.csv("Climate_trait_beta.csv")
pa1 <- ggplot(NECT_sample_climate,aes(x=NECT_sample_climate$MI..unitless., y=tg_chelsa_2006)) + geom_point(color = "blue")  + 
  geom_smooth(method = 'lm', formula = y ~ x) + 
  ylab(expression('Tg'*~'('*degree*'C'*')')) + 
  xlab(expression('MI'*~'('*'unitless'*')'))+
  theme_bw()+
  theme( axis.title.x=element_text(size=15),axis.title.y=element_text(size=13),
         axis.text=element_text(size=rel(1.5)),panel.grid.major=element_line(color="white"),
         panel.grid.minor=element_line(color="white")) +
  annotate("text", x=1, y=11.2, label="(a)",size=6.5)
pa1
pa2 <- ggplot(NECT_sample_climate,aes(x=MI..unitless., y=NECT_sample_climate$vpd_chelsa_2006)) + geom_point(color = "blue") + 
  geom_smooth(method = 'lm', formula = y ~ x) + 
  ylab(expression('VPD'*~'('*'Pa'*')')) + 
  xlab(expression('MI'*~'('*'unitless'*')'))+
  annotate("text", x=1, y=800, label="(b)",size=6.5)+
  theme_bw()+
  theme( axis.title.x=element_text(size=15),axis.title.y=element_text(size=15),
         axis.text=element_text(size=rel(1.5)),panel.grid.major=element_line(color="white"),
         panel.grid.minor=element_line(color="white")) 
pa2
pa3 <- ggplot(NECT_sample_climate,aes(x=MI..unitless., y=ppfd_chelsa_2006)) + geom_point(color = "blue") +
  geom_smooth(method = 'lm', formula = y ~ x) + 
  ylab(expression('PPFD'*~'('*'mol'*~m^-2*day^-1*')')) + 
  xlab(expression('MI'*~'('*'unitless'*')'))+
  theme_bw()+
  theme( axis.title.x=element_text(size=15),axis.title.y=element_text(size=15),
         axis.text=element_text(size=rel(1.5)),panel.grid.major=element_line(color="white"),
         panel.grid.minor=element_line(color="white")) +
  annotate("text", x=1, y=34, label="(c)",size=6.5)
pa3
plot_grid(pa1,pa2,pa3,label_x = 0.2,nrow = 2)
############################################################
lmtemp<-lm(tg_chelsa_2006~MI..unitless.,data=NECT_sample_climate)
summary(lmtemp)

lmvpd<-lm(vpd_chelsa_2006~MI..unitless.,data=NECT_sample_climate)
summary(lmvpd)

lmppfd<-lm(ppfd_chelsa_2006~MI..unitless.,data=NECT_sample_climate)
summary(lmppfd)

library(dplyr)
chicombine<-read.csv("/Users/vv934811@reading.ac.uk/Documents/work2024/ReworkALL/chelsa_2006/Chi_threemodel.csv")
NECT_mi <- NECT_sample_climate %>%
  dplyr::select(MAP..mm., MI..unitless.) %>%
  distinct(MAP..mm., .keep_all = TRUE)
chicombine2 <- chicombine %>%
  dplyr::left_join(NECT_mi, by = c("MAP_mm" = "MAP..mm."))
site_summary <- chicombine2 %>%
  group_by(Site_ID) %>%
  summarize(
    MI = mean(MI..unitless., na.rm = TRUE),  # Mean MAP for each site
    vpd_pa= mean(VPD_CHELSA_2006, na.rm = TRUE),
    Chi_Obs_Mean = mean(Chi_Obs, na.rm = TRUE),
    Chi_Obs_SD = sd(Chi_Obs, na.rm = TRUE),
    Chi_Constant_Mean = mean(Chi_Constant, na.rm = TRUE),
    Chi_Alienor_Mean = mean(Chi_Alienor, na.rm = TRUE),
    Chi_Own_Mean = mean(Chi_Own, na.rm = TRUE)
  )
# Load ggplot2
modelbeta<-lm(log(NECT_sample_climate$beta_own)~NECT_sample_climate$theta_m3m3,data=NECT_sample_climate)
summary(modelbeta)
############################################################
# 1. Fit model (you already did this)
############################################################
modelbeta <- lm(log(beta_own) ~ theta_m3m3, data = NECT_sample_climate)

############################################################
# 2. Prepare clean dataset
############################################################
beta_dat <- NECT_sample_climate %>%
  dplyr::select(beta_own, theta_m3m3) %>%
  na.omit()

# original prediction
beta_dat$pred_ln_beta <- predict(modelbeta, newdata = beta_dat)

############################################################
# 3. Monte Carlo simulation
############################################################
set.seed(123)

# extract coefficients and covariance
beta_hat <- coef(modelbeta)
V_beta   <- vcov(modelbeta)

# simulate coefficients (1000 draws)
nsim <- 1000
beta_sim <- MASS::mvrnorm(n = nsim, mu = beta_hat, Sigma = V_beta)

# design matrix
X <- model.matrix(~ theta_m3m3, data = beta_dat)

# predicted values (n × nsim)
pred_mat <- X %*% t(beta_sim)

############################################################
# 4. Summarize along theta (binning)
############################################################
theta <- beta_dat$theta_m3m3

# create bins
nbins <- 8
breaks <- seq(min(theta), max(theta), length.out = nbins + 1)
bin_id <- cut(theta, breaks = breaks, include.lowest = TRUE, labels = FALSE)

band_list <- vector("list", nbins)

for (b in 1:nbins) {
  idx <- which(bin_id == b)
  if (length(idx) == 0) next
  
  vals <- as.vector(pred_mat[idx, ])
  
  band_list[[b]] <- data.frame(
    theta_mid = mean(theta[idx]),
    pred_median = quantile(vals, 0.5),
    pred_low = quantile(vals, 0.025),
    pred_high = quantile(vals, 0.975)
  )
}

band_beta <- bind_rows(band_list)

############################################################
# 5. Plot
############################################################
p_beta_unc <- ggplot() +
  geom_point(
    data = beta_dat,
    aes(x = theta_m3m3, y = log(beta_own)),
    color = "blue", alpha = 0.7
  ) +
  geom_ribbon(
    data = band_beta,
    aes(x = theta_mid, ymin = pred_low, ymax = pred_high),
    fill = "red", alpha = 0.2
  ) +
  geom_line(
    data = band_beta,
    aes(x = theta_mid, y = pred_median),
    color = "red", linewidth = 1.2
  ) +
  xlab(expression(theta~'('*m^3~m^-3*')'))+
  ylab(expression(log(beta))) +
  theme_bw()+  theme( axis.title.x=element_text(size=15),axis.title.y=element_text(size=15),
                      axis.text=element_text(size=rel(1.5)),panel.grid.major=element_line(color="white"),
                      panel.grid.minor=element_line(color="white"))

p_beta_unc
##############################################Fig.S2, Same as the plot_all_figure.R, change MI
chi1<-ggplot(site_summary, aes(x = MI)) +
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
    width = 0.02
  ) +
  # Add smooth curves for other chi types
  geom_smooth(aes(y = Chi_Constant_Mean, color = "Chi_Constant"), method = "lm", se = TRUE, size = 1) +
  geom_smooth(aes(y = Chi_Alienor_Mean, color = "Chi_Alienor"), method = "lm", se = TRUE, size = 1) +
  geom_smooth(aes(y = Chi_Own_Mean, color = "Chi_Own"), method = "lm", se = TRUE, size = 1) +
  geom_smooth(aes(y = Chi_Obs_Mean, color = "Chi_Obs"), method = "lm", se = TRUE, size = 1) +
  # Customize labels and axis titles
  ylab(expression(~chi)) +
  xlab(expression('MI'*~'('*'unitless'*')')) +
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
vcmaxcombine1<-read.csv("/Users/vv934811@reading.ac.uk/Documents/work2024/ReworkALL/chelsa_2006/vcmax_threemodel_Ning2017_obsvcmax_newoption.csv")
data_combined <- data.frame(
  Narea = NECT_sample_climate$Narea..g.m2.,
  LMA = NECT_sample_climate$LMA..kg.m2. * 1000,  # Multiply by 1000 for correct units
  vcmax25_Own = vcmaxcombine1$vcmax25_Own
)
data_combined <- na.omit(data_combined)
lnneramodel <- lm(Narea ~ LMA + vcmax25_Own-1, data = data_combined)
summary(lnneramodel)

lnneramodel <- lm(Narea ~ LMA + vcmax25_Own, data = data_combined)
summary(lnneramodel)

Tg_K<-NECT_sample_climate$tg_chelsa_2006+273.15
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
site_summary2 <- site_summary2 %>%
  dplyr::left_join(NECT_mi, by = c("MAP_mm" = "MAP..mm."))


# Plot for Vcmax_Own
vcmax12 <- ggplot(site_summary2, aes(x = MI..unitless.)) +
  geom_point(aes(y = vcmax_Own_Mean, color = "vcmax_Own"), size = 1) +
  geom_smooth(aes(y = vcmax_Own_Mean, color = "vcmax_Own"), method = "lm", se = TRUE, size = 1) +
  geom_point(aes(y = Obs_vcmax, color = "Obs_vcmax")) +
  # Add smooth curves for other vcmax types
  geom_smooth(aes(y = Obs_vcmax, color = "Obs_vcmax"), method = "lm", se = TRUE, size = 1) +
  ylab(expression(V[cmax]~'('*mu*'mol'*m^-2*s^-1*')')) +
  xlab(expression('MI'*~'('*'unitless'*')')) +
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
  annotate("text", x = 1, y = 40, label = "(a)", size = 6.5)
vcmax12
vcmax22 <- ggplot(site_summary2, aes(x = MI..unitless.)) +
  geom_point(aes(y = vcmax25_Own_Mean, color = "vcmax25_Own_Mean")) +
  geom_smooth(aes(y = vcmax25_Own_Mean, color = "vcmax25_Own_Mean"), method = "lm", se = TRUE, size = 1) +
  geom_point(aes(y = vcmax25_obs_Mean, color = "vcmax25_obs_Mean")) +
  # Add smooth curves for other vcmax types
  geom_smooth(aes(y = vcmax25_obs_Mean, color = "vcmax25_obs_Mean"), method = "lm", se = TRUE, size = 1) +
  ylab(expression(V[cmax25]~'('*mu*'mol'*m^-2*s^-1*')')) +
  xlab(expression('MI'*~'('*'unitless'*')')) +
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
  annotate("text", x = 1, y = 200, label = "(b)", size = 6.5)
vcmax22
plot_grid(vcmax12, vcmax22, label_x = 0.2, nrow = 1)
write.csv(site_summary2,"/Users/vv934811@reading.ac.uk/Documents/work2024/ReworkALL/chelsa_2006/vcmax_vcmax25_site_obs_Ning2017_0.111IL.csv")
####the calculation of vcamx is from end of code /Users/vv934811@reading.ac.uk/Documents/work2024/ReworkALL/chelsa_2006/Pmodel.R
lmobvcmax<-lm(Obs_vcmax~MI..unitless.,data=site_summary2)
summary(lmobvcmax)

lmorevcmax<-lm(vcmax_Own_Mean~MI..unitless.,data=site_summary2)
summary(lmorevcmax)

lmorevcmax25<-lm(vcmax25_obs_Mean~MI..unitless.,data=site_summary2)
summary(lmorevcmax25)

lmprevcmax25<-lm(vcmax25_Own_Mean~MI..unitless.,data=site_summary2)
summary(lmprevcmax25)

#################################################################Fig2
# Annotate the plot
f1 <- ggplot(data=site_summary,aes(x=MI,y=Chi_Obs_Mean))+geom_point(color = "blue") + 
  geom_smooth(aes(y=Chi_Obs_Mean,x=MI),method = 'lm') +
  #geom_bar(stat="identity",position="dodge")+
  geom_errorbar(aes(ymax=Chi_Obs_Mean+Chi_Obs_SD,ymin=Chi_Obs_Mean - Chi_Obs_SD),position=position_dodge(0.9),width=0.01,color = "blue")+
  #scale_fill_brewer(palette = "Set1")+
  xlab(expression('MI'*~'('*'unitless'*')')) + ylab(expression(~chi)) +
  geom_point(aes(x=MI,y=Chi_Own_Mean),color="red") + 
  geom_smooth(aes(y=Chi_Own_Mean,x=MI, col='red'),method = 'lm') +
  theme_bw()+
  theme(legend.position="none",axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),
        axis.text=element_text(size=rel(1.5)),panel.grid.major=element_line(color="white"),panel.grid.minor=element_line(color="white")) +
  annotate("text", x=1, y=0.9, label="(a)",size=6.5)
f1

lmadata<-read.csv("predicted_LMA.csv")
lmadata <- lmadata %>%
  dplyr::left_join(NECT_mi, by = c("map" = "MAP..mm."))

f2<- ggplot(data=lmadata,aes(x=MI..unitless.,y=lmaobsmean))+geom_point(color = "blue") + 
  geom_smooth(aes(y=lmaobsmean,x=MI..unitless.),method = 'lm') +
  geom_errorbar(aes(ymax=lmaobsmean+lma.sd,ymin=lmaobsmean-lma.sd),position=position_dodge(0.9),width=0.01,color = "blue")+
  xlab(expression('MI'*~'('*'unitless'*')')) + ylab(expression('lnLMA'*~'('*g*'/'*m^2*')')) +
  geom_point(aes(x=MI..unitless.,y=lnprelma),color="red") + 
  geom_smooth(aes(y=lnprelma,x=MI..unitless., col='red'),method = 'lm') +
  theme_bw()+
  theme(legend.position="none",axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),
        axis.text=element_text(size=rel(1.5)),panel.grid.major=element_line(color="white"),panel.grid.minor=element_line(color="white"))+
  annotate("text", x=1, y=5, label="(c)",size=6.5)
f2

##################################Narea
datapre<-read.csv("predicted_Narea_NewVcmax25.csv")
datapre <- datapre %>%
  dplyr::left_join(NECT_mi, by = c("MAP_mm" = "MAP..mm."))

n4<- ggplot(data=datapre,aes(x=MI..unitless.,y=lnNarea_gm2_mean.))+geom_point(color = "blue") + 
  geom_smooth(aes(y=lnNarea_gm2_mean.,x=MI..unitless.),method = 'lm') +
  geom_errorbar(aes(ymax=lnNarea_gm2_mean.+lnNarea_sd,ymin=lnNarea_gm2_mean.-lnNarea_sd),position=position_dodge(0.9),width=0.01,color = "blue")+
  xlab(expression('MI'*~'('*'unitless'*')')) + ylab(expression('lnNarea'*~'('*g*'/'*m^2*')')) +
  geom_point(aes(x=MI..unitless.,y=ln_prenaera),color="red") + 
  geom_smooth(aes(y=ln_prenaera,x=MI..unitless., col='red'),method = 'lm') +
  theme_bw()+
  theme(legend.position="none",axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),
        axis.text=element_text(size=rel(1.5)),panel.grid.major=element_line(color="white"),panel.grid.minor=element_line(color="white"))+
  annotate("text", x=1, y=1.5, label="(d)",size=6.5)
n4
vcmax12again <- ggplot(site_summary2, aes(x = MI..unitless.)) +
  geom_point(aes(y = vcmax_Own_Mean, color = "vcmax_Own"), size = 1) +
  geom_smooth(aes(y = vcmax_Own_Mean, color = "vcmax_Own"), method = "lm", se = TRUE, size = 1) +
  geom_point(aes(y = Obs_vcmax, color = "Obs_vcmax")) +
  # Add smooth curves for other vcmax types
  geom_smooth(aes(y = Obs_vcmax, color = "Obs_vcmax"), method = "lm", se = TRUE, size = 1) +
  ylab(expression(V[cmax]~'('*mu*'mol'*m^-2*s^-1*')')) +
  xlab(expression('MI'*~'('*'unitless'*')')) +
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
  annotate("text", x = 1, y = 40, label = "(b)", size = 6.5)
vcmax12again
plot_grid(f1,vcmax12again,f2,n4,label_x = 0.2,nrow = 2)


data_combined$MI..unitless. <- NECT_sample_climate$MI..unitless.

# prepare data
narea_dat <- data_combined %>%
  dplyr::select(Narea, LMA, vcmax25_Own, MI..unitless.) %>%
  na.omit()

# fit model
lnneramodel <- lm(Narea ~ LMA + vcmax25_Own - 1, data = narea_dat)

# point prediction
narea_dat$pred_Narea <- predict(lnneramodel, newdata = narea_dat)
narea_dat$ln_obs <- log(narea_dat$Narea)
narea_dat$ln_pred <- log(narea_dat$pred_Narea)

# Monte Carlo
set.seed(123)
nsim <- 1000
beta_hat <- coef(lnneramodel)
V_beta <- vcov(lnneramodel)
beta_sim <- MASS::mvrnorm(n = nsim, mu = beta_hat, Sigma = V_beta)

X <- model.matrix(~ LMA + vcmax25_Own - 1, data = narea_dat)
pred_mat <- X %*% t(beta_sim)

# avoid log problems
pred_mat[pred_mat <= 0] <- NA
ln_pred_mat <- log(pred_mat)

# summarize along MI
mi <- narea_dat$MI..unitless.
nbins <- 8
breaks <- seq(min(mi, na.rm = TRUE), max(mi, na.rm = TRUE), length.out = nbins + 1)
bin_id <- cut(mi, breaks = breaks, include.lowest = TRUE, labels = FALSE)

band_list <- vector("list", nbins)

for (b in 1:nbins) {
  idx <- which(bin_id == b)
  if (length(idx) == 0) next
  
  vals <- as.vector(ln_pred_mat[idx, , drop = FALSE])
  
  band_list[[b]] <- data.frame(
    MI_mid = mean(mi[idx], na.rm = TRUE),
    pred_median = quantile(vals, 0.5, na.rm = TRUE),
    pred_low = quantile(vals, 0.025, na.rm = TRUE),
    pred_high = quantile(vals, 0.975, na.rm = TRUE)
  )
}

band_narea <- bind_rows(band_list)

# plot
p_narea_unc_MI <- ggplot() +
  geom_point(data = narea_dat,
             aes(x = MI..unitless., y = ln_obs),
             color = "blue", alpha = 0.7) +
  geom_ribbon(data = band_narea,
              aes(x = MI_mid, ymin = pred_low, ymax = pred_high),
              fill = "red", alpha = 0.2) +
  geom_line(data = band_narea,
            aes(x = MI_mid, y = pred_median),
            color = "red", linewidth = 1.2) +
  xlab(expression(MI~(unitless))) +
  ylab(expression(ln~N[area]~(g~m^-2))) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text = element_text(size = rel(1.5)),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_line(color = "white")
  )

p_narea_unc_MI
library(MASS)
library(dplyr)

mc_site_uncertainty <- function(model, data, site_var, nsim = 1000, seed = 123,
                                log_response = FALSE) {
  set.seed(seed)
  vars_needed <- all.vars(formula(model))
  vars_needed <- unique(c(vars_needed, site_var))
  dat <- data[, vars_needed, drop = FALSE]
  dat <- na.omit(dat)
  X <- model.matrix(delete.response(terms(model)), data = dat)
  beta_hat <- coef(model)
  V_beta   <- vcov(model)
  beta_sim <- MASS::mvrnorm(n = nsim, mu = beta_hat, Sigma = V_beta)
  pred_mat <- X %*% t(beta_sim)   # n_obs × nsim
  if (log_response) {
    pred_mat[pred_mat <= 0] <- NA
    pred_mat <- log(pred_mat)
  }
  pred_point <- predict(model, newdata = dat)
  if (log_response) {
    pred_point[pred_point <= 0] <- NA
    pred_point <- log(pred_point)
  }
  
  dat$pred_point <- pred_point
  dat[[site_var]] <- as.character(dat[[site_var]])
  
  site_ids <- unique(dat[[site_var]])
  
  out <- lapply(site_ids, function(s) {
    idx <- which(dat[[site_var]] == s)
    
    sim_means <- apply(pred_mat[idx, , drop = FALSE], 2, mean, na.rm = TRUE)
    
    data.frame(
      site = s,
      pred_mean = mean(dat$pred_point[idx], na.rm = TRUE),
      pred_low  = quantile(sim_means, 0.025, na.rm = TRUE),
      pred_high = quantile(sim_means, 0.975, na.rm = TRUE)
    )
  })
  
  bind_rows(out)
}
data_combined2<-data_combined
data_combined2$Site.ID<-NECT_sample_climate$Site.ID
narea_unc <- mc_site_uncertainty(
  model = lnneramodel,
  data = data_combined2,   
  site_var = "Site.ID",
  nsim = 1000,
  seed = 123,
  log_response = FALSE
)
narea_unc$ln_pred_mean <- log(narea_unc$pred_mean)
narea_unc$ln_pred_low  <- log(narea_unc$pred_low)
narea_unc$ln_pred_high <- log(narea_unc$pred_high)
datapre$Site.ID <- as.character(datapre$Site.ID)
narea_unc$site  <- as.character(narea_unc$site)
datapre2 <- datapre %>%
  left_join(narea_unc, by = c("Site.ID" = "site"))
n4_unc <- ggplot(data = datapre2, aes(x = MI..unitless.)) +
  geom_point(aes(y = lnNarea_gm2_mean.), color = "blue") +
  geom_smooth(aes(y = lnNarea_gm2_mean.), method = "lm", color = "blue") +
  geom_errorbar(aes(ymin = lnNarea_gm2_mean. - lnNarea_sd,
                    ymax = lnNarea_gm2_mean. + lnNarea_sd),
                width = 0.01, color = "blue") +
  
  geom_point(aes(y = ln_pred_mean), color = "red") +
  geom_smooth(aes(y = ln_pred_mean), method = "lm", color = "red") +
  geom_errorbar(aes(ymin = ln_pred_low, ymax = ln_pred_high),
                width = 0.01, color = "red") +
  
  xlab(expression(MI~(unitless))) +
  ylab(expression(ln~N[area]~(g~m^-2))) +
  theme_bw() + theme(
    legend.position = "none",
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text = element_text(size = rel(1.5)),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_line(color = "white")
  )
n4_unc
##########################################################
lmobservedx<-lm(Chi_Obs_Mean~MI,data=site_summary)
summary(lmobservedx)
lmpredicedx<-lm(Chi_Own_Mean~MI,data=site_summary)
summary(lmpredicedx)
lmpredicedxconstant<-lm(Chi_Constant_Mean~MI,data=site_summary)
summary(lmpredicedxconstant)
coefs <- summary(lmpredicedxconstant)$coefficients
round(coefs, 2)
lmpredicedxAlienor<-lm(Chi_Alienor_Mean~MI,data=site_summary)
summary(lmpredicedxAlienor)
lmlma<-lm(lmaobsmean~MI..unitless.,data=lmadata)
summary(lmlma)
lmlmapre<-lm(lnprelma~MI..unitless.,data=lmadata)
summary(lmlmapre)
lmobservednarea<-lm(lnNarea_gm2_mean.~MI..unitless.,data=datapre)
summary(lmobservednarea)
lmpredictednarea<-lm(ln_prenaera~MI..unitless.,data=datapre)
summary(lmpredictednarea)
#####################variation in the middle part of the transect (presumably related to temperature)
############################################################
# 1. Fit original MI relationships and calculate residuals
############################################################

# chi observed
lmobservedx <- lm(Chi_Obs_Mean ~ MI, data = site_summary)
site_summary$abs_resid_chi_obs <- abs(resid(lmobservedx))

# chi predicted
lmpredicedx <- lm(Chi_Own_Mean ~ MI, data = site_summary)
site_summary$abs_resid_chi_pred <- abs(resid(lmpredicedx))

# LMA observed
lmlma <- lm(lmaobsmean ~ MI..unitless., data = lmadata)
lmadata$abs_resid_lma_obs <- abs(resid(lmlma))

# LMA predicted
lmlmapre <- lm(lnprelma ~ MI..unitless., data = lmadata)
lmadata$abs_resid_lma_pred <- abs(resid(lmlmapre))

# Narea observed
lmobservednarea <- lm(lnNarea_gm2_mean. ~ MI..unitless., data = datapre)
datapre$abs_resid_narea_obs <- abs(resid(lmobservednarea))

# Narea predicted
lmpredictednarea <- lm(ln_prenaera ~ MI..unitless., data = datapre)
datapre$abs_resid_narea_pred <- abs(resid(lmpredictednarea))

# Vcmax observed
lmobvcmax <- lm(Obs_vcmax ~ MI..unitless., data = site_summary2)
site_summary2$abs_resid_vcmax_obs <- abs(resid(lmobvcmax))

# Vcmax predicted
lmorevcmax <- lm(vcmax_Own_Mean ~ MI..unitless., data = site_summary2)
site_summary2$abs_resid_vcmax_pred <- abs(resid(lmorevcmax))
# 1. make one site-level temperature table
temp_site <- NECT_sample_climate %>%
  dplyr::select(Site.ID, tg_chelsa_2006) %>%
  distinct() %>%
  group_by(Site.ID) %>%
  summarise(tg_chelsa_2006 = mean(tg_chelsa_2006, na.rm = TRUE), .groups = "drop")
site_summary <- site_summary %>%
  left_join(temp_site,  by = c("Site_ID" = "Site.ID"))
lmadata <- lmadata %>%
  left_join(temp_site, by = c("X" = "Site.ID"))
datapre$Site.ID <- as.character(datapre$Site.ID)
temp_site$Site.ID  <- as.character(temp_site$Site.ID)
datapre <- datapre %>%
  left_join(temp_site, by = c("Site.ID" = "Site.ID"))
site_summary2$Site_ID <- as.character(site_summary2$Site_ID)
site_summary2 <- site_summary2 %>%
  left_join(temp_site, by = c("Site_ID" = "Site.ID"))
############################################################
# 2. Define middle part of the transect
############################################################
# Adjust these if needed after checking your MI range.
# This version uses MI between 0.4 and 0.7.

mid_chi <- site_summary %>%
  filter(MI > 0.4, MI < 0.7)

mid_lma <- lmadata %>%
  filter(MI..unitless. > 0.4, MI..unitless. < 0.7)

mid_narea <- datapre %>%
  filter(MI..unitless. > 0.4, MI..unitless. < 0.7)

mid_vcmax <- site_summary2 %>%
  filter(MI..unitless. > 0.4, MI..unitless. < 0.7)

############################################################
# 3. Test whether residual spread is related to temperature
############################################################

# observed
mod_chi_obs_temp   <- lm(abs_resid_chi_obs   ~ tg_chelsa_2006,  data = mid_chi)
mod_lma_obs_temp   <- lm(abs_resid_lma_obs   ~ tg_chelsa_2006,       data = mid_lma)
mod_narea_obs_temp <- lm(abs_resid_narea_obs ~ tg_chelsa_2006,       data = mid_narea)
mod_vcmax_obs_temp <- lm(abs_resid_vcmax_obs ~ tg_chelsa_2006, data = mid_vcmax)

# predicted
mod_chi_pred_temp   <- lm(abs_resid_chi_pred   ~ tg_chelsa_2006,  data = mid_chi)
mod_lma_pred_temp   <- lm(abs_resid_lma_pred   ~ tg_chelsa_2006,       data = mid_lma)
mod_narea_pred_temp <- lm(abs_resid_narea_pred ~ tg_chelsa_2006,       data = mid_narea)
mod_vcmax_pred_temp <- lm(abs_resid_vcmax_pred ~ tg_chelsa_2006, data = mid_vcmax)

############################################################
# 4. Print model summaries
############################################################
summary(mod_chi_obs_temp)
summary(mod_lma_obs_temp)
summary(mod_narea_obs_temp)
summary(mod_vcmax_obs_temp)

summary(mod_chi_pred_temp)
summary(mod_lma_pred_temp)
summary(mod_narea_pred_temp)
summary(mod_vcmax_pred_temp)

############################################################
# 5. Make one combined dataframe for plotting
############################################################

plot_chi_obs <- data.frame(
  variable = "chi_obs",
  temperature = mid_chi$tg_chelsa_2006,
  abs_resid = mid_chi$abs_resid_chi_obs
)

plot_chi_pred <- data.frame(
  variable = "chi_pred",
  temperature = mid_chi$tg_chelsa_2006,
  abs_resid = mid_chi$abs_resid_chi_pred
)

plot_lma_obs <- data.frame(
  variable = "lma_obs",
  temperature = mid_lma$tg_chelsa_2006,
  abs_resid = mid_lma$abs_resid_lma_obs
)

plot_lma_pred <- data.frame(
  variable = "lma_pred",
  temperature = mid_lma$tg_chelsa_2006,
  abs_resid = mid_lma$abs_resid_lma_pred
)

plot_narea_obs <- data.frame(
  variable = "narea_obs",
  temperature = mid_narea$tg_chelsa_2006,
  abs_resid = mid_narea$abs_resid_narea_obs
)

plot_narea_pred <- data.frame(
  variable = "narea_pred",
  temperature = mid_narea$tg_chelsa_2006,
  abs_resid = mid_narea$abs_resid_narea_pred
)

plot_vcmax_obs <- data.frame(
  variable = "vcmax_obs",
  temperature = mid_vcmax$tg_chelsa_2006,
  abs_resid = mid_vcmax$abs_resid_vcmax_obs
)

plot_vcmax_pred <- data.frame(
  variable = "vcmax_pred",
  temperature = mid_vcmax$tg_chelsa_2006,
  abs_resid = mid_vcmax$abs_resid_vcmax_pred
)

plot_all <- bind_rows(
  plot_chi_obs, plot_chi_pred,
  plot_lma_obs, plot_lma_pred,
  plot_narea_obs, plot_narea_pred,
  plot_vcmax_obs, plot_vcmax_pred
)

############################################################
# 6. Plot: does spread increase with temperature?
############################################################

p_temp_resid <- ggplot(plot_all, aes(x = temperature, y = abs_resid)) +
  geom_point(color = "steelblue", alpha = 0.75, size = 2.5) +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  facet_wrap(~variable, scales = "free_y") +
  xlab("Temperature") +
  ylab("Absolute residual from MI relationship") +
  theme_bw() +
  theme(
    strip.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 11),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

print(p_temp_resid)

############################################################
# 8. Make a summary table of slopes and p-values
############################################################

extract_lm_stats <- function(model, varname) {
  s <- summary(model)
  data.frame(
    variable = varname,
    slope = s$coefficients[2, 1],
    SE = s$coefficients[2, 2],
    p_value = s$coefficients[2, 4],
    R2 = s$r.squared
  )
}

result_table <- bind_rows(
  extract_lm_stats(mod_chi_obs_temp,   "chi_obs"),
  extract_lm_stats(mod_chi_pred_temp,  "chi_pred"),
  extract_lm_stats(mod_lma_obs_temp,   "lma_obs"),
  extract_lm_stats(mod_lma_pred_temp,  "lma_pred"),
  extract_lm_stats(mod_narea_obs_temp, "narea_obs"),
  extract_lm_stats(mod_narea_pred_temp,"narea_pred"),
  extract_lm_stats(mod_vcmax_obs_temp, "vcmax_obs"),
  extract_lm_stats(mod_vcmax_pred_temp,"vcmax_pred")
)

print(result_table)

temp_site <- NECT_sample_climate %>%
  dplyr::select(Site.ID, tg_chelsa_2006,MI..unitless.,LMA..kg.m2.,Narea..g.m2.,d15N.14N) %>%
  distinct()
temp_sitemi <- temp_site %>%
  filter(MI..unitless. > 0.4, MI..unitless. < 0.7)
summary(lm(log(LMA..kg.m2.* 1000) ~ tg_chelsa_2006, data = temp_sitemi))
summary(lm(log(Narea..g.m2.) ~ tg_chelsa_2006, data = temp_sitemi))
summary(lm(d15N.14N ~ tg_chelsa_2006, data = temp_sitemi))
chicombine2mi <- chicombine2 %>%
  filter(MI..unitless. > 0.4, MI..unitless. < 0.7)
summary(lm(Chi_Obs ~ tg_chelsa_2006, data = chicombine2mi))
library(ggplot2)

p1 <- ggplot(temp_sitemi, aes(x = tg_chelsa_2006, y = log(LMA..kg.m2.*1000))) +
  geom_point(color = "blue", alpha = 0.6) +
  geom_smooth(method = "lm", color = "red") +
  annotate("text", x = Inf, y = Inf, label = "p < 0.001",
           hjust = 1.1, vjust = 1.5, size = 4) +ylab(expression(ln~LMA~(g~m^-2))) +xlab(expression('Tg'*~'('*degree*'C'*')')) + 
  theme_bw() + theme(
    legend.position = "none",
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text = element_text(size = rel(1.5)),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_line(color = "white")
  )+ annotate("text", x=6.7, y=5.2, label="(b)",size=6.5)

p2 <- ggplot(temp_sitemi, aes(x = tg_chelsa_2006, y = log(Narea..g.m2.))) +
  geom_point(color = "blue", alpha = 0.6) +
  geom_smooth(method = "lm", color = "red") +
  annotate("text", x = Inf, y = Inf, label = "p = 0.004",
           hjust = 1.1, vjust = 1.5, size = 4)  +ylab(expression(ln~N[area]~(g~m^-2))) +xlab(expression('Tg'*~'('*degree*'C'*')')) + 
  theme_bw() + theme(
    legend.position = "none",
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text = element_text(size = rel(1.5)),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_line(color = "white")
  )+ annotate("text", x=6.7, y=1.3, label="(c)",size=6.5)

p3 <- ggplot(temp_sitemi, aes(x = tg_chelsa_2006, y = d15N.14N)) +
  geom_point(color = "blue", alpha = 0.6) +
  geom_smooth(method = "lm", color = "red") +
  annotate("text", x = Inf, y = Inf, label = "p = 0.044",
           hjust = 1.1, vjust = 1.5, size = 4) + ylab(expression('Plant '*delta^15*'N')) +xlab(expression('Tg'*~'('*degree*'C'*')')) + 
  theme_bw() + theme(
    legend.position = "none",
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text = element_text(size = rel(1.5)),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_line(color = "white")
  )+ annotate("text", x=6.7, y=14, label="(d)",size=6.5)

p4 <- ggplot(chicombine2mi, aes(x = tg_chelsa_2006, y = Chi_Obs)) +
  geom_point(color = "blue", alpha = 0.6) +
  geom_smooth(method = "lm", color = "red") +
  annotate("text", x = Inf, y = Inf, label = "p < 0.001",
           hjust = 1.1, vjust = 1.5, size = 4)+ ylab(expression(~chi)) +xlab(expression('Tg'*~'('*degree*'C'*')')) + 
  theme_bw() + theme(
    legend.position = "none",
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text = element_text(size = rel(1.5)),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_line(color = "white")
  )+ annotate("text", x=6.7, y=0.75, label="(a)",size=6.5)

library(patchwork)
p4+p1 + p2 + p3

chi_dat <- chicombine2 %>%
  dplyr::select(Site_ID, Chi_Obs, MI..unitless., tg_chelsa_2006) %>%
  na.omit()

# sample-level MI model
mod_chi_mi <- lm(Chi_Obs ~ MI..unitless., data = chi_dat)

# residuals
chi_dat$resid_chi <- resid(mod_chi_mi)
chi_dat$abs_resid_chi <- abs(chi_dat$resid_chi)

# middle transect
chi_mid <- chi_dat %>%
  filter(MI..unitless. > 0.4, MI..unitless. < 0.7)

# test whether spread is related to temperature
mod_chi_temp <- lm(abs_resid_chi ~ tg_chelsa_2006, data = chi_mid)
summary(mod_chi_temp)##not siginifancet
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)     0.0519540  0.0135169   3.844 0.000185 ***
#   tg_chelsa_2006 -0.0008629  0.0014947  -0.577 0.564708    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.02978 on 136 degrees of freedom
# Multiple R-squared:  0.002444,	Adjusted R-squared:  -0.004891 
# F-statistic: 0.3332 on 1 and 136 DF,  p-value: 0.5647
###########################################################Fig 3
###################################seasonal fapar and LAI--Boya version
site_mean<-read.csv("/Users/vv934811@reading.ac.uk/Documents/work2024/faparmax_LAI/Pre_Obs_fapar_LAI_Boyaseasonal.csv")
site_mean <- site_mean %>%
  dplyr::left_join(NECT_mi, by = c("MAP" = "MAP..mm."))


pb1 <- ggplot(data=site_mean,aes(x=MI..unitless.,y=fapar_obs))+geom_point(color = "blue",size = 2.8) + 
  geom_smooth(aes(y=fapar_obs,x=MI..unitless.),method = 'lm') +
  xlab(expression('MI'*~'('*'unitless'*')')) + ylab(expression('fAPAR')) +
  geom_point(aes(x=MI..unitless.,y=fapar_pre),color="red",size = 2.8) + 
  geom_smooth(aes(y=fapar_pre,x=MI..unitless., col='red'),method = 'lm') +
  theme_bw()+
  theme(legend.position="none",axis.title.x=element_blank(),axis.title.y=element_text(size=29),
        axis.text=element_text(size=rel(2.5)),panel.grid.major=element_line(color="white"),panel.grid.minor=element_line(color="white"))
pb1
#######################################LAI
pb2 <- ggplot(data=site_mean,aes(x=MI..unitless.,y=LAI_obs))+geom_point(color = "blue",size = 2.8) + 
  geom_smooth(aes(y=LAI_obs,x=MI..unitless.),method = 'lm') +
  xlab(expression('MI'*~'('*'unitless'*')')) + ylab(expression('LAI')) +
  geom_point(aes(x=MI..unitless.,y=LAI_pre),color="red",,size = 2.8) + 
  geom_smooth(aes(y=LAI_pre,x=MI..unitless., col='red'),method = 'lm') +scale_y_continuous(breaks = c(0.3, 0.6, 0.9, 1.2))+
  theme_bw()+
  theme(legend.position="none",axis.title.x=element_blank(),axis.title.y=element_text(size=29),
        axis.text=element_text(size=rel(2.5)),panel.grid.major=element_line(color="white"),panel.grid.minor=element_line(color="white"))

pb2
lmfaparobs<-lm(fapar_obs~MI..unitless.,data=site_mean)
summary(lmfaparobs)
lmfaparpre<-lm(fapar_pre~MI..unitless.,data=site_mean)
summary(lmfaparpre)
lmlaiobs<-lm(LAI_obs~MI..unitless.,data=site_mean)
summary(lmlaiobs)
lmlai<-lm(LAI_pre~MI..unitless.,data=site_mean)
summary(lmlai)
#############################GPP gosif
model_season_df<-read.csv("gppobs_gCm2y_NECT_pmodel_FLUX_X_GOSIF_PML15.csv")

gpp_mean <- model_season_df %>%
  group_by(site) %>%
  summarize(across(where(is.numeric), mean, na.rm = TRUE))
gpp_mean <- gpp_mean %>%
  dplyr::left_join(NECT_mi, by = c("Precipitation" = "MAP..mm."))


pb3 <- ggplot(data=gpp_mean,aes(x=MI..unitless.,y=gpp))+geom_point(color="red",size = 2.8) + 
  geom_smooth(method = 'lm', formula = y ~ x, aes(color='red')) +
  xlab(expression('MI'*~'('*'unitless'*')')) + ylab(expression('GPP'*~'('*'gC'*~m^-1*d^-1*')')) +
  geom_point(aes(x=MI..unitless.,y=gpp_gosif),color="blue",size = 2.8) + 
  geom_smooth(aes(x=MI..unitless.,y=gpp_gosif),color="blue",method = 'lm') +
  theme_bw()+
  theme(legend.position="none",axis.title.x=element_blank(),axis.title.y=element_text(size=29),
        axis.text=element_text(size=rel(2.5)),panel.grid.major=element_line(color="white"),panel.grid.minor=element_line(color="white"))
pb3


p4 <- ggplot(data=gpp_mean,aes(x=MI..unitless.,y=gpp*BPE))+geom_point(color="red",size = 2.8) + 
  geom_smooth(method = 'lm', formula = y ~ x, aes(color='red')) +
  xlab(expression('MI'*~'('*'unitless'*')')) +  ylab(expression('NPP'*~'('*'gC'*~m^-1*d^-1*')'))  +
  geom_point(aes(x=MI..unitless.,y=gpp_gosif*BPE),color="blue",size = 2.8) + 
  geom_smooth(aes(x=MI..unitless.,y=gpp_gosif*BPE),color="blue",method = 'lm') +
  ylim(-1,2)+
  theme_bw()+
  theme(legend.position="none",axis.title.x=element_blank(),axis.title.y=element_text(size=29),
        axis.text=element_text(size=rel(2.5)),panel.grid.major=element_line(color="white"),panel.grid.minor=element_line(color="white"))
p4
lmgpp<-lm(gpp~MI..unitless.,data=gpp_mean)
summary(lmgpp) 
lmgppobx<-lm(gpp_gosif~MI..unitless.,data=gpp_mean)
summary(lmgppobx) 

lmnpp<-lm((gpp*BPE)~MI..unitless.,data=gpp_mean)
summary(lmnpp)
lmnppobs<-lm((gpp_gosif*BPE)~MI..unitless.,data=gpp_mean)
summary(lmnppobs)

RSimprove2<-read.csv("Climate_R_S_NECT_predictRS_biomass_BAADb0c0.csv")

pbRSgroup <- ggplot(RSimprove2, aes(x = MI..unitless.)) +
  geom_point(aes(y = log(R_S)), color = "blue",size = 2.8) +
  geom_point(aes(y = pred_lnRSnew), color = "red",size = 2.8) +
  geom_smooth(aes(y = log(R_S), group = sitegroup), method = 'lm',color = "blue", formula = y ~ x, se = TRUE) +
  geom_smooth(aes(y = pred_lnRSnew, group = sitegroup), method = 'lm',color = "red", formula = y ~ x, se = TRUE) +
  xlab(expression('MI'*~'('*'unitless'*')')) +
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

pbAGBgroup <- ggplot(RSimprove2, aes(x = MI..unitless.)) +
  geom_point(aes(y = log(AGB_gC.m.2)), color = "blue",size = 2.8) +
  geom_point(aes(y = log(preAGB_new)), color = "red",size = 2.8) +
  geom_smooth(aes(y = log(AGB_gC.m.2), group = sitegroup), method = 'lm',color = "blue", formula = y ~ x, se = TRUE) +
  geom_smooth(aes(y = log(preAGB_new), group = sitegroup), method = 'lm',color = "red", formula = y ~ x, se = TRUE) +
  
  # Labels and theme
  xlab(expression('MI'*~'('*'unitless'*')')) +
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


pbBBgroup <- ggplot(RSimprove2, aes(x = MI..unitless.)) +
  # Points: Black for observations, red for predictions
  geom_point(aes(y = log(BGB_gC.m.2)), color = "blue",size = 2.8) +
  geom_point(aes(y = log(prebelow_new)), color = "red",size = 2.8) +
  geom_smooth(aes(y = log(BGB_gC.m.2), group = sitegroup), method = 'lm',color = "blue", formula = y ~ x, se = TRUE) +
  geom_smooth(aes(y = log(prebelow_new), group = sitegroup), method = 'lm',color = "red", formula = y ~ x, se = TRUE) +
  xlab(expression('MI'*~'('*'unitless'*')')) +
  ylab(expression('lnBGB'*~'('*'gC'*~m^-2*')')) +scale_y_continuous(breaks = c(0, 2, 4, 6, 8))+
  theme_bw() +
  theme(
    legend.position = "right",  # Keep legend for sitegroup
    axis.title.x =element_blank(),
    axis.title.y = element_text(size = 29),
    axis.text = element_text(size = rel(2.5)),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_line(color = "white")
  ) 
pbBBgroup
lowBGB <- RSimprove2 %>% 
  filter(MI..unitless. < 0.6, !is.na(MI..unitless.))
BGBmodelobs_low<-lm(log(BGB_gC.m.2)~MI..unitless.,data=lowBGB)
summary(BGBmodelobs_low)
BGBmodelpre_low<-lm(log(prebelow_new)~MI..unitless.,data=lowBGB)
summary(BGBmodelpre_low)
AGBmodelobs_low<-lm(log(AGB_gC.m.2)~MI..unitless.,data=lowBGB)
summary(AGBmodelobs_low)
AGBmodelpre_low<-lm(log(preAGB_new)~MI..unitless.,data=lowBGB)
summary(AGBmodelpre_low)
RSmodelobs_low<-lm(log(R_S)~MI..unitless.,data=lowBGB)
summary(RSmodelobs_low)
RSmodelobs_lowpre<-lm(pred_lnRSnew~MI..unitless.,data=lowBGB)
summary(RSmodelobs_lowpre)
highBGB <- RSimprove2 %>% 
  filter(MI..unitless. > 0.6, !is.na(MI..unitless.))
BGBmodelobs_high<-lm(log(BGB_gC.m.2)~MI..unitless.,data=highBGB)
summary(BGBmodelobs_high)
BGBmodelpre_high<-lm(log(prebelow_new)~MI..unitless.,data=highBGB)
summary(BGBmodelpre_high)
AGBmodelobs_high<-lm(log(AGB_gC.m.2)~MI..unitless.,data=highBGB)
summary(AGBmodelobs_high)
AGBmodelpre_high<-lm(log(preAGB_new)~MI..unitless.,data=highBGB)
summary(AGBmodelpre_high)
RSmodelobs_high<-lm(log(R_S)~MI..unitless.,data=highBGB)
summary(RSmodelobs_high)
RSmodelobs_highpre<-lm(pred_lnRSnew~MI..unitless.,data=highBGB)
summary(RSmodelobs_highpre)
####################################################################
############################15N not change the original data and prediction
#####################################################################
ndata<-read.csv("Hard Traits15N(with_predictingN)_ftmean.csv")
nitrogendata<-read.csv("Hard Traits15N(repeatftemp_frunoff)_ftmean.csv")
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
ndata1 <- ndata1 %>%
  dplyr::left_join(NECT_mi, by = c("map" = "MAP..mm."))

istope1<- ggplot(data=ndata1,aes(x=MI..unitless.,y=nobsmean))+geom_point(color = "blue") + 
  geom_smooth(aes(y=nobsmean,x=MI..unitless.),method = 'lm') +
  geom_errorbar(aes(ymax=nobsmean+n.sd,ymin=nobsmean-n.sd),position=position_dodge(0.9),color = "blue",width=0.01)+
  xlab(expression('MI'*~'('*'unitless'*')')) +  ylab(expression('Plant '*delta^15*'N'))+
  geom_point(aes(x=MI..unitless.,y=lnpren),color="red") + 
  geom_smooth(aes(y=lnpren,x=MI..unitless., col='red'),method = 'lm') +
  theme_bw()+
  theme(legend.position="none",axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),
        axis.text=element_text(size=rel(1.5)),panel.grid.major=element_line(color="white"),panel.grid.minor=element_line(color="white"))+
  annotate("text", x=1, y=13, label="(a)",size=6.5)
istope1
lm15N<-lm(nobsmean~MI..unitless.,data=ndata1)
summary(lm15N)
lm15Npre<-lm(lnpren~MI..unitless.,data=ndata1)
summary(lm15Npre)
fTemp<-nitrogendata$ftemp_new
fR<-nitrogendata$frunoff_new
##nitrogenisotope_model<--7.6265 +(10.4682/(1+1.1405*(fR/fTemp)))
fgas<-1/(1+1.1405*(fR/fTemp))
fgasdata<-data.frame(cbind(fgas,ndata$MAP))
fgasdata <- na.omit(fgasdata)
fgas.sd<-aggregate(fgasdata$fgas,by=list(fgasdata$V2),FUN=sd)
fgas_mean<-aggregate(fgasdata$fgas,by=list(fgasdata$V2),FUN=mean)
relationshipdata1<-data.frame(cbind(fgas_mean,fgas.sd$x))
relationshipdata1[relationshipdata1<0]<-0
relationshipdata1 <- relationshipdata1 %>%
  dplyr::left_join(NECT_mi, by = c("Group.1" = "MAP..mm."))

mfas2<-ggplot(relationshipdata1, aes(x=MI..unitless., y=x)) + geom_point(color='red') + 
  geom_smooth(method = 'loess', formula = y ~ x, aes(color='red')) + 
  geom_errorbar(aes(ymax=x+fgas.sd$x,ymin=x-fgas.sd$x),position=position_dodge(0.9),color='red',width=0.01)+
  ylab(expression(f[gas])) +xlab(expression('MI'*'('*'unitless'*')'))+
  theme_bw()+theme(legend.position='none',axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),
                   axis.text=element_text(size=rel(1.5)),panel.grid.major=element_line(color="white"),panel.grid.minor=element_line(color="white"))+
  annotate("text", x=1, y=1, label="(b)",size=6.5)
mfas2
plot_grid(istope1,mfas2,label_x = 0.2,ncol = 2)
#################################################figS1
library(terra)
library(sf)
library(rnaturalearth)
library(tidyterra)
library(ggplot2)
library(viridis)

sta <- NECT_sample_climate %>%
  dplyr::select(Longitude, Latitude, MAP..mm.) %>%
  distinct(Longitude, Latitude, .keep_all = TRUE)

sta_33 <- sta[1:33, ]
sites_sf <- st_as_sf(sta_33, coords = c("Longitude","Latitude"), crs = 4326)

r_stack <- rast('/Users/vv934811@reading.ac.uk/Documents/JIaze_species_Data/Bioclimatic_variables/lnMI_01res_28.tif')             # SpatRaster with 12 layers
r_true <- exp(r_stack)
plot(r_true)
# (If you actually want ANNUAL TOTAL instead of mean, use this one line:)
# pr_mean <- app(r_stack, sum, na.rm = TRUE)

# 1) Bring in China’s boundary
china_sf <- ne_countries(country = "China", scale = "large", returnclass = "sf")

# 2) Match CRS to your raster and convert to terra vect
china_sf  <- st_transform(china_sf, crs(r_true))
china_vect <- vect(china_sf)

# 3) Clip/mask the raster to China → everything outside becomes NA (white)
pr_china <- r_true |> crop(china_vect) |> mask(china_vect)

# (Optional) 4) Crop to the same NE-China window as your example map
#    Adjust these bounds if you want it tighter/looser.
#    xmin, xmax, ymin, ymax  (lon/lat in degrees)
bb_ne <- ext(111, 140, 30, 55)
pr_ne <- crop(pr_china, bb_ne)

# 5) Plot (choose one object to plot: pr_china or pr_ne)

# Convert cropped raster to data frame
pr_df <- as.data.frame(pr_ne, xy = TRUE, na.rm = FALSE)
colnames(pr_df) <- c("x", "y", "MI")
# Plot
ggplot(pr_df) +
  geom_raster(aes(x = x, y = y, fill = MI)) +
  scale_fill_viridis_c(name = "MI (unitless)", na.value = "white")+
  geom_sf(data = sites_sf, shape = 21, size = 2.6, color = "black", fill = "black", stroke = 0.5) +
  xlab(expression('Longitude')) +  ylab(expression('Latitude'))  +
  coord_sf(expand = FALSE) +
  theme_bw()+
  theme(legend.position="right",axis.title.x=element_text(size=15),axis.title.y=element_text(size=15),
        axis.text=element_text(size=rel(1.5)),panel.grid.major=element_line(color="white"),panel.grid.minor=element_line(color="white"))

######################################################change scale color
# Crop to 40–50N latitude (and keep longitudes as before)
# Crop to 40–50N latitude
library(terra)
library(sf)
library(rnaturalearth)
library(ggplot2)
library(dplyr)

# MI raster
r_mi <- rast('/Users/vv934811@reading.ac.uk/Documents/JIaze_species_Data/Bioclimatic_variables/lnMI_01res_28.tif')
mi_land <- exp(r_mi)

# China boundary
china_sf <- ne_countries(country = "China", scale = "large", returnclass = "sf")
china_sf <- st_transform(china_sf, crs(mi_land))
china_vect <- vect(china_sf)

# same China crop/mask as old figure
mi_china <- mi_land |> crop(china_vect) |> mask(china_vect)

# same NE window as old figure
bb_ne <- ext(111, 140, 40, 50)
mi_ne <- crop(mi_china, bb_ne)

# convert to data frame
mi_df <- as.data.frame(mi_ne, xy = TRUE, na.rm = FALSE)
colnames(mi_df) <- c("x", "y", "MI")

# MI classes
mi_df$MI_class <- cut(
  mi_df$MI,
  breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1, 2, Inf),
  labels = c("0-0.2", "0.4", "0.6", "0.8", "1", "2", ">2"),
  right = TRUE,
  include.lowest = TRUE
)

legend_levels_mi <- c("0-0.2", "0.4", "0.6", "0.8", "1", "2", ">2")
mi_df$MI_class <- factor(mi_df$MI_class, levels = legend_levels_mi)

my_colors_mi <- colorRampPalette(c("yellow", "orchid", "blue"))(7)

# site data with MAP classes
sta <- NECT_sample_climate %>%
  dplyr::select(Longitude, Latitude, MAP..mm.) %>%
  distinct(Longitude, Latitude, .keep_all = TRUE)

sta_33 <- sta[1:33, ]

sta_33$MAP_class <- cut(
  sta_33$MAP..mm.,
  breaks = c(0, 200, 300, 400, 500, 600, 700, Inf),
  labels = c("0-200", "200-300", "300-400", "400-500", "500-600", "600-700", ">700"),
  right = TRUE,
  include.lowest = TRUE
)

legend_levels_map <- c("0-200", "200-300", "300-400", "400-500", "500-600", "600-700", ">700")
sta_33$MAP_class <- factor(sta_33$MAP_class, levels = legend_levels_map)

map_colors <- c(
  "0-200"   = "#440154",
  "200-300" = "#414487",
  "300-400" = "#2A788E",
  "400-500" = "#22A884",
  "500-600" = "#7AD151",
  "600-700" = "#B8DE29",
  ">700"    = "#FDE725"
)

ggplot(mi_df) +
  geom_raster(aes(x = x, y = y, fill = MI_class)) +
  geom_point(
    data = sta_33,
    aes(x = Longitude, y = Latitude, color = MAP_class),
    size = 4
  ) +
  scale_fill_manual(
    values = my_colors_mi,
    name = "MI (unitless)",
    limits = legend_levels_mi,
    breaks = legend_levels_mi,
    na.translate = FALSE
  ) +
  scale_color_manual(
    values = map_colors,
    name = "MAP (mm)",
    limits = legend_levels_map,
    breaks = legend_levels_map,
    drop = FALSE
  ) +
  scale_x_continuous(
    breaks = c(115, 120, 125, 130, 135),
    labels = c("115°E", "120°E", "125°E", "130°E", "135°E")
  ) +
  scale_y_continuous(
    breaks = c(42, 44, 46, 48, 50),
    labels = c("42°N", "44°N", "46°N", "48°N", "50°N")
  ) +
  xlab("Longitude") +
  ylab("Latitude") +
  coord_fixed(
    xlim = c(111, 135),
    ylim = c(40, 50),
    expand = FALSE
  ) +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme_bw() +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    legend.background = element_rect(fill = "white", colour = NA),
    legend.key = element_rect(fill = "white", colour = NA),
    legend.position = "right",
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text = element_text(size = rel(1.5)),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_line(color = "white")
  )
#################################################################Monte Carlo
library(MASS)
library(dplyr)
library(ggplot2)
RSimprove2 <- read.csv("Climate_R_S_NECT_predictRS_biomass.csv")

# woody / non-woody 子集
woodyRS   <- subset(RSimprove2, species == "woody")
nonwoodRS <- subset(RSimprove2, species == "non-woody")

# woody model
lmwoody <- lm(
  log(R_S) ~ tg_1987_2006 + log(AI) + log(RZwc) + log(GPP) + Sand +
    tg_1987_2006 * log(AI) + lnheightT + lnRRD + lnSLA +
    lnthickness + lnLDMC,
  data = woodyRS
)

# non-woody model
lmnonwoody <- lm(
  log(R_S) ~ tg_1987_2006 + log(AI) + log(GPP) + pH + Sand +
    lnheightT + lnRRD + lnSRL + lnthickness + lnLDMC,
  data = nonwoodRS
)

woodyRS$pred_lnRSnew   <- predict(lmwoody, newdata = woodyRS)
nonwoodRS$pred_lnRSnew <- predict(lmnonwoody, newdata = nonwoodRS)

RSimprove2$pred_lnRSnew <- NA
RSimprove2$pred_lnRSnew[RSimprove2$species == "woody"]     <- woodyRS$pred_lnRSnew
RSimprove2$pred_lnRSnew[RSimprove2$species == "non-woody"] <- nonwoodRS$pred_lnRSnew

############################################################
# 3. Monte Carlo
############################################################
mc_predict_sitegroup <- function(data_sub,
                                 lmwoody,
                                 lmnonwoody,
                                 mi_var = "MI..unitless.",
                                 nsim = 1000,
                                 nbins = 8,
                                 add_residual = FALSE,
                                 seed = 123) {
  
  set.seed(seed)
  
  woody_dat   <- subset(data_sub, species == "woody")
  nonwood_dat <- subset(data_sub, species == "non-woody")
  
  pred_list <- list()
  
  # woody
  if (nrow(woody_dat) > 0) {
    Xw <- model.matrix(delete.response(terms(lmwoody)), data = woody_dat)
    bw <- coef(lmwoody)
    Vw <- vcov(lmwoody)
    bw_sim <- MASS::mvrnorm(n = nsim, mu = bw, Sigma = Vw)
    pred_w <- Xw %*% t(bw_sim)
    
    if (add_residual) {
      pred_w <- pred_w + matrix(
        rnorm(nrow(woody_dat) * nsim, mean = 0, sd = sigma(lmwoody)),
        nrow = nrow(woody_dat), ncol = nsim
      )
    }
    
    pred_list$woody <- list(data = woody_dat, pred = pred_w)
  }
  
  # non-woody
  if (nrow(nonwood_dat) > 0) {
    Xn <- model.matrix(delete.response(terms(lmnonwoody)), data = nonwood_dat)
    bn <- coef(lmnonwoody)
    Vn <- vcov(lmnonwoody)
    bn_sim <- MASS::mvrnorm(n = nsim, mu = bn, Sigma = Vn)
    pred_n <- Xn %*% t(bn_sim)
    
    if (add_residual) {
      pred_n <- pred_n + matrix(
        rnorm(nrow(nonwood_dat) * nsim, mean = 0, sd = sigma(lmnonwoody)),
        nrow = nrow(nonwood_dat), ncol = nsim
      )
    }
    
    pred_list$nonwoody <- list(data = nonwood_dat, pred = pred_n)
  }
  
  all_dat <- bind_rows(
    if (!is.null(pred_list$woody)) pred_list$woody$data,
    if (!is.null(pred_list$nonwoody)) pred_list$nonwoody$data
  )
  
  pred_mat_all <- rbind(
    if (!is.null(pred_list$woody)) pred_list$woody$pred,
    if (!is.null(pred_list$nonwoody)) pred_list$nonwoody$pred
  )
  
  mi <- all_dat[[mi_var]]
  
  breaks <- seq(min(mi, na.rm = TRUE), max(mi, na.rm = TRUE), length.out = nbins + 1)
  bin_id <- cut(mi, breaks = breaks, include.lowest = TRUE, labels = FALSE)
  
  out_list <- vector("list", nbins)
  
  for (b in 1:nbins) {
    idx <- which(bin_id == b)
    if (length(idx) == 0) next
    
    vals <- as.vector(pred_mat_all[idx, , drop = FALSE])
    
    out_list[[b]] <- data.frame(
      MI_mid = mean(mi[idx], na.rm = TRUE),
      MI_min = min(mi[idx], na.rm = TRUE),
      MI_max = max(mi[idx], na.rm = TRUE),
      pred_median = unname(quantile(vals, 0.50, na.rm = TRUE)),
      pred_low    = unname(quantile(vals, 0.025, na.rm = TRUE)),
      pred_high   = unname(quantile(vals, 0.975, na.rm = TRUE)),
      n = length(idx)
    )
  }
  
  bind_rows(out_list)
}

############################################################
# 4.uncertainty analysis
############################################################
grass_dat <- subset(RSimprove2, sitegroup == "grasssite")
tree_dat  <- subset(RSimprove2, sitegroup == "treesite")

# grasssite
grass_band <- mc_predict_sitegroup(
  data_sub = grass_dat,
  lmwoody = lmwoody,
  lmnonwoody = lmnonwoody,
  mi_var = "MI..unitless.",
  nsim = 1000,
  nbins = 6,         
  add_residual = FALSE,
  seed = 123
)
grass_band$sitegroup <- "grasssite"

# treesite
tree_band <- mc_predict_sitegroup(
  data_sub = tree_dat,
  lmwoody = lmwoody,
  lmnonwoody = lmnonwoody,
  mi_var = "MI..unitless.",
  nsim = 1000,
  nbins = 4,           # treesite 点更少，bin更少
  add_residual = FALSE,
  seed = 123
)
tree_band$sitegroup <- "treesite"

band_site <- bind_rows(grass_band, tree_band)

############################################################
# 5. plot
############################################################
plot_site <- RSimprove2 %>%
  transmute(
    sitegroup = sitegroup,
    MI = MI..unitless.,
    obs_lnRS = log(R_S),
    pred_lnRS = pred_lnRSnew
  )
p_site_uncertainty <- ggplot() +
  geom_point(
    data = plot_site,
    aes(x = MI, y = obs_lnRS),
    color = "blue", size = 2.8, alpha = 0.75
  ) +
  geom_point(
    data = plot_site,
    aes(x = MI, y = pred_lnRS),
    color = "red", size = 2.8, alpha = 0.75
  ) +
  geom_ribbon(
    data = band_site,
    aes(x = MI_mid, ymin = pred_low, ymax = pred_high),
    fill = "red", alpha = 0.18
  ) +
  geom_line(
    data = band_site,
    aes(x = MI_mid, y = pred_median),
    color = "red", linewidth = 1.2
  ) +
  facet_wrap(~sitegroup, scales = "free_x") +
  xlab("MI (unitless)") +
  ylab("ln(R:S)") +
  theme_bw() +
  theme(
    strip.text = element_text(size = 22),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 17),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

print(p_site_uncertainty)

plot_site$sitegroup2 <- dplyr::recode(
  plot_site$sitegroup,
  "grasssite" = "Grass sites",
  "treesite"  = "Tree sites"
)

band_site$sitegroup2 <- dplyr::recode(
  band_site$sitegroup,
  "grasssite" = "Grass sites",
  "treesite"  = "Tree sites"
)


ggplot() +
  geom_point(
    data = plot_site,
    aes(x = MI, y = obs_lnRS),
    color = "blue", size = 2.6, alpha = 0.65
  ) +
  geom_point(
    data = plot_site,
    aes(x = MI, y = pred_lnRS),
    color = "red", size = 2.6, alpha = 0.55
  ) +
  geom_ribbon(
    data = band_site,
    aes(x = MI_mid, ymin = pred_low, ymax = pred_high),
    fill = "red", alpha = 0.16
  ) +
  geom_smooth(
    data = band_site,
    aes(x = MI_mid, y = pred_median),
    method = "loess", se = FALSE,
    span = 0.9, color = "red", linewidth = 1.2
  ) +
  facet_wrap(~sitegroup2, scales = "free_x") +
  xlab("MI (unitless)") +
  ylab("ln(R:S)") +
  theme_bw() +
  theme(
    strip.text = element_text(size = 22),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 17),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
###############################################################
load("3 clim.trait.RData")
climate.trait$fD <- ifelse(climate.trait$DecEv == "D", climate.trait$f, NA)
climate.trait$fE <- ifelse(climate.trait$DecEv == "E", climate.trait$f, NA)
wangCPT<-climate.trait %>% filter(grepl("CPT", source))

# Filter rows where the 'Site' column contains values from 1 to 33
wangCPT_filtered <- wangCPT %>%
  filter(site >= 1 & site <= 33)

wangCPT_filtered$fD <- ifelse(wangCPT_filtered$DecEv == "D", wangCPT_filtered$f, NA)
wangCPT_filtered$fE <- ifelse(wangCPT_filtered$DecEv == "E", wangCPT_filtered$f, NA)

library(MASS)
library(dplyr)
library(ggplot2)

############################################################
# 1. Fit deciduous LMA model
############################################################
modeldeciduous <- lm(
  log(LMA) ~ log(IPAR) + log(fD) + T_g_m + log(alpha_m),
  data = wangCPT_filtered
)

summary(modeldeciduous)

############################################################
# 2. Original fitted values
############################################################
wangCPT_filtered$pred_lnLMA_decid <- predict(modeldeciduous, newdata = wangCPT_filtered)
wangCPT_filtered$pred_LMA_decid   <- exp(wangCPT_filtered$pred_lnLMA_decid)
wangCPT_filtered <- wangCPT_filtered %>%
  dplyr::left_join(NECT_mi, by = c("Rain" = "MAP..mm."))

############################################################
# 3. Monte Carlo function for coefficient uncertainty
############################################################
mc_predict_binned_lma <- function(model,
                                  data,
                                  mi_var = "MI..unitless.",
                                  nsim = 1000,
                                  nbins = 8,
                                  add_residual = FALSE,
                                  seed = 123) {
  
  set.seed(seed)
  
  # Design matrix from original observations
  X <- model.matrix(delete.response(terms(model)), data = data)
  
  # Coefficient estimates and covariance matrix
  beta_hat <- coef(model)
  V_beta   <- vcov(model)
  
  # Draw coefficients from multivariate normal
  beta_sim <- MASS::mvrnorm(n = nsim, mu = beta_hat, Sigma = V_beta)
  # beta_sim: nsim x p
  
  # Simulated predictions
  pred_mat <- X %*% t(beta_sim)   # n_obs x nsim
  
  # Optional residual uncertainty
  if (add_residual) {
    pred_mat <- pred_mat + matrix(
      rnorm(nrow(data) * nsim, mean = 0, sd = sigma(model)),
      nrow = nrow(data), ncol = nsim
    )
  }
  
  # MI values
  mi <- data[[mi_var]]
  
  # Bin MI
  breaks <- seq(min(mi, na.rm = TRUE),
                max(mi, na.rm = TRUE),
                length.out = nbins + 1)
  
  bin_id <- cut(mi, breaks = breaks, include.lowest = TRUE, labels = FALSE)
  
  # Summarize within bins
  out_list <- vector("list", nbins)
  
  for (b in 1:nbins) {
    idx <- which(bin_id == b)
    if (length(idx) == 0) next
    
    vals <- as.vector(pred_mat[idx, , drop = FALSE])
    
    out_list[[b]] <- data.frame(
      MI_mid = mean(mi[idx], na.rm = TRUE),
      MI_min = min(mi[idx], na.rm = TRUE),
      MI_max = max(mi[idx], na.rm = TRUE),
      pred_median = unname(quantile(vals, 0.50, na.rm = TRUE)),
      pred_low    = unname(quantile(vals, 0.025, na.rm = TRUE)),
      pred_high   = unname(quantile(vals, 0.975, na.rm = TRUE)),
      n = length(idx)
    )
  }
  
  band_df <- bind_rows(out_list)
  
  return(list(
    band = band_df,
    pred_mat = pred_mat,
    beta_sim = beta_sim
  ))
}

############################################################
# 4. Keep only rows with complete data for the deciduous model
############################################################
decid_dat <- wangCPT_filtered %>%
  dplyr::select(LMA, IPAR, fD, T_g_m, alpha_m, MI..unitless.) %>%
  na.omit()

# Refit on the exact analysis dataset if you want full consistency
modeldeciduous2 <- lm(
  log(LMA) ~ log(IPAR) + log(fD) + T_g_m + log(alpha_m),
  data = decid_dat
)

decid_dat$pred_lnLMA <- predict(modeldeciduous2, newdata = decid_dat)

############################################################
# 5. Run Monte Carlo uncertainty analysis
############################################################
res_decid <- mc_predict_binned_lma(
  model = modeldeciduous2,
  data = decid_dat,
  mi_var = "MI..unitless.",
  nsim = 1000,
  nbins = 8,
  add_residual = FALSE,   # coefficient uncertainty only
  seed = 123
)

band_decid <- res_decid$band

############################################################
# 6. Plot observed and predicted lnLMA with uncertainty ribbon
############################################################
p_decid_unc <- ggplot() +
  geom_point(
    data = decid_dat,
    aes(x = MI..unitless., y = log(LMA)),
    color = "blue", alpha = 0.65, size = 2.5
  ) +
  geom_point(
    data = decid_dat,
    aes(x = MI..unitless., y = pred_lnLMA),
    color = "red", alpha = 0.55, size = 2.3
  ) +
  geom_ribbon(
    data = band_decid,
    aes(x = MI_mid, ymin = pred_low, ymax = pred_high),
    fill = "red", alpha = 0.18
  ) +
  geom_smooth(
    data = band_decid,
    aes(x = MI_mid, y = pred_median),
    method = "loess", se = FALSE,
    span = 0.9, color = "red", linewidth = 1.2
  ) +
  xlab("MI (unitless)") +
  ylab(expression('lnLMA'*~'('*g*'/'*m^2*')')) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

p_decid_unc

