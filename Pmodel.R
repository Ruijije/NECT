##################################################################################################### rerun the rpmodel
###############vpd, fapar, when predidct Vcmax,X, fapar and LAI we only use the data on 2006
###########################################################################
##data from chelsa 2.1 by first multiplying the raster values with the ‘scale’ value and then adding the ‘offset’ value.
###########################################################################
NECT_sample_climate<-read.csv("/Users/echo/E/eechange project/analysis/work2024/ReworkALL/NECT_Climate_trait.csv")
setwd("/Users/echo/E/eechange project/analysis/work2024/ReworkALL/chelsa_2006/tasmax")
tmaxfile <- list.files(pattern = "*.tif") 
tmaxfile_stack <- stack(tmaxfile)           # Load TIFF files into a RasterStack
site<-cbind(NECT_sample_climate$Longitude,NECT_sample_climate$Latitude)
tmax_NECT<-raster::extract(tmaxfile_stack,site,method = 'bilinear')
tmax_NECT<-data.frame(tmax_NECT)
tmax_C_NECT<- (tmax_NECT/10)-273.15##K/10 ,change unit from K to C
write.csv(tmax_C_NECT,"/Users/echo/E/eechange project/analysis/work2024/ReworkALL/chelsa_2006/tasmax/tmax_C_NECT.csv")

setwd("/Users/echo/E/eechange project/analysis/work2024/ReworkALL/chelsa_2006/tasmin")
tminfile <- list.files(pattern = "*.tif") 
tminfile_stack <- stack(tminfile)           # Load TIFF files into a RasterStack
site<-cbind(NECT_sample_climate$Longitude,NECT_sample_climate$Latitude)
tmin_NECT<-raster::extract(tminfile_stack,site,method = 'bilinear')
tmin_NECT<-data.frame(tmin_NECT)
tmin_C_NECT<- (tmin_NECT/10)-273.15##K/10 ,change unit from K to C
write.csv(tmin_C_NECT,"tmin_C_NECT.csv")

lat<-NECT_sample_climate$Latitude
tmx<-tmax_C_NECT
tmn<-tmin_C_NECT
s1 <- -23.1
s2 <- -17.3
s3 <- -8.0
s4 <- 4.1
s5 <- 14.8
s6 <- 21.9
s7 <- 23.2
s8 <- 18.3
s9 <- 8.6
s10 <- -2.8
s11 <- -14.1
s12 <- -21.6
s <- c(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12)
x <-data.frame(matrix(nrow = nrow(tmn),ncol=ncol(tmn))) ##create x dataframe: nrow=No. of sites, ncol=timestep
Tg_NECT <-data.frame(matrix(nrow = nrow(tmn),ncol=ncol(tmn)))
for (i in 1:12){
  
  x[,i]<- -tan(pi*lat/180)*tan(s[i]*pi/180)
  Tg_NECT[,i]<-tmx[,i]*(0.5+(1-x[,i]^2)^(0.5)/(2*acos(x[,i])))+ tmn[,i]*(0.5-(1-x[,i]^2)^(0.5)/(2*acos(x[,i])))
  
}
write.csv(Tg_NECT,"/Users/echo/E/eechange project/analysis/work2024/ReworkALL/chelsa_2006/Tg_chelsa_2006.csv")
#############################################################
setwd("/Users/echo/E/eechange project/analysis/work2024/ReworkALL/chelsa_2006/vpd")
vpdfile <- list.files(pattern = "*.tif") 
vpdfile_stack <- stack(vpdfile)           # Load TIFF files into a RasterStack
site<-cbind(NECT_sample_climate$Longitude,NECT_sample_climate$Latitude)
vpd_NECT<-raster::extract(vpdfile_stack,site,method = 'bilinear')
vpd_NECT<-data.frame(vpd_NECT)
vpd_Pa_NECT<-vpd_NECT*0.1 ##scale:0.1
write.csv(vpd_Pa_NECT,"vpd_Pa_NECT.csv")

#rsds Surface downwelling shortwave flux in air. MJ m-2 d-1  scale 0.001
setwd("/Users/echo/E/eechange project/analysis/work2024/ReworkALL/chelsa_2006/rsds")
rsdsfile <- list.files(pattern = "*.tif") 
rsdsfile_stack <- stack(rsdsfile)           # Load TIFF files into a RasterStack
plot(rsdsfile_stack)
site<-cbind(NECT_sample_climate$Longitude,NECT_sample_climate$Latitude)
rsds_NECT<-raster::extract(rsdsfile_stack,site,method = 'bilinear')
rsds_NECT<-data.frame(rsds_NECT)
sw<- rsds_NECT*11.574*0.001####MJ m-2 d-1 to w/m2; W m−2= MJ m−2d −1×1,000,000/86,400; scale for the data 0.001
write.csv(sw,"sw_w_m2_NECT.csv")

ppfd_NECT<-60*60*24*10^-6*sw[, 1:12]*2.04##molm-2d-1
ppfd_NECT<-data.frame(ppfd_NECT)
write.csv(ppfd_NECT,"/Users/echo/E/eechange project/analysis/work2024/ReworkALL/chelsa_2006/rsds/ppfd_molm-2d-1_NECT.csv")
#####################################fapar
fapar2006<-stack('/Users/echo/E/eechange project/analysis/work2024/ReworkALL/chelsa_2006/fapar_0.05/2006.nc')
plot(fapar2006)
fapar2006_NECT<-raster::extract(fapar2006,site,method = 'bilinear')
fapar2006_NECT<-data.frame(fapar2006_NECT)
write.csv(fapar2006_NECT,"/Users/echo/E/eechange project/analysis/work2024/ReworkALL/chelsa_2006/fapar_0.05/fapar2006_NECT.csv")
####################################################run pmodel when beta is constant
tc <- rowMeans(Tg_NECT[, 1:12], na.rm = TRUE) ##794
vpd<- rowMeans(vpd_Pa_NECT[, 1:12], na.rm = TRUE)
fapar<-rowMeans(fapar2006_NECT[, 1:12], na.rm = TRUE)
ppfd<-rowMeans(ppfd_NECT[, 1:12], na.rm = TRUE)
elv<-NECT_sample_climate$Elevation
beta<-146
pressure<-calc_patm(elv=NECT_sample_climate$Elevation, patm0 = 101325)
model_constant<-rpmodel(
  tc,
  vpd,
  co2=381.9,
  fapar,
  ppfd,
  patm=pressure,
  elv,
  kphio =   0.087182,
  beta ,
  soilm = 1,
  meanalpha = 1,
  apar_soilm = 0,
  bpar_soilm = 0.733,
  c4 = FALSE,
  method_jmaxlim = "wang17",
  do_ftemp_kphio = TRUE,
  do_soilmstress = FALSE,
  returnvar = NULL,
  verbose = FALSE
)

# Assuming model_constant is your list
# Convert the list to a dataframe
model_constant_df <- as.data.frame(model_constant)
# Inspect the dataframe
head(model_constant_df)
# Save the dataframe to a CSV file
write.csv(model_constant_df, "/Users/echo/E/eechange project/analysis/work2024/ReworkALL/chelsa_2006/Pmodel_constant_output.csv", row.names = FALSE)
#################################################################
##change beta, not use the constant 146. 
##Alienor paper lnβ <-1.73θ+4.55.  θ soil water content
##Monthly average soil water content (mm) from SPLASH 2.0 were converted into θ (m3/m3) over a 1 m soil depth.
##https://onlinelibrary.wiley.com/doi/full/10.1111/gcb.15364
################################################################# run rsplash 2.0

library(devtools)
##error install rsplash: topmodel’ is not available for package ‘rsplash’
#topmodle is moved from CRAN, install topmodel by hand 
install.packages("/Users/echo/E/eechange project/analysis/work2024/topmodel_0.7.5.tar.gz", repos = NULL, type = "source")
library(topmodel)
install_github("dsval/rsplash")
library(rsplash)
library(xts)

setwd("/Users/echo/E/eechange project/analysis/work2024/ReworkALL/chelsa_2006/pr")
pnfile <- list.files(pattern = "*.tif") 
pnfile_stack <- stack(pnfile)           # Load TIFF files into a RasterStack
site<-cbind(NECT_sample_climate$Longitude,NECT_sample_climate$Latitude)
pn_NECT<-raster::extract(pnfile_stack,site,method = 'bilinear')
pn_NECT<-data.frame(pn_NECT)
pn_mm_NECT<-pn_NECT/100 ##kg m-2 month1/100 kg/m²/month is numerically equivalent to mm/month
write.csv(pn_mm_NECT,"pn_mm_NECT.csv")

soildata<-read.csv("~/E/eechange project/data/splash_P_hydro/dailydata/soil_texture_from_splashtools.csv")
# Assuming `NECT_sample_climate` and `soildata` are your datasets

# Match the Site.ID in NECT_sample_climate with X in soildata
soildata_expanded <- soildata[match(NECT_sample_climate$Site.ID, soildata$X), ]
#soildata_expanded$Site.ID <- NECT_sample_climate$Site.ID #just check
Result_sw_in<-t(sw)
date <- seq(from = as.Date("2006-01-01"), 
            to = as.Date("2006-12-01"), 
            by = "month")
# Format the dates to "YYYY-MM-DD"
formatted_dates <- format(date, "%Y-%m-%d")
# Add the formatted dates as the first column in Result_sw_in
Result_sw_in <- data.frame(V1 = formatted_dates, Result_sw_in)
rownames(Result_sw_in)<-Result_sw_in[,1]
Result_sw_inasdata<-as.data.frame(Result_sw_in)
Result_sw_inlist<- as.list(Result_sw_inasdata)
Result_sw_inlist.ts <- xts(Result_sw_inlist$X1, order.by=as.POSIXct(Result_sw_inlist$V1))
Result_sw_inlist.ts

Result_tc<-t(Tg_NECT)
Result_tc <- data.frame(V1 = formatted_dates, Result_tc)
rownames(Result_tc)<-Result_tc[,1]
Result_tcasdata<-as.data.frame(Result_tc)
Result_tclist<- as.list(Result_tcasdata)
Result_tclist.ts <- xts(Result_tclist$X1, order.by=as.POSIXct(Result_tclist$V1))
Result_tclist.ts 

Result_pn<-t(pn_NECT)
days_in_month <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
# Expand monthly data to daily data
Result_pn <- sweep(Result_pn, 1, days_in_month, "/")
Result_pn <- data.frame(V1 = formatted_dates, Result_pn)
rownames(Result_pn)<-Result_pn[,1]
Result_pnasdata<-as.data.frame(Result_pn)
Result_pnlist<- as.list(Result_pnasdata)
Result_pnlist.ts <- xts(Result_pnlist$X1, order.by=as.POSIXct(Result_pnlist$V1))
Result_pnlist.ts 

soil<-soildata_expanded[,-1]
soil$depth <- as.numeric(soil$depth)
soil[] <- lapply(soil, as.numeric)
soil[1,]
##########################
soil_first_row <- as.numeric(unlist(soil[1, ]))
runNECT<-rsplash::splash.point(sw_in=Result_sw_inlist.ts,		# shortwave radiation W/m2
                               tc=Result_tclist.ts,		# air temperature C
                               pn=Result_pnlist.ts,		# precipitation mm
                               lat=NECT_sample_climate$Latitude[1],	# latitude deg
                               elev=NECT_sample_climate$Elevation[1],	# elevation masl
                               slop=0,	# slope deg
                               asp=0,	# aspect deg
                               soil_data=soil_first_row, 			# soil data: sand,clay,om,grvel in %, bulkdens g/cm3
                               Au=0,	# upslope area m2
                               resolution=250  			# resolution pixel dem used to get Au
)


# Initialize runNECT as an empty list to store results
runNECT1 <- list()

# Loop through each index i from 1 to 794
for (i in 2:795) {
  
  # Convert Result_sw_inlist, Result_pnlist, and Result_tclist to xts objects for the current iteration
  Result_sw_inlist[[i]] <- xts(Result_sw_inlist[[i]], order.by = as.POSIXct(Result_sw_inlist[[1]]))
  Result_pnlist[[i]] <- xts(Result_pnlist[[i]], order.by = as.POSIXct(Result_pnlist[[1]]))
  Result_tclist[[i]] <- xts(Result_tclist[[i]], order.by = as.POSIXct(Result_tclist[[1]]))
  
  # Extract the first row for the current soil (soil data for the i-th site)
  soil_first_row <- as.numeric(unlist(soil[i-1, ]))
  
  # Run the rsplash model for each iteration
  runNECT1[[i]] <- rsplash::splash.point(
    sw_in = Result_sw_inlist[[i]],         # Shortwave radiation W/m2
    tc = Result_tclist[[i]],               # Air temperature C
    pn = Result_pnlist[[i]],               # Precipitation mm
    lat = NECT_sample_climate$Latitude[i-1],   # Latitude in degrees
    elev = NECT_sample_climate$Elevation[i-1], # Elevation in meters
    slop = 0,                               # Slope in degrees (default 0)
    asp = 0,                                # Aspect in degrees (default 0)
    soil_data = soil_first_row,              # Soil data: sand, clay, OM, gravel, bulk density
    Au = 0,                                  # Upslope area in m2 (default 0)
    resolution = 250                         # Resolution of the DEM in pixels (default 250)
  )
  
  # Optionally, you can print progress to check the loop
  print(paste("Processed site", i))
}

daliy<- seq(from = as.Date("2006-01-01"), 
                   to = as.Date("2006-12-31"), 
                   by = "day")
all_wn_results <- do.call(cbind, lapply(runNECT1, function(x) x$wn))
# Convert the results to a data frame and include the date index if necessary
all_wn_results_df <- data.frame(Date = daliy, wn = all_wn_results)
# Save the results to a CSV file
write.csv(all_wn_results_df, "/Users/echo/E/eechange project/analysis/work2024/ReworkALL/chelsa_2006/splash_output/wn_mm_2006.csv", row.names = FALSE)
all_wn_results_m3[]<-((all_wn_results_df[,2:795]/1000)/2)##depth in soil data is 2, This results in a soil moisture content in volume percent
write.csv(all_wn_results_m3, "/Users/echo/E/eechange project/analysis/work2024/ReworkALL/chelsa_2006/splash_output/wn_m3m3_2006.csv", row.names = FALSE)

# Extract each variable from the runNECT1 results list
all_ro <- do.call(cbind, lapply(runNECT1, function(x) x$ro))
all_pet <- do.call(cbind, lapply(runNECT1, function(x) x$pet))
all_aet <- do.call(cbind, lapply(runNECT1, function(x) x$aet))
all_snow <- do.call(cbind, lapply(runNECT1, function(x) x$snow))
all_cond <- do.call(cbind, lapply(runNECT1, function(x) x$cond))
all_bflow <- do.call(cbind, lapply(runNECT1, function(x) x$bflow))
all_netr <- do.call(cbind, lapply(runNECT1, function(x) x$netr))
all_sm_lim <- do.call(cbind, lapply(runNECT1, function(x) x$sm_lim))

# Define a function to save each variable as a CSV file
save_variable_csv <- function(variable_data, variable_name) {
  # Create a data frame with Date index and the variable
  variable_df <- data.frame(Date = daliy, value = variable_data)
  
  # Define the file path
  file_path <- paste0("/Users/echo/E/eechange project/analysis/work2024/ReworkALL/chelsa_2006/splash_output/", variable_name, "_2006.csv")
  
  # Save the variable data to a CSV file
  write.csv(variable_df, file_path, row.names = FALSE)
  print(paste(variable_name, "saved to", file_path))
}

# Save each variable to a separate CSV
save_variable_csv(all_ro, "ro")
save_variable_csv(all_pet, "pet")
save_variable_csv(all_aet, "aet")
save_variable_csv(all_snow, "snow")
save_variable_csv(all_cond, "cond")
save_variable_csv(all_bflow, "bflow")
save_variable_csv(all_netr, "netr")
save_variable_csv(all_sm_lim, "sm_lim")
#####################################################################get the beta of Alinor paper
NECT_sample_climate<-read.csv("/Users/echo/E/eechange project/analysis/work2024/ReworkALL/NECT_Climate_trait.csv")
watercontent<-read.csv("/Users/echo/E/eechange project/analysis/work2024/ReworkALL/chelsa_2006/splash_output/wn_m3m3_2006.csv")
watercontent<-t(watercontent)
theta <- rowMeans(watercontent[, 1:365], na.rm = TRUE) ##794
NECT_sample_climate$theta_m3m3<-theta
lnbeta_Ali<-1.73*theta+4.55
NECT_sample_climate$beta_Ali<-exp(lnbeta_Ali)

NECT_sample_climate$tg_chelsa_2006 <- rowMeans(Tg_NECT[, 1:12], na.rm = TRUE) ##794
NECT_sample_climate$vpd_chelsa_2006<- rowMeans(vpd_Pa_NECT[, 1:12], na.rm = TRUE)
NECT_sample_climate$fapar_MODIS_2006<-rowMeans(fapar2006_NECT[, 1:12], na.rm = TRUE)
NECT_sample_climate$ppfd_chelsa_2006<-rowMeans(ppfd_NECT[, 1:12], na.rm = TRUE)
pressure<-calc_patm(elv=NECT_sample_climate$Elevation, patm0 = 101325)
model_Alienor<-rpmodel(
  tc=NECT_sample_climate$tg_chelsa_2006,
  vpd=NECT_sample_climate$vpd_chelsa_2006,
  co2=381.9,
  fapar=NECT_sample_climate$fapar_MODIS_2006,
  ppfd=NECT_sample_climate$ppfd_chelsa_2006,
  patm =pressure ,
  elv=NECT_sample_climate$Elevation,
  kphio =  0.087182,
  beta=NECT_sample_climate$beta_Ali ,
  soilm = 1,
  meanalpha = 1,
  apar_soilm = 0,
  bpar_soilm = 0.733,
  c4 = FALSE,
  method_jmaxlim = "wang17",
  do_ftemp_kphio = TRUE,
  do_soilmstress = FALSE,
  returnvar = NULL,
  verbose = FALSE
)


model_Alienor_df <- as.data.frame(model_Alienor)
head(model_Alienor_df)
# Save the dataframe to a CSV file
write.csv(model_Alienor_df, "/Users/echo/E/eechange project/analysis/work2024/ReworkALL/chelsa_2006/Pmodel_Alienor_output.csv", row.names = FALSE)

###################################################################calculate my beta
#################################first predict the observaition of chi based on carbon isotope
cisotope=NECT_sample_climate$d13C.12C
C13air <- data.frame(read.table("~/E/eechange project/data/co2c13_wlg_surface-flask_1_sil_event.txt", header=T))
colnames(C13air)[2] <- "year"
C13air2006 <- C13air[which(C13air$year%in%2006),]
C13air2006data <- C13air2006[,12]
cisotopeair<-mean(C13air2006data)
deta=((cisotopeair-cisotope)/(1+(cisotope/1000)))  # deta = (13Cair- 13Cleaf) / (1+(13Cleaf/1000))
NECT_sample_climate$deta_obs<-deta
####################################calculate gammastar
##equation from rpmodel calc_patm(elv, patm0 = 101325) 
pressure<-calc_patm(elv=NECT_sample_climate$Elevation, patm0 = 101325)
NECT_sample_climate$pressure<-pressure
##equation from rpmodel pdf.calc_gammastar(tc, patm), 
NECT_sample_climate$gammastar<-calc_gammastar(tc=NECT_sample_climate$tg_chelsa_2006, patm=NECT_sample_climate$pressure)
NECT_sample_climate$ca<-co2_to_ca(co2=381.9, patm=NECT_sample_climate$pressure)
NECT_sample_climate$kmm<-calc_kmm(tc=NECT_sample_climate$tg_chelsa_2006, patm=NECT_sample_climate$pressure)

NECT_sample_climate$chi_obs=((NECT_sample_climate$deta-4.4+16*NECT_sample_climate$gammastar/NECT_sample_climate$ca)/(30-4.4)) ##
NECT_sample_climate$chi_obs[
  is.na(NECT_sample_climate$deta) | 
    is.na(NECT_sample_climate$gammastar) | 
    is.na(NECT_sample_climate$ca)
] <- NA
NECT_sample_climate$chi_obs[NECT_sample_climate$chi_obs<0.3] <- NA
###equation from https://onlinelibrary.wiley.com/doi/full/10.1111/gcb.15364 equation (2)


##########################calclulation x by equation, just test
Tg_K<-NECT_sample_climate$tg_chelsa_2006+273.15     # in K
Gstar <- pressure*4.275*exp((37830/8.31447)*(1/298.15 - 1/(Tg_K)))/101325    #Gstar in Pa [Bernacchi et al. (2001)], KC, KO and ??* equal to 79430, 36380 and 37830 J/mol# standard values of KC, KO and ??* at 25?? equals to 39.37, 27480 and 4.22 Pa. O
##1.34-2.03
Kc <- 40.49*exp((79430/8.31447)*(1/298.15 - 1/(Tg_K)))                      #Kc in Pa -> 1
Oc <-0.21*pressure #Oc in Pa -> month/daily
#Po = 21000*exp(-0.114*ele*0.001)
Ko <-27840*exp((36380/8.31447)*(1/298.15 - 1/(Tg_K)))                       #Ko in Pa -> month/daily
K  <- Kc*(1+(Oc/Ko)) #O in Pa -> month/daily
##13.36-22.21



NECT_sample_climate$ns_star <- exp(580*(1/(Tg_K-138)-1/160))
##equation(3)
NECT_sample_climate$beta_own<-1.6*NECT_sample_climate$ns_star*(NECT_sample_climate$vpd_chelsa_2006)*((NECT_sample_climate$chi_obs-NECT_sample_climate$gammastar/NECT_sample_climate$ca)^2/((1-NECT_sample_climate$chi_obs)^2*(NECT_sample_climate$kmm+NECT_sample_climate$gammastar)))
NECT_sample_climate$beta_own[ is.na(NECT_sample_climate$chi_obs)] <- NA
#vpd=Pa here
##in Alinor paper lnbepa is around 3-7 so I filter my data from e3-e7, because some are unreasonable
lower_bound <- exp(3)  # e^3
upper_bound <- exp(7)  # e^7

# Replace values outside the range with NA
NECT_sample_climate$beta_own <- ifelse(
  NECT_sample_climate$beta_own >= lower_bound & NECT_sample_climate$beta_own <= upper_bound,
  NECT_sample_climate$beta_own,
  NA
)
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
model_own<-rpmodel(
  tc=NECT_sample_climate$tg_chelsa_2006,
  vpd=NECT_sample_climate$vpd_chelsa_2006,
  co2=381.9,
  fapar=NECT_sample_climate$fapar_MODIS_2006,
  ppfd=NECT_sample_climate$ppfd_chelsa_2006,
  patm =pressure ,
  elv=NECT_sample_climate$Elevation,
  kphio =   0.087182,
  beta=NECT_sample_climate$beta_own ,
  soilm = 1,
  meanalpha = 1,
  apar_soilm = 0,
  bpar_soilm = 0.733,
  c4 = FALSE,
  method_jmaxlim = "wang17",
  do_ftemp_kphio = TRUE,
  do_soilmstress = FALSE,
  returnvar = NULL,
  verbose = FALSE
)

model_own_df <- as.data.frame(model_own)
head(model_own_df)
# Save the dataframe to a CSV file
write.csv(model_own_df, "/Users/echo/E/eechange project/analysis/work2024/ReworkALL/chelsa_2006/Pmodel_ownbeta_output.csv", row.names = FALSE)

write.csv(NECT_sample_climate,"/Users/echo/E/eechange project/analysis/work2024/ReworkALL/chelsa_2006/Climate_trait_beta.csv")
#############################################################################
##################plot figure
#############################################################################
library(cowplot)
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

plot_grid(pa1,pa2,pa3,label_x = 0.2,nrow = 2)
ggsave(filename = "Fig1_environmental.jpg", path = "/Users/echo/E/eechange project/analysis/work2024/ReworkALL/chelsa_2006/plot")
############################################three chi compare
chicombine<-data.frame(cbind(NECT_sample_climate$Site.ID.1,NECT_sample_climate$MAP..mm.,NECT_sample_climate$vpd_chelsa_2006,NECT_sample_climate$chi_obs,model_constant_df$chi,model_Alienor_df$chi,model_own_df$chi))

colnames(chicombine) <- c(
  "Site_ID", 
  "MAP_mm", 
  "VPD_CHELSA_2006", 
  "Chi_Obs", 
  "Chi_Constant", 
  "Chi_Alienor", 
  "Chi_Own"
)
write.csv(chicombine,"/Users/echo/E/eechange project/analysis/work2024/ReworkALL/chelsa_2006/Chi_threemodel.csv")
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
  ) +
  annotate("text", x=693, y=0.9, label="(a)",size=6.5)
  # Annotate the plot

chi2<-ggplot(site_summary, aes(x = vpd_pa)) +
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
  xlab(expression('VPD'*~'('*'Pa'*')')) +
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
  ) +
  annotate("text", x=740, y=0.9, label="(b)",size=6.5)
# Annotate the plot
chi2

plot_grid(chi1,chi2,label_x = 0.2,nrow = 1)
ggsave(filename = "FigS2_X_map_vpd.jpg", path = "/Users/echo/E/eechange project/analysis/work2024/ReworkALL/chelsa_2006/plot")

############################################################ three vcmax compare
##change PPFD data unit from mol/m2/d to umol/m2/s to get vcmax and vcmax 25
#PPFD (mol m−2 d−1) is derived from shortwave downwelling radiation as , where kEC=2.04 µmol J−1
#umol/m2/s must
##these are not right, I will calculate by equation.

Tg_K<-NECT_sample_climate$tg_chelsa_2006+273.15

Gstar <- pressure*4.275*exp((37830/8.31447)*(1/298.15 - 1/(Tg_K)))/101325    #Gstar in Pa [Bernacchi et al. (2001)], KC, KO and ??* equal to 79430, 36380 and 37830 J/mol# standard values of KC, KO and ??* at 25?? equals to 39.37, 27480 and 4.22 Pa. O
Ca_constant<-model_constant_df$ca
Ci_constant<-model_constant_df$ci
m_constant <- (Ci_constant-Gstar)/(Ci_constant+2*Gstar)

Kc <- 40.49*exp((79430/8.31447)*(1/298.15 - 1/(Tg_K)))                      #Kc in Pa -> 1
Oc <-0.21*pressure #Oc in Pa -> month/daily
Ko <-27840*exp((36380/8.31447)*(1/298.15 - 1/(Tg_K)))                       #Ko in Pa -> month/daily
K  <- Kc*(1+(Oc/Ko)) 

PPFD<-rowMeans(sw[, 1:12], na.rm = TRUE)*0.45/0.22

Vcmax_constant <- ((0.352+0.021*(Tg_K-273.15)-3.4*0.0001*(Tg_K-273.15)^2)/8)*PPFD*((Ci_constant+K)/(Ci_constant+2*Gstar))*sqrt(1-(0.41/m_constant)^(2/3))
Vcmax25_constant <- Vcmax*exp((65330/8.31447)*((1/Tg_K)-(1/298.15)))
####
Gstar <- pressure*4.275*exp((37830/8.31447)*(1/298.15 - 1/(Tg_K)))/101325    #Gstar in Pa [Bernacchi et al. (2001)], KC, KO and ??* equal to 79430, 36380 and 37830 J/mol# standard values of KC, KO and ??* at 25?? equals to 39.37, 27480 and 4.22 Pa. O
Ca_Alienor<-model_Alienor_df$ca
Ci_Alienor<-model_Alienor_df$ci
m_Alienor <- (Ci_Alienor-Gstar)/(Ci_Alienor+2*Gstar)

Kc <- 40.49*exp((79430/8.31447)*(1/298.15 - 1/(Tg_K)))                      #Kc in Pa -> 1
Oc <-0.21*pressure #Oc in Pa -> month/daily
Ko <-27840*exp((36380/8.31447)*(1/298.15 - 1/(Tg_K)))                       #Ko in Pa -> month/daily
K  <- Kc*(1+(Oc/Ko)) 

PPFD<-rowMeans(sw[, 1:12], na.rm = TRUE)*0.45/0.22

Vcmax_Alienor <- ((0.352+0.021*(Tg_K-273.15)-3.4*0.0001*(Tg_K-273.15)^2)/8)*PPFD*((Ci_Alienor+K)/(Ci_Alienor+2*Gstar))*sqrt(1-(0.41/m_Alienor)^(2/3))
Vcmax25_Alienor <- Vcmax*exp((65330/8.31447)*((1/Tg_K)-(1/298.15)))
#####
Gstar <- pressure*4.275*exp((37830/8.31447)*(1/298.15 - 1/(Tg_K)))/101325    #Gstar in Pa [Bernacchi et al. (2001)], KC, KO and ??* equal to 79430, 36380 and 37830 J/mol# standard values of KC, KO and ??* at 25?? equals to 39.37, 27480 and 4.22 Pa. O
Ca_own<-model_own_df$ca
Ci_own<-model_own_df$ci
m_own <- (Ci_own-Gstar)/(Ci_own+2*Gstar)

Kc <- 40.49*exp((79430/8.31447)*(1/298.15 - 1/(Tg_K)))                      #Kc in Pa -> 1
Oc <-0.21*pressure #Oc in Pa -> month/daily
Ko <-27840*exp((36380/8.31447)*(1/298.15 - 1/(Tg_K)))                       #Ko in Pa -> month/daily
K  <- Kc*(1+(Oc/Ko)) 

PPFD<-rowMeans(sw[, 1:12], na.rm = TRUE)*0.45/0.22

Vcmax_own <- ((0.352+0.021*(Tg_K-273.15)-3.4*0.0001*(Tg_K-273.15)^2)/8)*PPFD*((Ci_own+K)/(Ci_own+2*Gstar))*sqrt(1-(0.41/m_own)^(2/3))
Vcmax25_own <- Vcmax*exp((65330/8.31447)*((1/Tg_K)-(1/298.15)))

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
write.csv(vcmaxcombine,"/Users/echo/E/eechange project/analysis/work2024/ReworkALL/chelsa_2006/vcmax_threemodel.csv")
site_summary <- vcmaxcombine %>%
  group_by(Site_ID) %>%
  summarize(
    MAP_mm = mean(MAP_mm, na.rm = TRUE),  # Mean MAP for each site
    vcmax_Constant_Mean = mean(vcmax_Constant, na.rm = TRUE),
    vcmax_Alienor_Mean = mean(vcmax_Alienor, na.rm = TRUE),
    vcmax_Own_Mean = mean(vcmax_Own, na.rm = TRUE),
    vcmax25_Constant_Mean = mean(vcmax25_Constant, na.rm = TRUE),
    vcmax25_Alienor_Mean = mean(vcmax25_Alienor, na.rm = TRUE),
    vcmax25_Own_Mean = mean(vcmax25_Own, na.rm = TRUE)
  )

vcmax1<-ggplot(site_summary, aes(x = MAP_mm)) +
  # Add points and error bars for vcmax_Obs
  geom_point(aes(y = vcmax_Constant_Mean, color = "vcmax_Constant")) +
  geom_point(aes(y = vcmax_Alienor_Mean, color = "vcmax_Alienor")) +
  geom_point(aes(y = vcmax_Own_Mean, color = "vcmax_Own")) +
  # Add smooth curves for other vcmax types
  geom_smooth(aes(y = vcmax_Constant_Mean, color = "vcmax_Constant"), method = "lm", se = TRUE, size = 1) +
  geom_smooth(aes(y = vcmax_Alienor_Mean, color = "vcmax_Alienor"), method = "lm", se = TRUE, size = 1) +
  geom_smooth(aes(y = vcmax_Own_Mean, color = "vcmax_Own"), method = "lm", se = TRUE, size = 1) +
  # Customize labels and axis titles
  ylab(expression('Vcmax'*~'('*mu*'mol'*m^-2*s^-1*')')) +
  xlab(expression('MAP'*~'('*'mm'*')')) +
  # Use a minimal theme and customize appearance
  theme_bw()+
  scale_color_manual(
    values = c(
      "vcmax_Constant" = "purple",
      "vcmax_Alienor" = "#65bab7",
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
  annotate("text", x=693, y=30, label="(a)",size=6.5)
# Annotate the plot

vcmax2<-ggplot(site_summary, aes(x = MAP_mm)) +
  # Add points and error bars for vcmax_Obs
  geom_point(aes(y = vcmax25_Constant_Mean, color = "vcmax_Constant")) +
  geom_point(aes(y = vcmax25_Alienor_Mean, color = "vcmax_Alienor")) +
  geom_point(aes(y = vcmax25_Own_Mean, color = "vcmax_Own")) +
  # Add smooth curves for other vcmax types
  geom_smooth(aes(y = vcmax25_Constant_Mean, color = "vcmax_Constant"), method = "lm", se = TRUE, size = 1) +
  geom_smooth(aes(y = vcmax25_Alienor_Mean, color = "vcmax_Alienor"), method = "lm", se = TRUE, size = 1) +
  geom_smooth(aes(y = vcmax25_Own_Mean, color = "vcmax_Own"), method = "lm", se = TRUE, size = 1) +
  # Customize labels and axis titles
  ylab(expression('Vcmax25'*~'('*mu*'mol'*m^-2*s^-1*')')) +
  xlab(expression('MAP'*~'('*'mm'*')')) +
  # Use a minimal theme and customize appearance
  theme_bw()+
  scale_color_manual(
    values = c(
      "vcmax_Constant" = "purple",
      "vcmax_Alienor" = "#65bab7",
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
  annotate("text", x=740, y=120, label="(b)",size=6.5)
# Annotate the plot
vcmax2

plot_grid(vcmax1,vcmax2,label_x = 0.2,nrow = 1)
ggsave(filename = "FigS23_Vcmax_Vcmax25_map.jpg", path = "/Users/echo/E/eechange project/analysis/work2024/ReworkALL/chelsa_2006/plot")

# Summarize only necessary columns
site_summary1 <- vcmaxcombine %>%
  group_by(Site_ID) %>%
  summarize(
    MAP_mm = mean(MAP_mm, na.rm = TRUE),  # Mean MAP for each site
    vcmax_Own_Mean = mean(vcmax_Own, na.rm = TRUE),
    vcmax25_Own_Mean = mean(vcmax25_Own, na.rm = TRUE)
  )

# Plot for Vcmax_Own
vcmax11 <- ggplot(site_summary1, aes(x = MAP_mm)) +
  geom_point(aes(y = vcmax_Own_Mean, color = "vcmax_Own"), size = 1) +
  geom_smooth(aes(y = vcmax_Own_Mean, color = "vcmax_Own"), method = "lm", se = TRUE, size = 1) +
  ylab(expression('Vcmax'*~'('*mu*'mol'*m^-2*s^-1*')')) +
  xlab(expression('MAP'*~'('*'mm'*')')) +
  theme_bw() +
  scale_color_manual(values = c("vcmax_Own" = "red")) +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text = element_text(size = rel(1.5)),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_line(color = "white")
  ) +
  annotate("text", x = 693, y = 30, label = "(a)", size = 6.5)

# Plot for Vcmax25_Own
vcmax21 <- ggplot(site_summary1, aes(x = MAP_mm)) +
  geom_point(aes(y = vcmax25_Own_Mean, color = "vcmax_Own"), size = 1) +
  geom_smooth(aes(y = vcmax25_Own_Mean, color = "vcmax_Own"), method = "lm", se = TRUE, size = 1) +
  ylab(expression('Vcmax25'*~'('*mu*'mol'*m^-2*s^-1*')')) +
  xlab(expression('MAP'*~'('*'mm'*')')) +
  theme_bw() +
  scale_color_manual(values = c("vcmax_Own" = "red")) +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text = element_text(size = rel(1.5)),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_line(color = "white")
  ) +
  annotate("text", x = 740, y = 120, label = "(b)", size = 6.5)

# Combine plots
plot_grid(vcmax11, vcmax21, label_x = 0.2, nrow = 1)

# Save the plots
ggsave(
  filename = "FigS23_Vcmax_Vcmax25_map_one model.jpg",
  path = "/Users/echo/E/eechange project/analysis/work2024/ReworkALL/chelsa_2006/plot"
)
#######compare with other vcmax
##############################vcmax observation from https://zenodo.org/records/6466968
############################################################
vcmaxG<-stack("/Users/echo/E/eechange project/analysis/work2024/vcmax/GOME2_VcmaxTg_05deg.tif")
plot(vcmaxG)##unit:µmol m−2 s−1 Vcmax from GOME-2 SIF

site<-cbind(NECT_sample_climate$Longitude,NECT_sample_climate$Latitude)
obsvcmax<-raster::extract(vcmax,site,method = 'bilinear')
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
write.csv(vcmaxcombine1,"/Users/echo/E/eechange project/analysis/work2024/ReworkALL/chelsa_2006/vcmax_threemodel_obsvcmax.csv")
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
  geom_point(aes(y = vcmax_Own_Mean, color = "vcmax_Own")) +
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
    vcmax25_Own_Mean = mean(vcmax25_Own, na.rm = TRUE)
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

# Plot for Vcmax25_Own
vcmax22 <- ggplot(site_summary2, aes(x = MAP_mm)) +
  geom_point(aes(y = vcmax25_Own_Mean, color = "vcmax_Own"), size = 1) +
  geom_smooth(aes(y = vcmax25_Own_Mean, color = "vcmax_Own"), method = "lm", se = TRUE, size = 1) +
  ylab(expression('Vcmax25'*~'('*mu*'mol'*m^-2*s^-1*')')) +
  xlab(expression('MAP'*~'('*'mm'*')')) +
  theme_bw() +
  scale_color_manual(values = c("vcmax_Own" = "red")) +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text = element_text(size = rel(1.5)),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_line(color = "white")
  ) +
  annotate("text", x = 740, y = 120, label = "(b)", size = 6.5)

# Combine plots
plot_grid(vcmax12, vcmax22, label_x = 0.2, nrow = 1)

# Save the plots
ggsave(
  filename = "FigS23_Vcmax_Vcmax25_map_one model_obs.jpg",
  path = "/Users/echo/E/eechange project/analysis/work2024/ReworkALL/chelsa_2006/plot"
)




vcmax1<-ggplot(site_summary, aes(x = MAP_mm)) +
  # Add points and error bars for vcmax_Obs
  geom_point(aes(y = vcmax_Constant_Mean, color = "vcmax_Constant")) +
  geom_point(aes(y = vcmax_Alienor_Mean, color = "vcmax_Alienor")) +
  geom_point(aes(y = vcmax_Own_Mean, color = "vcmax_Own")) +
  # Add smooth curves for other vcmax types
  geom_smooth(aes(y = vcmax_Constant_Mean, color = "vcmax_Constant"), method = "lm", se = TRUE, size = 1) +
  geom_smooth(aes(y = vcmax_Alienor_Mean, color = "vcmax_Alienor"), method = "lm", se = TRUE, size = 1) +
  geom_smooth(aes(y = vcmax_Own_Mean, color = "vcmax_Own"), method = "lm", se = TRUE, size = 1) +
  # Customize labels and axis titles
  ylab(expression('Vcmax'*~'('*mu*'mol'*m^-2*s^-1*')')) +
  xlab(expression('MAP'*~'('*'mm'*')')) +
  # Use a minimal theme and customize appearance
  theme_bw()+
  scale_color_manual(
    values = c(
      "vcmax_Constant" = "purple",
      "vcmax_Alienor" = "#65bab7",
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
  annotate("text", x=693, y=30, label="(a)",size=6.5)
# Annotate the plot

vcmax2<-ggplot(site_summary, aes(x = MAP_mm)) +
  # Add points and error bars for vcmax_Obs
  geom_point(aes(y = vcmax25_Constant_Mean, color = "vcmax_Constant")) +
  geom_point(aes(y = vcmax25_Alienor_Mean, color = "vcmax_Alienor")) +
  geom_point(aes(y = vcmax25_Own_Mean, color = "vcmax_Own")) +
  # Add smooth curves for other vcmax types
  geom_smooth(aes(y = vcmax25_Constant_Mean, color = "vcmax_Constant"), method = "lm", se = TRUE, size = 1) +
  geom_smooth(aes(y = vcmax25_Alienor_Mean, color = "vcmax_Alienor"), method = "lm", se = TRUE, size = 1) +
  geom_smooth(aes(y = vcmax25_Own_Mean, color = "vcmax_Own"), method = "lm", se = TRUE, size = 1) +
  # Customize labels and axis titles
  ylab(expression('Vcmax25'*~'('*mu*'mol'*m^-2*s^-1*')')) +
  xlab(expression('MAP'*~'('*'mm'*')')) +
  # Use a minimal theme and customize appearance
  theme_bw()+
  scale_color_manual(
    values = c(
      "vcmax_Constant" = "purple",
      "vcmax_Alienor" = "#65bab7",
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
  annotate("text", x=740, y=120, label="(b)",size=6.5)
# Annotate the plot
vcmax2

plot_grid(vcmax1,vcmax2,label_x = 0.2,nrow = 1)
ggsave(filename = "FigS23_Vcmax_Vcmax25_map.jpg", path = "/Users/echo/E/eechange project/analysis/work2024/ReworkALL/chelsa_2006/plot")

# Summarize only necessary columns
site_summary1 <- vcmaxcombine %>%
  group_by(Site_ID) %>%
  summarize(
    MAP_mm = mean(MAP_mm, na.rm = TRUE),  # Mean MAP for each site
    vcmax_Own_Mean = mean(vcmax_Own, na.rm = TRUE),
    vcmax25_Own_Mean = mean(vcmax25_Own, na.rm = TRUE)
  )

# Plot for Vcmax_Own
vcmax11 <- ggplot(site_summary1, aes(x = MAP_mm)) +
  geom_point(aes(y = vcmax_Own_Mean, color = "vcmax_Own"), size = 3) +
  geom_smooth(aes(y = vcmax_Own_Mean, color = "vcmax_Own"), method = "lm", se = TRUE, size = 1) +
  ylab(expression('Vcmax'*~'('*mu*'mol'*m^-2*s^-1*')')) +
  xlab(expression('MAP'*~'('*'mm'*')')) +
  theme_bw() +
  scale_color_manual(values = c("vcmax_Own" = "red")) +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text = element_text(size = rel(1.5)),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_line(color = "white")
  ) +
  annotate("text", x = 693, y = 30, label = "(a)", size = 6.5)

# Plot for Vcmax25_Own
vcmax21 <- ggplot(site_summary1, aes(x = MAP_mm)) +
  geom_point(aes(y = vcmax25_Own_Mean, color = "vcmax_Own"), size = 3) +
  geom_smooth(aes(y = vcmax25_Own_Mean, color = "vcmax_Own"), method = "lm", se = TRUE, size = 1) +
  ylab(expression('Vcmax25'*~'('*mu*'mol'*m^-2*s^-1*')')) +
  xlab(expression('MAP'*~'('*'mm'*')')) +
  theme_bw() +
  scale_color_manual(values = c("vcmax_Own" = "red")) +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text = element_text(size = rel(1.5)),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_line(color = "white")
  ) +
  annotate("text", x = 740, y = 120, label = "(b)", size = 6.5)

# Combine plots
plot_grid(vcmax11, vcmax21, label_x = 0.2, nrow = 1)

# Save the plots
ggsave(
  filename = "FigS23_Vcmax_Vcmax25_map_one model.jpg",
  path = "/Users/echo/E/eechange project/analysis/work2024/ReworkALL/chelsa_2006/plot"
)

vcmax13<-ggplot(site_summary2, aes(x = MAP_mm)) +
  # Add points and error bars for vcmax_Obs
  geom_point(aes(y = vcmax_Constant_Mean, color = "vcmax_Constant")) +
  geom_point(aes(y = vcmax_Alienor_Mean, color = "vcmax_Alienor")) +
  geom_point(aes(y = vcmax_Own_Mean, color = "vcmax_Own")) +
  geom_smooth(aes(y = Obs_vcmax, color = "Obs_vcmax"), method = "lm", se = TRUE, size = 1) +
  geom_point(aes(y = Obs_vcmax, color = "Obs_vcmax")) +
  # Add smooth curves for other vcmax types
  geom_smooth(aes(y = vcmax_Constant_Mean, color = "vcmax_Constant"), method = "lm", se = TRUE, size = 1) +
  geom_smooth(aes(y = vcmax_Alienor_Mean, color = "vcmax_Alienor"), method = "lm", se = TRUE, size = 1) +
  geom_smooth(aes(y = vcmax_Own_Mean, color = "vcmax_Own"), method = "lm", se = TRUE, size = 1) +
  # Customize labels and axis titles
  ylab(expression('Vcmax'*~'('*mu*'mol'*m^-2*s^-1*')')) +
  xlab(expression('MAP'*~'('*'mm'*')')) +
  # Use a minimal theme and customize appearance
  theme_bw()+
  scale_color_manual(
    values = c(
      "Obs_vcmax"="blue",
      "vcmax_Constant" = "purple",
      "vcmax_Alienor" = "#65bab7",
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
  ) 
vcmax13

ggsave(
  filename = "FigS23_Vcmax_threemodel_obs.jpg",
  path = "/Users/echo/E/eechange project/analysis/work2024/ReworkALL/chelsa_2006/plot"
)
###########################################################
####test Vcmax and Vcmax25 in another way
###################In smith 2019, we used the equation about  φ 0=(0.352+0.021*(Tg_K-273.15)-3.4*0.0001*(Tg_K-273.15)^2)/8)
###################But in Wang 2017 it is  φ 0 is the intrinsic quantum yield (1.02 g C mol–1), change unit: 1.02gCmol −1÷12g/mol=0.085molC/mol
######φ0 is the intrinsic quantum yield with a maximum value usually taken to be 0.125
NECT_sample_climate<-read.csv("/Users/echo/E/eechange project/analysis/work2024/ReworkALL/chelsa_2006/Climate_trait_beta.csv")
model_constant_df<-read.csv("/Users/echo/E/eechange project/analysis/work2024/ReworkALL/chelsa_2006/Pmodel_constant_output.csv")
sw<-read.csv("/Users/echo/E/eechange project/analysis/work2024/ReworkALL/chelsa_2006/rsds/sw_w_m2_NECT.csv")
Tg_K<-NECT_sample_climate$tg_chelsa_2006+273.15

Gstar <- pressure*4.275*exp((37830/8.31447)*(1/298.15 - 1/(Tg_K)))/101325    #Gstar in Pa [Bernacchi et al. (2001)], KC, KO and ??* equal to 79430, 36380 and 37830 J/mol# standard values of KC, KO and ??* at 25?? equals to 39.37, 27480 and 4.22 Pa. O
Ca_constant<-model_constant_df$ca
Ci_constant<-model_constant_df$ci
m_constant <- (Ci_constant-Gstar)/(Ci_constant+2*Gstar)

Kc <- 40.49*exp((79430/8.31447)*(1/298.15 - 1/(Tg_K)))                      #Kc in Pa -> 1
Oc <-0.21*pressure #Oc in Pa -> month/daily
Ko <-27840*exp((36380/8.31447)*(1/298.15 - 1/(Tg_K)))                       #Ko in Pa -> month/daily
K  <- Kc*(1+(Oc/Ko)) 

PPFD<-rowMeans(sw[, 2:13], na.rm = TRUE)*0.45/0.22

Vcmax_constant <-(1.02/12)*PPFD*((Ci_constant+K)/(Ci_constant+2*Gstar))*sqrt(1-(0.41/m_constant)^(2/3))
Vcmax25_constant <- Vcmax_constant*exp((65330/8.31447)*((1/Tg_K)-(1/298.15)))
####
model_Alienor_df<-read.csv("/Users/echo/E/eechange project/analysis/work2024/ReworkALL/chelsa_2006/Pmodel_Alienor_output.csv")
Gstar <- pressure*4.275*exp((37830/8.31447)*(1/298.15 - 1/(Tg_K)))/101325    #Gstar in Pa [Bernacchi et al. (2001)], KC, KO and ??* equal to 79430, 36380 and 37830 J/mol# standard values of KC, KO and ??* at 25?? equals to 39.37, 27480 and 4.22 Pa. O
Ca_Alienor<-model_Alienor_df$ca
Ci_Alienor<-model_Alienor_df$ci
m_Alienor <- (Ci_Alienor-Gstar)/(Ci_Alienor+2*Gstar)

Kc <- 40.49*exp((79430/8.31447)*(1/298.15 - 1/(Tg_K)))                      #Kc in Pa -> 1
Oc <-0.21*pressure #Oc in Pa -> month/daily
Ko <-27840*exp((36380/8.31447)*(1/298.15 - 1/(Tg_K)))                       #Ko in Pa -> month/daily
K  <- Kc*(1+(Oc/Ko)) 


Vcmax_Alienor <- (1.02/12)*PPFD*((Ci_Alienor+K)/(Ci_Alienor+2*Gstar))*sqrt(1-(0.41/m_Alienor)^(2/3))
Vcmax25_Alienor <- Vcmax_Alienor*exp((65330/8.31447)*((1/Tg_K)-(1/298.15)))
#####
model_own_df<-read.csv("/Users/echo/E/eechange project/analysis/work2024/ReworkALL/chelsa_2006/Pmodel_ownbeta_output.csv")
Gstar <- pressure*4.275*exp((37830/8.31447)*(1/298.15 - 1/(Tg_K)))/101325    #Gstar in Pa [Bernacchi et al. (2001)], KC, KO and ??* equal to 79430, 36380 and 37830 J/mol# standard values of KC, KO and ??* at 25?? equals to 39.37, 27480 and 4.22 Pa. O
Ca_own<-model_own_df$ca
Ci_own<-model_own_df$ci
m_own <- (Ci_own-Gstar)/(Ci_own+2*Gstar)

Kc <- 40.49*exp((79430/8.31447)*(1/298.15 - 1/(Tg_K)))                      #Kc in Pa -> 1
Oc <-0.21*pressure #Oc in Pa -> month/daily
Ko <-27840*exp((36380/8.31447)*(1/298.15 - 1/(Tg_K)))                       #Ko in Pa -> month/daily
K  <- Kc*(1+(Oc/Ko)) 

Vcmax_own <- (1.02/12)*PPFD*((Ci_own+K)/(Ci_own+2*Gstar))*sqrt(1-(0.41/m_own)^(2/3))
Vcmax25_own <- Vcmax_own*exp((65330/8.31447)*((1/Tg_K)-(1/298.15)))

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
write.csv(vcmaxcombine,"/Users/echo/E/eechange project/analysis/work2024/ReworkALL/chelsa_2006/vcmax_threemodel_f0constant.csv")

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
write.csv(vcmaxcombine1,"/Users/echo/E/eechange project/analysis/work2024/ReworkALL/chelsa_2006/vcmax_threemodel_f0_constant_obsvcmax.csv")
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

# Save the plots
ggsave(
  filename = "FigS23_Vcmax_Vcmax25_map_one model_obs_foconstant.jpg",
  path = "/Users/echo/E/eechange project/analysis/work2024/ReworkALL/chelsa_2006/plot"
)
write.csv(site_summary2,"/Users/echo/E/eechange project/analysis/work2024/ReworkALL/chelsa_2006/vcmax_vcmax25_site_obs_f0_constant_.csv")

####################################################################
###################The light incident on that leaf would (on average) be less than the light incident on the canopy. 
###################The average reduction depends on the LAI. Look at the simple correction for this presented in Dong et al. (2017) Biogoesciences.
##IL ≈ I0.1 − e−kL/=L ≈ I0kfv/ln[1/1 − fv]; fv is fapar
##the difference is the calclulation of Iabs 
##IL ( mol m s ). Vcmax ≈ f0IL(ci + K/ci + 2r*)
##I0 is the incident PAR above the canopy
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

Vcmax_constant <-0.093*I*((Ci_constant+K)/(Ci_constant+2*Gstar))
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

Vcmax_own <- 0.093*I*((Ci_own+K)/(Ci_own+2*Gstar))
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
write.csv(vcmaxcombine,"/Users/echo/E/eechange project/analysis/work2024/ReworkALL/chelsa_2006/vcmax_threemodel_Ning2017.csv")

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
write.csv(vcmaxcombine1,"/Users/echo/E/eechange project/analysis/work2024/ReworkALL/chelsa_2006/vcmax_threemodel_Ning2017_obsvcmax.csv")
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

# Save the plots
ggsave(
  filename = "FigS23_Vcmax_Vcmax25_map_one model_obs_Ning2017.jpg",
  path = "/Users/echo/E/eechange project/analysis/work2024/ReworkALL/chelsa_2006/plot"
)
write.csv(site_summary2,"/Users/echo/E/eechange project/analysis/work2024/ReworkALL/chelsa_2006/vcmax_vcmax25_site_obs_Ning2017.csv")
lm1<-lm(Obs_vcmax~MAP_mm,data = site_summary2)
summary(lm1)
#MAP_mm      -0.031172  Multiple R-squared:  0.4403,
lm2<-lm(vcmax_Own_Mean~MAP_mm,data = site_summary2)
summary(lm2)
#MAP_mm      -0.022731 Multiple R-squared:  0.8803,
lm3<-lm(vcmax25_obs_Mean~MAP_mm,data = site_summary2)
summary(lm3)
##MAP_mm       -0.13594. Multiple R-squared:  0.4979,
lm4<-lm(vcmax25_Own_Mean~MAP_mm,data = site_summary2)
summary(lm4)
##MAP_mm       -0.09863 Multiple R-squared:  0.5554
##################################access another data
# Install the R.matlab package (if not already installed)
if (!require("R.matlab")) {
  install.packages("R.matlab")
}
library(R.matlab)
data <- readMat("/Users/echo/E/eechange project/analysis/work2024/vcmax/TROPOMI_Vmax_Tg_mean.mat")

library(reshape2)
library(ggplot2)

# Extract matrix data
matrix_data <- data$trop.vcTg.mean.05deg

# Generate correct latitude and longitude
lat <- seq(90, -90, length.out = 360)  # Reverse latitude
lon <- seq(-180, 180, length.out = 720)

# Melt the matrix into a long format
df <- melt(matrix_data, varnames = c("lat_index", "lon_index"), value.name = "value")

# Assign actual latitude and longitude
df$latitude <- lat[df$lat_index]
df$longitude <- lon[df$lon_index]

# Remove NA or NaN values
df <- na.omit(df)

# Plot the corrected map
ggplot(df, aes(x = longitude, y = latitude, fill = value)) +
  geom_tile() +
  scale_fill_viridis_c(na.value = "white", name = "Vcmax") +
  coord_equal() +
  labs(title = "Spatial Distribution of Vcmax",
       x = "Longitude", y = "Latitude") +
  theme_minimal()

#############################this data is not good
library(raster)
library(rnaturalearth)
library(sp)

# Load the matrix data
matrix_data <- data$trop.vcTg.mean.05deg  # Replace with your matrix name

# Flip the matrix vertically (reverse rows)
matrix_data_flipped <- matrix_data[nrow(matrix_data):1, ]

# Create the raster object from the flipped matrix
vcmax_raster <- raster(matrix_data_flipped,  
                       xmn = -180, xmx = 180,   # Longitude range
                       ymn = -90, ymx = 90,     # Latitude range
                       crs = "+proj=longlat +datum=WGS84")

# Plot the corrected raster
plot(vcmax_raster, 
     main = "Corrected Spatial Distribution of Vcmax", 
     xlab = "Longitude", ylab = "Latitude", 
     col = terrain.colors(100), 
     axes = TRUE)
vcmax_raster_flipped <- flip(vcmax_raster, direction = 2)  #  'y' or 'x'; or 1 (=x) or 2 (=y)

# Plot the corrected raster after flipping
plot(vcmax_raster_flipped, 
     main = "Corrected Spatial Distribution of Vcmax", 
     xlab = "Longitude", ylab = "Latitude", 
     col = terrain.colors(100))

obsvcmax_tro<-raster::extract(vcmax_raster_flipped,site,method = 'bilinear')
vcmaxcombine1<-data.frame(cbind(NECT_sample_climate$Site.ID.1,NECT_sample_climate$MAP..mm.,Vcmax_constant,Vcmax_Alienor,Vcmax_own,obsvcmax_tro,Vcmax25_constant,Vcmax25_Alienor,Vcmax25_own))
colnames(vcmaxcombine1) <- c(
  "Site_ID", 
  "MAP_mm", 
  "vcmax_Constant", 
  "vcmax_Alienor", 
  "vcmax_Own",
  "obsvcmax_tro",
  "vcmax25_Constant", 
  "vcmax25_Alienor", 
  "vcmax25_Own"
)
#write.csv(vcmaxcombine1,"/Users/echo/E/eechange project/analysis/work2024/ReworkALL/chelsa_2006/vcmax_threemodel_Ning2017_obsvcmax.csv")
vcmaxcombine1$Vcmax25_obs_tro <- vcmaxcombine1$obsvcmax_tro*exp((65330/8.31447)*((1/Tg_K)-(1/298.15)))
# Summarize only necessary columns
site_summary2 <- vcmaxcombine1 %>%
  group_by(Site_ID) %>%
  summarize(
    MAP_mm = mean(MAP_mm, na.rm = TRUE),  # Mean MAP for each site
    vcmax_Constant_Mean = mean(vcmax_Constant, na.rm = TRUE),
    vcmax_Alienor_Mean = mean(vcmax_Alienor, na.rm = TRUE),
    vcmax_Own_Mean = mean(vcmax_Own, na.rm = TRUE),
    Obs_vcmax = mean(obsvcmax_tro, na.rm = TRUE),
    vcmax25_Constant_Mean = mean(vcmax25_Constant, na.rm = TRUE),
    vcmax25_Alienor_Mean = mean(vcmax25_Alienor, na.rm = TRUE),
    vcmax25_Own_Mean = mean(vcmax25_Own, na.rm = TRUE),
    vcmax25_obs_Mean = mean(Vcmax25_obs_tro, na.rm = TRUE)
  )


# Plot for Vcmax_Own
vcmaxtro <- ggplot(site_summary2, aes(x = MAP_mm)) +
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
vcmaxtro
vcmax25tro <- ggplot(site_summary2, aes(x = MAP_mm)) +
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
vcmax25tro
# Combine plots
plot_grid(vcmaxtro, vcmax25tro, label_x = 0.2, nrow = 1)
#################################################################another data is LCC
data <- readMat("/Users/echo/E/eechange project/analysis/work2024/vcmax/LCC_Vcmax_Tg_mean.mat")
matrix_data <- data$Vcmax.LCC.GW.Tg
lat <- seq(90, -90, length.out = 360)  # Reverse latitude
lon <- seq(-180, 180, length.out = 720)
df <- melt(matrix_data, varnames = c("lat_index", "lon_index"), value.name = "value")
df$latitude <- lat[df$lat_index]
df$longitude <- lon[df$lon_index]

df <- na.omit(df)

ggplot(df, aes(x = longitude, y = latitude, fill = value)) +
  geom_tile() +
  scale_fill_viridis_c(na.value = "white", name = "Vcmax") +
  coord_equal() +
  labs(title = "Spatial Distribution of Vcmax",
       x = "Longitude", y = "Latitude") +
  theme_minimal()

#############################this data is not good
library(raster)
library(rnaturalearth)
library(sp)

# Load the matrix data
matrix_data <- data$Vcmax.LCC.GW.Tg  # Replace with your matrix name

# Flip the matrix vertically (reverse rows)
matrix_data_flipped <- matrix_data[nrow(matrix_data):1, ]

# Create the raster object from the flipped matrix
vcmax_raster <- raster(matrix_data_flipped,  
                       xmn = -180, xmx = 180,   # Longitude range
                       ymn = -90, ymx = 90,     # Latitude range
                       crs = "+proj=longlat +datum=WGS84")

# Plot the corrected raster
plot(vcmax_raster, 
     main = "Corrected Spatial Distribution of Vcmax", 
     xlab = "Longitude", ylab = "Latitude", 
     col = terrain.colors(100), 
     axes = TRUE)
vcmax_raster_flipped <- flip(vcmax_raster, direction = 2)  #  'y' or 'x'; or 1 (=x) or 2 (=y)

# Plot the corrected raster after flipping
plot(vcmax_raster_flipped, 
     main = "Corrected Spatial Distribution of Vcmax", 
     xlab = "Longitude", ylab = "Latitude", 
     col = terrain.colors(100))

obsvcmax_LCC<-raster::extract(vcmax_raster_flipped,site,method = 'bilinear')
vcmaxcombine1<-data.frame(cbind(NECT_sample_climate$Site.ID.1,NECT_sample_climate$MAP..mm.,Vcmax_constant,Vcmax_Alienor,Vcmax_own,obsvcmax_LCC,Vcmax25_constant,Vcmax25_Alienor,Vcmax25_own))
colnames(vcmaxcombine1) <- c(
  "Site_ID", 
  "MAP_mm", 
  "vcmax_Constant", 
  "vcmax_Alienor", 
  "vcmax_Own",
  "obsvcmax_LCC",
  "vcmax25_Constant", 
  "vcmax25_Alienor", 
  "vcmax25_Own"
)
#write.csv(vcmaxcombine1,"/Users/echo/E/eechange project/analysis/work2024/ReworkALL/chelsa_2006/vcmax_threemodel_Ning2017_obsvcmax.csv")
vcmaxcombine1$Vcmax25_obs_LCC <- vcmaxcombine1$obsvcmax_LCC*exp((65330/8.31447)*((1/Tg_K)-(1/298.15)))
# Summarize only necessary columns
site_summary2 <- vcmaxcombine1 %>%
  group_by(Site_ID) %>%
  summarize(
    MAP_mm = mean(MAP_mm, na.rm = TRUE),  # Mean MAP for each site
    vcmax_Constant_Mean = mean(vcmax_Constant, na.rm = TRUE),
    vcmax_Alienor_Mean = mean(vcmax_Alienor, na.rm = TRUE),
    vcmax_Own_Mean = mean(vcmax_Own, na.rm = TRUE),
    Obs_vcmax = mean(obsvcmax_LCC, na.rm = TRUE),
    vcmax25_Constant_Mean = mean(vcmax25_Constant, na.rm = TRUE),
    vcmax25_Alienor_Mean = mean(vcmax25_Alienor, na.rm = TRUE),
    vcmax25_Own_Mean = mean(vcmax25_Own, na.rm = TRUE),
    vcmax25_obs_Mean = mean(Vcmax25_obs_LCC, na.rm = TRUE)
  )


# Plot for Vcmax_Own
vcmaxLCC <- ggplot(site_summary2, aes(x = MAP_mm)) +
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
vcmaxLCC
vcmax25LCC <- ggplot(site_summary2, aes(x = MAP_mm)) +
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
vcmax25LCC
# Combine plots
plot_grid(vcmaxLCC, vcmax25LCC, label_x = 0.2, nrow = 1)
#############################################################
##################################
###############################################################
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
###########################################################
##more option from Colin suggestion (email title: Many ways to calculate theoretical Vcmax)
#####if it's constant, it could be 0.085., or 0.111, or 0.125.I used 0.111, and the equation is same as line 1225-1392
####showed as P model use sample version
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

Vcmax_constant <-0.111*I*((Ci_constant+K)/(Ci_constant+2*Gstar))
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


Vcmax_Alienor <- 0.111*I*((Ci_Alienor+K)/(Ci_Alienor+2*Gstar))
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

Vcmax_own <- 0.111*I*((Ci_own+K)/(Ci_own+2*Gstar))
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
write.csv(vcmaxcombine1,"/Users/echo/E/eechange project/analysis/work2024/ReworkALL/chelsa_2006/vcmax_threemodel_Ning2017_obsvcmax_newoption.csv")
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
plot_grid(vcmax12, vcmax22, label_x = 0.2, nrow = 1)
