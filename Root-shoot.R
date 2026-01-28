#############################################root:shoot
library(raster)
NECT_sample_climate<-read.csv("/Volumes/LaCie/data_chelsa/NECT_Climate_trait.csv")

setwd("/Volumes/LaCie/data_chelsa/tmx1987_2006")
tmaxfile <- list.files(pattern = "*.tif") 
tmaxfile_stack <- stack(tmaxfile)           # Load TIFF files into a RasterStack
site<-cbind(NECT_sample_climate$Longitude,NECT_sample_climate$Latitude)
tmax_NECT<-raster::extract(tmaxfile_stack,site,method = 'bilinear')
tmax_NECT<-data.frame(tmax_NECT)
tmax_C_NECT<- (tmax_NECT/10)-273.15##K/10 ,change unit from K to C
write.csv(tmax_C_NECT,"/Volumes/LaCie/data_chelsa/tmx1987_2006/tmax_C_NECT_1987_2006.csv")

setwd("/Volumes/LaCie/data_chelsa/tmin1987_2006")
tminfile <- list.files(pattern = "*.tif") 
tminfile_stack <- stack(tminfile)           # Load TIFF files into a RasterStack
tmin_NECT<-raster::extract(tminfile_stack,site,method = 'bilinear')
tmin_NECT<-data.frame(tmin_NECT)
tmin_C_NECT<- (tmin_NECT/10)-273.15##K/10 ,change unit from K to C
write.csv(tmin_C_NECT,"tmin_C_NECT__1987_2006.csv")

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
write.csv(Tg_NECT,"/Volumes/LaCie/data_chelsa/Tg_chelsa_1987_2006.csv")
#############################################################
setwd("/Volumes/LaCie/data_chelsa/vpd1987_2006")
vpdfile <- list.files(pattern = "*.tif") 
vpdfile_stack <- stack(vpdfile)           # Load TIFF files into a RasterStack
site<-cbind(NECT_sample_climate$Longitude,NECT_sample_climate$Latitude)
vpd_NECT<-raster::extract(vpdfile_stack,site,method = 'bilinear')
vpd_NECT<-data.frame(vpd_NECT)
vpd_Pa_NECT<-vpd_NECT*0.1 ##scale:0.1
write.csv(vpd_Pa_NECT,"vpd_Pa_NECT_1987_2006.csv")

#rsds Surface downwelling shortwave flux in air. MJ m-2 d-1  scale 0.001
setwd("/Volumes/LaCie/data_chelsa/rsds1987_2006")
rsdsfile <- list.files(pattern = "*.tif") 
rsdsfile_stack <- stack(rsdsfile)           # Load TIFF files into a RasterStack
plot(rsdsfile_stack)
rsds_NECT<-raster::extract(rsdsfile_stack,site,method = 'bilinear')
rsds_NECT<-data.frame(rsds_NECT)
sw<- rsds_NECT*11.574*0.001####MJ m-2 d-1 to w/m2; W m−2= MJ m−2d −1×1,000,000/86,400; scale for the data 0.001
write.csv(sw,"sw_w_m2_NECT_1987_2006.csv")

ppfd_NECT<-60*60*24*10^-6*sw[, 1:12]*2.04##molm-2d-1
ppfd_NECT<-data.frame(ppfd_NECT)
write.csv(ppfd_NECT,"ppfd_molm-2d-1_NECT_1987_2006.csv")
#####################################fapar
fapar<-stack('/Volumes/LaCie/data_chelsa/0.05fapar/add_time/fapar_1987_2006.nc')
plot(fapar)
fapar_NECT<-raster::extract(fapar,site,method = 'bilinear')
fapar_NECT<-data.frame(fapar_NECT)
write.csv(fapar_NECT,"/Volumes/LaCie/data_chelsa/0.05fapar/fapar1987_2006_NECT.csv")
################################################################# run rsplash 2.0

library(devtools)
##error install rsplash: topmodel’ is not available for package ‘rsplash’
#topmodle is moved from CRAN, install topmodel by hand 
#install.packages("/Users/echo/E/eechange project/analysis/work2024/topmodel_0.7.5.tar.gz", repos = NULL, type = "source")
library(topmodel)
#install_github("dsval/rsplash")
library(rsplash)
library(xts)

setwd("/Volumes/LaCie/data_chelsa/pr1987_2006")
pnfile <- list.files(pattern = "*.tif") 
pnfile_stack <- stack(pnfile)           # Load TIFF files into a RasterStack
site<-cbind(NECT_sample_climate$Longitude,NECT_sample_climate$Latitude)
pn_NECT<-raster::extract(pnfile_stack,site,method = 'bilinear')
pn_NECT<-data.frame(pn_NECT)
pn_mm_NECT<-pn_NECT/100 ##kg m-2 month1/100 kg/m²/month is numerically equivalent to mm/month
write.csv(pn_mm_NECT,"pn_mm_NECT.csv")

setwd("/Users/echo/E/eechange project/analysis/work2024/ReworkALL/chelsa_1987_2006/")
Tg_NECT<-read.csv("Tg_chelsa_1987_2006.csv")
vpd_Pa_NECT<-read.csv("vpd_Pa_NECT_1987_2006.csv")
ppfd_NECT<-read.csv("ppfd_molm-2d-1_NECT_1987_2006.csv")
fapar_NECT<-read.csv("fapar1987_2006_NECT.csv")
pn_NECT<-read.csv("pn_mm_NECT.csv")
sw<-read.csv("sw_w_m2_NECT_1987_2006.csv")
soildata<-read.csv("~/E/eechange project/data/splash_P_hydro/dailydata/soil_texture_from_splashtools.csv")
# Assuming `NECT_sample_climate` and `soildata` are your datasets

# Match the Site.ID in NECT_sample_climate with X in soildata
soildata_expanded <- soildata[match(NECT_sample_climate$Site.ID, soildata$X), ]
#soildata_expanded$Site.ID <- NECT_sample_climate$Site.ID #just check
Result_sw_in<-t(sw)
Result_sw_in<-Result_sw_in[-1,]
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
Result_tc<- Result_tc[-1,]
Result_tc <- data.frame(V1 = formatted_dates, Result_tc)
rownames(Result_tc)<-Result_tc[,1]
Result_tcasdata<-as.data.frame(Result_tc)
Result_tclist<- as.list(Result_tcasdata)
Result_tclist.ts <- xts(Result_tclist$X1, order.by=as.POSIXct(Result_tclist$V1))
Result_tclist.ts 

Result_pn<-t(pn_NECT)
Result_pn<-Result_pn[-1,]
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
write.csv(all_wn_results_df, "splash_output/wn_mm_1987_2006.csv", row.names = FALSE)
all_wn_results_m3<-all_wn_results_df
all_wn_results_m3[,2:795]<-((all_wn_results_df[,2:795]/1000)/2)##depth in soil data is 2, This results in a soil moisture content in volume percent
write.csv(all_wn_results_m3, "splash_output/wn_m3m3_1987_2006.csv", row.names = FALSE)

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
  file_path <- paste0("splash_output/", variable_name, "_1987_2006.csv")
  
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
#############################################################################
sm_lm<-read.csv("/Users/echo/E/eechange project/analysis/work2024/ReworkALL/chelsa_1987_2006/splash_output/sm_lim_1987_2006.csv")
library(lubridate)

sm_lm <- sm_lm %>%
  mutate(Date = as.Date(Date))

# Add Year-Month column
sm_lm_data <- sm_lm %>%
  mutate(YearMonth = floor_date(Date, "month")) %>%
  group_by(YearMonth) %>%
  summarise(across(starts_with("value.sm_lim"), sum, na.rm = TRUE))
sm_lm_data[,2:795][sm_lm_data[,2:795]>1]<-1

sm_lm_data1<-t(sm_lm_data)
sm_lm_data1<-sm_lm_data1[-1,]
sm_lm_data1<-data.frame(sm_lm_data1)

sm_lm_data1 <- sm_lm_data1 %>%
  mutate(across(1:12, as.numeric))


NECT_sample_climate$tg_1987_2006 <- rowMeans(Tg_NECT[, 2:13], na.rm = TRUE) ##794
NECT_sample_climate$vpd_1987_2006<- rowMeans(vpd_Pa_NECT[, 2:13], na.rm = TRUE)
NECT_sample_climate$fapar_1987_2006<-rowMeans(fapar_NECT[, 2:13], na.rm = TRUE)
NECT_sample_climate$ppfd_1987_2006<-rowMeans(ppfd_NECT[, 2:13], na.rm = TRUE)
NECT_sample_climate$sm_1987_2006 <- rowMeans(sm_lm_data1[, 1:12], na.rm = TRUE) ##794
beta<-146
library(rpmodel)
pressure<-calc_patm(elv=NECT_sample_climate$Elevation, patm0 = 101325)
################################################################run P model with soil moisture limitation
model_sm<-rpmodel(
  tc=NECT_sample_climate$tg_1987_2006,
  vpd=NECT_sample_climate$vpd_1987_2006,
  co2=381.9,
  fapar=NECT_sample_climate$fapar_1987_2006,
  ppfd=NECT_sample_climate$ppfd_1987_2006,
  patm =pressure ,
  elv=NECT_sample_climate$Elevation,
  kphio =   0.049977,
  beta=146 ,
  soilm = NECT_sample_climate$sm_1987_2006,
  meanalpha = 1,
  apar_soilm = 0,
  bpar_soilm = 0.733,
  c4 = FALSE,
  method_jmaxlim = "wang17",
  do_ftemp_kphio = TRUE,
  do_soilmstress = TRUE,
  returnvar = NULL,
  verbose = FALSE
)

model_sm_df <- as.data.frame(model_sm)
head(model_sm_df)
# Save the dataframe to a CSV file
write.csv(model_sm_df, "/Users/echo/E/eechange project/analysis/work2024/ReworkALL/chelsa_1987_2006/Pmodel_sm_constanbeta.csv", row.names = FALSE)

##AI annual potential evapotranspiration to annual precipitation
pet<-read.csv("/Users/echo/E/eechange project/analysis/work2024/ReworkALL/chelsa_1987_2006/splash_output/pet_1987_2006.csv")

pet <- pet %>%
  mutate(Date = as.Date(Date))

# Add Year-Month column
pet_data <- pet %>%
  mutate(Year = year(Date)) %>% # Extract the year from the Date column
  group_by(Year) %>% # Group by year
  summarise(across(starts_with("value.pet"), sum, na.rm = TRUE)) # Summarize by summing values

pet_data1<-t(pet_data)
pet_data1<-pet_data1[-1,]
pet_data1<-data.frame(pet_data1)

NECT_sample_climate$pet<-pet_data1$pet_data1
########use MAP as the original data
NECT_sample_climate$AI<-NECT_sample_climate$pet/NECT_sample_climate$MAP..mm.

NECT_sample_climate$pr_annual_mm<-rowSums(pn_NECT[, 2:13], na.rm = TRUE)

NECT_sample_climate$pr_mean_1987_2006<-rowMeans(pn_NECT[, 2:13], na.rm = TRUE)
################################################
########################root zone water capacity
#data foe GEE
siteupload<-data.frame(cbind(NECT_sample_climate$SAMPLE.ID,NECT_sample_climate$Longitude,NECT_sample_climate$Latitude))
names(siteupload)<-c("SAMPLE.ID","Longitude","Latitude")
write.csv(siteupload,"/Users/echo/E/eechange project/analysis/work2024/ReworkALL/chelsa_1987_2006/Site_NECT_GEE.csv")
##MOD16A2.061: Terra Net Evapotranspiration 8-Day Global 500m only from 2001 year
##I use 2001 to 2006
########################################RZWC
######china dataset
##calculate the rzwc
##ET from google earth engine
et_other1<-read.csv("/Users/echo/E/eechange project/analysis/work2024/ReworkALL/chelsa_1987_2006/ET_NECT_mod16_01_06.csv",header=T)
et_other1<-et_other1[order(et_other1$SAMPLEID),]

etdata<-data.frame(cbind(et_other1$SAMPLEID,et_other1[,2:276]))
etdata1<-t(etdata)
colnames(etdata1)<-etdata1[1,]
etdata1<-etdata1[-1,]
etdata1<-data.frame(etdata1)

ET1<-(etdata1[,1:794])/8

rownames(ET1) <- as.Date(gsub("_ET", "", gsub("X", "", rownames(ET1))), format = "%Y_%m_%d")

# Verify the row names
head(rownames(ET1))

# If you also want the row names as a new column in the data:
ET1$Date <- rownames(ET1)

# Reorder the columns if needed (move 'Date' to the first column)
ET1 <- ET1[, c(ncol(ET1), 1:(ncol(ET1) - 1))]

######################################################################the way from Rodolfo to make 8 day data into everyday,then monthly
library(lubridate)
library(dplyr)
library(zoo)

# Now, let's create a daily continuous time series fo dates (matching the dataset period)
dates_daily <- as.Date(seq(from = as.Date("2001-01-01"), to = as.Date("2006-12-31"), by = "1 day"))

# Transforming the date column into a date object
ET1$Date <- as.Date(ET1$Date)

# Let's merge the 8-day dataset and the daily time series, the dates that have no data will be automatically "NA"
new_ET <- merge(dates_daily, ET1, by.x = "x",by.y = "Date", all = TRUE) 

##new_ET[2:1908] <- na.locf(new_ET[,2:1908], fromLast = FALSE) the reult is NA,I do not why

num<-NULL
num<-rep(1:46,each=8)
data1<-read.csv("day_for_year.csv")
number<-NULL
for (i in 1:dim(data1)[1]) {
  number<-c(number,num[1:data1[i,2]])
  
}
number<-data.frame(number)
daydata<-data.frame(cbind(new_ET,number))

#sapply(dayjoin, class)  
library(dplyr)
library(tidyr)

daydata<-daydata %>% 
  group_by(number) %>%
  fill(X1:X794, .direction = 'down') %>%
  ungroup()

write.csv(daydata,"ET_NECT_mod16_)_01_06_daily.csv")
################################
daydata$yearmonth <- strftime(daydata$x, "%Y-%m")
daydata<-daydata[,-796]
#daydata$yearmonth <- as.Date(daydata$yearmonth)
monthlyET<-aggregate(daydata[,2:795], by=list(daydata$yearmonth), sum,na.rm = TRUE)

##new_ET_monthly <- aggregate(x = daydata[,2:17550], by = list(month=daydata$yearmonth), FUN = sum)
colnames(monthlyET)[1] <- "date"
write.csv(monthlyET,"ET_NECT_mod16_01_06_monthly.csv")
###################################################################

prec_sites <- read.csv("/Users/echo/E/eechange project/analysis/work2024/ReworkALL/chelsa_1987_2006/pn_mm_NECT_2001_2006.csv",header=TRUE,row.names = 1)
ET <- read.csv("ET_NECT_mod16_01_06_monthly.csv",header=TRUE,row.names = 1)
prec_sites<-t(prec_sites)
prec_sites<-as.data.frame(prec_sites)
prec_sites <- prec_sites %>%
  mutate(across(2:795, as.numeric))

# Example: Assume ET1 is your data frame with row names
# Extract the row names
rownames_data <- rownames(prec_sites)
# Extract the year and month from the row names
date_column <- gsub("CHELSA_pr_(\\d{2})_(\\d{4})_.*", "\\2-\\1", rownames_data)
# Add the extracted date column as the first column in the data frame
prec_sites <- cbind(date = date_column, prec_sites)

prec_sites<-as.data.frame(prec_sites)
prec_sites <- prec_sites %>%
  mutate(across(2:795, as.numeric))

#install.packages("dplyr")
library(dplyr)
get_month <- function(series) {
  # this performs
  # decomposition
  i <- ts(series, frequency = 12) # we want to capture a 12-month pattern
  d <- decompose(i, type = c("additive"), filter = NULL) # additive approach is always a good practice
  seasonality <- as.vector(d$figure)
  seasonality <- replace(seasonality, seasonality >= 0, 1) # this highlights the rainy season
  seasonality <- replace(seasonality, seasonality < 0, 0) # this highlights the dry season
  seasonality <- with(as.data.frame(seasonality), ifelse(seasonality == 1 & dplyr::lag(seasonality, n = 1, default = seasonality[12]) == 0, 1, 0)) # this masks the start of the rainy season
  rainy_month <- max(match(1,seasonality))
  return(rainy_month)
}

get_rainy_month <- function(series) {
  r <- list()
  for (i in colnames(series)[2:length(colnames(series))]) { 
    d <- subset(series, select = c('date', i))
    r[[i]] <- get_month(d)
  }
  r
}

##  Water deficit calculation 
prec_sites_monthly <- prec_sites
subset_et <- ET
rainy_month <- get_rainy_month(prec_sites_monthly)

rainy_month_reset <- rainy_month
do_hydrological_calc <- TRUE
use_runoff_regression <- FALSE 
use_runoff_data <- FALSE

sr <- matrix(ncol = (ncol(prec_sites_monthly) - 1), nrow = nrow(prec_sites_monthly)) 
sr <- data.frame(sr)
names(sr) <- names(subset_et[2:ncol(subset_et)])
for (i in 2:ncol(prec_sites_monthly)) {
  for (j in 1:nrow(prec_sites_monthly)) {
    if (use_runoff_regression) {
      runoff_site <- prec_sites_monthly[j,i] * runoff_results$runoff_slope[(i - 1)] # Estimating runoff based on adhoc regression
      if (runoff_site < 0) { runoff_site <- 0 }
    } else {
      if (use_runoff_data) {
        runoff_site <- subset_runoff[j,i]
      } else {
        runoff_site <- 0
      }
    }
    
    if (rainy_month[i - 1] == 1 | !(do_hydrological_calc)) {
      day_reset <- 13 # 13 is the first month of the year
      if ((j %% day_reset) == 0) { # The day_reset (= 13) means to reset the water balance at the beginning of each year
        sr[j,(i - 1)] <- max(subset_et[j,i] + runoff_site - prec_sites_monthly[j,i] , 0) # + subset_runoff[j,i]
      } else {
        sr[j,(i - 1)] <- max(subset_et[j,i] + runoff_site - prec_sites_monthly[j,i]   + sr[j - 1,(i - 1)], 0) # + subset_runoff[j,i]
      }
    } else {
      day_reset <- as.numeric(rainy_month_reset[i - 1])
      if ((j == day_reset | j %% (day_reset + 12) == 0)) { # The day_reset here comes from the previous loop calculation
        sr[j,(i - 1)] <- max(subset_et[j,i] + runoff_site - prec_sites_monthly[j,i] , 0) # + subset_runoff[j,i]
      } else {
        sr[j,(i - 1)] <- max(subset_et[j,i] + runoff_site - prec_sites_monthly[j,i]   + sr[j - 1,(i - 1)], 0) # + subset_runoff[j,i]
      }
    }
  }
}
write.csv(sr,"sr_2001_2006__monthly.csv")
# A quick test to check if precipitation and ET values are reliable
# We are going to run through a selection of points that Ruijie has done
# In general, as rule of thumb, we expect the ET in dry sites to be around 500 mm (+/- 300 mm)
# and in the Wet areas around 1000 mm (+/- 300 mm) - Of course there are many exceptions
# Climate, Site Nr, Lon, Lat, Country
rzwc_site <- data.frame(apply(sr,2,max))

write.csv(rzwc_site,"sr_2001_2006_right_root_zone_water_capacity.csv")

NECT_sample_climate$RZwc<-rzwc_site$apply.sr..2..max.

###################sand an pH

NECT_sample_climate$Sand<-soildata_expanded$sand
PH<-stack("/Volumes/rd619/home/soil_fertility/PHH2O1.tif")

PH_NECT<-raster::extract(PH,site,method = 'bilinear')
PH_NECT<-data.frame(PH_NECT)
PH_NECT<-PH_NECT/10
NECT_sample_climate$pH<-PH_NECT$PH_NECT

PFTdata<-read.csv("/Users/echo/E/eechange project/analysis/work2024/ecy2091-sup-0002-datas1/PFT data.csv")
PFTdataNECT <- subset(PFTdata, SAMPLE.ID >= 1 & SAMPLE.ID <= 794) ##794
NECT_sample_climate$life.form<-PFTdataNECT$Life.form
NECT_sample_climate <- NECT_sample_climate %>%
  mutate(
    species = ifelse(grepl("shrub|tree", life.form, ignore.case = TRUE), "woody", "non-woody")
  )

biomassdata<-read.csv("/Users/echo/E/eechange project/data/biomass1.14/siteAGBperarea1.21 with total biomass.csv")

dataselect<-biomassdata[match(NECT_sample_climate$Site.ID, biomassdata$Site), ]
wholedata<-data.frame(cbind(NECT_sample_climate,dataselect))
wholedata$R_S<-wholedata$BGB_MgC.ha/wholedata$AGB_MgC.ha
wholedata$GPP<-model_sm_df$gpp
lmtwotypesnoin<-lm(log(R_S)~tg_1987_2006+log(AI)+log(RZwc)+log(GPP)+pH+Sand+tg_1987_2006*log(AI)+species,data=wholedata)
summary(lmtwotypesnoin)
table(wholedata$species)
# non-woody     woody 
# 611       183
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)          -0.253999   0.571718  -0.444  0.65697    
# tg_1987_2006          0.137095   0.066348   2.066  0.03913 *  
#   log(AI)               1.931206   0.591117   3.267  0.00113 ** 
#   log(RZwc)             1.006881   0.101708   9.900  < 2e-16 ***
#   log(GPP)             -1.004996   0.325602  -3.087  0.00210 ** 
#   pH                   -0.184787   0.059178  -3.123  0.00186 ** 
#   Sand                 -0.164950   0.004197 -39.304  < 2e-16 ***
#   specieswoody          0.016320   0.060682   0.269  0.78805    
# tg_1987_2006:log(AI) -0.363330   0.059381  -6.119 1.49e-09 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.6786 on 785 degrees of freedom
# Multiple R-squared:  0.8894,	Adjusted R-squared:  0.8882 
# F-statistic: 788.8 on 8 and 785 DF,  p-value: < 2.2e-16
wholedata$prelnR_S<-predict(lmtwotypesnoin, data = wholedata)
wholedata$preR_S <- exp(wholedata$prelnR_S)

pbln2<-ggplot(wholedata,aes(x=MAP..mm., y=R_S)) + geom_point() + 
  geom_smooth(method = 'lm', formula = y ~ x)+
  xlab(expression('MAP'*~'('*'mm'*')')) + ylab(expression('R:S'))+
  geom_point(aes(x=MAP..mm.,y=preR_S),color="red") + 
  geom_smooth(aes(y=preR_S,x=MAP..mm., col='red'),method = 'lm') +
  theme_bw()+
  theme(legend.position="none",axis.title.x=element_text(size=15),axis.title.y=element_text(size=15),
        axis.text=element_text(size=rel(1.5)),panel.grid.major=element_line(color="white"),panel.grid.minor=element_line(color="white"))+
  annotate("text", x=655, y=1, label="(e)",size=6.5)
pbln2

rsmodelobs<-lm(log(R_S)~MAP..mm.,data=wholedata)
summary(rsmodelobs)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -7.6217003  0.1430755  -53.27   <2e-16 ***
#   MAP..mm.     0.0092345  0.0003088   29.90   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 1.392 on 792 degrees of freedom
# Multiple R-squared:  0.5303,	Adjusted R-squared:  0.5297 
# F-statistic: 894.2 on 1 and 792 DF,  p-value: < 2.2e-16
rsmodelpre<-lm(log(preR_S)~MAP..mm.,data=wholedata)
summary(rsmodelpre)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -7.6058255  0.1258215  -60.45   <2e-16 ***
#   MAP..mm.     0.0091980  0.0002716   33.87   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 1.224 on 792 degrees of freedom
# Multiple R-squared:  0.5916,	Adjusted R-squared:  0.5911 
# F-statistic:  1147 on 1 and 792 DF,  p-value: < 2.2e-16
write.csv(wholedata,"/Users/echo/E/eechange project/analysis/work2024/ReworkALL/chelsa_1987_2006/Climate_R_S_NECT.csv")

########################################visreg of plot
wholedata<-read.csv("/Users/echo/E/eechange project/analysis/work2024/ReworkALL/chelsa_1987_2006/Climate_R_S_NECT.csv")
wholedata$log_AI <- log(wholedata$AI)
wholedata$log_RZwc <- log(wholedata$RZwc)
wholedata$log_GPP <- log(wholedata$GPP)

# Refit the model with the precomputed variables
lmtwotypesnoin <- lm(log(R_S) ~ tg_1987_2006 + log_AI + log_RZwc + log_GPP + pH + Sand + tg_1987_2006 * log_AI + species, 
                     data = wholedata)
summary(lmtwotypesnoin)
f1t_types1<-visreg(lmtwotypesnoin,"tg_1987_2006",ylab="log(R_S)",line=c(col="darkred"))

f1t_typesplot1<-ggplot(data=f1t_types1$res, aes(x=tg_1987_2006, y=visregRes))+
  geom_point(size = 1.5) +
  #scale_color_viridis(option = 'inferno', direction = -1, breaks = c(100, 500, 1000, 1500,2000,3000,4000,5000), begin = 0.25,limits=c(floor(0), ceiling(5000)))+
  geom_ribbon(data=f1t_types1$fit,aes(x=tg_1987_2006, y=visregFit,ymin = visregLwr, ymax = visregUpr), alpha = 0.5)+
  geom_smooth(data=f1t_types1$fit, aes(x=tg_1987_2006, y=visregFit),col="dodgerblue3",method='lm',se=TRUE,size=2)+
  xlab(expression('Tg'*" ("*paste(degree,C)*")")) +
  theme_bw() +theme( legend.position = "none",axis.title.x=element_text(size=25),axis.title.y=element_blank(),
                     axis.text=element_text(size=rel(1.5)),panel.grid.major=element_line(color="white"),panel.grid.minor=element_line(color="white"))+
  annotate("text", x=12, y=0, label="(a)",size=10)
f1t_typesplot1

f2t_types1<-visreg(lmtwotypesnoin,"log_RZwc",ylab="ln(Root:Shoot)",line=c(col="deeppink"))

f2t_typesplot1<- ggplot(data=f2t_types1$res, aes(x=log_RZwc, y=visregRes))+
  geom_point(size = 1.5) +
  #scale_color_viridis(option = 'inferno', direction = -1, breaks = c(100, 500, 1000, 1500,2000,3000,4000,5000),begin = 0.25,limits=c(floor(0), ceiling(5000)))+
  geom_ribbon(data=f2t_types1$fit,aes(x=log_RZwc, y=visregFit,ymin = visregLwr, ymax = visregUpr), alpha = 0.5)+
  geom_smooth(data=f2t_types1$fit, aes(x=log_RZwc, y=visregFit),col="dodgerblue3",method='lm',se=TRUE,size=2)+
  xlab(expression('lnRZ'[wc] *~'(mm)')) +
  theme_bw() +theme( legend.position = "none",axis.title.x=element_text(size=25),axis.title.y=element_blank(),
                     axis.text=element_text(size=rel(1.5)),panel.grid.major=element_line(color="white"),panel.grid.minor=element_line(color="white"))+
  annotate("text", x=5.5, y=0, label="(b)",size=10)
f2t_typesplot1

f3t_types1<-visreg(lmtwotypesnoin,"log_GPP",ylab="ln(Root:Shoot)",line=c(col="deeppink"))
f3t_typesplot1<- ggplot(data=f3t_types1$res, aes(x=log_GPP, y=visregRes))+ 
  geom_point(size = 1.5) +
  #scale_color_viridis(option = 'inferno', direction = -1,breaks = c(100, 500, 1000, 1500,2000,3000,4000,5000), begin = 0.25,limits=c(floor(0), ceiling(5000)))+
  geom_ribbon(data=f3t_types1$fit,aes(x=log_GPP, y=visregFit,ymin = visregLwr, ymax = visregUpr), alpha = 0.5)+
  geom_smooth(data=f3t_types1$fit, aes(x=log_GPP, y=visregFit),col="dodgerblue3",method='lm',se=TRUE,size=2)+
  xlab(expression('lnGPP (gC'*~'m'^-2*~'d'^-1*')')) +
  theme_bw() +theme( legend.position = "none",axis.title.x=element_text(size=25),axis.title.y=element_blank(),
                     axis.text=element_text(size=rel(1.5)),panel.grid.major=element_line(color="white"),panel.grid.minor=element_line(color="white"))+
  annotate("text", x=0.5, y=0, label="(c)",size=10)
f3t_typesplot1

f4t_types1<-visreg(lmtwotypesnoin,"pH",ylab="ln(Root:Shoot)",line=c(col="deeppink"))

f4t_typesplot1<- ggplot(data=f4t_types1$res, aes(x=pH, y=visregRes))+ 
  geom_point(size = 1.5) +
  #scale_color_viridis(option = 'inferno', direction = -1, breaks = c( 100, 500, 1000, 1500,2000,3000,4000,5000), begin = 0.25,limits=c(floor(0), ceiling(5000)))+
  geom_ribbon(data=f4t_types1$fit,aes(x=pH, y=visregFit,ymin = visregLwr, ymax = visregUpr), alpha = 0.5)+
  geom_smooth(data=f4t_types1$fit, aes(x=pH, y=visregFit),col="dodgerblue3",method='lm',se=TRUE,size=2)+
  xlab(expression('Soil pH')) +
  theme_bw() +theme( legend.position = "none",axis.title.x=element_text(size=25),axis.title.y=element_blank(),
                     axis.text=element_text(size=rel(1.5)),panel.grid.major=element_line(color="white"),panel.grid.minor=element_line(color="white"))+
  annotate("text", x=8.5, y=0, label="(d)",size=10)
f4t_typesplot1

f5t_types1<-visreg(lmtwotypesnoin,"Sand",ylab="ln(Root:Shoot)",line=c(col="deeppink"))
f5t_typesplot1<- ggplot(data=f5t_types1$res, aes(x=Sand, y=visregRes))+ 
  geom_point(size = 1.5) +
  #scale_color_viridis(option = 'inferno', direction = -1, breaks = c( 100, 500, 1000, 1500,2000,3000,4000,5000), begin = 0.25,limits=c(floor(0), ceiling(5000)))+
  geom_ribbon(data=f5t_types1$fit,aes(x=Sand, y=visregFit,ymin = visregLwr, ymax = visregUpr), alpha = 0.5)+
  geom_smooth(data=f5t_types1$fit, aes(x=Sand, y=visregFit),col="dodgerblue3",method='lm',se=TRUE,size=2)+
  xlab(expression('Sand (% of weight)')) +
  theme_bw() +theme( legend.position = "none",axis.title.x=element_text(size=25),axis.title.y=element_blank(),
                     axis.text=element_text(size=rel(1.5)),panel.grid.major=element_line(color="white"),panel.grid.minor=element_line(color="white"))+
  annotate("text", x=54, y=0, label="(e)",size=10)
f5t_typesplot1


f6t_types1<-visreg(lmtwotypesnoin,"species",ylab="ln(Root:Shoot)",line=c(col="deeppink"))
f6t_types1$res$x_numeric<-as.numeric(f6t_types1$res$species)
f6t_types1$fit$x_numeric<-as.numeric(f6t_types1$fit$species)
#f6t_types1$res<-na.omit(f6t_types1$res)

f6t_typesplot1<- ggplot(data=f6t_types1$res, aes(x=x_numeric, y=visregRes))+ 
  geom_point(size = 1.5,position = "jitter") +
  #scale_color_viridis(option = 'inferno', direction = -1, breaks = c( 1000,2000,3000,6000,9000,13000), begin = 0.25,limits=c(floor(0), ceiling(13000)))+
  #geom_ribbon(data=f6t_types1$fit,aes(x=as.factor(Types), y=visregFit,ymin = visregLwr, ymax = visregUpr), alpha = 0.8)+
  #geom_smooth(data=f6t_types1$fit, aes(x=x_numeric, y=visregFit, group = as.factor(Types)),col="deeppink",method='lm',se=TRUE,size=1)+
  stat_summary(fun= mean, geom = "crossbar", width = .5, color = "dodgerblue3")+
  #scale_x_discrete(labels = c("1" = "Woody", "2" = "Herbaceous")) +
  xlab(expression('Types')) +
  theme_bw() +theme( legend.position = "none",axis.title.x=element_text(size=25),axis.title.y=element_blank(),axis.text.x = element_blank(),
                     axis.text=element_text(size=rel(1.8)),panel.grid.major=element_line(color="white"),panel.grid.minor=element_line(color="white"))+
  annotate("text", x=2.3, y=0, label="(g)",size=10)+
  #facet_wrap(~ Types, ncol = 1, scales = "free_x") +
  geom_smooth(aes(group = species), method = 'lm', se = TRUE, size = 2, col = "dodgerblue3")
f6t_typesplot1

f7t_types1<-visreg(lmtwotypesnoin,"log_AI",ylab="ln(Root:Shoot)",line=c(col="darkred"))

f7t_typesplot1<-ggplot(data=f7t_types1$res, aes(x=log_AI, y=visregRes))+
  geom_point(size = 1.5) +
  #scale_color_viridis(option = 'inferno', direction = -1, breaks = c(100, 500, 1000, 1500,2000,3000,4000,5000), begin = 0.25,limits=c(floor(0), ceiling(5000)))+
  geom_ribbon(data=f7t_types1$fit,aes(x=log_AI, y=visregFit,ymin = visregLwr, ymax = visregUpr), alpha = 0.5)+
  geom_smooth(data=f7t_types1$fit, aes(x=log_AI, y=visregFit),col="dodgerblue3",method='lm',se=TRUE,size=2)+
  xlab(expression('lnAI')) +
  theme_bw() +theme( legend.position = "none",axis.title.x=element_text(size=25),axis.title.y=element_blank(),
                     axis.text=element_text(size=rel(1.5)),panel.grid.major=element_line(color="white"),panel.grid.minor=element_line(color="white"))+
  annotate("text", x=2, y=0, label="(f)",size=10)
f7t_typesplot1


test1<-visreg(lmtwotypesnoin,"tg_1987_2006",ylab="ln(Root:Shoot)",by="log_AI",layout=c(3,1))

test1plot<-ggplot(data=test1$res, aes(x=tg_1987_2006, y=visregRes))+
  geom_point(size = 1.5) +
  facet_wrap(~log_AI,labeller = labeller(log_AI = 
                                           c("0.375694792420571" = "ln AI: 0.38",
                                             "0.973676501472261" = "ln AI: 0.97",
                                             "1.68462843670931" = "ln AI: 1.68")))+
  #scale_color_viridis(option = 'inferno', direction = -1,breaks = c( 100, 500, 1000, 1500,2000,3000,4000,5000), begin = 0.25,limits=c(floor(0), ceiling(5000)))+
  geom_ribbon(data=test1$fit,aes(x=tg_1987_2006, y=visregFit,ymin = visregLwr, ymax = visregUpr), alpha = 0.5)+
  geom_smooth(data=test1$fit, aes(x=tg_1987_2006, y=visregFit),col="dodgerblue3",method='lm',se=TRUE,size=2)+
  xlab(expression('Tg'*" ("*paste(degree,C)*")")) +
  theme_bw() +theme( legend.position = "none",axis.title.x=element_text(size=25),axis.title.y=element_blank(),
                     axis.text=element_text(size=rel(1.5)),panel.grid.major=element_line(color="white"),panel.grid.minor=element_line(color="white"),
                     strip.text = element_text(size = 12)) #+
#annotate("text", x=30, y=3.8, label="(a)",size=5)
test1plot

plot_grid(
  test1plot, f2t_typesplot1, f3t_typesplot1, 
  f4t_typesplot1,f5t_typesplot1, f7t_typesplot1, f6t_typesplot1,
  ncol = 4, nrow = 2
)

#######################predict BGB
wholedata<-read.csv("/Users/echo/E/eechange project/analysis/work2024/ReworkALL/chelsa_1987_2006/Climate_R_S_NECT.csv")
wholedata$preRMF<-wholedata$preR_S/(wholedata$preR_S+1)
wholedata$prebelow<-wholedata$preRMF*wholedata$total.biomass_gC.m.2

pbBB<-ggplot(wholedata,aes(x=MAP..mm., y=BGB_gC.m.2)) + geom_point() + 
  geom_smooth(aes(x=MAP..mm., y=BGB_gC.m.2),method = 'lm', formula = y ~ x)+
  xlab(expression('MAP'*~'('*'mm'*')')) + ylab(expression('BGB'*~'('*'gC'*~m^-2*')'))+
  geom_point(aes(x=MAP..mm.,y=prebelow),color="red") + 
  geom_smooth(aes(y=prebelow,x=MAP..mm., col='red'),method = 'lm') +
  theme_bw()+
  theme(legend.position="none",axis.title.x=element_text(size=15),axis.title.y=element_text(size=15),
        axis.text=element_text(size=rel(1.5)),panel.grid.major=element_line(color="white"),panel.grid.minor=element_line(color="white"))+
  annotate("text", x=655, y=8000, label="(e)",size=6.5)
pbBB

wholedata$sitegroup <- ifelse(wholedata$BGB.tree_gC.m.2 > 0, "treesite", "grasssite")

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

write.csv(wholedata,"/Users/echo/E/eechange project/analysis/work2024/ReworkALL/chelsa_1987_2006/Climate_R_S_NECT_2.0.csv")
###############################change the R:S prediction, because we change the model include the trait data
RSimprove<-read.csv("/Users/echo/E/eechange project/analysis/work2024/ReworkALL/chelsa_1987_2006/Climate_R_S_NECT_2.0.csv")
heightT<-rast("/Users/echo/PhD_third/Taje_trait_map/X3106_mean_Shrub_Tree_Grass_1km.tif")
plot(heightT)
crs(heightT)
height_raster <- project(heightT, "EPSG:4326")
plot(height_raster)
site_point<-cbind(RSimprove$Longitude,RSimprove$Latitude)
heightTsite<-extract(height_raster,site_point,method='bilinear')
heightTsite<-data.frame(heightTsite)
heightTsite$Plant.height..veg....mean.[heightTsite$Plant.height..veg....mean.<0]<-NA
RSimprove$heightT <- heightTsite$Plant.height..veg....mean.
RSimprove$lnheightT<-log(RSimprove$heightT)

SLA<-rast("/Users/echo/PhD_third/Taje_trait_map/X3117_mean_Shrub_Tree_Grass_1km.tif")
SLA <- project(SLA, "EPSG:4326")
plot(SLA)
SLAsite<-raster::extract(SLA,site_point,method='bilinear')
SLAsite<-data.frame(SLAsite)
RSimprove$SLA <- SLAsite$SLA..mean.
RSimprove$lnSLA<-log(RSimprove$SLA)

SRL<-rast("/Users/echo/PhD_third/Taje_trait_map/X1080_mean_Shrub_Tree_Grass_1km.tif")
SRL <- project(SRL, "EPSG:4326")
plot(SRL)
SRLsite<-raster::extract(SRL,site_point,method='bilinear')
SRLsite<-data.frame(SRLsite)
RSimprove$SRL <- SRLsite$SRL..mean.
RSimprove$lnSRL<-log(RSimprove$SRL)

RRD<-rast("/Users/echo/PhD_third/Taje_trait_map/X6_mean_Shrub_Tree_Grass_1km.tif")
RRD <- project(RRD, "EPSG:4326")
plot(RRD)
RRDsite<-raster::extract(RRD,site_point,method='bilinear')
RRDsite<-data.frame(RRDsite)
RSimprove$RRD <- RRDsite$RRD..mean.
RSimprove$lnRRD<-log(RSimprove$RRD)

thickness<-rast("/Users/echo/PhD_third/Taje_trait_map/X46_mean_Shrub_Tree_Grass_1km.tif")
thickness <- project(thickness, "EPSG:4326")
plot(thickness)
thicknesssite<-raster::extract(thickness,site_point,method='bilinear')
thicknesssite<-data.frame(thicknesssite)
RSimprove$thickness <- thicknesssite$Leaf.thickness..mean.
RSimprove$lnthickness<-log(RSimprove$thickness)

LDMC<-rast("/Users/echo/PhD_third/Taje_trait_map/X47_mean_Shrub_Tree_Grass_1km.tif")
LDMC <- project(LDMC, "EPSG:4326")
plot(LDMC)
LDMCsite<-raster::extract(LDMC,site_point,method='bilinear')
LDMCsite<-data.frame(LDMCsite)
RSimprove$LDMC <- LDMCsite$LDMC..mean.
RSimprove$lnLDMC<-log(RSimprove$LDMC)

###for woody
woodyRS<-subset(RSimprove,RSimprove$species=="woody") 
lmwoody<-lm(log(R_S)~tg_1987_2006+log(AI)+log(RZwc)+log(GPP)+Sand+tg_1987_2006*log(AI)+lnheightT+lnRRD+lnSLA+lnthickness+lnLDMC,data=woodyRS)
summary(lmwoody)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)          27.903434  13.651656   2.044 0.042491 *  
#   tg_1987_2006          0.080279   0.105927   0.758 0.449571    
# log(AI)              -0.810657   1.306261  -0.621 0.535693    
# log(RZwc)             1.158424   0.136413   8.492 9.49e-15 ***
#   log(GPP)             -1.669358   0.497746  -3.354 0.000982 ***
#   Sand                 -0.186174   0.013406 -13.887  < 2e-16 ***
#   lnheightT             0.008138   0.166580   0.049 0.961094    
# lnRRD                -0.929691   0.411129  -2.261 0.025000 *  
#   lnSLA                -1.018472   0.936533  -1.087 0.278350    
# lnthickness           4.198963   1.112095   3.776 0.000220 ***
#   lnLDMC                8.057084   5.981458   1.347 0.179759    
# tg_1987_2006:log(AI) -0.284079   0.105175  -2.701 0.007609 ** 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.4202 on 171 degrees of freedom
# Multiple R-squared:  0.9435,	Adjusted R-squared:  0.9399 
# F-statistic: 259.7 on 11 and 171 DF,  p-value: < 2.2e-16
nonwoodRS<-subset(RSimprove,RSimprove$species=="non-woody")
lmnonwoody<-lm(log(R_S)~tg_1987_2006+log(AI)+log(GPP)+pH+Sand+lnheightT+lnRRD+lnSRL+lnthickness+lnLDMC,data=nonwoodRS)
summary(lmnonwoody)
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  -2.802644   7.041744  -0.398 0.690769    
# tg_1987_2006 -0.212455   0.022449  -9.464  < 2e-16 ***
#   log(AI)       0.165110   0.309011   0.534 0.593320    
# log(GPP)      0.906878   0.342638   2.647 0.008340 ** 
#   pH           -0.240342   0.077846  -3.087 0.002112 ** 
#   Sand         -0.165805   0.006284 -26.385  < 2e-16 ***
#   lnheightT     0.837852   0.192454   4.354 1.58e-05 ***
#   lnRRD        -1.171537   0.331310  -3.536 0.000437 ***
#   lnSRL         2.532292   0.322058   7.863 1.75e-14 ***
#   lnthickness   2.241561   0.801108   2.798 0.005306 ** 
#   lnLDMC        3.171943   2.854762   1.111 0.266968    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.7144 on 600 degrees of freedom
# Multiple R-squared:  0.8753,	Adjusted R-squared:  0.8732 
# F-statistic: 421.2 on 10 and 600 DF,  p-value: < 2.2e-16
# Get predictions from woody model
woodyRS$pred_lnRS <- predict(lmwoody, newdata = woodyRS)

# Get predictions from non-woody model
nonwoodRS$pred_lnRS <- predict(lmnonwoody, newdata = nonwoodRS)

# Combine the predictions back into original dataframe
RSimprove$pred_lnRSnew <- NA  # Initialize column
RSimprove$pred_lnRSnew[RSimprove$species == "woody"] <- woodyRS$pred_lnRS
RSimprove$pred_lnRSnew[RSimprove$species == "non-woody"] <- nonwoodRS$pred_lnRS

# Convert back from log scale if needed
RSimprove$pred_RS_new <- exp(RSimprove$pred_lnRSnew)

pbRSgroup <- ggplot(RSimprove, aes(x = MAP..mm.)) +
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
  annotate("text", x = 655, y = 10, label = "(e)", size = 6.5)

pbRSgroup
########################################################3
##predict aboveground biomass
##total biomass
##turnover time from turnover rate
cbinddata2<-read.csv("/Users/echo/E/eechange project/data/biomass1.14/treeAGBperarea_with shrubs.csv")
###assume as aginosperm
##dbh cm, width m, H, m
##width (冠幅, CW) = (longest diameter + perpendicular diameter) / 2
##Crown Area=π×(a/2)×(b/2)
# Split the column into major (a) and minor (b) axes
split_sizes <- strsplit(cbinddata2$crown.diameter, "\\*")  # Split by "*"

# Convert to numeric and store as separate columns
cbinddata2$a <- sapply(split_sizes, function(x) as.numeric(x[1]))  # Major axis (a)
cbinddata2$b <- sapply(split_sizes, function(x) as.numeric(x[2]))  # Minor axis (b)

cbinddata2$crown_area_m2 <- pi * (cbinddata2$a / 2) * (cbinddata2$b / 2)
cbinddata2$DH_m2<-(cbinddata2$DBH/100)*cbinddata2$tree.height

crownareamodel <- lm(crown_area_m2 ~ DH_m2 - 1, data =cbinddata2)
summary(crownareamodel) ##c0=4.17
# Call:
#   lm(formula = crown_area_m2 ~ DH_m2 - 1, data = cbinddata2)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -18.280  -0.849   1.513   4.420  74.388 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# DH_m2   4.1732     0.2378   17.55   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 8.281 on 207 degrees of freedom
# (1469 observations deleted due to missingness)
# Multiple R-squared:  0.5981,	Adjusted R-squared:  0.5961 
# F-statistic:   308 on 1 and 207 DF,  p-value: < 2.2e-16
rmse<-sqrt(mean(crownareamodel$residuals^2)) ##8.2606

crownareamodelplot <- ggplot(data = cbinddata2, aes(x = DH_m2, y = crown_area_m2)) +
  geom_point(size = 2,alpha = 0.4) +
  xlab(expression(D * H ~ (m^2))) +
  ylab(expression(Tree~crown~area ~ (m^2))) +
  theme_bw() +
  theme(
    legend.position = "none",
    legend.title=element_text(size=13), 
    legend.text=element_text(size=15),
    axis.title.x = element_text(size = 25),
    axis.title.y = element_text(size = 25),
    axis.text = element_text(size = rel(2)),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_line(color = "white")
  )+
  geom_smooth( method = 'lm', se = TRUE, size = 2)

print(crownareamodelplot)

cbinddata2$D2H_m3<-(cbinddata2$DBH/100)^2*cbinddata2$tree.height
cbinddata2$treebiomss_Kg<-cbinddata2$treeagb_area*100/1000 ##Kg
biomassmodel<-lm((treebiomss_Kg)*0.47~D2H_m3-1, data = cbinddata2) ##Kg to KgC
summary(biomassmodel) ##b0=156.79
# Call:
#   lm(formula = (treebiomss_Kg) * 0.47 ~ D2H_m3 - 1, data = cbinddata2)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -75.988  -5.612  -0.166   0.005  50.808 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# D2H_m3  156.789      2.032   77.15   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 15.18 on 178 degrees of freedom
# (1498 observations deleted due to missingness)
# Multiple R-squared:  0.971,	Adjusted R-squared:  0.9708 
# F-statistic:  5952 on 1 and 178 DF,  p-value: < 2.2e-16
rmse<-sqrt(mean(biomassmodel$residuals^2)) ##15.1393

biomassmodelplot <- ggplot(data = cbinddata2, aes(x=D2H_m3, y=treebiomss_Kg*0.47)) +
  geom_point(size = 2,alpha = 0.4) +
  xlab(expression(D^2~H~(m^3)))+ylab(expression('Tree biomass (Kg C)'))+
  theme_bw() +
  theme(
    legend.position = "none",
    legend.title=element_text(size=13), 
    legend.text=element_text(size=15),
    axis.title.x = element_text(size = 25),
    axis.title.y = element_text(size = 25),
    axis.text = element_text(size = rel(2)),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_line(color = "white")
  )+
  geom_smooth( method = 'lm', se = TRUE, size = 2)

print(biomassmodelplot)
###Wtot  =  n W  = (b0/c0) D.
b0=156.79
c0=4.17
cbinddata2$pretreebiomass_KgCm2<-(b0/c0)*(cbinddata2$DBH/100)
cbinddata2$pretreebiomass_gCm2<-((b0/c0)*(cbinddata2$DBH/100))*1000

mean_biomass <- cbinddata2 %>%
  group_by(Site) %>%
  summarise(mean_pretreebiomass_gCm2 = mean(pretreebiomass_gCm2, na.rm = TRUE))

RSimprove1 <- RSimprove %>%
  left_join(mean_biomass, by = c("Site" = "Site"))

RSimprove1 <- RSimprove1 %>%
  mutate(pre_totalbiomass = ifelse(
    is.nan(mean_pretreebiomass_gCm2),  # Condition: Check for NaN
    total.biomass_gC.m.2,              # Use this if NaN
    mean_pretreebiomass_gCm2           # Keep original value otherwise
  ))
woodbiomass<- RSimprove1 %>% 
  filter(!is.na(mean_pretreebiomass_gCm2))  # Removes NA and NaN
pbtotolbiomass<-ggplot(woodbiomass,aes(x=MAP..mm., y=log(total.biomass_gC.m.2))) + geom_point() + 
  geom_smooth(aes(x=MAP..mm., y=log(total.biomass_gC.m.2)),method = 'lm', formula = y ~ x)+
  xlab(expression('MAP'*~'('*'mm'*')')) + ylab(expression('lnTotal biomass'*~'('*'gC'*~m^-2*')'))+
  geom_point(aes(x=MAP..mm.,y=log(pre_totalbiomass)),color="red") + 
  geom_smooth(aes(y=log(pre_totalbiomass),x=MAP..mm., col='red'),method = 'lm') +
  theme_bw()+
  theme(legend.position="none",axis.title.x=element_text(size=15),axis.title.y=element_text(size=15),
        axis.text=element_text(size=rel(1.5)),panel.grid.major=element_line(color="white"),panel.grid.minor=element_line(color="white"))
pbtotolbiomass


RSimprove1$preRMF_new<-RSimprove1$pred_RS_new/(RSimprove1$pred_RS_new+1)
RSimprove1$prebelow_new<-RSimprove1$preRMF_new*RSimprove1$total.biomass_gC.m.2

pbBB<-ggplot(RSimprove1,aes(x=MAP..mm., y=BGB_gC.m.2)) + geom_point() + 
  geom_smooth(aes(x=MAP..mm., y=BGB_gC.m.2),method = 'lm', formula = y ~ x)+
  xlab(expression('MAP'*~'('*'mm'*')')) + ylab(expression('BGB'*~'('*'gC'*~m^-2*')'))+
  geom_point(aes(x=MAP..mm.,y=prebelow_new),color="red") + 
  geom_smooth(aes(y=prebelow_new,x=MAP..mm., col='red'),method = 'lm') +
  theme_bw()+
  theme(legend.position="none",axis.title.x=element_text(size=15),axis.title.y=element_text(size=15),
        axis.text=element_text(size=rel(1.5)),panel.grid.major=element_line(color="white"),panel.grid.minor=element_line(color="white"))+
  annotate("text", x=655, y=8000, label="(e)",size=6.5)
pbBB

pbBBgroup <- ggplot(RSimprove1, aes(x = MAP..mm.)) +
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
  annotate("text", x = 655, y = 10, label = "(f)", size = 6.5)

pbBBgroup

RSimprove1$preAGB_new<-RSimprove1$total.biomass_gC.m.2-RSimprove1$prebelow_new
pbAGBgroup <- ggplot(RSimprove1, aes(x = MAP..mm.)) +
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
  annotate("text", x = 655, y = 10, label = "(g)", size = 6.5)

pbAGBgroup


write.csv(RSimprove1,"/Users/echo/E/eechange project/analysis/work2024/ReworkALL/chelsa_1987_2006/Climate_R_S_NECT_predictRS_biomass.csv")
