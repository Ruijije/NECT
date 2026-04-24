load("3 clim.trait.RData")
##################################################
wangCPT<-climate.trait %>% filter(grepl("CPT", source))

# Filter rows where the 'Site' column contains values from 1 to 33
wangCPT_filtered <- wangCPT %>%
  filter(site >= 1 & site <= 33)

wangCPT_filtered$fD <- ifelse(wangCPT_filtered$DecEv == "D", wangCPT_filtered$f, NA)
wangCPT_filtered$fE <- ifelse(wangCPT_filtered$DecEv == "E", wangCPT_filtered$f, NA)
############################################################re-fitted the LMA model evergreen use Wang2023, decedeous use self data fitted
NECT_sample_climate<-read.csv("NECT_Climate_trait.csv")
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
select<-data.frame(ClimatecombineNECT$Site.ID,ClimatecombineNECT$MAP..mm.)
select_filtered <- select %>%
  filter(!ClimatecombineNECT.Site.ID %in% c(11, 12, 13, 14, 19))

lmadata<-data.frame(cbind(select_filtered$ClimatecombineNECT.MAP..mm.,lnpre$x,lma.sd$x,lmaobsmean$x))
names(lmadata) <- c("map", "lnprelma", "lma.sd","lmaobsmean")
##################################Predict Narea
vamaxdata<-read.csv("vcmax_threemodel_obsvcmax.csv")
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


