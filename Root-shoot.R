#############################################Predict total biomass
cbinddata2<-read.csv("treeAGBperarea_with shrubs.csv")
library(dplyr)

gymnosperm_genera <- c(
  "Abies", "Araucaria", "Callitris", "Cathaya", "Cedrus", "Cephalotaxus",
  "Chamaecyparis", "Cryptomeria", "Cunninghamia", "Cupressus", "Cycas",
  "Dacrycarpus", "Dacrydium", "Ephedra", "Fitzroya", "Ginkgo", "Juniperus",
  "Keteleeria", "Larix", "Metasequoia", "Microbiota", "Nageia", "Picea",
  "Pinus", "Platycladus", "Podocarpus", "Pseudolarix", "Pseudotsuga",
  "Sciadopitys", "Sequoia", "Sequoiadendron", "Taxodium", "Taxus", "Thuja",
  "Thujopsis", "Torreya", "Tsuga", "Wollemia", "Zamia"
)

angiosperm_genera <- c(
  "Acer", "Aesculus", "Alnus", "Betula", "Carpinus", "Castanea", "Celtis",
  "Cerasus", "Cinnamomum", "Dipterocarpus", "Eucalyptus", "Fagus",
  "Fraxinus", "Juglans", "Liquidambar", "Liriodendron", "Magnolia",
  "Mallotus", "Nyssa", "Olea", "Ostrya", "Paulownia", "Phellodendron",
  "Platanus", "Populus", "Prunus", "Pyrus", "Quercus", "Rhamnus",
  "Robinia", "Salix", "Sorbus", "Styrax", "Tilia", "Ulmus", "Zelkova"
)

df2 <- cbinddata2 %>%
  filter(!is.na(tree.height), !is.na(DBH)) %>%
  mutate(
    plant_group = case_when(
      Genus %in% gymnosperm_genera ~ "gymnosperm",
      Genus %in% angiosperm_genera ~ "angiosperm",
      TRUE ~ NA_character_
    )
  )


split_sizes <- strsplit(df2$crown.diameter, "\\*")  # Split by "*"

RSimprove<-read.csv("Climate_R_S_NECT_predictRS_biomass.csv")

###Wtot  =  n W  = (b0/c0) D.
library(dplyr)

df2 <- df2 %>%
  mutate(
    b0 = case_when(
      plant_group == "angiosperm" ~ 177.74,
      plant_group == "gymnosperm" ~ 92.09,
      TRUE ~ NA_real_
    ),
    c0 = case_when(
      plant_group == "angiosperm" ~ 8.31,
      plant_group == "gymnosperm" ~ 1.82,
      TRUE ~ NA_real_
    ),
    pretreebiomass_KgCm2 = (b0 / c0) * (DBH / 100),
    pretreebiomass_gCm2  = pretreebiomass_KgCm2 * 1000
  )

mean_biomass <- df2 %>%
  group_by(Site) %>%
  summarise(mean_pretreebiomass_gCm2 = mean(pretreebiomass_gCm2, na.rm = TRUE))
RSimprove1 <- RSimprove %>%
  left_join(mean_biomass, by = c("Site" = "Site"))

RSimprove1 <- RSimprove1 %>%
  left_join(
    df2 %>% dplyr::select(Site, Genus, Species, tree.height),
    by = c(
      "Site.ID" = "Site",
      "ACCEPTED.GENUS" = "Genus",
      "ACCEPTED.SPECIES" = "Species"
    )
  )
RSimprove1 <- RSimprove1 %>%
  mutate(pre_totalbiomass = ifelse(
    is.nan(mean_pretreebiomass_gCm2.y),  # Condition: Check for NaN
    total.biomass_gC.m.2,              # Use this if NaN
    mean_pretreebiomass_gCm2.y           # Keep original value otherwise
  ))
dataheight<-read_xlsx("/Volumes/LaCie/E/eechange project/data/NECT2006_quadrats_latin_height_cleaned.xlsx")

library(purrr)

set.seed(123)

# 1. prepare RSimprove2
RSimprove2 <- RSimprove1 %>%
  mutate(
    row_id = row_number(),
    Site.ID = trimws(as.character(Site.ID)),
    ACCEPTED.GENUS = iconv(as.character(ACCEPTED.GENUS), from = "", to = "UTF-8", sub = ""),
    ACCEPTED.SPECIES = iconv(as.character(ACCEPTED.SPECIES), from = "", to = "UTF-8", sub = ""),
    species_name = trimws(paste(ACCEPTED.GENUS, ACCEPTED.SPECIES)),
    tree.height_fill = tree.height
  )

# 2. prepare dataheight
dataheight2 <- dataheight %>%
  mutate(
    SiteID = trimws(as.character(SiteID)),
    species_name = iconv(as.character(species_name), from = "", to = "UTF-8", sub = ""),
    species_name = trimws(species_name)
  ) %>%
  filter(!is.na(height_m))

# 3. pool height values by SiteID + species_name
height_pool <- dataheight2 %>%
  group_by(SiteID, species_name) %>%
  summarise(height_values = list(height_m), .groups = "drop")

# 4. join the pool to RSimprove2
RSimprove2 <- RSimprove2 %>%
  left_join(
    height_pool,
    by = c("Site.ID" = "SiteID", "species_name" = "species_name")
  )

# 5. only fill missing tree.height_fill, keep original tree.height unchange

RSimprove2 <- RSimprove2 %>%
  group_by(Site.ID, species_name) %>%
  group_modify(~{
    df_group <- .x
    vals <- df_group$height_values[[1]]
    
    if (is.null(vals) || length(vals) == 0 || all(is.na(vals))) {
      return(df_group)
    }
    
    miss_idx <- which(is.na(df_group$tree.height_fill))
    n_miss <- length(miss_idx)
    
    if (n_miss == 0) {
      return(df_group)
    }
    
    vals_sorted <- sort(vals, decreasing = TRUE)
    
    if (length(vals_sorted) >= n_miss) {
      sampled_vals <- vals_sorted[1:n_miss]
    } else {
      sampled_vals <- vals_sorted
    }
    
    df_group$tree.height_fill[miss_idx[1:length(sampled_vals)]] <- sampled_vals
    df_group
  }) %>%
  ungroup() %>%
  dplyr::select(-height_values)
woodyRS<-subset(RSimprove2,RSimprove2$species=="woody") ####281

lmwoody <- lm(
  log(R_S) ~ tg_1987_2006 + log(AI) + log(RZwc) + log(GPP) + Sand +
    tg_1987_2006 * log(AI) + log(tree.height_fill) + lnRRD +
    log(SLA..m2.kg.) + lnthickness + log(LDMC..mg.g./1000),
  data = woodyRS %>%
    filter(
      !is.na(R_S), R_S > 0,
      !is.na(AI), AI > 0,
      !is.na(RZwc), RZwc > 0,
      !is.na(GPP), GPP > 0,
      !is.na(tree.height_fill), tree.height_fill > 0,
      !is.na(SLA..m2.kg.), SLA..m2.kg. > 0,
      !is.na(LDMC..mg.g.), LDMC..mg.g. > 0,
      !is.na(tg_1987_2006),
      !is.na(Sand),
      !is.na(lnRRD),
      !is.na(lnthickness)
    )
)

summary(lmwoody)
nonwoodRS<-subset(RSimprove2,RSimprove2$species=="non-woody")##611
lmnonwoody<-lm(log(R_S)~tg_1987_2006+log(AI)+log(GPP)+pH+Sand+lnheightT+lnRRD+lnSRL+lnthickness+log(LDMC..mg.g./1000),data=nonwoodRS)
summary(lmnonwoody)
woodyRS$pred_lnRS <- predict(lmwoody, newdata = woodyRS)

# Get predictions from non-woody model
nonwoodRS$pred_lnRS <- predict(lmnonwoody, newdata = nonwoodRS)

# Combine the predictions back into original dataframe
RSimprove2$pred_lnRSnew <- NA  # Initialize column
RSimprove2$pred_lnRSnew[RSimprove2$species == "woody"] <- woodyRS$pred_lnRS
RSimprove2$pred_lnRSnew[RSimprove2$species == "non-woody"] <- nonwoodRS$pred_lnRS

# Convert back from log scale if needed
RSimprove2$pred_RS_new <- exp(RSimprove2$pred_lnRSnew)

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
  annotate("text", x = 655, y = 0, label = "(e)", size = 6.5)

pbRSgroup
########
RSimprove2$preRMF_new<-RSimprove2$pred_RS_new/(RSimprove2$pred_RS_new+1)
RSimprove2$prebelow_new<-RSimprove2$preRMF_new*RSimprove2$total.biomass_gC.m.2

pbBB<-ggplot(RSimprove2,aes(x=MAP..mm., y=BGB_gC.m.2)) + geom_point() + 
  geom_smooth(aes(x=MAP..mm., y=BGB_gC.m.2),method = 'lm', formula = y ~ x)+
  xlab(expression('MAP'*~'('*'mm'*')')) + ylab(expression('BGB'*~'('*'gC'*~m^-2*')'))+
  geom_point(aes(x=MAP..mm.,y=prebelow_new),color="red") + 
  geom_smooth(aes(y=prebelow_new,x=MAP..mm., col='red'),method = 'lm') +
  theme_bw()+
  theme(legend.position="none",axis.title.x=element_text(size=15),axis.title.y=element_text(size=15),
        axis.text=element_text(size=rel(1.5)),panel.grid.major=element_line(color="white"),panel.grid.minor=element_line(color="white"))+
  annotate("text", x=655, y=8000, label="(e)",size=6.5)
pbBB

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
  annotate("text", x = 655, y = 10, label = "(f)", size = 6.5)

pbBBgroup

RSimprove2$preAGB_new<-RSimprove2$total.biomass_gC.m.2-RSimprove2$prebelow_new
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
  annotate("text", x = 655, y = 10, label = "(g)", size = 6.5)

pbAGBgroup
write.csv(RSimprove2,"Climate_R_S_NECT_predictRS_biomass_BAADb0c0.csv")
