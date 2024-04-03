# File: NRSA_compile_variables.R
# Purpose: Compile NRSA 0809 dataset with variables of interest. Add in StreamCat variables of interest, to match NRSA sites.
# Date: December 12, 2022
# Beginning of code aided by Donald Benkendorf at EPA.

# All data and metadata are available from: https://www.epa.gov/national-aquatic-resource-surveys/data-national-aquatic-resource-surveys
library(tidyverse)
library(dplyr) 
library(ggpubr)

# See https://rstudio-pubs-static.s3.amazonaws.com/646006_232ccb16922444ff9dacb2441f38d9c3.html for instructions on downloading and using StreamCatTools.
library(StreamCatTools)
library(sf)
library(nhdplusTools)

### Siteinfo ---- 08/09 -----

# Read in siteinfo from NARS website that hosts all NARS data
NRSA0809_siteinfo <- read.csv('https://www.epa.gov/sites/default/files/2015-09/siteinfo_0.csv')

# View variables
dplyr::glimpse(NRSA0809_siteinfo)
# US_L3CODE useful to filter down to Level III Ecoregions of interest
# US_L3CODE2015 in NRSA0809

# Select variables of interest
NRSA0809_siteinfo_df <- NRSA0809_siteinfo %>% dplyr::select(UID, SITE_ID, VISIT_NO, DATE_COL, US_L3CODE_2015, LAT_DD83, LON_DD83, RT_NRSA, AGGR_ECO9_2015)
# METADATA: RT_NRSA	= Reference/So-so/Trashed designations for sites used in NRSA	

# Chemistry ----

# Read in chem data.
NRSA0809_chem <- read.csv('https://www.epa.gov/sites/default/files/2015-09/chem.csv')
#pH in 0809 is in this document, under "PHLAB"
# View variables
dplyr::glimpse(NRSA0809_chem)

# Select variables of interest
NRSA0809_chem_df <- NRSA0809_chem %>% dplyr::select(UID, COND, K, SODIUM, SO4, MG, CL, CA, AL, PTL, NTL, TURB, PHLAB)
# METADATA: ALERT	= "R-Value below reporting limit for parameter, S-Sample shipping time exceeded two day limit, H-Holding time of sample exceeded parameter-specific limit"	
# METADATA: NTL	= total nitrogen
# METADATA: PTL	= total phosphorus


# Physical habitat ----
# Read in phab data
NRSA0809_phab <- read.csv('https://www.epa.gov/sites/default/files/2015-09/phablow.csv')

# View variables
dplyr::glimpse(NRSA0809_phab)

# Select variables of interest
NRSA0809_phab_df <- NRSA0809_phab %>% dplyr::select(UID, PCT_FN, LSUB_DMM, XSLOPE, XEMBED)
# METADATA: PCT_FN = "Bed Surface % Fines <0.06mm"

# Landscape metrics ----

# This file has basin slope variables and PRISM temperature variables.
# Read in landscape data
NRSA0809_land <- read.csv('https://www.epa.gov/sites/default/files/2015-09/land.csv')
# View variables
dplyr::glimpse(NRSA0809_land)

# Select variables of interest
NRSA0809_land_df <- NRSA0809_land %>% dplyr::select(UID, NHDWAT_SLOPE)

# Merge data ----

# Merge 0809 NRSA siteinfo, chemistry, phab, and landscape variables into one dataset. 
merge1 <- dplyr::left_join(NRSA0809_siteinfo, NRSA0809_chem_df, by = c("UID"))
merge2 <- dplyr::left_join(merge1, NRSA0809_phab_df, by = c("UID"))
merge3 <- dplyr::left_join(merge2, NRSA0809_land_df, by = c("UID"))

# Get 0809 COMIDs to use to obtain StreamCat variables.
# Promote data frame to sf spatial points data frame
nrsa_sf <- st_as_sf(merge3, coords = c("LON_DD83", "LAT_DD83"), crs = 4269)
# COMID is already in 2013/2014, col 16 of merge3

# Get COMIDs using nhdplusTools package
nrsa_sf$COMID<- NA
for (i in 1:nrow(nrsa_sf)){
  print (i)
  nrsa_sf[i,'COMID'] <- discover_nhdplus_id(nrsa_sf[i,c('geometry')])
}

# Temperature data from StreamCat 
# Create a vector of the COMIDs in 0809, pass into the StreamCat function to get MAST, MSST, MWST for 08
COMID <- as.character(nrsa_sf$COMID)
COMID <- paste(COMID, collapse = ',')
sc_get_params(param='comid') #tells us which predictors I can get w/ comid, also ensures I'm matching how predictor is formatted
SCtemp08 <- sc_get_data(metric = 'mast_2008, msst_2008, mwst_2008', aoi='other', comid = COMID) 

# Merge nrsa_sf and SCtemp08 to get full NRSA 08 df
NRSA_08 <- dplyr::left_join(nrsa_sf, SCtemp08, by = c("COMID"))
rm(NRSA0809_chem, NRSA0809_chem_df, NRSA0809_land, NRSA0809_land_df, NRSA0809_phab, NRSA0809_phab_df, NRSA0809_siteinfo, 
   NRSA0809_siteinfo_df, nrsa_sf, SCtemp08)
#write.csv(NRSA_08, 'NRSA0809_compiled.csv')


### Siteinfo 13/14 -----
# Read in siteinfo from NARS website that hosts all NARS data
NRSA1314_siteinfo <- read.csv('https://www.epa.gov/sites/default/files/2019-04/nrsa1314_siteinformation_wide_04292019.csv')

# View variables
dplyr::glimpse(NRSA1314_siteinfo)
# US_L3CODE useful to filter down to Level III Ecoregions of interest

# Select variables of interest
NRSA1314_siteinfo_df <- NRSA1314_siteinfo %>% dplyr::select(UID, SITE_ID, VISIT_NO, DATE_COL, US_L3CODE, LAT_DD83, LON_DD83, RT_NRSA, AG_ECO9)
# METADATA: RT_NRSA	= Reference/So-so/Trashed designations for sites used in NRSA	

# Chemistry ----

# Read in chem data.
NRSA1314_chem <- read.csv('https://www.epa.gov/sites/default/files/2019-04/nrsa1314_widechem_04232019.csv')

# View variables
dplyr::glimpse(NRSA1314_chem)

# Select variables of interest
NRSA1314_chem_df <- NRSA1314_chem %>% dplyr::select(UID, COND_RESULT, POTASSIUM_RESULT, SODIUM_RESULT, SULFATE_RESULT, 
                                                    MAGNESIUM_RESULT, CHLORIDE_RESULT, CALCIUM_RESULT, PTL_RESULT, NTL_RESULT, 
                                                    TURB_RESULT, ANC_RESULT, PH_RESULT, SILICA_RESULT, TSS_RESULT)
# no Al in 13/14?
# METADATA: ALERT	= "R-Value below reporting limit for parameter, S-Sample shipping time exceeded two day limit, H-Holding time of sample exceeded parameter-specific limit"	
# METADATA: NTL	= total nitrogen
# METADATA: PTL	= total phosphorus


# Physical habitat ----

# Read in phab data
NRSA1314_phab <- read.csv('https://www.epa.gov/sites/default/files/2019-04/nrsa1314_phabmed_04232019.csv')

# View variables
dplyr::glimpse(NRSA1314_phab)

# Select variables of interest
NRSA1314_phab_df <- NRSA1314_phab %>% dplyr::select(UID, PCT_FN, LSUB_DMM, XSLOPE, XEMBED)
# METADATA: PCT_FN = "Bed Surface % Fines <0.06mm"

# Landscape metrics ----
# This file has basin slope variables and PRISM temperature variables.

# Read in landscape data
NRSA1314_land <- read.csv('https://www.epa.gov/sites/default/files/2019-04/nrsa1314_landmet_02132019.csv')

# View variables
dplyr::glimpse(NRSA1314_land)

# Select variables of interest
NRSA1314_land_df <- NRSA1314_land %>% dplyr::select(UID, RDDENS_WS, ELEV_WS_MAX)
#no NHDWAT_SLOPE in 13/14

# Merge data ----
# Merge 13/14 NRSA siteinfo, chemistry, phab, and landscape variables into one dataset. 
merge1 <- dplyr::left_join(NRSA1314_siteinfo, NRSA1314_chem_df, by = c("UID"))
merge2 <- dplyr::left_join(merge1, NRSA1314_phab_df, by = c("UID"))
merge3 <- dplyr::left_join(merge2, NRSA1314_land_df, by = c("UID"))

# COMIDs are included in the NRSA 13/14 data 

# Temperature data from StreamCat 
# Create a vector of the COMIDs in 0809, pass into the StreamCat function to get MAST, MSST, MWST for 08
COMID <- as.character(merge3$COMID)
COMID <- paste(COMID, collapse = ',')
sc_get_params(param='comid') #tells us which predictors I can get w/ comid, also ensures I'm matching how predictor is formatted
SCtemp13 <- sc_get_data(metric = 'mast_2013, msst_2013, mwst_2013', aoi='other', comid = COMID) 

# Merge nrsa_sf and SCtemp08 to get full NRSA 08 df
NRSA_13 <- dplyr::left_join(merge3, SCtemp13, by = c("COMID"))
rm(NRSA1314_chem, NRSA1314_chem_df, NRSA1314_land, NRSA1314_land_df, NRSA1314_phab, NRSA1314_phab_df, NRSA1314_siteinfo, 
   NRSA1314_siteinfo_df, SCtemp13, merge1, merge2, merge3, COMID, i)

#write.csv(NRSA_13, 'NRSA1314_compiled.csv')

### -------
### Can skip data compiling steps above if I load in the files I saved, to save time
#NRSA_08 <- read.csv('NRSA0809_compiled.csv')
#NRSA_13 <- read.csv('NRSA1314_compiled.csv')


#### ----------
# Subset NRSA data to WV data using StreamCat

# Downloading data from EPA's StreamCat data portal by state. While WV69 & WV70 ecoregions were used, they cut off at WV state line.

# Read in WV MSST data from StreamCat website
url <- 'https://gaftp.epa.gov/epadatacommons/ORD/NHDPlusLandscapeAttributes/StreamCat/States/RefStreamTempPred_WV.zip'
download.file(url, 'MSST.zip')
unzip('MSST.zip')
MSST_StreamCatWV <- read.csv('RefStreamTempPred_WV.csv')

# Reduce nrsa_sf nrsa dataframe to include only the COMIDs for the region of interest, by joining with MSST_StreamCatWV above.
WV_full_08 <- left_join(MSST_StreamCatWV, NRSA_08, by="COMID") #this is not working for some reason
WV_full_13 <- left_join(MSST_StreamCatWV, NRSA_13, by="COMID")
# Delete rows that have blank "UID" 
WV_full_08 <- WV_full_08[!is.na(WV_full_08$UID),]
WV_full_13 <- WV_full_13[!is.na(WV_full_13$UID),]
# Only 46 sites remain (08/09) and 34 (13/14) - is this right?

# Reduce to unique SITE_ID
WV_full_08 <- WV_full_08[!duplicated(WV_full_08$COMID),] #40 remain 08/09
WV_full_13 <- WV_full_13[!duplicated(WV_full_13$COMID),] #29 remain 13/14

#merge the 2 WV files, but need to reduce variables and rename variables to match first
#keep: COMID, MAST_2008, MAST_2013, 
WV_vars08 <- WV_full_08[,c(1,2,6,10,15,46:52, 54:61)]
WV_vars13 <- WV_full_13[,c(1,4,8,12,15,95:104, 106, 109:112)] #USE STREAMCAT MAST2013
#rename to match
WV_vars13 <- WV_vars13 %>% rename(COND=COND_RESULT,
                                  K = POTASSIUM_RESULT,
                                  'NA' = SODIUM_RESULT,
                                  SO4 = SULFATE_RESULT,
                                  MG = MAGNESIUM_RESULT,
                                  CL = CHLORIDE_RESULT,
                                  CA = CALCIUM_RESULT,
                                  PTL = PTL_RESULT,
                                  NTL = NTL_RESULT,
                                  TURB = TURB_RESULT,
                                  PH = PH_RESULT,
                                  MAST = MAST_2013.x,
                                  MSST = MSST_2013.x,
                                  MWST = MWST_2013.x)
WV_vars08 <- WV_vars08 %>% rename(MAST=MAST_2008.x,
                                  MSST = MSST_2008.x,
                                  MWST = MWST_2008.x,
                                  PH= PHLAB,
                                  'NA' = SODIUM)

#there are some repeat comids - make a year column and make a counter, then combine the two by just adding them together, then
#go through and pick out duplicates and divide duplicates evenly amongst 2008 and 2013, only keeping one row per duplicate COMID
WV_vars08$year <- rep(2008, 40)
WV_vars13$year <- rep(2013, 29)

WV_vars <- rbind(WV_vars08, WV_vars13)
# Go through and manually select 1 entry for each duplicate COMID
rownames(WV_vars) = NULL
# get rid of row 68 (2013), 39 (2008), 45 (2013), 7 (2008), 43 (2013), 25 (2008), 65 (2013), 30 (2008), 54 (2013), 14 (2008)

WV_vars <- WV_vars[-c(68,39,45,7,43,25,65,30,54,14),]

#Very few XEMBED values, going to remove that column
#WV_vars <- WV_vars[,-c(19:21)]

#write.csv(WV_vars, 'NRSA_WVsubset_yearscombined.csv')
rm(MSST_StreamCatWV, url)

### Can skip data compiling steps above by loading the results from previous runs below
#WV_vars <- read.csv('NRSA_WVsubset_yearscombined.csv')

#### ---
# NRSA SPLOMS --- (skip to line ~400 for more WV)

# use NRSA_13 for code for now

# Plot to see correlated variables in the data ---
NRSA_vars <- NRSA_13[,c(83:92, 94, 96:98, 100, 103:105)] #16 for COMID, 63 for reference designations

NRSA_vars <- NRSA_vars %>% rename(COND=COND_RESULT,
                                  K = POTASSIUM_RESULT,
                                  'NA' = SODIUM_RESULT,
                                  SO4 = SULFATE_RESULT,
                                  MG = MAGNESIUM_RESULT,
                                  CL = CHLORIDE_RESULT,
                                  CA = CALCIUM_RESULT,
                                  PTL = PTL_RESULT,
                                  NTL = NTL_RESULT,
                                  TURB = TURB_RESULT,
                                  PH = PH_RESULT,
                                  TSS = TSS_RESULT,
                                  MSST = MSST_2013,
                                  MAST = MAST_2013,
                                  MWST = MWST_2013)
# Change WV_vars to NA for sodium as well
WV_vars <- WV_vars %>% rename('NA' = SODIUM)

#delete lines w/ NAs
#NRSA_corr <- na.omit(NRSA_corr) #this deletes a lot - maybe reduce variables w/ lots of NAs instead

# Find skewed variables.
library(dlookr)
find_skewness(NRSA_vars, index=FALSE, value=TRUE, thres=0.5)
# Set above skewness threshold to 0.5.
# All skews are positively skewed except MSST. Log transform the positively skewed.
# Columns w/ 0's should be log+1 transformed. (PCT_FN) but this is extreme skew so 1/x+1 transform instead?

# Log transform skewed variables, then do Pearson's

#cols with negative values, try log1p() or cube root (negative skewness w/ negative values)
#cols with negative skewness, try square transformation or log10(max(x+1) - x) or 1/(max(x+1) - x) 
#severe skew: try inverse
NRSA_trans <- NRSA_vars %>%
  mutate(COND = log10(COND),
         K = log10(K),
         'NA' = log10(NRSA_vars$`NA`),
         SO4 = log10(SO4),
         MG = log10(MG),
         CL = log10(CL),
         CA = log10(CA),
         PTL = log1p(PTL),
         NTL = 1/(NTL+1), #inverse for severe skew, +1 because of 0's
         TURB = log10(TURB),
         TSS = log1p(TSS),
         PCT_FN = log1p(PCT_FN),
         #XEMBED = (XEMBED)^2, 
         MWST = log10(MWST))

#Check skewness post-transformation
find_skewness(NRSA_trans, index=FALSE, value=TRUE, thres=0.5) # A lot better

library(GGally) #for ggpairs() SPLOM function
# Use Spearman's for untransformed
#NRSA_vars %>% ggpairs(upper = list(continuous = wrap("cor", method = "spearman")))

# SPLOM w/ transformed data
ggpairs(NRSA_trans)


### Reference sites in NRSA surveys ONLY ----

#RT_NRSA = reference designation?
#Subset NRSA_full to just "R"
NRSA_ref <- NRSA_13[NRSA_13$RT_NRSA=="R",]
NRSA_ref <- NRSA_ref[,c(83:92, 94, 96:98, 100, 103:105)]
NRSA_ref <- NRSA_ref %>% rename(COND=COND_RESULT,
                                K = POTASSIUM_RESULT,
                                'NA' = SODIUM_RESULT,
                                SO4 = SULFATE_RESULT,
                                MG = MAGNESIUM_RESULT,
                                CL = CHLORIDE_RESULT,
                                CA = CALCIUM_RESULT,
                                PTL = PTL_RESULT,
                                NTL = NTL_RESULT,
                                TURB = TURB_RESULT,
                                PH = PH_RESULT,
                                TSS = TSS_RESULT)
#check for skewness
find_skewness(NRSA_ref, index=FALSE, value=TRUE, thres=0.5)

#transform
NRSA_transref <- NRSA_ref %>%
  mutate(
    COND = log10(COND),
    K = log10(K),
    'NA' = log10("NA"),
    SO4 = log10(SO4),
    MG = log10(MG),
    CL = log10(CL),
    CA = log10(CA),
    PTL = log1p(PTL),
    NTL = 1/(NTL+1), #inverse for severe skew, +1 because of 0's
    TURB = log10(TURB),
    TSS = log10(TSS),
    PCT_FN = log1p(PCT_FN),
    MAST_2013 = log10(MAST_2013),
    MWST_2013 = log10(MWST_2013))

#check skewness again
find_skewness(NRSA_transref, index=FALSE, value=TRUE, thres = 0.5)

#SPLOM w/ untransformed data
NRSA_ref %>% ggpairs(upper = list(continuous = wrap("cor", method = "spearman")))

#SPLOM w/ transformed data
ggpairs(NRSA_transref) #SPLOM for reference only, transformed data


### SPLOM for manuscript - all NRSA 13, transformed
options(scipen=0)

my_fn1 <- function(data, mapping, ...){
  p <- ggplot(data = data, mapping=mapping) + 
    stat_cor(aes(label=paste(..r.label.., ..p.label.., sep = "~")), 
             method="pearson", size=2.5, label.y.npc=0.1, label.x.npc = 0.1, hjust=0.5, digits=2)
  p
}  

### To get plots w/ r^2 and p values
p1 <- ggpairs(NRSA_trans, columnLabels = colnames(NRSA_trans),
              #upper=list(continuous = wrap("cor", method = "pearson"), size=2),
              upper=list(continuous = my_fn1),
              lower=list(continuous = "points"))+ 
  theme_bw() + theme(axis.text.x=(element_text(size=rel(0.75), angle=0)),
                     axis.text.y=(element_text(size=rel(0.75), angle=0)), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), panel.border = element_rect(fill = NA,colour = "grey35")) 
p1




### -----
# Work with the WV_full subset, after combining 0809 with 1314
#not transformed, Spearman's
WV_vars <- WV_vars[,-c(1, 5, 19, 21)] #get rid of COMID, SLOPE, year
#WV_vars %>% ggpairs(upper = list(continuous = wrap("cor", method = "spearman")))

# Check for skewness
find_skewness(WV_vars, index=FALSE, value=TRUE, thres=0.5)

WV_trans <- WV_vars %>% mutate(
  COND = log10(COND),
  K = log10(K),
  'NA' = log10(WV_vars$`NA`),
  SO4 = log10(SO4),
  MG = log10(MG),
  CL = log10(CL),
  CA = sqrt(CA),
  PTL = log10(PTL),
  NTL = log10(NTL),
  TURB = log10(TURB),
  PCT_FN = log1p(PCT_FN))

#reorder to be the same order as NRSA_trans SPLOM plot: COND, K, NA, SO4, MG, CL, CA, PTL, NTL, TURB, PH, TSS,
# PCT_FN, LSUB_DMM, XEMBED, MSST, MWST, MAST
NRSA_trans2 <- NRSA_trans[,c(1:11,13:18)]

WV_trans <- WV_trans[names(NRSA_trans2)]

#check skewness
find_skewness(WV_trans, index=FALSE, value=TRUE, thres=0.5)

#transformed, Pearson's
ggpairs(WV_trans)

### To get plots w/ r^2 and p values
p2 <- ggpairs(WV_trans, columnLabels = colnames(WV_trans),
              #upper=list(continuous = wrap("cor", method = "pearson"), size=2),
              upper=list(continuous = my_fn1),
              lower=list(continuous = "points"))+ 
  theme_bw() + theme(axis.text.x=(element_text(size=rel(0.75), angle=0)),
                     axis.text.y=(element_text(size=rel(0.75), angle=0)), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), panel.border = element_rect(fill = NA,colour = "grey35")) 
p2


### ---------
### For all 3 datasets (NRSA full, NRSA reference sites, and NRSA subset to WV), build a random forest model in which
### conductivity is the response and use % embeddedness, pH, % fines, LSUB_DMM, and MSST as predictors.

library(randomForest)

# with just response transformed
NRSA.rf.condt.df <- NRSA_vars[,c(1,11,13,14,15,16)]
#transform conductivity
NRSA.rf.condt.df$COND <- log10(NRSA.rf.condt.df$COND)
NRSA.rf.condt.df <- na.omit(NRSA.rf.condt.df)
NRSA.rf.condt <- randomForest(data=NRSA.rf.condt.df, COND ~ PCT_FN + PH + XEMBED + LSUB_DMM + MSST_2013, importance = TRUE)
print(NRSA.rf.condt)
varImpPlot(NRSA.rf.condt, sort = T, main="NRSA",n.var=5,type=1)
importance(NRSA.rf.condt)


#generate PartialDependencePlots 

#NRSA.rf.condt
a <- as.data.frame(partialPlot(NRSA.rf.condt, NRSA.rf.condt.df, PCT_FN,
                               lwd=2, cex.axis=1.25, cex.lab=1.25, las=1, cex.main=2, which.class=1, ylab="COND"))
b <- as.data.frame(partialPlot(NRSA.rf.condt, NRSA.rf.condt.df, PH,
                               lwd=2, cex.axis=1.25, cex.lab=1.25, las=1, cex.main=2, which.class=1, ylab="COND"))
c <- as.data.frame(partialPlot(NRSA.rf.condt, NRSA.rf.condt.df, XEMBED,
                               lwd=2, cex.axis=1.25, cex.lab=1.25, las=1, cex.main=2, which.class=1, ylab="COND"))
d <- as.data.frame(partialPlot(NRSA.rf.condt, NRSA.rf.condt.df, LSUB_DMM,
                               lwd=2, cex.axis=1.25, cex.lab=1.25, las=1, cex.main=2, which.class=1, ylab="COND"))
e <- as.data.frame(partialPlot(NRSA.rf.condt, NRSA.rf.condt.df, MSST_2013,
                               lwd=2, cex.axis=1.25, cex.lab=1.25, las=1, cex.main=2, which.class=1, ylab="COND"))

#that doesn't allow me to put them all into one figure, patchwork doesn't work here.
#it does put all x and y values into a list
z <- as.data.frame(quantile(NRSA.rf.condt.df$MSST_2013, probs = seq(.1, .9, by = .1)))
z <- z %>% rename(quants = 1)

l <- ggplot(e, aes(x=x, y=y)) +
  geom_line(size=1) +
  geom_rug(data=z, aes(x=quants), sides="b", inherit.aes=FALSE) +
  theme_classic() +
  labs(x="Mean summer stream temperature (°C)") +
  theme(axis.title = element_text(size = 14),
        axis.title.y = element_blank(),
        plot.title = element_text(size= 14),
        #axis.text.y=element_blank(),
        axis.text = element_text(size = 12, color='black')) +
  ylim(c(1.9,2.8)) 

#do for all. store, put in patchwork plot
library(patchwork)
p4 <- ggplot(data.frame(l = "log10(Conductivity)", x = 1, y = 1)) +
  geom_text(aes(x, y, label = l), angle = 90, size=5.5) + 
  theme_void() +
  coord_cartesian(clip = "off")

p4+((h+i)/(j+k)/(l+plot_spacer())) + plot_layout(widths=c(1,25))

#would be best to rewrite lines ~452-480 to be looped through, instead did manually for time's sake

### RF for WV subset
WV.rf.df <- WV_vars[,c(2,4,13,14,15,16)]
WV.rf.df <- na.omit(WV.rf.df)
#transform just the response variable
WV.rf.df$COND <- log10(WV.rf.df$COND)
WV.rf <- randomForest(data=WV.rf.df, COND ~ PCT_FN + PH + LSUB_DMM + MSST + TURB, importance = TRUE)
print(WV.rf)
importance(WV.rf)

varImpPlot(WV.rf, sort = T, main="WV subset of NRSA",n.var=5,type=1)

a <- as.data.frame(partialPlot(WV.rf, WV.rf.df, PCT_FN,
                               lwd=2, cex.axis=1.25, cex.lab=1.25, las=1, cex.main=2, which.class=1, ylab="COND"))
b <- as.data.frame(partialPlot(WV.rf, WV.rf.df, PH,
                               lwd=2, cex.axis=1.25, cex.lab=1.25, las=1, cex.main=2, which.class=1, ylab="COND"))
c <- as.data.frame(partialPlot(WV.rf, WV.rf.df, TURB,
                               lwd=2, cex.axis=1.25, cex.lab=1.25, las=1, cex.main=2, which.class=1, ylab="COND"))
d <- as.data.frame(partialPlot(WV.rf, WV.rf.df, LSUB_DMM,
                               lwd=2, cex.axis=1.25, cex.lab=1.25, las=1, cex.main=2, which.class=1, ylab="COND"))
e <- as.data.frame(partialPlot(WV.rf, WV.rf.df, MSST,
                               lwd=2, cex.axis=1.25, cex.lab=1.25, las=1, cex.main=2, which.class=1, ylab="COND"))

z <- as.data.frame(quantile(WV.rf.df$MSST, probs = seq(.1, .9, by = .1)))
z <- z %>% rename(quants = 1)

q <- ggplot(e, aes(x=x, y=y)) +
  geom_line(size=1) +
  geom_rug(data=z, aes(x=quants), sides="b", inherit.aes=FALSE) +
  theme_classic() +
  labs(x="Mean summer stream temperature (°C)") +
  theme(axis.title = element_text(size = 14),
        axis.title.y = element_blank(),
        plot.title = element_text(size=14),
        #axis.text.y=element_blank(),
        axis.text = element_text(size = 12, color='black')) +
  ylim(c(2.1, 2.5))

#do for all. store, put in patchwork plot

p4+((m+n)/(o+p)/(q+plot_spacer())) + plot_layout(widths=c(1,25))
#saved manually

### -----
### NMDS on NRSA data
### Use NRSA_vars for NMDS, add in macros data
library(vegan)
library(goeveg) #for dimcheck

#load in EPA NRSA 2013/2014 abundance data
#NRSA_bent <- read.csv('https://www.epa.gov/sites/default/files/2019-04/nrsa1314_bentcnts_04232019.csv') ###This code not working anymore
NRSA_bent <- read.csv('../../data/nrsa1314/nrsa1314_BentCnts_04232019.csv')

#delete unnecessary columns
NRSA_bent <- NRSA_bent[,c(2,7,9)]
NRSA_bent <- NRSA_bent %>%
  spread(TAXA_ID, TOTAL)  #spread to get into wide format for use in rankabundance()
NRSA_bent[is.na(NRSA_bent)] = 0 #replace NA values with 0
#dimcheckMDS(na.omit(NRSA_bent), distance = "bray")

# Get rid of LSUB for now because of negative numbers, bray's dissimilarity cannot handle
NRSA_env <- NRSA_13[,c(2, 83:92, 94, 96, 97, 100, 103:105)] #Add in col2 for UID, remove LSUB col
#subset NRSA_env to the UIDs that are represented in the NRSA_bent
#also rename columns for plotting purposes
NRSA_env <- NRSA_env %>% filter(UID %in% NRSA_bent$UID) %>%
  rename(COND=COND_RESULT,
         K = POTASSIUM_RESULT,
         'NA' = SODIUM_RESULT,
         SO4 = SULFATE_RESULT,
         MG = MAGNESIUM_RESULT,
         CL = CHLORIDE_RESULT,
         CA = CALCIUM_RESULT,
         PTL = PTL_RESULT,
         NTL = NTL_RESULT,
         TURB = TURB_RESULT,
         PH = PH_RESULT,
         TSS = TSS_RESULT)
#remove UID?
NRSA_env <- NRSA_env[,-1]
NRSA_bent <- NRSA_bent[,-1]

NRSA_nmds <- metaMDS(na.omit(NRSA_bent), distance = "bray", noshare = TRUE, k=3) #should I be using autotransform?
plot(NRSA_nmds)

ordiplot(NRSA_nmds, display = "sites")
en <- envfit(NRSA_nmds, NRSA_env, add=TRUE, na.rm=TRUE)
plot(en)

### use transformed data
NRSA_env <- NRSA_13[,c(2, 83:92, 94, 96, 97, 100, 103:105)]

NRSA_bent <- read.csv('../../data/nrsa1314/nrsa1314_BentCnts_04232019.csv')

#delete unnecessary columns
NRSA_bent <- NRSA_bent[,c(2,7,9)]
NRSA_bent <- NRSA_bent %>%
  spread(TAXA_ID, TOTAL)  #spread to get into wide format for 

NRSA_trans_env <- NRSA_env %>% filter(UID %in% NRSA_bent$UID) %>%
  rename(COND=COND_RESULT,
         K = POTASSIUM_RESULT,
         'NA' = SODIUM_RESULT,
         SO4 = SULFATE_RESULT,
         MG = MAGNESIUM_RESULT,
         CL = CHLORIDE_RESULT,
         CA = CALCIUM_RESULT,
         PTL = PTL_RESULT,
         NTL = NTL_RESULT,
         TURB = TURB_RESULT,
         PH = PH_RESULT,
         TSS = TSS_RESULT,
         MSST = MSST_2013,
         MAST = MAST_2013,
         MWST = MWST_2013) %>%
  mutate(COND = log10(COND),
         K = log10(K),
         'NA' = log10(`NA`),
         SO4 = log10(SO4),
         MG = log10(MG),
         CL = log10(CL),
         CA = log10(CA),
         PTL = log1p(PTL),
         NTL = 1/(NTL+1), #inverse for severe skew, +1 because of 0's
         TURB = log10(TURB),
         TSS = log1p(TSS),
         PCT_FN = log1p(PCT_FN),
         #XEMBED = (XEMBED)^2, 
         MWST = log10(MWST))

NRSA_trans_env <- NRSA_trans_env[,-1]
NRSA_bent <- NRSA_bent[,-1]

entrans <- envfit(NRSA_nmds, NRSA_trans_env, add=TRUE, na.rm=TRUE)
plot(entrans)

#make ggplots instead of baseR plots
#extract NMDS scores (x and y coordinates) for sites from newer versions of vegan package
data.scores.nrsa = as.data.frame(scores(NRSA_nmds)$sites)

#extract envfit vectors
entrans_coord = as.data.frame(scores(entrans, "vectors")) * ordiArrowMul(entrans)

library(ggrepel)
#ggplot
a <-ggplot(data.scores.nrsa, aes(x=NMDS1, y=NMDS2)) +
  geom_point( size=2, alpha=0.75) +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = entrans_coord, size = 1, alpha = 0.75, colour = 'black') +
  geom_text_repel(data = entrans_coord, aes(x = NMDS1, y = NMDS2), colour = "black", 
            label = row.names(entrans_coord), size=5.5, force=0.1, min.segment.length = 0.25) + 
  ggtitle('a') +
  theme(axis.title = element_text(size = 14, colour = "black"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "black"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
        legend.text = element_text(size = 12, colour = "black"),
        plot.title = element_text(size= 20, vjust = -9, hjust = 0.025)) + xlim(-3,2)


### ----------
### NMDS using benthic data, but just for WV subset, but use both years of WV data
# How to add benthic data to WV_vars? 
# maybe start w/ the 'WV_full_08' and 'WV_full_13', match benthic data to COMID, and then join and pick 1 COMID per year?
# No COMIDs in the benthic file, going to have to join by UID, add UID to the full files way above

NRSA_bent08 <- read.csv('https://www.epa.gov/sites/default/files/2016-11/nrsa0809bentcts.csv')

#reduce variables down to taxon name, uid, total
NRSA_bent08 <- NRSA_bent08[,c(2,9,11)]
#spread by site (UID)
NRSA_bent08 <- NRSA_bent08 %>%
  pivot_wider(names_from=TAXA_ID, values_from = TOTAL, values_fill = 0, values_fn = sum) #this might not have worked correctly


NRSA_bent13 <- read.csv('https://www.epa.gov/sites/default/files/2019-04/nrsa1314_bentcnts_04232019.csv')
#reduce variables down to taxon name, uid, total
NRSA_bent13 <- NRSA_bent13[,c(2,7,9)]
#spread by site (UID)
NRSA_bent13 <- NRSA_bent13 %>%
  pivot_wider(names_from=TAXA_ID, values_from = TOTAL, values_fill = 0, values_fn = sum) #this might not have worked correctly

# Add the benthic data to each WV_vars file
WV_all08 <- left_join(WV_vars08, NRSA_bent08, by="UID") #this does not work because UIDs are different?
WV_all13 <- left_join(WV_vars13, NRSA_bent13, by='UID') #this one works

WV_all08 <- WV_all08 %>% rename(`NA` = `NA.x`)

#join files back together, like before
library(gtools)
WV_all <- smartbind(WV_all08, WV_all13)
# Go through and manually select 1 entry for each duplicate COMID
rownames(WV_all) = NULL
# get rid of row 68 (2013), 39 (2008), 45 (2013), 7 (2008), 43 (2013), 25 (2008), 65 (2013), 30 (2008), 54 (2013), 14 (2008)

WV_all <- WV_all[-c(68,39,45,7,43,25,65,30,54,14),]

#Very few XEMBED values, going to remove that column
WV_all <- WV_all[,-c(5,19:21)]
#remove rows 9, 26, 27 because no benthic data
WV_all <- WV_all[-c(9,26,27),]
WV_bent <- WV_all[,c(18:ncol(WV_all))] 

#make all NAs 0, remove colSums that = 0
WV_bent[is.na(WV_bent)] = 0
WV_bent <- WV_bent[,colSums(WV_bent) != 0] #removes species that don't have presencse in the WV subset
WV_bent <- WV_bent[rowSums(WV_bent) != 0, ] #okay, there were still rows that had no bugs? will need to rejoin to subset WV_all

#filter WV_env to pH >=6
WV_env <- WV_all[,c(2:17)]
WV_env <- WV_env[WV_env$PH>=6,]
rows <- rownames(WV_env)
#filter WV_bent to be the sites in WV_env
WV_bent <- WV_bent[rownames(WV_bent) %in% rows,]

#transform
WV_env <- WV_env %>% mutate(
  COND = log10(COND),
  K = log10(K),
  `NA` = log10(`NA`),
  SO4 = log10(SO4),
  MG = log10(MG),
  CL = log10(CL),
  CA = sqrt(CA),
  PTL = log10(PTL),
  NTL = log10(NTL),
  TURB = log10(TURB),
  PCT_FN = log1p(PCT_FN))

rows <- rownames(WV_bent)
WV_env <- WV_env[rownames(WV_env) %in% rows,]

#dimcheckMDS(na.omit(WV_bent), distance = "bray")
WV_nmds <- metaMDS(na.omit(WV_bent), distance = "bray", noshare = TRUE, k=3)
ordiplot(WV_nmds, display = "sites")

en <- envfit(WV_nmds, WV_env, add=TRUE, na.rm=TRUE)
plot(en)

#make ggplots instead of baseR plots
#extract NMDS scores (x and y coordinates) for sites from newer versions of vegan package
data.scores = as.data.frame(scores(WV_nmds)$sites)

#extract envfit vectors
en_coord = as.data.frame(scores(en, "vectors")) * ordiArrowMul(en)

#ggplot
b <- ggplot(data.scores, aes(x=NMDS1, y=NMDS2)) +
  geom_point(size=2, alpha=0.75) +
  #theme_bw() +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = en_coord, size = 1, alpha = 0.75, colour = "black") +
  geom_text(data = en_coord, aes(x = NMDS1, y = NMDS2), 
            label = row.names(en_coord), size=5.5, nudge_x = 0.07) + 
  ggtitle("b") +
  theme(axis.title = element_text(size = 14), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 12), 
        legend.text = element_text(size = 12),
        plot.title = element_text(size= 20, vjust = -9, hjust = 0.025))

library(patchwork)
a/b

#### Generate what # of NRSA stream sites would be above each HC05 value from the simulation
# First, load in hc05_vals.csv
hc05_vals <- read.csv("hc05_results.csv")
#NRSA_13 <- read.csv('NRSA1314_compiled.csv')
#WGT_EXT_SP variable = "Length of sampled streams represented by site in km"
#column 80 in NRSA_13, also keep column 64 for reference designation, col 83 = COND_RESULT
NRSA_lengthdf <- NRSA_13[,c(3,64,80,83)]
#filter to just ref
#NRSA_lengthdf <- NRSA_lengthdf[NRSA_lengthdf$RT_NRSA=="R",]

#create a loop that creates a df w/ the following columns:
#subsample size, hc05 of that subsample, # streams above that hc05, stream length above that hc05
lengthcond <- data.frame(matrix(vector(), 0, ncol = 6,
                                dimnames = list(c(), c("subsample_size", "meanbsHC05", "no_above", "perc_above", "length_above_km", "perc_length_above"))))

#cols 3 and 4 can't work with NAs

for(i in 1:nrow(hc05_vals)){
  lengthcond[i,1] <- hc05_vals[i,2]
  lengthcond[i,2] <- hc05_vals[i,3]
  lengthcond[i,3] <- sum(NRSA_lengthdf$COND_RESULT > hc05_vals[i,3])
  lengthcond[i,4] <- (lengthcond[i,3]/nrow(NRSA_lengthdf))*100
  #there are NAs in the length column, need to na.omit for this column only
  lengthcond[i,5] <- sum(na.omit(NRSA_lengthdf$WGT_EXT_SP[which(NRSA_lengthdf$COND_RESULT>=hc05_vals[i,3])]))
  lengthcond[i,6] <- (lengthcond[i,5]/sum(na.omit(NRSA_lengthdf$WGT_EXT_SP)))*100
}

write.csv(lengthcond, file = "NRSAref_HC05_length.csv")

sum(is.na(NRSA_lengthdf$WGT_EXT_SP))

#for some reason, cols 3 and 4 not working in lengthcond, delete them
lengthcond <- lengthcond[,c(1,2,5,6)]


### Suter and Cormier (2013) filter so that all sites w/ pH <6 are excluded. Do that and rerun processes above.

#filter out pH < 6
NRSA_trans_pH <- NRSA_trans[NRSA_trans$PH>=6,]

p3 <- ggpairs(NRSA_trans_pH, columnLabels = colnames(NRSA_trans_pH),
              #upper=list(continuous = wrap("cor", method = "pearson"), size=2),
              upper=list(continuous = my_fn1), #my_fn1 is in script above
              lower=list(continuous = "points"))+ 
  theme_bw() + theme(axis.text.x=(element_text(size=rel(0.75), angle=0)),
                     axis.text.y=(element_text(size=rel(0.75), angle=0)), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), panel.border = element_rect(fill = NA,colour = "grey35")) 
p3

#### --------------
#### use pH-subset data in RF and NMDS

# just response transformed, used transformed data in pH filter above, need to refilter non-transformed data and then transform conductivity
NRSAph.rf.condt.df <- NRSA_vars[NRSA_vars$PH >= 6, c(1,11,13,14,15,16)]

#transform conductivity
NRSAph.rf.condt.df$COND <- log10(NRSAph.rf.condt.df$COND)
NRSAph.rf.condt.df <- na.omit(NRSAph.rf.condt.df)
NRSAph.rf.condt <- randomForest(data=NRSAph.rf.condt.df, COND ~ PCT_FN + PH + XEMBED + LSUB_DMM + MSST_2013, importance = TRUE)
print(NRSAph.rf.condt)
varImpPlot(NRSAph.rf.condt, sort = T, main="NRSA",n.var=5,type=1)
importance(NRSAph.rf.condt)


#generate PartialDependencePlots 

#NRSA.rf.condt
a <- as.data.frame(partialPlot(NRSAph.rf.condt, NRSAph.rf.condt.df, PCT_FN,
                               lwd=2, cex.axis=1.25, cex.lab=1.25, las=1, cex.main=2, which.class=1, ylab="COND"))
b <- as.data.frame(partialPlot(NRSAph.rf.condt, NRSAph.rf.condt.df, PH,
                               lwd=2, cex.axis=1.25, cex.lab=1.25, las=1, cex.main=2, which.class=1, ylab="COND"))
c <- as.data.frame(partialPlot(NRSAph.rf.condt, NRSAph.rf.condt.df, XEMBED,
                               lwd=2, cex.axis=1.25, cex.lab=1.25, las=1, cex.main=2, which.class=1, ylab="COND"))
d <- as.data.frame(partialPlot(NRSAph.rf.condt, NRSAph.rf.condt.df, LSUB_DMM,
                               lwd=2, cex.axis=1.25, cex.lab=1.25, las=1, cex.main=2, which.class=1, ylab="COND"))
e <- as.data.frame(partialPlot(NRSAph.rf.condt, NRSAph.rf.condt.df, MSST_2013,
                               lwd=2, cex.axis=1.25, cex.lab=1.25, las=1, cex.main=2, which.class=1, ylab="COND"))

#that doesn't allow me to put them all into one figure, patchwork doesn't work here.
#it does put all x and y values into a list
z <- as.data.frame(quantile(NRSAph.rf.condt.df$MSST_2013, probs = seq(.1, .9, by = .1))) #for the rug
z <- z %>% rename(quants = 1)

# run through this code 5 times to get each PDP in ggplot format for publication
l <- ggplot(e, aes(x=x, y=y)) +
  geom_line(size=1) +
  geom_rug(data=z, aes(x=quants), sides="b", inherit.aes=FALSE) +
  theme_classic() +
  labs(x="Mean summer stream temperature (°C)") +
  theme(axis.title = element_text(size = 14),
        axis.title.y = element_blank(),
        plot.title = element_text(size= 14),
        #axis.text.y=element_blank(),
        axis.text = element_text(size = 12, color='black')) +
  ylim(c(1.9,2.85)) 

#do for all. store, put in patchwork plot
library(patchwork)
p4 <- ggplot(data.frame(l = "log10(Conductivity)", x = 1, y = 1)) +
  geom_text(aes(x, y, label = l), angle = 90, size=5.5) + 
  theme_void() +
  coord_cartesian(clip = "off")

p4+((h+i)/(j+k)/(l+plot_spacer())) + plot_layout(widths=c(1,25))


### WV filtered to ph >= 6

WV_trans_pH <- WV_trans[WV_trans$PH>=6,] # only 1 site filtered out

p4 <- ggpairs(WV_trans_pH, columnLabels = colnames(WV_trans_pH),
              #upper=list(continuous = wrap("cor", method = "pearson"), size=2),
              upper=list(continuous = my_fn1), #my_fn1 is in script above
              lower=list(continuous = "points"))+ 
  theme_bw() + theme(axis.text.x=(element_text(size=rel(0.75), angle=0)),
                     axis.text.y=(element_text(size=rel(0.75), angle=0)), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), panel.border = element_rect(fill = NA,colour = "grey35")) 
p4

### RF for WV subset
WVph.rf.df <- WV_vars[,c(2,4,13,14,15,16)]
WVph.rf.df <- WVph.rf.df[WVph.rf.df$PH >= 6, ]
WVph.rf.df <- na.omit(WVph.rf.df)
#transform just the response variable
WVph.rf.df$COND <- log10(WVph.rf.df$COND)
WVph.rf <- randomForest(data=WVph.rf.df, COND ~ PCT_FN + PH + LSUB_DMM + MSST + TURB, importance = TRUE)
print(WVph.rf)
importance(WVph.rf)

varImpPlot(WVph.rf, sort = T, main="WV subset of NRSA",n.var=5,type=1)

a <- as.data.frame(partialPlot(WVph.rf, WVph.rf.df, PCT_FN,
                               lwd=2, cex.axis=1.25, cex.lab=1.25, las=1, cex.main=2, which.class=1, ylab="COND"))
b <- as.data.frame(partialPlot(WVph.rf, WVph.rf.df, PH,
                               lwd=2, cex.axis=1.25, cex.lab=1.25, las=1, cex.main=2, which.class=1, ylab="COND"))
c <- as.data.frame(partialPlot(WVph.rf, WVph.rf.df, TURB,
                               lwd=2, cex.axis=1.25, cex.lab=1.25, las=1, cex.main=2, which.class=1, ylab="COND"))
d <- as.data.frame(partialPlot(WVph.rf, WVph.rf.df, LSUB_DMM,
                               lwd=2, cex.axis=1.25, cex.lab=1.25, las=1, cex.main=2, which.class=1, ylab="COND"))
e <- as.data.frame(partialPlot(WVph.rf, WVph.rf.df, MSST,
                               lwd=2, cex.axis=1.25, cex.lab=1.25, las=1, cex.main=2, which.class=1, ylab="COND"))

z <- as.data.frame(quantile(WVph.rf.df$LSUB_DMM, probs = seq(.1, .9, by = .1)))
z <- z %>% rename(quants = 1)

p <- ggplot(d, aes(x=x, y=y)) +
  geom_line(size=1) +
  geom_rug(data=z, aes(x=quants), sides="b", inherit.aes=FALSE) +
  theme_classic() +
  labs(x="Log geometric mean particle diameter (mm)") +
  theme(axis.title = element_text(size = 14),
        axis.title.y = element_blank(),
        plot.title = element_text(size= 14),
        axis.text.y=element_blank(),
        axis.text = element_text(size = 12, color='black')) +
  ylim(c(2.0, 2.6))

#do for all. store, put in patchwork plot

p4+((m+n)/(o+p)/(q+plot_spacer())) + plot_layout(widths=c(1,25))



### ------------------------
### Make NRSA species richness figure to match the one from simulation_SSD.R 
nrsa <- read.csv("../../data/nrsa1314_bentcnts.csv")
#delete unnecessary columns
nrsa <- nrsa[,c(2,7,9)]
nrsa_chem <- read.csv("../../data/nrsa1314_widechem.csv")
#keep only UID and conductivity data
nrsa_chem <- nrsa_chem[, c(2, 119)]
nrsa <- left_join(nrsa_chem, nrsa, by = "UID")

#to get presence absence for species richness, we want wide format
nrsa <- nrsa %>%
  rename(site = UID, conductivity = COND_RESULT, species = TAXA_ID, abunds = TOTAL) %>%
  pivot_wider(names_from = species, values_from = abunds)
nrsa[is.na(nrsa)] = 0

nrsa_pa <- nrsa[, 3:ncol(nrsa)] 
cond_sites <- nrsa[, 2]

nrsa_pa <- nrsa_pa[, colSums(nrsa_pa == 0) <= (nrow(nrsa_pa)-25)] #only species w/ >=25 presences, 2261 total sites
nrsa_pa[nrsa_pa>=1] = 1 #getting rid of abundances, changing to pres/abs
total <- rowSums(nrsa_pa) 
nrsa_pa <- cbind(nrsa_pa, total)
PA <- cbind(total, cond_sites)
PA <- as.data.frame(PA)

p5 <- ggplot(PA, aes(x = conductivity, y = total)) +
  geom_jitter() +
  geom_quantile(quantiles = 0.95, size = 1.5, col = "grey30") + 
  xlab("Conductivity (µS/cm)") +
  ylab("") +
  xlim(0,2000) + 
  ylim(0,110) +
  theme_classic()+
  ggtitle("b") +
  theme(legend.position = "none",
        axis.title = element_text(size = 14),
        plot.title = element_text(size= 20, vjust = -9, hjust = 0.025),
        axis.text = element_text(size = 12, color='black'))

p4 <- ggplot(data.frame(l = "Species richness", x = 1, y = 1)) +
  geom_text(aes(x, y, label = l), angle = 90, size=5.5) + 
  theme_void() +
  coord_cartesian(clip = "off")

#when spec_rich_plot loaded from simulation_SSD R file, can use below code to make manuscript plot
#library(patchwork)
#richp <- spec_rich_plot/plot_spacer()/p5 + plot_layout(heights = c(5.2, -1.01 ,5.2))
richp <- spec_rich_plot/p5
species_rich <- p4+richp+plot_layout(widths=c(0.5,25))
#ggsave("../../../manuscript/plots/species_rich", plot = species_rich, width = 7, height = 8, units = "in", device = "tiff")
# manually save, that looks bad