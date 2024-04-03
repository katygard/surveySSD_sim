# File: simulation_SSD.R
# Purpose: Simulate community, subsample from community, create SSDs and bootstrap, compare results. Cleaned up version for publication.
# Date: cleaned 2024, created 2021-2023

options(scipen = 999999)
library(coenoflex)
library(xlsx)
library(plyr)
library(tidyverse) #need dplyr, ggplot2, stringr (for str_sort()), tidyr
library(ggforce) #for facet_zoom()
library(patchwork) 
#library(BiodiversityR) #crashes R as of March '24, used only in rankabundance figures (which we do not use for manuscript) and in Table 2 calculations

# Run rarify script (credit: John Van Sickle, USEPA, 2005) before running code below, needed to subsample.

#environment data saved after all steps = post_bootstrap.R.RData

### generate community using coenoflex -----
set.seed(102)
sim_comm <- coenoflex(numgrd = 1, numplt = 1000, numspc = 400, #numspc needs to be higher than actual target
                      grdtyp = 'e', #one environmental gradient over 1000 sites
                      grdlen = 2000, width = 1000, variab = 75, grdprd = -100, #gradient from 0-2000, species amplitude width = 1,000
                      #with 75% variability (250-1750 width)
                      #grdprd = change in productivity from low to high x-axis value
                      alphad = 1.5, pdist = 'g', sdist = 'r', skew = 3.5, 
                      #skew >1 to make abundances more highly concentrated at low conductivity values
                      aacorr = 0, cmpasy = 1, maxtot = 5000000, noise = 10, 
                      slack = 0.2, #slack = % of time they're absent from suitable sites
                      autlin = 'ave(1)') #'doesn't matter' for 1 gradient

sites <- sim_comm$site
abundances <- floor(sim_comm$veg)
spec_rich <- apply(floor(sim_comm$veg)>0,1,sum) #get data into pres/abs format, sum by row to get richness at each site
pres_abs <- cbind(spec_rich, sites)
rm(spec_rich)
pres_abs <- as.data.frame(pres_abs)
pres_abs <- pres_abs %>% rename(sites = V2)

# plot to see species richness pattern
spec_rich_plot <- ggplot(pres_abs, aes(x = sites, y = spec_rich)) +
  geom_point() +
  xlab("") +
  ylab("") +
  theme_classic() +
  ylim(0,250) +
  ggtitle("a") +
  theme(legend.position = "none",
        axis.title = element_text(size = 14),
        plot.title = element_text(size= 20, vjust = -9, hjust = 0.025),
        axis.text = element_text(size = 12, color='black'),
        axis.text.x = element_blank())

#----
# need the dataframe in 3 column format (site, species, abundance) w/ only non-zero datapoints
sites <- floor(sites)
abundances <- cbind(abundances, sites)
abundances <- as.data.frame(abundances)
#write.csv(abundances, file = "../../data/simulated_total_abundances.csv") #export before gather
abundances <- abundances %>% 
  rename(site_cond = rev(names(abundances))[1]) %>% #added site column renamed to site_cond to distinguish it from species
  pivot_longer(names_to="species", values_to="abunds", cols=!site_cond) #gather to put data in rarify() format
#keep only observations with non-zero values
abundances<-abundances[abundances$abunds>0,]

samp_sizes <- c(50, 75, 100, 150, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 2000, 5000, 10000)

#create empty dataframes for subsampling
for(i in samp_sizes){
  assign(paste("outbug", i, sep = "_"), 
         data.frame(matrix(vector(), 0, 0))) #creates blank df of name outbug_i for each samp_size
}

#put empty dataframes into a list
samp_list <- mget(ls(pattern = "outbug_"))
rm(list=ls(pattern="outbug_")) #useful code to remove many items from environment at once, matching pattern
nameorder <- c(names(samp_list))
samp_list <- samp_list[str_sort(nameorder, numeric = TRUE)]

#created a list of empty dataframes, now need to populate them using rarify to subsample (loop takes >5 hrs on laptop, ~3 on desktop)
### Must run rarify script before this step, if rarify isn't already loaded into environment
set.seed(123)
for (i in 1:length(samp_sizes)){
  samp_list[[i]] <- rarify(inbug = as.data.frame(abundances), sample.ID = 'site_cond', abund = 'abunds', subsiz = samp_sizes[i])
}

#next step: sort, remove absences, put into wide format
for (i in 1:length(samp_sizes)) {
  samp_list[[i]] <- samp_list[[i]][
    order(samp_list[[i]][,1]),] #orders list by site conductivity
  samp_list[[i]] <- samp_list[[i]][samp_list[[i]][3] > 0, ] #removes all rows with 0 abundance
  samp_list[[i]] <- samp_list[[i]] %>%
    spread(species, abunds)#spreads dataframe back out so that each species is its own column (will have NAs instead of 0s now)
  samp_list[[i]] <- samp_list[[i]][, colSums(is.na(samp_list[[i]])) <= 975] 
  # ^ only keeps species w/ at least 25 presences (deletes other columns)
}

# after subsampling, need to determine conductivity upper limit/XC95 value for each species
# Use the ewcdf function in the spatstat library.
library(spatstat)

#loop to create blank dataframes for the xc95 vals to be stored in
for(i in samp_sizes){
  assign(paste("xc_vals", i, sep = "_"), 
         data.frame(matrix(vector(), 0, 2,
                           dimnames=list(c(), c("Species_ID", "XC95"))),
                    stringsAsFactors=F)) #creates blank df of name xc_vals_i
}

#put blank dataframes into a list
xc95_list <- mget(ls(pattern = "xc_vals_"))
rm(list=ls(pattern="xc_vals_")) 
nameorder <- c(names(xc95_list))
xc95_list <- xc95_list[str_sort(nameorder, numeric = TRUE)]

# loop to populate dataframes with XC95 values, calculated from ewcdf() 
for(i in 1:length(samp_sizes)){
  for(j in 2:ncol(samp_list[[i]])){
    pres_only <- na.omit(samp_list[[i]][,c(1,j)])
    cdfunc <- ewcdf(pres_only[,1], pres_only[,2])
    pres_only$ewcdf_vals <- cdfunc(pres_only[,1])
    sort_xy <- sortedXyData(expression(site_cond), expression(ewcdf_vals), pres_only)
    (XC95 <- NLSstClosestX(sort_xy, yval = 0.95))
    #store above into new line on blank dataframe
    xc95_list[[i]][j-1,] <- c(names(pres_only[2]), XC95)
  }
}
rm(pres_only, sort_xy, i, j, XC95)

### SSD - do manually instead of using SSDtools. ssdtools distfits were bad
# manual calculation of frequency distribution for SSD curve
# generalized additive model with 3 degrees of freedom) to model the 
# likelihood of a taxon being observed with increasing conductivity

## We will want to store the HC05 values for each list name
# create empty dataframe, one column for names of df in main list, one column for corresponding HC05 value
hc05_vals <- data.frame(matrix(vector(), 16, 2,
                               dimnames=list(c(), c("subsample_size", "HC05_val"))),
                        stringsAsFactors=F) #creates blank df
hc05_vals$subsample_size <- samp_sizes #populate sample size values

#create storage list for SSD plots for each subsample size iteration
plot_list = list()

#loop to manually calculate SSD and corresponding HC05 for each subsample df, and generate and store SSD ggplot
for (i in 1:length(samp_sizes)) {
  xc95_list[[i]] <- xc95_list[[i]] %>%
   rename(Conc = XC95)
  xc95_list[[i]][,2] <- as.numeric(xc95_list[[i]][,2])
  xc95_list[[i]] <- xc95_list[[i]][
    order(xc95_list[[i]][,2]),]
  # cumulative proportion for each genus P is calculated as P=R/(N+1), 
  # where R is the rank of the genus and N is the number of genera
  xc95_list[[i]]$rank <- rank(xc95_list[[i]]$Conc, ties.method = "first")
  xc95_list[[i]]$p <- (xc95_list[[i]][,3])/(nrow(xc95_list[[i]])+1) 
  # make a sortedXyData object, to be used in inverse interpolation to get HC05 value ----
  xysort <- sortedXyData(expression(Conc), expression(p), xc95_list[[i]])
  ### inverse interpolation 
  hc05_vals[i,2] <- (HC05 <- NLSstClosestX(xysort, yval = 0.05))  
  #plot (stored)
  z <- ggplot(xc95_list[[i]], aes(x = Conc, y = p)) +
    geom_point() +
    geom_hline(yintercept = 0.05, linetype = "dashed") +
    stat_ecdf(pad=FALSE) +
    labs(y = "Proportion of species", x = "Conductivity µS/cm") +
    ggtitle(paste0(samp_sizes[i])) +
    #geom_vline(xintercept = HC05) + 
    theme_bw()
  plot_list[[i]] <- z
}
rm(HC05, cdfunc, i, xysort, z)

library(gridExtra)
#export plots - pdf does all 16, choose 1 and do .tiff for single plots
ggsave(
  filename = "plots.pdf", 
  plot = marrangeGrob(plot_list, nrow=2, ncol=2),
  width = 8, height = 11, units = 'in'
)

# Export HC05 results
#write.xlsx(hc05_vals, "../../data/hc05results.xlsx")

# Calculate uncertainty (upper and lower 95% CIs) for all XC95 values, and then their corresponding HC05 ----
# First: bootstrapped XC95 for each species in each of the 16 rarify subsets of data (outbug_#)
# goes back to samp_list data

# Empty dataframes to store resampled values, can reuse for each species in loop because we don't want to store re-samples (I think)
for(i in 1:1000){
  assign(paste("resample", i, sep = "_"), 
         data.frame(matrix(vector(), 0, 2,
                           dimnames=list(c(), c("site_cond", "abundance"))),
                    stringsAsFactors=F))}

# Put empty dataframes into a list
rs_list <- mget(ls(pattern = "resample_")) #rs for resample
rm(list=ls(pattern="resample_")) #useful code to remove many items from environment at once, matching pattern
nameorder <- c(names(rs_list))
rs_list <- rs_list[str_sort(nameorder, numeric = TRUE)]

# Create list of 16 dataframes to store bootstrapped xc95 values for each species in each outbug
# dataframe details: ncol=ncol(samp_list[[i]]), nrow = 1000
# each column is named after the species, each row is an XC95 value for that species
for(i in 1:length(samp_sizes)){
  assign(paste("bs_xc95", names(samp_list[i]), sep = "_"), 
         data.frame(matrix(vector(), 0, ncol = ncol(samp_list[[i]])-1,
                           dimnames = list(c(), c(colnames(samp_list[[i]][,2:ncol(samp_list[[i]])])))),
                    stringsAsFactors=F)) #creates blank df of name xc_vals_i
}
bs_xc95s <- mget(ls(pattern = "bs_xc95_outbug"))
rm(list=ls(pattern="bs_xc95_outbug_"))
nameorder <- c(names(bs_xc95s))
bs_xc95s <- bs_xc95s[str_sort(nameorder, numeric = TRUE)]

# WARNING: loop takes ~5-6 hrs to run on laptop
for(h in 1:length(samp_list)){ #first loop specifies that we're going to do the next two loops 16 different times, 1 for each rarify
  for(i in 2:length(samp_list[[h]])){ #second loop isolates species to re-sample
    sp_i <- samp_list[[h]][,c(1, i)]
    sp_i <- sp_i[complete.cases(sp_i),] #only want to sample from sites that have presences
    for(j in 1:length(rs_list)){ 
      #third loop resamples 1000 times, calculates values for cdf, 
      x <- sample(sp_i$site_cond, size = nrow(sp_i), replace = TRUE) #this samples sites w/ replacement, and replicates that n times
      x <- as.data.frame(x)
      x <- x %>% 
        rename(site_cond = x) 
      rs_list[[j]] <- rs_list[[j]][0,] #need to clear current dataframe of past species data before rbind() and left_join() will work properly
      rs_list[[j]] <- rbind(rs_list[[j]][,1], x)
      rs_list[[j]][,1] <- as.numeric(rs_list[[j]][,1])
      #need to reassign abundance values to site values in resampled iterations, do through left_join
      rs_list[[j]] <- left_join(rs_list[[j]], sp_i, by = "site_cond")
      cdfunc <- ewcdf(rs_list[[j]][,1], rs_list[[j]][,2])
      rs_list[[j]]$ewcdf_vals <- cdfunc(rs_list[[j]][,1]) 
      sort_xy <- sortedXyData(expression(site_cond), expression(ewcdf_vals), rs_list[[j]])
      XC95 <- NLSstClosestX(sort_xy, yval = 0.95)
      #store above into new line on blank dataframe
      bs_xc95s[[h]][j,i-1] <- XC95
    }
  }
}

rm(sp_i, x)
#it worked !!

#now that I have 1 list of 16 dfs of 1,000 XC95 values for each species:
#c&s 2013: after they bootstrapped each genus 1000 times to determine XC95 (like I did above), then calculated a 95% CI for each
#genus' XC95 value
#after this, they "generated an HC05 from each of the 1,000 sets of bootstrapped XC95 estimates. 
#"The distribution of 1,000 HC05 values is used to generate two-tailed 95% confidence bounds on these bootstrap-derived values"

#add "real" XC95s from above as an additional row onto each of these???

#create function to calculate 95% CI
normConfInt <- function(x, alpha = 0.05)
  mean(x) + qt(1 - alpha / 2, length(x) - 1) * sd(x) / sqrt(length(x)) * c(-1, 1) #credit: julius vainora on stackoverflow

#create loop to calculate mean XC95 and 95% CI for each species in each dataframe of the bs_xc95s list
#should be stored as: 1 df per rarify "outbug", with columns for: species, mean XC95, upper limit, lower limit

#start by creating 16 dfs to store this
for(i in 1:length(samp_sizes)){
  assign(paste("ci_xc95", names(samp_list[i]), sep = "_"), 
         data.frame(matrix(vector(), nrow = ncol(samp_list[[i]])-1, ncol = 4,
                           dimnames = list(c(), 
                                           c("Species_ID", "mean_XC95", "low_95ci", "up_95ci"))), 
                    #if using quantile line of code instead of normConfInt (below), columns 3 & 4 aren't actually ci's
                    stringsAsFactors=F))
}

ci_xc95 <- mget(ls(pattern = "ci_xc95_outbug"))
rm(list=ls(pattern="ci_xc95_outbug_"))
nameorder <- c(names(ci_xc95))
ci_xc95 <- ci_xc95[str_sort(nameorder, numeric = TRUE)]
rm(nameorder)

#create forloop to populate those dfs with species ids, mean xc95s, upper and lower ci bounds
for(k in 1:length(ci_xc95)){
  for(i in 1:ncol(bs_xc95s[[k]])){
    ci_xc95[[k]][i,1] <- colnames(bs_xc95s[[k]])[i]
    ci_xc95[[k]][i,2] <- mean(bs_xc95s[[k]][,i])
    #ci_xc95[[k]][i,3:4] <- normConfInt(bs_xc95s[[k]][,i]) #this computes 95% CI
    ci_xc95[[k]][i,3:4] <- quantile(bs_xc95s[[k]][,i], probs = c(0.025, 0.975)) #this computes 95% distribution bounds
  }
}

#now, need to create an HC05 value for each of the 1,000 resampled versions of the community in each rarify output "outbug"
#this process generates the uncertainty bounds for the HC05

#need to create an SSD for each row of each of the bs_xc95s dataframes
#will use 1 SSD df that will be reused throughout the loop, not storing the SSD data
#will store output HC05 values in 1 dataframe, with 16 columns (1 for each rarify)

#create bootstrapped HC05 storage dataframe
bs_hc05 <- data.frame(matrix(vector(), 0, ncol = 16,
                             dimnames = list(c(), c(names(samp_list)))))

#create loop to populate bs_hc05 from bs_xc95s list 
for(l in 1:length(bs_xc95s)){
  ssd_df <- data.frame(matrix(vector(), nrow = ncol(bs_xc95s[[l]]), ncol = 4, 
                              dimnames = list(c(), c("Species_ID", "Conc", "rank", "p"))))
  for(m in 1:nrow(bs_xc95s[[l]])){
    ssd_df[,1] <- colnames(bs_xc95s[[l]])
    ssd_df[,2] <- as.numeric(bs_xc95s[[l]][m,]) #need to convert to numeric or else it doesn't work
    ssd_df[,3] <- rank(ssd_df[,2], ties.method = "first")
    ssd_df[,4] <- ssd_df[,3]/(nrow(ssd_df)+1) 
    xysort <- sortedXyData(expression(Conc), expression(p), ssd_df)
    bs_hc05[m,l] <- NLSstClosestX(xysort, yval = 0.05)
  }
}

#create dataframe to store mean, upper and lower CI bounds for each of 16 rarify outputs "outbug"
bs_hc05_ci <- data.frame(matrix(vector(), nrow = length(samp_list), ncol = 5,
                                dimnames = list(c(), c("sample size", "mean_bs_hc05", "low_95ci", "up_95ci", "sd"))), 
                         #if using quantile line of code instead of normConfInt (below), columns 3 & 4 aren't actually ci's
                         stringsAsFactors=F)

#create loop to calculate mean, upper and lower CI bounds for hc05 for each samplesize/rarify output 
for(n in 1:ncol(bs_hc05)){
  bs_hc05_ci[n,1] <- samp_sizes[n]
  bs_hc05_ci[n,2] <- mean(bs_hc05[,n])
  #bs_hc05_ci[n,3:4] <- normConfInt(bs_hc05[,n]) #this computes 95% CI
  bs_hc05_ci[n,3:4] <- quantile(bs_hc05[,n], probs = c(0.025, 0.975)) #this computes 95% distribution bounds
  #also compute SD or some measure of variance
  bs_hc05_ci[n,5] <- sd(bs_hc05[,n])
}

#export bootstrapped HC05 results
#write.xlsx(bs_hc05_ci, "../../data/bs_hc05.xlsx")

#create 16 SSD plots similar to c&s 2013 figure 5
#start with 1, if you want all 16 change index in first row - unimportant to manuscript
ggplot(xc95_list[[1]], aes(x = Conc, y = p)) +
  geom_point() +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  stat_ecdf(pad=FALSE) +
  labs(y = "Proportion of species", x = "Conductivity µS/cm") +
  ggtitle(paste0(samp_sizes[i])) + #I think this paste line stopped working feb 2024, prints NA now
  geom_vline(xintercept = hc05_vals[1,2])

#want to add the XC95 points from the 95% distributions for each species

#determine "true" HC05 from entire species assemblage ("abundances" df)
#first, put abundances back in wide format
abundances <- abundances %>% spread(species, abunds)
#keep only species w/ at least 25 presences
abund_25 <- abundances[, colSums(is.na(abundances)) <= 975] 
#create xc95 storage
realXC95s <- data.frame(matrix(vector(), 0, ncol = 2,
                               dimnames = list(c(), c("Species_ID", "XC95"))))

for(i in 2:ncol(abund_25)){          #want it to start at 2 to skip the sites column
  pres_only <- na.omit(as.data.frame(abund_25[,c(1,i)]))
  cdfunc <- ewcdf(pres_only[,1], pres_only[,2])
  pres_only$ewcdf_vals <- cdfunc(pres_only[,1])
  sort_xy <- sortedXyData(expression(site_cond), expression(ewcdf_vals), pres_only)
  XC95 <- NLSstClosestX(sort_xy, yval = 0.95)
  #store above into new line on blank dataframe
  realXC95s[i-1,] <- c(names(pres_only[2]), XC95)
}
rm(pres_only, sort_xy, i, j, h, k, l, m, n, cdfunc, XC95)

### SSD
realSSD <- realXC95s %>%
  rename(Conc = XC95)
realSSD$Conc <- as.numeric(realSSD$Conc)
realSSD <- realSSD[order(realSSD$Conc),]

# cumulative proportion for each genus P is calculated as P=R/(N+1), 
# where R is the rank of the genus and N is the number of genera
realSSD$rank <- rank(realSSD$Conc, ties.method = "first") #decide how to handle ties?
realSSD$p <- realSSD$rank/(nrow(realSSD)+1) 
# make a sortedXyData object, to be used in inverse interpolation to get HC05 value ----
xysort <- sortedXyData(expression(Conc), expression(p), realSSD)
### inverse interpolation ----
realHC05 <- NLSstClosestX(xysort, yval = 0.05)

#could be interesting to see how many taxa would be affected from "true" SSD by the subsampled HC05 values
#ex: if subsample_50 estimates HC05 at 432 but true HC05 is 196, how many taxa (%) are negatively affected by a 432 benchmark?
#try switching the xysort (x becomes y, y becomes x) and then using inverse interpolation
#also do this for bootstrapped mean HC05 vals
xysort <- xysort %>% rename("x"="y", "y"="x")
#create dataframe to store values
overest_effects <- data.frame(matrix(vector(), 0, ncol = 5,
                                     dimnames = list(c(), c("subsample_size", "HC05_val", "true_perc_affected", "bs_mean_HC05",
                                                            "true_perc_affected_bs"))))
#create loop to store the interpolated points
for(i in 1:length(samp_sizes)){
  overest_effects[i,1] <- samp_sizes[i]
  overest_effects[i,2] <- hc05_vals[i,2]
  overest_effects[i,3] <- (NLSstClosestX(xysort, yval = hc05_vals[i,2])*100)
  overest_effects[i,4] <- bs_hc05_ci[i,2]
  overest_effects[i,5] <- (NLSstClosestX(xysort, yval = bs_hc05_ci[i,2])*100)
}

#save to excel sheet
#write.xlsx(overest_effects, "../../data/percent_off_trueHC05.xlsx")

#create lm to add 2nd y-axis
y1 <-(c(max(overest_effects$bs_mean_HC05), min(overest_effects$bs_mean_HC05)))

#y2 is the Minimum and maximum of the secondary axis data.
y2<-(c(max(overest_effects$true_perc_affected_bs), min(overest_effects$true_perc_affected_bs)))

#axis combines y1 and y2 into a dataframe used for regressions.
axis<-cbind(y1,y2)
axis<-data.frame(axis)

#Regression of hc05 to percents:
perc2hc<-lm(formula = y2 ~ y1, data = axis)
p2hc_summary <- summary(lm(formula = y2 ~ y1, data = axis))
p2hc_summary   

#Identifies the intercept and slope of regressing hc05 to percents:
p2hcInt<-p2hc_summary$coefficients[1, 1] 
p2hcSlope<-p2hc_summary$coefficients[2, 1] 

#plot true % affected
library(yarrr)
pal <- piratepal(palette = "nemo")
p1 <- ggplot(overest_effects, aes(x=log(subsample_size), y=bs_mean_HC05)) + #col=log(subsample_size) add in for color points
  geom_point(size=2) +
  theme_classic() +
  geom_hline(yintercept = realHC05, linetype="dashed") + #for HC05 plot
  #geom_hline(yintercept = 5, linetype='dashed') + #for percent plot
  scale_x_continuous(breaks = log(samp_sizes), labels = samp_sizes) +
  scale_y_continuous("HC05 value (µS/cm)", sec.axis = sec_axis(~.*p2hcSlope + p2hcInt, name = "% species extirpated at HC05 value")) +
  #scale_color_gradientn(colours = pal) + #for color plots for presentations
  xlab("Subsample size") + 
  #ggtitle("Subsample size effect on bootstrapped mean HC05") +
  theme(axis.title = element_text(size = 14),
        #plot.title = element_text(size=22, face = 'bold'),
        axis.text = element_text(size = 12, color = 'black'),
        axis.text.x = element_text(angle = 90, hjust = 0.9, vjust=0.5, size=11),
        legend.position = "none") #last line makes value labels fit better by tilting on x-axis
#second axis fit is really close - not sure it fits *exactly* around the 8% mark?

ggsave("../../../manuscript/plots/subsamp-effect-hc05-bw", plot = p1, width = 7, height = 5, units = "in", device = "tiff")


pal <- piratepal(palette = "nemo", trans = 0.5)
#plot true SSD
p2 <- ggplot(realSSD, aes(x = Conc, y = p)) +
  stat_ecdf(pad=FALSE) +
  geom_point(size=1) +
  geom_hline(yintercept = 0.05, linetype = "dashed", linewidth = 1) +
  labs(y = "Proportion of species", x = "Conductivity uS/cm") +
  #ggtitle("'True' SSD - Coenoflex Full Community") +
  geom_vline(xintercept = realHC05, linewidth = 1) +
  #scale_color_gradientn(colours = pal) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12, color = 'black')) +
  geom_label(label="HC05 = 196 µS/cm", x=800, y=.11,
             label.padding = unit(0.55, "lines"), # Rectangle size around label
             label.size = 0.5, size=5)
ggsave("../../../manuscript/plots/truessd", plot = p2, width = 7, height = 5, units = "in", device = "tiff")
rm(xysort)

#plot hc05 pattern between all subsampling iterations
bs_hc05_ci <- bs_hc05_ci %>% rename(subsample_size = sample.size)
hc05_results <- left_join(bs_hc05_ci, hc05_vals, by = "subsample_size")
hc05_results <- hc05_results %>% rename(fs_hc05 = HC05_val) #fs = full sample HC05... need to better rephrase
#write.xlsx(hc05_results, "../../data/hc05_results.xlsx")

ggplot_hc <- hc05_results[,c(1:4)]
ggplot_hc <- ggplot_hc %>% pivot_longer(!subsample_size, names_to = "hc05_type", values_to = "hc05_val")
ggplot(hc05_results, aes(y=fs_hc05, x=subsample_size)) +
  scale_x_continuous(trans='log10') +
  geom_point() +
  geom_smooth(se=FALSE, fullrange = TRUE) +
  theme_bw() +
  xlab("Subsample size") +
  ylab("HC05") +
  geom_hline(yintercept = realHC05, linetype = "dashed") + #would be cool to add bootstrapped range around this realHC05 value?
  geom_ribbon(aes(ymin=low_95ci, ymax=up_95ci), alpha=0.25) +
  theme(text = element_text(size = 13),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13))


###standard dev calcs and plots -----
#standard dev should be col 5 in hc05_results, if not, rerun because wrong env saved.

#or, if standard dev not in pre-saved env, could recalculate from bs_hc05 df
#each column is a subsample size
hc05_results$sd <- rep(NA, 16)

for(k in 1:ncol(bs_hc05)){
  hc05_results[k,6] <- sd(bs_hc05[,k])
}

ggsave("../../../manuscript/plots/subsamp-effect-hc05", plot = p1, width = 7, height = 5, units = "in", device = "tiff")


# Could put scatter behind each subsample size to show variability?
bp_df <- pivot_longer(bs_hc05, cols=starts_with("outbug"), names_to = "sample_count", values_to = "bs_hc05", names_prefix = "outbug_" )
bp_df$sample_count <- as.numeric(bp_df$sample_count)

ggplot(bp_df, aes(x=sample_count, y=bs_hc05, group=sample_count)) +
  geom_jitter(alpha=0.25, aes(col=log(sample_count))) +
  geom_boxplot(outlier.shape=NA, alpha=0.75) +
  scale_color_gradientn(colours = pal) +
  scale_x_continuous(trans='log10') +
  theme_bw() +
  labs(x = "Subsample Size", y = "HC05") +
  theme(legend.position = "none",
        axis.title = element_text(size = 16, face = 'bold'),
        plot.title = element_text(size=18, face = 'bold'),
        axis.text = element_text(face='bold', size = 12))

# export hc05_results and use in a NRSA script to see what % of streams would be declared imperiled by each HC05
#write.csv(hc05_results, file = "hc05_results.csv")


##############################
# PROBLEM FEB 2024: somehow the same species did not get bootstrapped/sampled when re-running this code, 
# which means the below code no longer works 
# (because species must be present in all subsample sizes to work in this code)

# likely due to tidyverse or R updates since first running the code
# for just the example plot below, load in the saved bootstrap results, post_bootstrap.R.RData


### Pick a few species from the simulation that differ in abundance, range, etc. and plot out how their XC95s change
### for each subsample iteration. We need tangible examples of how sample count is affecting XC95s to explain 
### the differences we are seeing in HC05s.

# to do this: get XC95s from "xc95_list" for 20 species?
# pick species IDs and put them in a vector for the forloop, also use as col headers

spvec <- c('V152', 'V7', 'V165', 'V177', 'V184', 'V12', 'V179', 'V201', 'V140', 'V170', 'V191', 'V89', 'V13', 'V84', 'V49',
           'V224', 'V252', 'V265', 'V376', 'V59')

# write as a loop, store in one df
xc95_effects <- data.frame(matrix(vector(), 0, ncol = length(spvec)+1,
                                  dimnames = list(c(), c("subsample_size", 'V152','V7', 'V165', 'V177', 'V184', 'V12', 'V179',
                                                         'V201', 'V140', 'V170', 'V191', 'V89', 'V13', 'V84', 'V49', 'V224', 
                                                         'V252', 'V265', 'V376', 'V59'))))

for(k in 2:ncol(xc95_effects)){
  for(y in 1:length(samp_sizes)){
    xc95_effects[y,1] <- samp_sizes[y]  
    xc95_effects[y,k] <- mean(bs_xc95s[[y]][,spvec[k-1]]) 
  }}

#cool, plot (have to plot species by species?)
xc95_effects <- xc95_effects %>% pivot_longer(!subsample_size, names_to = "species_ID", values_to = "mean_xc95_val")
ggplot(subset(xc95_effects, species_ID==spvec[20]), aes(x=subsample_size, y=mean_xc95_val)) +
  geom_point()+
  theme_bw()+
  scale_x_continuous(trans='log10') 

### Can I create a figure that shows change in XC95 (+-) from true XC95 (from full dataset) in %?
# add true XC95 to the xc95_effects df

xc95_effects <- xc95_effects %>% pivot_wider(names_from = subsample_size, values_from = mean_xc95_val, id_cols = species_ID) %>%
  rename(Species_ID = species_ID)
xc95_effects <- left_join(xc95_effects, realXC95s, by = "Species_ID")
xc95_effects <- xc95_effects %>% rename(trueXC95 = XC95) %>%
  pivot_longer(!c(Species_ID, trueXC95), names_to = "subsample_size", values_to = "mean_xc95_val") 
xc95_effects$trueXC95 <- as.numeric(xc95_effects$trueXC95)

xc95_effects <- xc95_effects %>% mutate(diff_from_true = ((mean_xc95_val/trueXC95)-1)*100) #puts it into a percent

xc95_effects$subsample_size <- factor(xc95_effects$subsample_size, levels=samp_sizes) #puts x-axis in correct order for plot below

p3 <- ggplot(xc95_effects, aes(x=subsample_size, y=diff_from_true)) + 
  geom_boxplot() +
  theme_classic() +
  labs(x = "Subsample size", y = "% difference from true XC95 value") +
  theme(legend.position = "none",
        axis.title = element_text(size = 14),
        axis.text.y = element_text(size=12, colour="black"),
        axis.text.x = element_text(angle = 90, hjust = 0.9, vjust=0.5, size=12, colour="black"))

ggsave("../../../manuscript/plots/xc95change", plot = p3, width = 7, height = 5, units = "in", device = "tiff")

### Calculate the range (max – min) of the true values from the simulation as a measure of niche width. 
### We can rank species by minimum conductivity value and plot ‘stacked’ ranges across the 2000 uS gradient. 
### That should allow us to visualize if there are more ‘specialists’ at low conductivities and more generalists at high conductivities.

# use abund_25, remove # from front of one of the arrange options below depending on order preference

long_abund_25 <- abund_25 %>% pivot_longer(!site_cond, names_to = 'Species_ID', values_to = 'abund')
long_abund_25 <- na.omit(long_abund_25)
nichewidth <- long_abund_25 %>% group_by(Species_ID) %>% 
  summarise(count = n(), min_cond=min(site_cond), max_cond=max(site_cond), range = max_cond-min_cond, tot_abund = sum(abund)) 
#  %>% arrange(min_cond, max_cond)
#nichewidth$order <- rev(seq(1:nrow(nichewidth)))

# Add in the true XC95 values
nichewidth <- left_join(nichewidth, realXC95s, by = "Species_ID")
nichewidth$XC95 <- round(as.numeric(nichewidth$XC95, digits = 5))
nichewidth <- nichewidth %>% arrange(XC95)
nichewidth$order <- rev(seq(1:nrow(nichewidth)))

#write.csv(nichewidth, file = "sim_nichewidths.csv")

### ----------------------------------------------------------------------------------------------------------------
### MANUSCRIPT PLOTS

ggplot(nichewidth, aes(x=min_cond, y=order)) +
  geom_linerange(aes(xmin=min_cond, xmax=max_cond, alpha=desc(range))) + 
  geom_point(aes(x=XC95), shape=17, col="red")+
  theme_bw() +
  xlim(c(0,2000)) +
  theme(legend.position = "none") +
  labs(y="",x="Conductivity uS/cm")
#subset this to every ~5th species, to improve visualization?

ggplot(nichewidth[seq(1, nrow(nichewidth), 5),], aes(x=min_cond, y=order)) + #every 5th species
  geom_linerange(aes(xmin=min_cond, xmax=max_cond), size=1) + 
  geom_point(aes(x=XC95), shape=23, col="black", fill="red", alpha=0.8, size=3)+
  theme_bw() +
  xlim(c(0,2000)) +
  labs(y="",x="Conductivity uS/cm") +
  theme(legend.position = "none",
        axis.title = element_text(size = 16, face = 'bold'),
        plot.title = element_text(size=18, face = 'bold'),
        axis.text = element_text(face='bold', size = 12))

a <- ggplot(nichewidth, aes(y=range, x=tot_abund)) +
  geom_point() +
  theme_classic() +
  labs(y = "Niche width (range)", x = "Millions of individuals") +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.x = element_text(size=12, colour="black"),
        axis.title.x = element_text(size = 14),
        axis.text.y=element_text(size=12, colour="black"))+
  scale_x_continuous(labels = function(x) format(x/1000000))

b <- ggplot(nichewidth, aes(y=range, x=count)) +
  geom_point() +
  theme_classic() +
  labs(y = "Niche width (range)", x = "Number of sites present") +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.x = element_text(size=12, colour="black"),
        axis.title.x = element_text(size = 14),
        axis.text.y=element_text(size=12, colour="black"))

### need to create a plot for a shared y-axis to add to patchwork plot (below)
p4 <- ggplot(data.frame(l = "Niche width (μS/cm)", x = 1, y = 1)) +
  geom_text(aes(x, y, label = l), angle = 90, size=5) + 
  theme_void() +
  coord_cartesian(clip = "off")

c <- p4 + (a/b) + plot_layout(widths=c(1,25))

ggsave("../../../manuscript/plots/nichewidthvsbothflipped", plot = c, width = 7, height = 8, units = "in", device = "tiff")

d <- ggplot(nichewidth, aes(x=tot_abund, y=XC95)) +
  geom_point() +
  theme_classic() +
  labs(x = "Millions of individuals", y = "XC95 value (μS/cm)") +
  theme(legend.position = "none",
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12, colour="black")) +
  scale_x_continuous(labels = function(x) format(x/1000000)) #to force x-axis into scientific notation

ggsave("../../../manuscript/plots/abundvsXC95", plot = d, width = 7, height = 5, units = "in", device = "tiff")


## Plot a frequency distribution of niche width
ggplot(nichewidth, aes(range)) +
  geom_histogram() +
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_text(size = 16),
        plot.title = element_text(size=18),
        axis.text = element_text(size = 12, colour="black")) +
  labs(x = 'Niche Width', y = 'Frequency')


### -

#----
### Make rank abundance plots for single sites - 3-5 along gradient 

#make a vector of row #s of sites you're interested in
sites <- c(26, 126, 251, 501, 750, 901)
conds <- c(50, 250, 500, 1001, 1499, 1801)
samps <- c(1,4,6,9,15,17) # the 5 subsample index places in samp_list[[]]
sampsizeRA <- c(50, 150, 300, 600, 5000, 'all')

#create a nested list to store results in, main list shows lists for each sample count, within each sample count list is dataframes by site
#create each list of sites for each sample count in a loop, then combine them all
for(i in 1:length(sites)){
  assign(paste("site", sites[i], sep = "_"), 
         data.frame(matrix(vector(), 0, 0))) 
}
#put dataframes from loop into a list
sitelist <- mget(ls(pattern = "site_"))
rm(list=ls(pattern="site_")) 
nameorder <- c(names(sitelist))
sitelist <- sitelist[str_sort(nameorder, numeric = TRUE)]
#ok now we want each site list to be named a subsample size, and then put all the subsample size named lists into a list

for(g in 1:length(samps)){
  assign(paste("sampsize", sampsizeRA[g], sep = "_"), 
         data.frame(matrix(vector(), 0, 0))) 
}
RAnest <- mget(ls(pattern = "sampsize_"))
rm(list=ls(pattern="sampsize_")) 
index <- c(2,3,1,5,4,6)
RAnest  <- RAnest[order(index)] ####THIS ONLY WORKS SOMETIMES, TRIPLE CHECK ORDER IS HOW YOU WANT

for(g in 1:length(samps)){
  RAnest[[g]] <- sitelist
}

#Nested list done. Now to populate.

for(k in 1:length(samps)){ #run this 6 times, for our 6 main lists
  for(g in 1:length(sites)){ #run this for each site for each of the above lists
    df1 <- as.data.frame(samp_list[[samps[k]]][sites[g], c(2:ncol(samp_list[[samps[k]]]))])
    df1 <- df1[, colSums(is.na(df1))==0]
    df1 <- as.data.frame(rankabundance(df1))
    #df1$proportion <- log10(df1$proportion)
    RAnest[[k]][[g]] <- df1 #ok now I need to deposit results into a lil df in a lil list BUT I don't know how to reference dif lists 
  }
} 

#all didn't fill in bc that data is in abund_25, not samp_list, new loop to fill in last bit
for(g in 1:length(sites)){ #run this for each site for each of the above lists
  df1 <- as.data.frame(abund_25[sites[g], c(2:ncol(abund_25))])
  df1 <- df1[, colSums(is.na(df1))==0]
  df1 <- as.data.frame(rankabundance(df1))
  #df1$proportion <- log10(df1$proportion)
  RAnest[[6]][[g]] <- df1 #ok now I need to deposit results into a lil df in a lil list BUT I don't know how to reference dif lists 
}

#RAnest is populated
rm(df1)
write.csv(RAnest[[5]][[2]], "ex_RA3.csv")

#loop through to create 6 ggplots, 1 for each site
library(viridis)
pal=viridis(n=6)

plots <- list()

for(f in 1:length(sites)){
  plots[[f]] <- ggplot(RAnest[[1]][[f]], aes(y=proportion, x=rank), col=pal[1]) +
    geom_point() +
    geom_point(data=RAnest[[2]][[f]], aes(y=proportion, x=rank), col=pal[2]) +
    geom_point(data=RAnest[[3]][[f]], aes(y=proportion, x=rank), col=pal[3]) +
    geom_point(data=RAnest[[4]][[f]], aes(y=proportion, x=rank), col=pal[4]) +
    geom_point(data=RAnest[[5]][[f]], aes(y=proportion, x=rank), col=pal[5]) +
    geom_point(data=RAnest[[6]][[f]], aes(y=proportion, x=rank), col=pal[6]) +
    theme_bw() +
    labs(x = "", y = "") +
    ggtitle(f) +
    theme(axis.title = element_text(size = 14),
          plot.title = element_text(size=14),
          axis.text = element_text(size = 12))
}

#ggsave("../../../manuscript/plots/RAbysites", plot = z, width = 7.5, height = 7, units = "in", device = "tiff")


##### --------
### For manuscript. create a table w/ the following info:
### -gradient location
### subsample count
### % of individuals in the top 5 most common taxa
### No. of taxa w/ <=1% of individuals
### % of taxa w/ <=1% of individuals 

### To create this figure, going to need to calculate the last 2 points from the dfs in the RAnest nested lists

#create a df to store the table info, w/ 6 sites and 6 sample counts (50, 150, 300, 600, 5000, all)

RA_table <- data.frame(matrix(vector(), nrow=36, ncol = 7,
                              dimnames = list(c(), c('target_count', 'real_count', 'site_cond', 'perc_ind_top5', 'no.taxa_1perc', 'perc.taxa_1perc', 'no.taxa'))))
RA_table[,1] <- c(rep(50,6), rep(150,6), rep(300,6), rep(600,6), rep(5000,6), rep('all',6))
RA_table[,3] <- rep(conds, 6)

#to get col2, colSum col2 in each RAnest df (summing abundance)
#ex:
#RA_table[1,2] <- sum(RAnest[[1]][[1]][,2]) #first 1 in index = sampsize, second 1 = site

#to get col4, will want to sum values in col3 for only the first 5 rows of RAnest
#ex:
#RA_table[1,4] <- sum(RAnest[[1]][[1]][c(1:5),3])

#to get col5, count # of taxa that have a value <=1 in col3 of RANest
#ex:
#RA_table[1,5] <- nrow(RAnest[[1]][[1]][(RAnest[[1]][[1]][,3]<=1),])

k=1 #to properly index row in loop below
# ok, loops
for(i in 1:length(samps)){
  for(h in 1:length(sites)){
    RA_table[k,2] <- sum(RAnest[[i]][[h]][,2]) 
    RA_table[k,4] <- sum(RAnest[[i]][[h]][c(1:5),3])
    RA_table[k,5] <- nrow(RAnest[[i]][[h]][(RAnest[[i]][[h]][,3]<=1),]) #this being 1 excludes any taxa from 50 - may change to 2.5?
    RA_table[k,6] <- (RA_table[k,5]/nrow(RAnest[[i]][[h]])*100)
    RA_table[k,7] <- nrow(RAnest[[i]][[h]])
    k= 1+k
  }
}

write.csv(RA_table, "table2.csv")