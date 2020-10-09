library(dplyr)
library(ggcorrplot)
library(lsr)
library(ppcor)

setwd("/Users/ltozzi/Desktop/Oak/ltozzi/hcpdes/dmn_subsystems_project/")

clin<-read.csv('clinical/hcpdes_nohist.csv')
dmn<-read.csv('dmn_subsystems/DMNsubsys_fc_hcpdes.csv')

# Remove subjects with nans
dmn<-na.omit(dmn)
merged<-merge(clin, dmn, by='ID')

# Create composite diagnosis columns
merged$currmdd=factor(merged$mini7_mdd_current==1)
merged$curranx=factor(merged$mini7_panic_current | merged$mini7_gad_current | merged$mini7_agoraphobia_current |  merged$mini7_sad )
merged$currmdd_only=factor(merged$currmdd==TRUE & merged$curranx==FALSE)
merged$curranx_only=factor(merged$currmdd==FALSE & merged$curranx==TRUE)
merged$curranxmdd=factor(merged$currmdd==TRUE & merged$curranx==TRUE)

# merged<-merged[!merged$ID=='sub-CONN147', ] # remove potential outlier

# Run Spearman partial correlations between connectivity and ruminations given age and sex
corrdf<-merged[, c('Core', 'rrs_total', 'bio_sex', 'demo_age')]
corr_core<-spcor(corrdf, method="spearman")
corrdf<-merged[, c('DMPFC', 'rrs_total', 'bio_sex', 'demo_age')]
corr_dmpfc<-spcor(corrdf, method="spearman")
corrdf<-merged[, c('Core_DMPFC', 'rrs_total', 'bio_sex', 'demo_age')]
corr_coredmpfc<-spcor(corrdf, method="spearman")

# Regress age and sex from data at each site for each network
merged_agesexreg<-merged
nws=c('Core', 'DMPFC', 'MTL', 'Core_DMPFC')
for (nw in nws){
  agesex.mdl<-lm(paste(nw,'~demo_age+bio_sex'),data=merged_agesexreg)
  merged_agesexreg[, nw]<-agesex.mdl$residuals
}

# T-test between any MDD/anxiety or not 
t_core_diag<-t.test(merged_agesexreg$Core ~ merged_agesexreg$curranxmdd, alternative='greater')
t_dmpfc_diag<-t.test(merged_agesexreg$DMPFC ~ merged_agesexreg$curranxmdd, alternative='greater')
t_coredmpfc_diag<-t.test(merged_agesexreg$Core_DMPFC ~ merged_agesexreg$curranxmdd, alternative='greater')

# Convert to Cohen's d
d_core_diag<-cohensD(merged_agesexreg$Core ~ merged_agesexreg$curranxmdd)
d_dmpfc_diag<-cohensD(merged_agesexreg$DMPFC ~ merged_agesexreg$curranxmdd)
d_coredmpfc_diag<-cohensD(merged_agesexreg$Core_DMPFC ~ merged_agesexreg$curranxmdd)

# T-test between current MDD or not 
t_core_mdd<-t.test(merged_agesexreg$Core ~ merged_agesexreg$mini7_mdd_current, alternative='greater')
t_dmpfc_mdd<-t.test(merged_agesexreg$DMPFC ~ merged_agesexreg$mini7_mdd_current, alternative='greater')
t_coredmpfc_mdd<-t.test(merged_agesexreg$Core_DMPFC ~ merged_agesexreg$mini7_mdd_current, alternative='greater')

# Convert to Cohen's d
d_core_mdd<-cohensD(merged_agesexreg$Core ~ merged_agesexreg$mini7_mdd_current)
d_dmpfc_mdd<-cohensD(merged_agesexreg$DMPFC ~ merged_agesexreg$mini7_mdd_current)
d_coredmpfc_mdd<-cohensD(merged_agesexreg$Core_DMPFC ~ merged_agesexreg$mini7_mdd_current)


# Plot connectivity difference in subnetwork between Anx+MDD and No MDD for Core
merged_agesexreg_summary <- merged_agesexreg %>% # the names of the new data frame and the data frame to be summarised
  group_by(curranxmdd) %>%   # the grouping variable
  summarise(mean_PL = mean(Core, na.rm = TRUE),  # calculates the mean of each group
            sd_PL = sd(Core, na.rm = TRUE), # calculates the standard deviation of each group
            n_PL = n(),  # calculates the sample size per group
            SE_PL = sd(Core, na.rm = TRUE)/sqrt(n())) # calculates the standard error of each group

png('visualizations/core_anxmdd.png',width=10,height=10,units='in',res=300)
merged_agesexreg_Plot <- ggplot(merged_agesexreg_summary, aes(curranxmdd, mean_PL)) + 
  geom_col(fill = "grey", colour='black') +  
  geom_errorbar(aes(ymin = mean_PL - SE_PL, ymax = mean_PL + SE_PL), width=0.2)+
  scale_x_discrete(labels=c("No Anx or MDD","Anx or MDD"))
merged_agesexreg_Plot + labs(y="Core functional connectivity (adjusted) ± SE", x='') + theme_bw(base_size = 26)
dev.off()

# Plot connectivity difference in subnetwork between Anx+MDD and No MDD for DMPFC
merged_agesexreg_summary <- merged_agesexreg %>% # the names of the new data frame and the data frame to be summarised
  group_by(curranxmdd) %>%   # the grouping variable
  summarise(mean_PL = mean(DMPFC, na.rm = TRUE),  # calculates the mean of each group
            sd_PL = sd(DMPFC, na.rm = TRUE), # calculates the standard deviation of each group
            n_PL = n(),  # calculates the sample size per group
            SE_PL = sd(DMPFC, na.rm = TRUE)/sqrt(n())) # calculates the standard error of each group

png('visualizations/dmpfc_anxmdd.png',width=10,height=10,units='in',res=300)
merged_agesexreg_Plot <- ggplot(merged_agesexreg_summary, aes(curranxmdd, mean_PL)) + 
  geom_col(fill = "grey", colour='black') +  
  geom_errorbar(aes(ymin = mean_PL - SE_PL, ymax = mean_PL + SE_PL), width=0.2)+
  scale_x_discrete(labels=c("No Anx or MDD","Anx or MDD"))
merged_agesexreg_Plot + labs(y="DMPFC functional connectivity (adjusted) ± SE", x='') + theme_bw(base_size = 26)
dev.off()

# Plot connectivity difference in subnetwork between Anx+MDD and No MDD for Core-DMPFC
merged_agesexreg_summary <- merged_agesexreg %>% # the names of the new data frame and the data frame to be summarised
  group_by(curranxmdd) %>%   # the grouping variable
  summarise(mean_PL = mean(Core_DMPFC, na.rm = TRUE),  # calculates the mean of each group
            sd_PL = sd(Core_DMPFC, na.rm = TRUE), # calculates the standard deviation of each group
            n_PL = n(),  # calculates the sample size per group
            SE_PL = sd(Core_DMPFC, na.rm = TRUE)/sqrt(n())) # calculates the standard error of each group

png('visualizations/coredmpfc_anxmdd.png',width=10,height=10,units='in',res=300)
merged_agesexreg_Plot <- ggplot(merged_agesexreg_summary, aes(curranxmdd, mean_PL)) + 
  geom_col(fill = "grey", colour='black') +  
  geom_errorbar(aes(ymin = mean_PL - SE_PL, ymax = mean_PL + SE_PL), width=0.2)+
  scale_x_discrete(labels=c("No Anx or MDD","Anx or MDD"))
merged_agesexreg_Plot + labs(y="Core-DMPFC functional connectivity (adjusted) ± SE", x='') + theme_bw(base_size = 26)
dev.off()


# Plot connectivity difference in subnetwork between MDD and No MDD for Core
merged_agesexreg_summary <- merged_agesexreg %>% # the names of the new data frame and the data frame to be summarised
  group_by(mini7_mdd_current) %>%   # the grouping variable
  summarise(mean_PL = mean(Core, na.rm = TRUE),  # calculates the mean of each group
            sd_PL = sd(Core, na.rm = TRUE), # calculates the standard deviation of each group
            n_PL = n(),  # calculates the sample size per group
            SE_PL = sd(Core, na.rm = TRUE)/sqrt(n())) # calculates the standard error of each group

png('visualizations/core_mdd.png',width=10,height=10,units='in',res=300)
merged_agesexreg_Plot <- ggplot(merged_agesexreg_summary, aes(mini7_mdd_current, mean_PL)) + 
  geom_col(fill = "grey", colour='black') +  
  geom_errorbar(aes(ymin = mean_PL - SE_PL, ymax = mean_PL + SE_PL), width=0.2)+
  scale_x_continuous(breaks=0:1,labels=c("No MDD","MDD"))
merged_agesexreg_Plot + labs(y="Core functional connectivity (adjusted) ± SE", x='') + theme_bw(base_size = 22)
dev.off()

# Plot connectivity difference in subnetwork between MDD and No MDD for DMPFC
merged_agesexreg_summary <- merged_agesexreg %>% # the names of the new data frame and the data frame to be summarised
  group_by(mini7_mdd_current) %>%   # the grouping variable
  summarise(mean_PL = mean(DMPFC, na.rm = TRUE),  # calculates the mean of each group
            sd_PL = sd(DMPFC, na.rm = TRUE), # calculates the standard deviation of each group
            n_PL = n(),  # calculates the sample size per group
            SE_PL = sd(DMPFC, na.rm = TRUE)/sqrt(n())) # calculates the standard error of each group

png('visualizations/dmpfc_mdd.png',width=10,height=10,units='in',res=300)
merged_agesexreg_Plot <- ggplot(merged_agesexreg_summary, aes(mini7_mdd_current, mean_PL)) + 
  geom_col(fill = "grey", colour='black') +  
  geom_errorbar(aes(ymin = mean_PL - SE_PL, ymax = mean_PL + SE_PL), width=0.2)+
  scale_x_continuous(breaks=0:1,labels=c("No MDD","MDD"))
merged_agesexreg_Plot + labs(y="DMPFC functional connectivity (adjusted) ± SE", x='') + theme_bw(base_size = 22)
dev.off()

# Plot connectivity difference in subnetwork between MDD and No MDD for Core-DMPFC
merged_agesexreg_summary <- merged_agesexreg %>% # the names of the new data frame and the data frame to be summarised
  group_by(mini7_mdd_current) %>%   # the grouping variable
  summarise(mean_PL = mean(Core_DMPFC, na.rm = TRUE),  # calculates the mean of each group
            sd_PL = sd(Core_DMPFC, na.rm = TRUE), # calculates the standard deviation of each group
            n_PL = n(),  # calculates the sample size per group
            SE_PL = sd(Core_DMPFC, na.rm = TRUE)/sqrt(n())) # calculates the standard error of each group

png('visualizations/coredmpfc_mdd.png',width=10,height=10,units='in',res=300)
merged_agesexreg_Plot <- ggplot(merged_agesexreg_summary, aes(mini7_mdd_current, mean_PL)) + 
  geom_col(fill = "grey", colour='black') +  
  geom_errorbar(aes(ymin = mean_PL - SE_PL, ymax = mean_PL + SE_PL), width=0.2)+
  scale_x_continuous(breaks=0:1,labels=c("No MDD","MDD"))
merged_agesexreg_Plot + labs(y="Core-DMPFC functional connectivity (adjusted) ± SE", x='') + theme_bw(base_size = 22)
dev.off()


