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

# Total RRS
corrdf<-merged[, c('Core', 'rrs_total', 'bio_sex', 'demo_age')]
corr_core<-spcor(corrdf, method="spearman")
corrdf<-merged[, c('DMPFC', 'rrs_total', 'bio_sex', 'demo_age')]
corr_dmpfc<-spcor(corrdf, method="spearman")
corrdf<-merged[, c('Core_DMPFC', 'rrs_total', 'bio_sex', 'demo_age')]
corr_coredmpfc<-spcor(corrdf, method="spearman")

# RRS reflection
corrdf<-merged[, c('Core', 'reflection_total', 'bio_sex', 'demo_age')]
corr_core<-spcor(corrdf, method="spearman")
corrdf<-merged[, c('DMPFC', 'reflection_total', 'bio_sex', 'demo_age')]
corr_dmpfc<-spcor(corrdf, method="spearman")
corrdf<-merged[, c('Core_DMPFC', 'reflection_total', 'bio_sex', 'demo_age')]
corr_coredmpfc<-spcor(corrdf, method="spearman")

# RRS brooding
corrdf<-merged[, c('Core', 'brooding_total', 'bio_sex', 'demo_age')]
corr_core<-spcor(corrdf, method="spearman")
corrdf<-merged[, c('DMPFC', 'brooding_total', 'bio_sex', 'demo_age')]
corr_dmpfc<-spcor(corrdf, method="spearman")
corrdf<-merged[, c('Core_DMPFC', 'brooding_total', 'bio_sex', 'demo_age')]
corr_coredmpfc<-spcor(corrdf, method="spearman")

# RRS deprelated
corrdf<-merged[, c('Core', 'deprelated_total', 'bio_sex', 'demo_age')]
corr_core<-spcor(corrdf, method="spearman")
corrdf<-merged[, c('DMPFC', 'deprelated_total', 'bio_sex', 'demo_age')]
corr_dmpfc<-spcor(corrdf, method="spearman")
corrdf<-merged[, c('Core_DMPFC', 'deprelated_total', 'bio_sex', 'demo_age')]
corr_coredmpfc<-spcor(corrdf, method="spearman")

# Test for an interaction effect between group and total RRS 

# MDD/anxiety or not
core_inter_lm=lm(Core ~ rrs_total+curranxmdd+rrs_total*curranxmdd+bio_sex+demo_age, data=merged) 
dmpfc_inter_lm=lm(DMPFC ~ rrs_total+curranxmdd+rrs_total*curranxmdd+bio_sex+demo_age, data=merged) 
coredmpfc_inter_lm=lm(Core_DMPFC ~ rrs_total+curranxmdd+rrs_total*curranxmdd+bio_sex+demo_age, data=merged) 
# current MDD or not 
core_inter_lm=lm(Core ~ rrs_total+mini7_mdd_current+rrs_total*mini7_mdd_current+bio_sex+demo_age, data=merged) 
dmpfc_inter_lm=lm(DMPFC ~ rrs_total+mini7_mdd_current+rrs_total*mini7_mdd_current+bio_sex+demo_age, data=merged) 
coredmpfc_inter_lm=lm(Core_DMPFC ~ rrs_total+mini7_mdd_current+rrs_total*mini7_mdd_current+bio_sex+demo_age, data=merged) 


# Test for an interaction effect between group and RRS reflection 

# MDD/anxiety or not
core_inter_lm=lm(Core ~ reflection_total+curranxmdd+reflection_total*curranxmdd+bio_sex+demo_age, data=merged) 
dmpfc_inter_lm=lm(DMPFC ~ reflection_total+curranxmdd+reflection_total*curranxmdd+bio_sex+demo_age, data=merged) 
coredmpfc_inter_lm=lm(Core_DMPFC ~ reflection_total+curranxmdd+reflection_total*curranxmdd+bio_sex+demo_age, data=merged) 
# current MDD or not 
core_inter_lm=lm(Core ~ reflection_total+mini7_mdd_current+reflection_total*mini7_mdd_current+bio_sex+demo_age, data=merged) 
dmpfc_inter_lm=lm(DMPFC ~ reflection_total+mini7_mdd_current+reflection_total*mini7_mdd_current+bio_sex+demo_age, data=merged) 
coredmpfc_inter_lm=lm(Core_DMPFC ~ reflection_total+mini7_mdd_current+reflection_total*mini7_mdd_current+bio_sex+demo_age, data=merged) 

# Test for an interaction effect between group and RRS brooding 

# MDD/anxiety or not
core_inter_lm=lm(Core ~ brooding_total+curranxmdd+brooding_total*curranxmdd+bio_sex+demo_age, data=merged) 
dmpfc_inter_lm=lm(DMPFC ~ brooding_total+curranxmdd+brooding_total*curranxmdd+bio_sex+demo_age, data=merged) 
coredmpfc_inter_lm=lm(Core_DMPFC ~ brooding_total+curranxmdd+brooding_total*curranxmdd+bio_sex+demo_age, data=merged) 
# current MDD or not 
core_inter_lm=lm(Core ~ brooding_total+mini7_mdd_current+brooding_total*mini7_mdd_current+bio_sex+demo_age, data=merged) 
dmpfc_inter_lm=lm(DMPFC ~ brooding_total+mini7_mdd_current+brooding_total*mini7_mdd_current+bio_sex+demo_age, data=merged) 
coredmpfc_inter_lm=lm(Core_DMPFC ~ brooding_total+mini7_mdd_current+brooding_total*mini7_mdd_current+bio_sex+demo_age, data=merged) 

# Test for an interaction effect between group and RRS deprelated 

# MDD/anxiety or not
core_inter_lm=lm(Core ~ deprelated_total+curranxmdd+deprelated_total*curranxmdd+bio_sex+demo_age, data=merged) 
dmpfc_inter_lm=lm(DMPFC ~ deprelated_total+curranxmdd+deprelated_total*curranxmdd+bio_sex+demo_age, data=merged) 
coredmpfc_inter_lm=lm(Core_DMPFC ~ deprelated_total+curranxmdd+deprelated_total*curranxmdd+bio_sex+demo_age, data=merged) 
# current MDD or not 
core_inter_lm=lm(Core ~ deprelated_total+mini7_mdd_current+deprelated_total*mini7_mdd_current+bio_sex+demo_age, data=merged) 
dmpfc_inter_lm=lm(DMPFC ~ deprelated_total+mini7_mdd_current+deprelated_total*mini7_mdd_current+bio_sex+demo_age, data=merged) 
coredmpfc_inter_lm=lm(Core_DMPFC ~ deprelated_total+mini7_mdd_current+deprelated_total*mini7_mdd_current+bio_sex+demo_age, data=merged) 

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

# Helper function for plots
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

# Plot connectivity difference in subnetwork between Anx+MDD and No MDD for Core
png('visualizations/core_anxmdd.png',width=10,height=10,units='in',res=300)
ggplot(merged_agesexreg, aes(x=curranxmdd, y=Core, fill=curranxmdd)) + 
  geom_violin() +  
  scale_x_discrete(labels=c("No Anx or MDD","Anx or MDD")) +
  theme_minimal(base_size=22) +
  labs(y="Core functional connectivity (adjusted)", x='') + 
  theme(legend.position = "none") + 
  stat_summary(fun.data=data_summary)
dev.off()

# Plot connectivity difference in subnetwork between Anx+MDD and No MDD for DMPFC
png('visualizations/dmpfc_anxmdd.png',width=10,height=10,units='in',res=300)
ggplot(merged_agesexreg, aes(x=curranxmdd, y=DMPFC, fill=curranxmdd)) + 
  geom_violin() +  
  scale_x_discrete(labels=c("No Anx or MDD","Anx or MDD")) +
  theme_minimal(base_size=22) +
  labs(y="DMPFC functional connectivity (adjusted)", x='') + 
  theme(legend.position = "none") + 
  stat_summary(fun.data=data_summary)
dev.off()

# Plot connectivity difference in subnetwork between Anx+MDD and No MDD for Core-DMPFC
png('visualizations/coredmpfc_anxmdd.png',width=10,height=10,units='in',res=300)
ggplot(merged_agesexreg, aes(x=curranxmdd, y=Core_DMPFC, fill=curranxmdd)) + 
  geom_violin() +  
  scale_x_discrete(labels=c("No Anx or MDD","Anx or MDD")) +
  theme_minimal(base_size=22) +
  labs(y="Core-DMPFC functional connectivity (adjusted)", x='') + 
  theme(legend.position = "none") + 
  stat_summary(fun.data=data_summary)
dev.off()

# Plot connectivity difference in subnetwork between MDD and No MDD for Core
png('visualizations/core_mdd.png',width=10,height=10,units='in',res=300)
ggplot(merged_agesexreg, aes(x=currmdd, y=Core, fill=currmdd)) + 
  geom_violin() +  
  scale_x_discrete(labels=c("No MDD","MDD")) +
  theme_minimal(base_size=22) +
  labs(y="Core functional connectivity (adjusted)", x='') + 
  theme(legend.position = "none") + 
  stat_summary(fun.data=data_summary)
dev.off()

# Plot connectivity difference in subnetwork between MDD and No MDD for DMPFC
png('visualizations/dmpfc_mdd.png',width=10,height=10,units='in',res=300)
ggplot(merged_agesexreg, aes(x=currmdd, y=DMPFC, fill=currmdd)) + 
  geom_violin() +  
  scale_x_discrete(labels=c("No MDD","MDD")) +
  theme_minimal(base_size=22) +
  labs(y="DMPFC functional connectivity (adjusted)", x='') + 
  theme(legend.position = "none") + 
  stat_summary(fun.data=data_summary)
dev.off()

# Plot connectivity difference in subnetwork between MDD and No MDD for Core-DMPFC
png('visualizations/coredmpfc_mdd.png',width=10,height=10,units='in',res=300)
ggplot(merged_agesexreg, aes(x=currmdd, y=Core_DMPFC, fill=currmdd)) + 
  geom_violin() +  
  scale_x_discrete(labels=c("No MDD","MDD")) +
  theme_minimal(base_size=22) +
  labs(y="Core-DMPFC functional connectivity (adjusted)", x='') + 
  theme(legend.position = "none") + 
  stat_summary(fun.data=data_summary)
dev.off()




