setwd("/Users/ltozzi/Desktop/Oak/ltozzi/hcpdes/dmn_subsystems_project/")

library(meta)
library(tidyverse)
library(grid)
library(plyr)

source('scripts/meta_analysis_functions.R')
dir.create('meta')
dir.create('meta/tables')
dir.create('meta/plots')


###### WITHIN NETWORK CONNECTIVITY

# Create merged dataframe of clinical and nws
clinical_data=read_csv('clinical/metamdd_CG_LT.csv') # clinical data only for people included in Yan (2019)
nw_data=read_csv('dmn_subsystems/within_DMNsubsys_fc_metamdd.csv') # average network data
merged=merge(clinical_data, nw_data)

# Labeling networks and sites
nws=sort(c('Core','DMPFC', 'MTL'))
sites=sort(unique(merged$Site))

# Label sex as -1 and 1 
merged[merged$Sex==2,'Sex']=-1

# Consider removing S09 per reviewer suggestion 
# merged=merged[merged$Site!='S09',]

################ ANALYSIS

#### AGE AND SEX REGRESSION

# Regress age and sex from data at each site for each network
merged_agesexreg<-merged
for (site in sites){
  sitedata=merged[merged$Site==site,]
  for (nw in nws){
  agesex.mdl<-lm(paste(nw,'~Age+Sex'),data=sitedata)
  merged_agesexreg[merged_agesexreg$Site==site, nw]<-agesex.mdl$residuals
    }
}

# Remove sites with less than 20 people per group
factor_count<-by(merged_agesexreg$Site, merged_agesexreg$Group, count)
sites_keep1<-factor_count[['-1']][factor_count[["-1"]]$freq>=20, 'x']
sites_keep2<-factor_count[['1']][factor_count[["1"]]$freq>=20, 'x']
sites_keep<-intersect(sites_keep1,sites_keep2)

# Run meta-analysis for all networks for comparison between groups
MA_list_hcmdd_agesexreg = list()
for (nw in nws){
  nwdf<-create_ES_df_binintvar(df=merged_agesexreg, sites = sites_keep, interestvar = 'Group', nw=nw)
  MA_list_hcmdd_agesexreg[[nw]]<-run_btw_group_meta(nwdf)
}

# Remove sites with less than 20 values of HAMD
factor_count<-by(merged_agesexreg$Site, !is.na(merged_agesexreg$HAMD), count)
sites_keep<-factor_count$'TRUE'[factor_count[["TRUE"]]$freq>=20,'x']

# Run meta-analysis for all networks for correlation with HAMD
MA_list_hamdcorr_agesexreg = list()
for (nw in nws){
  nwdf<-create_ES_df_corr(df=merged_agesexreg, sites = sites_keep, interestvar = 'HAMD', nw=nw)
  MA_list_hamdcorr_agesexreg[[nw]]<-run_corr_meta(nwdf)
}

# Remove sites with less than 20 values of HAMA
factor_count<-by(merged_agesexreg$Site, !is.na(merged_agesexreg$HAMA), count)
sites_keep<-factor_count$'TRUE'[factor_count[["TRUE"]]$freq>=20,'x']

# Run meta-analysis for all networks for correlation with HAMA
MA_list_hamacorr_agesexreg = list()
for (nw in nws){
  nwdf<-create_ES_df_corr(df=merged_agesexreg, sites = sites_keep, interestvar = 'HAMA', nw=nw)
  MA_list_hamacorr_agesexreg[[nw]]<-run_corr_meta(nwdf)
}

# Remove sites with less than 20 values of Illness duration
factor_count<-by(merged_agesexreg$Site, !is.na(merged_agesexreg$'IllnessDuration_months_'), count)
sites_keep<-factor_count$'TRUE'[factor_count[["TRUE"]]$freq>=20,'x']

# Run meta-analysis for all networks for correlation with Illness duration
MA_list_duration_agesexreg = list()
for (nw in nws){
  nwdf<-create_ES_df_corr(df=merged_agesexreg, sites = sites_keep, interestvar = 'IllnessDuration_months_', nw=nw)
  MA_list_duration_agesexreg[[nw]]<-run_corr_meta(nwdf)
}

# Remove sites with less than 20 values of fdvals
factor_count<-by(merged_agesexreg$Site, !is.na(merged_agesexreg$fdvols), count)
sites_keep<-factor_count$'TRUE'[factor_count[["TRUE"]]$freq>=20,'x']

# Run meta-analysis for all networks for correlation with fdvals
MA_list_fdvols_agesexreg = list()
for (nw in nws){
  nwdf<-create_ES_df_corr(df=merged_agesexreg, sites = sites_keep, interestvar = 'fdvols', nw=nw)
  MA_list_fdvols_agesexreg[[nw]]<-run_corr_meta(nwdf)
}

# Remove sites with less than 20 values of Age
factor_count<-by(merged$Site, !is.na(merged$Age), count)
sites_keep<-factor_count$'TRUE'[factor_count[["TRUE"]]$freq>=20,'x']

# Run meta-analysis for all networks for correlation with age
MA_list_age_agesexreg = list()
for (nw in nws){
  nwdf<-create_ES_df_corr(df=merged, sites = sites_keep, interestvar = 'Age', nw=nw)
  MA_list_age_agesexreg[[nw]]<-run_corr_meta(nwdf)
}

# Remove sites with less than 20 values of sex
factor_count<-by(merged$Site, !is.na(merged$Sex), count)
sites_keep<-factor_count$'TRUE'[factor_count[["TRUE"]]$freq>=20,'x']

# Run meta-analysis for all networks for comparison between sexes
MA_list_sex_agesexreg = list()
for (nw in nws){
  nwdf<-create_ES_df_binintvar(df=merged, sites = sites_keep, interestvar = 'Sex', nw=nw)
  MA_list_sex_agesexreg[[nw]]<-run_btw_group_meta(nwdf)
}


###### EXPORT RESULTS

# Draw and save forest plots
for (nw in nws){
  drawforplot(MA_list_hcmdd_agesexreg[[nw]], paste(nw, '_hcmdd_agesexreg_dmnsubsys', sep=''))
  drawforplot(MA_list_hamdcorr_agesexreg[[nw]], paste(nw, '_hamdcorr_agesexreg_dmnsubsys', sep=''))
  drawforplot(MA_list_hamacorr_agesexreg[[nw]], paste(nw, '_hamacorr_agesexreg_dmnsubsys', sep='')) 
  drawforplot(MA_list_duration_agesexreg[[nw]], paste(nw, '_duration_agesexreg_dmnsubsys', sep='')) 
  drawforplot(MA_list_fdvols_agesexreg[[nw]], paste(nw, '_fdvols_agesexreg_dmnsubsys', sep=''))  
  drawforplot(MA_list_age_agesexreg[[nw]], paste(nw, '_age_agesexreg_dmnsubsys', sep=''))  
  drawforplot(MA_list_sex_agesexreg[[nw]], paste(nw, '_sex_agesexreg_dmnsubsys', sep=''))  
}

# Export group meta-analyses
meta_export_hcmdd_agesexreg <- data.frame(matrix(ncol = 13, nrow = 0))
colnames(meta_export_hcmdd_agesexreg)<-colnames(export_meta_group(meta_export_hcmdd_agesexreg[[nw]], nw))
meta_export_sex_agesexreg <- data.frame(matrix(ncol = 13, nrow = 0))
colnames(meta_export_sex_agesexreg)<-colnames(export_meta_group(meta_export_sex_agesexreg[[nw]], nw))

for (nw in nws){
  meta_export_hcmdd_agesexreg<-rbind(meta_export_hcmdd_agesexreg, export_meta_group(MA_list_hcmdd_agesexreg[[nw]], nw))
  meta_export_sex_agesexreg<-rbind(meta_export_sex_agesexreg, export_meta_group(MA_list_sex_agesexreg[[nw]], nw))
}

# Compute FDR correction
meta_export_hcmdd_agesexreg$effect_pfdr<-p.adjust(meta_export_hcmdd_agesexreg$effect_p, method = "fdr")
meta_export_hcmdd_agesexreg$qpfdr<-p.adjust(meta_export_hcmdd_agesexreg$qp, method = "fdr")
meta_export_sex_agesexreg$effect_pfdr<-p.adjust(meta_export_sex_agesexreg$effect_p, method = "fdr")
meta_export_sex_agesexreg$qpfdr<-p.adjust(meta_export_sex_agesexreg$qp, method = "fdr")

# Export correlation meta-analyses
meta_export_hamd_agesexreg <- data.frame(matrix(ncol = 12, nrow = 0))
colnames(meta_export_hamd_agesexreg)<-colnames(export_meta_corr(MA_list_hamdcorr_agesexreg[[nw]], nw))
meta_export_hama_agesexreg <- data.frame(matrix(ncol = 12, nrow = 0))
colnames(meta_export_hama_agesexreg)<-colnames(export_meta_corr(MA_list_hamacorr_agesexreg[[nw]], nw))
meta_export_duration_agesexreg <- data.frame(matrix(ncol = 12, nrow = 0))
colnames(meta_export_duration_agesexreg)<-colnames(export_meta_corr(MA_list_duration_agesexreg[[nw]], nw))
meta_export_fdvols_agesexreg <- data.frame(matrix(ncol = 12, nrow = 0))
colnames(meta_export_fdvols_agesexreg)<-colnames(export_meta_corr(MA_list_fdvols_agesexreg[[nw]], nw))
meta_export_age_agesexreg <- data.frame(matrix(ncol = 12, nrow = 0))
colnames(meta_export_age_agesexreg)<-colnames(export_meta_corr(MA_list_age_agesexreg[[nw]], nw))

for (nw in nws){
  meta_export_hamd_agesexreg<-rbind(meta_export_hamd_agesexreg, export_meta_corr(MA_list_hamdcorr_agesexreg[[nw]], nw))
  meta_export_hama_agesexreg<-rbind(meta_export_hama_agesexreg, export_meta_corr(MA_list_hamacorr_agesexreg[[nw]], nw))
  meta_export_duration_agesexreg<-rbind(meta_export_duration_agesexreg, export_meta_corr(MA_list_duration_agesexreg[[nw]], nw))
  meta_export_fdvols_agesexreg<-rbind(meta_export_fdvols_agesexreg, export_meta_corr(MA_list_fdvols_agesexreg[[nw]], nw))
  meta_export_age_agesexreg<-rbind(meta_export_age_agesexreg, export_meta_corr(MA_list_age_agesexreg[[nw]], nw))
}

# Compute FDR correction
meta_export_hamd_agesexreg$effect_pfdr<-p.adjust(meta_export_hamd_agesexreg$effect_p, method = "fdr")
meta_export_hama_agesexreg$effect_pfdr<-p.adjust(meta_export_hama_agesexreg$effect_p, method = "fdr")
meta_export_duration_agesexreg$effect_pfdr<-p.adjust(meta_export_duration_agesexreg$effect_p, method = "fdr")
meta_export_fdvols_agesexreg$effect_pfdr<-p.adjust(meta_export_fdvols_agesexreg$effect_p, method = "fdr")
meta_export_age_agesexreg$effect_pfdr<-p.adjust(meta_export_age_agesexreg$effect_p, method = "fdr")

meta_export_hamd_agesexreg$qpfdr<-p.adjust(meta_export_hamd_agesexreg$qp, method = "fdr")
meta_export_hama_agesexreg$qpfdr<-p.adjust(meta_export_hama_agesexreg$qp, method = "fdr")
meta_export_duration_agesexreg$qpfdr<-p.adjust(meta_export_duration_agesexreg$qp, method = "fdr")
meta_export_fdvols_agesexreg$qpfdr<-p.adjust(meta_export_fdvols_agesexreg$qp, method = "fdr")
meta_export_age_agesexreg$qpfdr<-p.adjust(meta_export_age_agesexreg$qp, method = "fdr")

# Save CSVs
write_csv(x = meta_export_hcmdd_agesexreg, path ='meta/tables/summary_hcmdd_agesexreg_dmnsubsys.csv')
write_csv(x = meta_export_hamd_agesexreg, path ='meta/tables/summary_hamd_agesexreg_dmnsubsys.csv')
write_csv(x = meta_export_hama_agesexreg, path ='meta/tables/summary_hama_agesexreg_dmnsubsys.csv')
write_csv(x = meta_export_duration_agesexreg, path ='meta/tables/summary_duration_agesexreg_dmnsubsys.csv')
write_csv(x = meta_export_fdvols_agesexreg, path ='meta/tables/summary_fdvols_agesexreg_dmnsubsys.csv')
write_csv(x = meta_export_age_agesexreg, path ='meta/tables/summary_age_agesexreg_dmnsubsys.csv')
write_csv(x = meta_export_sex_agesexreg, path ='meta/tables/summary_sex_agesexreg_dmnsubsys.csv')

# Run meta-analysis for all networks for comparison of motion between groups
fddf<-create_ES_df_binintvar(df=merged_agesexreg, sites = sites, interestvar = 'Group', nw='fdvols')
MA_fd_hcmdd_agesexreg<-run_btw_group_meta(fddf)

# Compare first episode vs non first episode
factor_count<-by(merged_agesexreg$Site, merged_agesexreg$IfFirstEpisode, count)
sites_keep1<-factor_count[['-1']][factor_count[["-1"]]$freq>=20, 'x']
sites_keep2<-factor_count[['1']][factor_count[["1"]]$freq>=20, 'x']
sites_keep<-intersect(sites_keep1,sites_keep2)

# Compare recurrent and first episode in S20
s20<-merged_agesexreg[merged_agesexreg$Site=='S20',]
ttest_core<-t.test(s20$Core~s20$IfFirstEpisode)
ttest_dmpfc<-t.test(s20$DMPFC~s20$IfFirstEpisode)

# Compare medicated vs non nedicated
factor_count<-by(merged_agesexreg$Site, merged_agesexreg$OnMedication, count)
sites_keep1<-factor_count[['-1']][factor_count[["-1"]]$freq>=20, 'x']
sites_keep2<-factor_count[['1']][factor_count[["1"]]$freq>=20, 'x']
sites_keep<-intersect(sites_keep1,sites_keep2)

# Compare medicated and unmedicated in S20
s20<-merged_agesexreg[merged_agesexreg$Site=='S20',]
ttest_core<-t.test(s20$Core~s20$OnMedication)
ttest_dmpfc<-t.test(s20$DMPFC~s20$OnMedication)

# Compare motion between HC and MDD
wilcox.test(merged_agesexreg$fdvols~merged_agesexreg$Group, na.rm=TRUE) 
