setwd("/Users/ltozzi/Desktop/Oak/ltozzi/hcpdes/dmn_subsystems_project/")

library(meta)
library(tidyverse)
library(grid)
library(plyr)

source('scripts/meta_analysis_functions.R')
dir.create('meta')
dir.create('meta/tables')
dir.create('meta/plots')

###### BETWEEN NETWORK CONNECTIVITY

# Create merged dataframe of clinical and nw_pairs
clinical_data=read_csv('clinical/metamdd_CG_LT.csv') # clinical data only for people included in Yan (2019)
nw_data=read_csv('dmn_subsystems/between_DMNsubsys_fc_metamdd.csv') # average network data
merged=merge(clinical_data, nw_data)

# Labeling networks and sites
nw_pairs=c('Core_DMPFC', 'Core_MTL', 'DMPFC_MTL')
sites=sort(unique(merged$Site))

################ ANALYSIS

#### AGE AND SEX REGRESSION

# Regress age and sex from data at each site for each network
merged_agesexreg<-merged
for (site in sites){
  sitedata=merged[merged$Site==site,]
  for (nw in nw_pairs){
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
for (nw in nw_pairs){
  nwdf<-create_ES_df_binintvar(df=merged_agesexreg, sites = sites_keep, interestvar = 'Group', nw=nw)
  MA_list_hcmdd_agesexreg[[nw]]<-run_btw_group_meta(nwdf)
}

# Remove sites with less than 20 values of HAMD
factor_count<-by(merged_agesexreg$Site, !is.na(merged_agesexreg$HAMD), count)
sites_keep<-factor_count$'TRUE'[factor_count[["TRUE"]]$freq>=20,'x']

# Run meta-analysis for all networks for correlation with HAMD
MA_list_hamdcorr_agesexreg = list()
for (nw in nw_pairs){
  nwdf<-create_ES_df_corr(df=merged_agesexreg, sites = sites_keep, interestvar = 'HAMD', nw=nw)
  MA_list_hamdcorr_agesexreg[[nw]]<-run_corr_meta(nwdf)
}

# Remove sites with less than 20 values of HAMA
factor_count<-by(merged_agesexreg$Site, !is.na(merged_agesexreg$HAMA), count)
sites_keep<-factor_count$'TRUE'[factor_count[["TRUE"]]$freq>=20,'x']

# Run meta-analysis for all networks for correlation with HAMA
MA_list_hamacorr_agesexreg = list()
for (nw in nw_pairs){
  nwdf<-create_ES_df_corr(df=merged_agesexreg, sites = sites_keep, interestvar = 'HAMA', nw=nw)
  MA_list_hamacorr_agesexreg[[nw]]<-run_corr_meta(nwdf)
}

# Remove sites with less than 20 values of fdvals
factor_count<-by(merged_agesexreg$Site, !is.na(merged_agesexreg$fdvols), count)
sites_keep<-factor_count$'TRUE'[factor_count[["TRUE"]]$freq>=20,'x']

# Run meta-analysis for all networks for correlation with fdvals
MA_list_fdvols_agesexreg = list()
for (nw in nw_pairs){
  nwdf<-create_ES_df_corr(df=merged_agesexreg, sites = sites_keep, interestvar = 'fdvols', nw=nw)
  MA_list_fdvols_agesexreg[[nw]]<-run_corr_meta(nwdf)
}

###### EXPORT RESULTS

# Draw and save forest plots
for (nw in nw_pairs){
  drawforplot(MA_list_hcmdd_agesexreg[[nw]], paste(nw, '_hcmdd_agesexreg_dmnsubsys', sep=''))
  drawforplot(MA_list_hamdcorr_agesexreg[[nw]], paste(nw, '_hamdcorr_agesexreg_dmnsubsys', sep=''))
  drawforplot(MA_list_hamacorr_agesexreg[[nw]], paste(nw, '_hamacorr_agesexreg_dmnsubsys', sep='')) 
  drawforplot(MA_list_fdvols_agesexreg[[nw]], paste(nw, '_fdvols_agesexreg_dmnsubsys', sep=''))  
}

# Export group meta-analyses
meta_export_hcmdd_agesexreg <- data.frame(matrix(ncol = 13, nrow = 0))
colnames(meta_export_hcmdd_agesexreg)<-colnames(export_meta_group(meta_export_hcmdd_agesexreg[[nw]], nw))

for (nw in nw_pairs){
  meta_export_hcmdd_agesexreg<-rbind(meta_export_hcmdd_agesexreg, export_meta_group(MA_list_hcmdd_agesexreg[[nw]], nw))
}

# Compute FDR correction
meta_export_hcmdd_agesexreg$effect_pfdr<-p.adjust(meta_export_hcmdd_agesexreg$effect_p, method = "fdr")
meta_export_hcmdd_agesexreg$qpfdr<-p.adjust(meta_export_hcmdd_agesexreg$qp, method = "fdr")

# Export correlation meta-analyses
meta_export_hamd_agesexreg <- data.frame(matrix(ncol = 12, nrow = 0))
colnames(meta_export_hamd_agesexreg)<-colnames(export_meta_corr(MA_list_hamdcorr_agesexreg[[nw]], nw))
meta_export_hama_agesexreg <- data.frame(matrix(ncol = 12, nrow = 0))
colnames(meta_export_hama_agesexreg)<-colnames(export_meta_corr(MA_list_hamacorr_agesexreg[[nw]], nw))
meta_export_fdvols_agesexreg <- data.frame(matrix(ncol = 12, nrow = 0))
colnames(meta_export_fdvols_agesexreg)<-colnames(export_meta_corr(MA_list_fdvols_agesexreg[[nw]], nw))

for (nw in nw_pairs){
  meta_export_hamd_agesexreg<-rbind(meta_export_hamd_agesexreg, export_meta_corr(MA_list_hamdcorr_agesexreg[[nw]], nw))
  meta_export_hama_agesexreg<-rbind(meta_export_hama_agesexreg, export_meta_corr(MA_list_hamacorr_agesexreg[[nw]], nw))
  meta_export_fdvols_agesexreg<-rbind(meta_export_fdvols_agesexreg, export_meta_corr(MA_list_fdvols_agesexreg[[nw]], nw))
}

# Compute FDR correction
meta_export_hamd_agesexreg$effect_pfdr<-p.adjust(meta_export_hamd_agesexreg$effect_p, method = "fdr")
meta_export_hama_agesexreg$effect_pfdr<-p.adjust(meta_export_hama_agesexreg$effect_p, method = "fdr")
meta_export_fdvols_agesexreg$effect_pfdr<-p.adjust(meta_export_fdvols_agesexreg$effect_p, method = "fdr")

meta_export_hamd_agesexreg$qpfdr<-p.adjust(meta_export_hamd_agesexreg$qp, method = "fdr")
meta_export_hama_agesexreg$qpfdr<-p.adjust(meta_export_hama_agesexreg$qp, method = "fdr")
meta_export_fdvols_agesexreg$qpfdr<-p.adjust(meta_export_fdvols_agesexreg$qp, method = "fdr")

# Save CSVs
write_csv(x = meta_export_hcmdd_agesexreg, path ='meta/tables/btwFC_summary_hcmdd_agesexreg_dmnsubsys.csv')
write_csv(x = meta_export_hamd_agesexreg, path ='meta/tables/btwFC_summary_hamd_agesexreg_dmnsubsys.csv')
write_csv(x = meta_export_hama_agesexreg, path ='meta/tables/btwFC_summary_hama_agesexreg_dmnsubsys.csv')
write_csv(x = meta_export_fdvols_agesexreg, path ='meta/tables/btwFC_summary_fdvols_agesexreg_dmnsubsys.csv')

# Compare recurrent and first episode in S20
s20<-merged_agesexreg[merged_agesexreg$Site=='S20',]
ttest_coredmpfc<-t.test(s20$Core_DMPFC~s20$IfFirstEpisode)


