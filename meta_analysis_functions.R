###### DEFINE FUNCTIONS FOR RUNNING THE META-ANALYSES

# Function to create group mean and std dataframe for a network for a specific variable of interest
create_ES_df_binintvar<-function(df, sites, interestvar, nw){
  
  vars=c('g1_n', 'g1_mean', 'g1_std', 'g2_n', 'g2_mean', 'g2_std')
  nw_df=data.frame(matrix(nrow = length(sites), ncol = length(vars)+1))
  colnames(nw_df) <- c('Site', vars)
  nw_df$Site=sites
  interest1=unique(df[,interestvar])[1]
  interest2=unique(df[,interestvar])[2]
  
  for (site in sites)
  {
    nw_df[nw_df$Site==site,vars[1]]<-length(df[df[,interestvar]==interest1 & df$Site==site,nw])
    nw_df[nw_df$Site==site,vars[2]]<-mean(df[df[,interestvar]==interest1& df$Site==site,nw])
    nw_df[nw_df$Site==site,vars[3]]<-sd(df[df[,interestvar]==interest1 & df$Site==site,nw])
    nw_df[nw_df$Site==site,vars[4]]<-length(df[df[,interestvar]==interest2 & df$Site==site,nw])
    nw_df[nw_df$Site==site,vars[5]]<-mean(df[df[,interestvar]==interest2 & df$Site==site,nw])
    nw_df[nw_df$Site==site,vars[6]]<-sd(df[df[,interestvar]==interest2 & df$Site==site,nw])
  }
  return(nw_df)
}

# Function to create correlation dataframe for a network for a specific variable of interest
create_ES_df_corr<-function(df, sites, nw, interestvar){
  
  vars=c('n', 'corr')
  nw_df=data.frame(matrix(nrow = length(sites), ncol = length(vars)+1))
  colnames(nw_df) <- c('Site', vars)
  nw_df$Site=sites
  
  for (site in sites)
  {
    nw_df[nw_df$Site==site,vars[1]]<-length(df[df$Site==site,nw])
    nw_df[nw_df$Site==site,vars[2]]<-cor(df[df$Site==site,nw], df[df$Site==site,interestvar],method = 'spearman', use='complete.obs')
  }
  
  return(nw_df)
}

# Function to run between groups meta-analysis for a network given the dataframes returned from the functions above
run_btw_group_meta <- function(dataframe) {
  m.hksj <- metacont(
    g1_n,
    g1_mean,
    g1_std,
    g2_n,
    g2_mean,
    g2_std,
    studlab = dataframe$Site,
    data = dataframe,
    comb.fixed = FALSE,
    comb.random = TRUE,
    method.tau = "REML",
    #choosing REML according to https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4950030/
    hakn = TRUE,
    prediction = TRUE,
    sm = "SMD", 
    control=list(maxiter=10000, stepadj=0.5)
  )
  return(m.hksj)
}

# Function to run between correlation meta-analysis for a network given the dataframes returned from the functions above
run_corr_meta <- function(dataframe){
  m.cor<-metacor(cor = corr,
                 n = n,
                 studlab = dataframe$Site,
                 data=dataframe,
                 method.tau = "REML",
                 hakn = TRUE,
                 comb.fixed = FALSE,
                 comb.random = TRUE,
                 prediction = TRUE,
                 control=list(maxiter=10000, stepadj=0.5),
                 sm = "ZCOR")
  return(m.cor)
}

# Function to generate forest plots and add title
drawforplot<-function(metaandata, title)
{
  postscript(file=paste('meta/plots/', title,'.ps',sep=''))
  forest(metaandata, main=title,layout = "RevMan5",
         digits.sd = 2)
  grid.text(nw, .5, .80, gp=gpar(cex=2))
  dev.off()
}

# Function to export meta-analysis details for groups
export_meta_group <-function(metaandata, nw){
  df=data.frame(network = nw)
  df$ngroup1<-sum(metaandata$n.e)
  df$ngroup2<-sum(metaandata$n.c)
  df$effect<-metaandata$TE.random
  df$effect_l<-metaandata$lower.random
  df$effect_u<-metaandata$upper.random
  df$effect_p<-metaandata$pval.random
  df$predict_l<-metaandata$lower.predict
  df$predict_u<-metaandata$upper.predict
  df$tau2<-metaandata$tau2
  df$q<-metaandata$Q
  df$qp<-metaandata$pval.Q
  df$I2<-metaandata$I2
  return(df)
}

# Function to export meta-analysis details for correlations
export_meta_corr <-function(metaandata, nw){
  df=data.frame(network = nw)
  df$ngroup<-sum(metaandata$n)
  df$effect<-metaandata$TE.random
  df$effect_l<-metaandata$lower.random
  df$effect_u<-metaandata$upper.random
  df$effect_p<-metaandata$pval.random
  df$predict_l<-metaandata$lower.predict
  df$predict_u<-metaandata$upper.predict
  df$tau2<-metaandata$tau2
  df$q<-metaandata$Q
  df$qp<-metaandata$pval.Q
  df$I2<-metaandata$I2
  return(df)
}

