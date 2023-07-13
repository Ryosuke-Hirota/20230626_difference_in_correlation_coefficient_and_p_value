# This script is to compare correlation coefficient and p.value between analysis using RPKM or TPM
# 2023/06/26 made

# activate package
library(ggplot2)

# make new directory
setwd("C:/Rdata")
dir.create("20230626_difference_in_correlation_coefficient_and_p_value_between_analysis_using_RPKM_data_and_TPM_data")

# import summary of correlation analysis using TPM data
# this summary is located at "https://github.com/Ryosuke-Hirota/20230627_ROC_curve_of_pvalue_after_cutoff_using_TPM_data"
setwd("C:/Rdata/20230517_CCLE_correlation_between_residual_and_RBPDB_RBP_RSEM_TPM")
TPM.sm <-read.table("summary_of_correlation_between_residual_and_RBPDB_RBP_RSEM_TPM.txt",sep="\t",header=T,stringsAsFactors = F)

# remove unnecessary rows and arrange dataframe
TPM.sm <-subset(TPM.sm,!is.na(TPM.sm[,5]))
TPM.sm[,1] <-paste0(TPM.sm[,1],"_",TPM.sm[,2],"_",TPM.sm[,3],"_vs_",TPM.sm[,4])
TPM.sm <-TPM.sm[,c(1,5,6,7)]
colnames(TPM.sm)[1:3] <-c("combination","TPM.rho","TPM.p.value")

# import summary of correlation analysis using RPKM data
# this summary is located at "https://github.com/Ryosuke-Hirota/20230621_CCLE_correlation_between_residual_and_RBPDB_RBP_RPKM"
setwd("C:/Rdata/20230621_CCLE_correlation_between_residual_and_RBPDB_RBP_RPKM")
RPKM.sm <-read.table("summary_of_correlation_between_residual_and_RBPDB_RBP_RPKM.txt",sep="\t",header=T,stringsAsFactors = F)

#  remove unnecessary rows and arrange dataframe
RPKM.sm <-subset(RPKM.sm,!is.na(RPKM.sm[,5]))
RPKM.sm[,1] <-paste0(RPKM.sm[,1],"_",RPKM.sm[,2],"_",RPKM.sm[,3],"_vs_",RPKM.sm[,4])
RPKM.sm <-RPKM.sm[,c(1,5,6,7)]
colnames(RPKM.sm)[1:3] <-c("combination","RPKM.rho","RPKM.p.value")

# merge summaries of TPM and RPKM
TPM.RPKM <-merge(TPM.sm,RPKM.sm,by=c("combination","number_of_cell"))
TPM.RPKM <-TPM.RPKM[,c(1,3:6,2)]

# -log10(p.value)
TPM.RPKM[,3] <-log10(TPM.RPKM[,3])*-1
TPM.RPKM[,5] <-log10(TPM.RPKM[,5])*-1

setwd("C:/Rdata/20230626_difference_in_correlation_coefficient_and_p_value_between_analysis_using_RPKM_data_and_TPM_data")

# draw density plot about difference of correlation coefficient between analysis using TPM data and RPKM data
p <-ggplot(data=TPM.RPKM,aes(x=TPM.rho,y=RPKM.rho))+
    geom_bin2d(bins=200)+
    geom_abline(intercept = 0, slope = 1,col="black",linewidth=1)+
    scale_fill_gradient(low='#c80000', high='#FFFF00')+
    coord_cartesian(xlim = c(-1, 1), ylim = c(-1, 1))+
    theme_classic()

ggsave("density_plot_of_correlation_coefficient_between_analysis_using_TPM_data_and_analysis_using_RPKM_data.pdf",plot = p)

# set cutoff cell lines † 50, draw plot
TPM.RPKM50 <-TPM.RPKM[TPM.RPKM[,6]>=50,]
p50 <-ggplot(data=TPM.RPKM50,aes(x=TPM.rho,y=RPKM.rho))+
  geom_bin2d(bins=200)+
  geom_abline(intercept = 0, slope = 1,col="black",linewidth=1)+
  scale_fill_gradient(low='#c80000', high='#FFFF00')+
  coord_cartesian(xlim = c(-1, 1), ylim = c(-1, 1))+
  theme_classic()

ggsave("density_plot_of_correlation_coefficient_between_analysis_using_TPM_data_and_analysis_using_RPKM_data_cutoff_50.pdf",plot = p50)

# draw density plot about difference of p.value between analysis using TPM data and RPKM data
pp <-ggplot(data=TPM.RPKM,aes(x=TPM.p.value,y=RPKM.p.value))+
  geom_bin2d(bins=200)+
  geom_abline(intercept = 0, slope = 1,col="black",linewidth=1)+
  scale_fill_gradient(low='#c80000', high='#FFFF00')+
  xlab("-log10(TPM.p.value)")+
  ylab("-log10(RPKM.p.value)")+
  theme_classic()

ggsave("density_plot_of_p_value_between_analysis_using_TPM_data_and_analysis_using_RPKM_data.pdf",plot = pp)

# set cutoff cell lines † 50, draw plot
pp50 <-ggplot(data=TPM.RPKM50,aes(x=TPM.p.value,y=RPKM.p.value))+
  geom_bin2d(bins=200)+
  geom_abline(intercept = 0, slope = 1,col="black",linewidth=1)+
  scale_fill_gradient(low='#c80000', high='#FFFF00')+
  xlab("-log10(TPM.p.value)")+
  ylab("-log10(RPKM.p.value)")+
  theme_classic()

ggsave("density_plot_of_p_value_between_analysis_using_TPM_data_and_analysis_using_RPKM_data_cutoff_50.pdf",plot = pp50)
