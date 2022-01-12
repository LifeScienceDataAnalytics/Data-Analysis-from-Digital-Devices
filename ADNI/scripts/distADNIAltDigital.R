rm(list=ls())
library(tidyverse)
library(plyr)
library(grid)
library(gridExtra)
library(Matching)
library(LaplacesDemon)
compare_dists_dig<-function(data,typecol,folder){
  #data = all_nm
  typecol = "type"
  gnames<-unique(data[,typecol])
  vnames<-colnames(data[,!grepl('type',colnames(data))])
  
  lv<-levels(factor(data[,typecol], levels = c("altoida", "adni")))
  group.1<-subset(data,eval(parse(text = typecol))==lv[1])
  group.2<-subset(data,eval(parse(text = typecol))==lv[2])
  
  print(paste0('1: ',lv[1],' 2: ',lv[2]))
  
  pvals<-NULL
  for (col in vnames){
    #print(col)
    #col = "BIT_AR_accVarianceX"
    dat<-data[,grepl(paste0(col,'$|type'),colnames(data))]
    if (is.numeric(data[,col])){
      dv1<-group.1[,col]
      dv2<-group.2[,col]
      set.seed(123)
      ##decoded real and decoded virtual
      #px =  dnorm(group.1[,col])
      #py =  dnorm(group.2[,col])
      #kld = KLD(px, py)
      #klValue = round(kld$sum.KLD.py.px, 4)
      stats<-matrix(c(lv[1],round(mean(dv1,na.rm=T),2),round(sd(dv1,na.rm=T),2),round(quantile(dv1,na.rm=T)['25%'],2),round(median(dv1,na.rm=T),2),round(quantile(dv1,na.rm=T)['75%'],2),
                      lv[2],round(mean(dv2,na.rm=T),2),round(sd(dv2,na.rm=T),2),round(quantile(dv2,na.rm=T)['25%'],2),round(median(dv2,na.rm=T),2),round(quantile(dv2,na.rm=T)['75%'],2)
      ),ncol=6,byrow=TRUE)
      colnames(stats) <- c("Type","Mean","SD","25%","Median","75%")
      plot1<-ggplot(dat, aes(x=type,y=eval(parse(text = col)),fill=type))+
        geom_violin(position=position_dodge(1)) + 
        geom_boxplot(position=position_dodge(1),width=0.05) +# ggtitle(paste('kldivergence:', klValue))+
        xlab('type')+ylab(col)+theme(legend.position="None",axis.title.x=element_blank())
      gr<-grey.colors(3)
      tt3 <- ttheme_minimal(
        core=list(bg_params = list(fill = c(gr[2],gr[3],gr[2]), col=NA),
                  fg_params=list(fontface=3)))
      plot<-grid.arrange(plot1,tableGrob(stats,theme=tt3),ncol=1, 
                         widths=unit(15,"cm"),
                         heights=unit(c(8,3), c("cm", "cm")))
    }else{
      gd <- dat %>% group_by(type) %>% count
      plot<-ggplot(gd, aes(x=eval(parse(text = col)),y=freq,fill=type))+geom_bar(position="dodge", stat="identity")+scale_fill_brewer(palette="Dark2")+xlab(col)+ylab('Count')#+ggtitle(paste('Permutation Pval:',pval,'Alpha:',signif(alpha,digits = 3),'Sig:',sig))
    }
    ggsave(paste0(folder,col,'.png'), plot, device = "png", width=7.26, height = 4.35)
    ggsave(paste0(folder,col,'.eps'), plot, device = "eps", width=7.26, height = 4.35)
  }
}

reconRP_Altoida <- read_csv("~/Documents/Documents_IT/paper/Altoida_VAMBN_paper/data/HI-VAE/reconRP.csv")

reconRP_ADNI_Dig <- read_csv("data/HI-VAE/reconRP_all_Dig.csv")
reconRP_ADNI_Dig$SUBJID = NULL
reconRP_Altoida_Dig = reconRP_Altoida[,colnames(reconRP_ADNI_Dig)]
reconRP_Altoida_Dig$type<-'altoida'
reconRP_ADNI_Dig$type<-'adni'

reconRP_Altoida_Dig = as.data.frame(reconRP_Altoida_Dig)
reconRP_ADNI_Dig = as.data.frame(reconRP_ADNI_Dig)

all<-rbind(reconRP_Altoida_Dig, reconRP_ADNI_Dig)

compare_dists_dig(all,'type','data/data_out/dist_violin_joint_dig/')

##missing data 
reconRP_Altoida_nm <- read_csv("~/Documents/Documents_IT/paper/Altoida_VAMBN_paper/data/HI-VAE/reconRP_nm.csv")
reconRP_ADNI_nm_Dig <- read_csv("data/HI-VAE/reconRP_all_Dig_nm.csv")
reconRP_ADNI_nm_Dig$SUBJID = NULL
reconRP_Altoida_nm_Dig = reconRP_Altoida_nm[,colnames(reconRP_ADNI_nm_Dig)]
reconRP_Altoida_nm_Dig$type<-'altoida'
reconRP_ADNI_nm_Dig$type<-'adni'

reconRP_Altoida_nm_Dig = as.data.frame(reconRP_Altoida_nm_Dig)
reconRP_ADNI_nm_Dig = as.data.frame(reconRP_ADNI_nm_Dig)


all_nm<-rbind(reconRP_Altoida_nm_Dig, reconRP_ADNI_nm_Dig)

compare_dists_dig(all_nm,'type','data/data_out/dist_violin_joint_dig/')


######load the last rows in the data#
reconRP_ADNI_Dig_last178 <- read_csv("data/HI-VAE/reconRP_all_Dig_last178.csv")
reconRP_ADNI_Dig_last178$SUBJID = NULL
reconRP_Altoida_Dig = reconRP_Altoida[,colnames(reconRP_ADNI_Dig_last178)]
reconRP_Altoida_Dig$type<-'altoida'
reconRP_ADNI_Dig_last178$type<-'adni'


reconRP_Altoida_Dig = as.data.frame(reconRP_Altoida_Dig)
reconRP_ADNI_Dig_last178 = as.data.frame(reconRP_ADNI_Dig_last178)

all_last_178<-rbind(reconRP_ADNI_Dig_last178, reconRP_Altoida_Dig)

compare_dists_dig(all_last_178,'type','data/data_out/dist_violin_joint_dig_last178/')



######load the last rows in the data#
reconRP_ADNI_Dig_last_nm_178 <- read_csv("data/HI-VAE/reconRP_all_Dig_las_nm_178.csv")
reconRP_ADNI_Dig_last_nm_178$SUBJID = NULL
reconRP_Altoida_nm_Dig$type<-'altoida'
reconRP_ADNI_Dig_last_nm_178$type<-'adni'


reconRP_Altoida_nm_Dig = as.data.frame(reconRP_Altoida_nm_Dig)
reconRP_ADNI_Dig_last_nm_178 = as.data.frame(reconRP_ADNI_Dig_last_nm_178)

all_nm_last_178<-rbind(reconRP_ADNI_Dig_last_nm_178, reconRP_Altoida_nm_Dig)

compare_dists_dig(all_21,'type','data/data_out/dist_violin_joint_dig_nm_last178/')


