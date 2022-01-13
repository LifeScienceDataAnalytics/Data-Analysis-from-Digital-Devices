##change December 2021
#load("~/Documents/Documents_IT/paper/Altoida_VAMBN_paper/classifiersgl_upd.RData")
load("~/Documents/Documents_IT/paper/Altoida_VAMBN_paper/auc_upd.RData")
##change December 2021
load("~/Documents/Documents_IT/paper/Altoida_VAMBN_paper/classifiersglInd_upd.RData")
rm(list=setdiff(ls(),  c("rocRESMMSE", "rocRESCog", "rocDig", "rocCogMMSE", "rocDigMMSE","rocCogMMSEAUC","rocDigMMSEAUC", "rocDigAUC")))
#####
boxplotDf = cbind.data.frame(unlist(rocRESMMSE[[1]][1:48]),unlist(rocRESCog[[1]][1:48]), 
                             unlist(rocCogMMSEAUC), unlist(rocDigAUC), unlist(rocDigMMSEAUC))
fill <- "#00BFC4"
line <- "#1F3552"

names(boxplotDf) = c("MMSE", "Digital cognitive domains", "MMSE and Digital cognitive domains", "Digital tasks", "MMSE and Digital tasks")
library(reshape)
boxplotDf = reshape2::melt(boxplotDf)
plt <- ggplot(data = boxplotDf, aes(x = variable, y = value))
png("sgl_R/plots_revised/AUC_boxplot.png", width = 150, height = 150, units='mm',res = 400) 
plt + geom_boxplot(fill = fill, colour = line) + labs(x = "Features", y = "AUC") + ylim(0.5,1) + theme(axis.text = element_text(size = 10)) +theme(text = element_text(size = 14))+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) + theme_classic(base_size = 14)
dev.off()

#####
##f1 boxplot
boxplotDfF1 = cbind.data.frame(unlist(rocRESMMSE[[2]][1:48]),unlist(rocRESCog[[2]][1:48]), 
                             unlist(rocCogMMSE[4]), unlist(rocDig[4]), unlist(rocDigMMSE[4]))
fill <- "#00BFC4"
line <- "#1F3552"

names(boxplotDfF1) = c("MMSE", "Digital cognitive domains", "MMSE and Digital cognitive domains", "Digital tasks", "MMSE and Digital tasks")
library(reshape)
boxplotDfF1 = reshape2::melt(boxplotDfF1)
plt <- ggplot(data = boxplotDfF1, aes(x = variable, y = value))
png("sgl_R/plots_revised/F1_boxplot.png", width = 150, height = 150, units='mm',res = 400) 
plt + geom_boxplot(fill = fill, colour = line) + labs(x = "Features", y = "F1 score") + ylim(0.5,1) + theme(axis.text = element_text(size = 14)) +theme(text = element_text(size = 18))+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) + theme_classic(base_size = 14)
dev.off()

##checking ROC
mmseWithCog = wilcox.test(unlist(rocRESMMSE[[1]][1:48]), unlist(rocCogMMSE[[1]]), paired = TRUE, alternative = "two.sided")
mmseWithDig = wilcox.test(unlist(rocRESMMSE[[1]][1:48]), unlist(rocDigMMSEAUC), paired = TRUE, alternative = "two.sided")
mmseAndDig = wilcox.test(unlist(rocRESMMSE[[1]][1:48]), unlist(rocDigAUC), paired = TRUE, alternative = "two.sided")

##checking F1
mmseWithCogF1 = wilcox.test(unlist(rocRESMMSE[[2]][1:48]), unlist(rocCogMMSE[[4]]), paired = TRUE, alternative = "two.sided")
mmseWithDigF1 = wilcox.test(unlist(rocRESMMSE[[2]][1:48]), unlist(rocDigMMSE[[4]]), paired = TRUE, alternative = "two.sided")
mmseAndDigF1 = wilcox.test(unlist(rocRESMMSE[[2]][1:48]), unlist(rocDig[[4]]), paired = TRUE, alternative = "two.sided")


##FOR POSTER
# boxplotDf = cbind.data.frame(unlist(rocRESMMSE), unlist(rocDig))
# fill <- "#00BFC4"
# line <- "#1F3552"
# 
# names(boxplotDf) = c("MMSE",  "Digital tasks")
# library(reshape)
# boxplotDf = reshape2::melt(boxplotDf)
# plt <- ggplot(data = boxplotDf, aes(x = variable, y = value))
# png("sgl_R/plots/AUC_boxplot_POSTER.png", width = 100, height = 100, units='mm',res = 300) 
# plt + geom_boxplot(fill = fill, colour = line) + labs(x = "Features", y = "AUC") + ylim(0.5,1) + theme(axis.text = element_text(size = 14)) +theme(text = element_text(size = 18))+
#   scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) + theme_classic(base_size = 15)
# dev.off()
