rm(list=ls())
  

# load libraries
library(corrplot)
library(pROC)
library(tidyverse)
library(mixOmics)
library(NMF)
library(RColorBrewer)
library(cowplot)
library(scales)
library(knitr)
library(UpSetR)
library(grid)
library(xlsx)
library(reshape2)
library(Hmisc)
library(pheatmap)
library(igraph)

#create wd
setwd('R:\\Para-CIH\\CIH Group member folders\\Simon\\Results\\Nasal allergen challenge mixomics\\Output')


dir.create('Two datasets')
setwd('Two datasets/')

#make session info
sink('session info.txt')
sessionInfo()
sink()

## load data
data <- read.xlsx('../../Input/Raw_data_20221004.xlsx', 1)
groups <- read.xlsx('../../Input/Raw_data_20221004.xlsx', 2)


#make in correct format
rownames(data) <- data$Clustername #add rownmaes
colnames(data) <- gsub('X', '', colnames(data))

data_t <- t(data[,3:44]) #transpose data

#get the metadata sorted
info <- data.frame('ID' = rownames(data_t))
info$Time <- 't0'
info$Time[grepl('After', info$ID)] <- 't1'
info$Patient <- do.call(c, lapply(strsplit(info$ID, '_'), function(i)i[1]))
info2 <- merge(info, groups, by.x='Patient', by.y='Donor.Id')
info2$Group2 <- paste(info2$Group, info2$Time, sep = '_')

#check the order is the same
identical(info2$ID,rownames(data_t))

X <- list(CyTOF=data_t[,39:74], scRNA=data_t[,1:38]) #create the data list
lapply(X, dim) #check the number of people and measurements in the normalized dataset


#make matrix for the time corrected samples
Cov <- data.frame(sample = info2$Patient, time = info2$Time)
Cov[,1] <- as.character(Cov[,1])
Cov[,1] <- as.numeric(factor(Cov[,1]))


#################################################################################################
#do for HDs and patients seperately
dir.create('withinvariation/Sep groups')
setwd('withinvariation/Sep groups')

X.ctrl <- lapply(X, function(i) i[info2$Group == 'Ctrl',])
lapply(X.ctrl, dim)
lapply(X.ctrl, rownames)

#make new metadata
info2.ctrl <- info2[info2$Group == 'Ctrl',]
Cov.ctrl <- data.frame(sample = info2.ctrl$Patient, time = info2.ctrl$Time)
Cov.ctrl[,1] <- as.character(Cov.ctrl[,1])
Cov.ctrl[,1] <- as.numeric(factor(Cov.ctrl[,1]))
sample.ctrl = Cov.ctrl[,1]


#withinvariation for controls
A.ctrl = lapply(X.ctrl, function(i) suppressMessages(withinVariation(X = i, design = Cov.ctrl))) #normalize per individual
lapply(A.ctrl, dim) #check the number of people and measurements in the normalized dataset
Y.ctrl = info2.ctrl$Group2


#tune the model with 3 components
design <- matrix(1, nrow = length(X), ncol = length(X)) #create a design matrix with max correlation between datasets
keepX = list(CyTOF = 2:6, scRNA = 2:6)
ncomp=2
tune.Actrl <- tune.block.splsda(A.ctrl, Y.ctrl, test.keepX = keepX, ncomp = ncomp, design=design, validation = 'loo')
plot(tune.Actrl)
dev.off()




#tune manually with Leave person using weighted vote
  predict2 <- matrix(nrow=0,ncol=10)
  colnames(predict2) <- c( 'CYTOF1', 'CYTOF2', 'RNA1', 'RNA2', 'Accuracy0',
                           'Accuracy1', 'Accuracy2', 'Accuracy3', 'Accuracy4',
                           'Accuracy5')
  
  for(k in 2:6){
    for(l in 2:6){
      for(m in 2:6){
        for(n in 2:6){
              
              list.keepX <- list(CyTOF = c(k,l), scRNA = c(m,n))
              predict <- matrix(ncol=12,nrow=0)
              for(rm in 1 : length(unique(sample.ctrl))){
                sequence <- which(rm == sample.ctrl)
                trainX <- lapply(A.ctrl, function(i) i[-sequence, ])
                trainY <- Y.ctrl[-sequence]
                xtest <- lapply(A.ctrl, function(i) i[sequence, ])
                
                diablo = block.splsda(X = trainX, Y = trainY, ncomp = 2, keepX = list.keepX, design = design)
                test = predict(diablo, xtest, dist  = "centroids.dist")
                predict <- rbind(predict, c(
                                 test$AveragedPredict.class$max.dist[,2],
                                 test$AveragedPredict.class$max.dist[,2],
                                 test$WeightedVote$centroids.dist[,2],
                                 test$AveragedPredict.class$max.dist[,1],
                                 test$AveragedPredict.class$max.dist[,1],
                                 test$WeightedVote$centroids.dist[,1]
                                 ))
              }
              accuracy0 <- (sum(predict[,1] == 'Ctrl_t0')+sum(predict[,2] == 'Ctrl_t1'))/length(sample.ctrl)
              accuracy1 <- (sum(predict[,3] == 'Ctrl_t0')+sum(predict[,4] == 'Ctrl_t1'))/length(sample.ctrl)
              accuracy2 <- (sum(predict[,5] == 'Ctrl_t0')+sum(predict[,6] == 'Ctrl_t1'))/length(sample.ctrl)
              accuracy3 <- (sum(predict[,7] == 'Ctrl_t0')+sum(predict[,8] == 'Ctrl_t1'))/length(sample.ctrl)
              accuracy4 <- (sum(predict[,9] == 'Ctrl_t0')+sum(predict[,10] == 'Ctrl_t1'))/length(sample.ctrl)
              accuracy5 <- (sum(predict[,11] == 'Ctrl_t0')+sum(predict[,12] == 'Ctrl_t1'))/length(sample.ctrl)
              
              
              result <- c(k,l,m,n, accuracy0,accuracy1, accuracy2,  accuracy3,accuracy4, accuracy5)
              predict2 <- rbind(predict2, result)
            }
          }
        }
      }

predict.ctrl <- predict2

#plot tuning results
pdf('ctrls_tuning_A.pdf', width = 6, height = 4)
plot(tune.Actrl)
dev.off()

#make model with best tuning parameters
list.keepActrl = tune.Actrl$choice.keepX

#maximum numbers of features
list.keepActrl$CyTOF = c(6,3)
list.keepActrl$scRNA = c(6,3)

#minimum numbers of features
list.keepActrl$CyTOF = c(5,2)
list.keepActrl$scRNA = c(4,5)


final.diablo.model.A.ctrl = block.splsda(X = A.ctrl, Y = Y.ctrl, ncomp = ncomp, 
                                    keepX = list.keepActrl, design = design)


perf.ctrl <- perf(final.diablo.model.A.ctrl, validation = 'loo')
plot(perf.ctrl)
ctr_cons <- perf.ctrl$features


#make the plots to look at the models
pdf('ctrls_arrowplot comp1_2_A.pdf', width = 6, height = 6)
plotArrow(final.diablo.model.A.ctrl, comp = c(1,2), ind.names = FALSE, legend = TRUE, 
          title = 'DIABLO')
dev.off()


pdf('ctrls_loadings_A.pdf', width = 8, height = 5)
plotLoadings(final.diablo.model.A.ctrl, comp = 1, contrib = 'max', method = 'median')
plotLoadings(final.diablo.model.A.ctrl, comp = 2, contrib = 'max', method = 'median')
dev.off()



dat <- as.data.frame(Reduce("+", final.diablo.model.A.ctrl$variates)/length(X)) %>% 
  mutate(time = Y.ctrl, sample = factor(sample.ctrl))
ggplot(dat, aes(x = `comp 1`, y = `comp 2`, group = time, color = time)) +
  geom_point(size = 4) +
  stat_ellipse() +
  geom_line(aes(group = sample), color = "gray")+
  xlab("Consensus Component 1") +
  ylab("Consensus Component 2") +
  ggtitle("Non-allergic") +
  theme(legend.position = c(0.1, 0.9)) + 
  theme_bw() + 
  theme(aspect.ratio = 1) +
  scale_color_manual(values=c("#388ECC", "darkgrey"))
ggsave('Non_allergic_indiv.pdf', width = 10, height = 10, units = 'cm')


plotIndiv(final.diablo.model.A.ctrl)
plotDiablo(final.diablo.model.A.ctrl)


#################################################################################################
#do only for AR
X.AR <- lapply(X, function(i) i[info2$Group == 'AR',])
lapply(X.AR, dim)
lapply(X.AR, rownames)

#make new metadata
info2.AR <- info2[info2$Group == 'AR',]
Cov.AR <- data.frame(sample = info2.AR$Patient, time = info2.AR$Time)
Cov.AR[,1] <- as.character(Cov.AR[,1])
Cov.AR[,1] <- as.numeric(factor(Cov.AR[,1]))
sample.AR = Cov.AR[,1]


#withinvariation for controls
A.AR = lapply(X.AR, function(i) suppressMessages(withinVariation(X = i, design = Cov.AR))) #normalize per individual
lapply(A.AR, dim) #check the number of people and measurements in the normalized dataset
Y.AR = info2.AR$Group2


#tune manually with Leave person using weighted vote
predict2 <- matrix(nrow=0,ncol=10)
colnames(predict2) <- c( 'CYTOF1', 'CYTOF2', 'RNA1', 'RNA2', 'Accuracy0',
                         'Accuracy1', 'Accuracy2', 'Accuracy3', 'Accuracy4',
                         'Accuracy5')

for(k in 2:4){
  for(l in 2:4){
    for(m in 2:4){
      for(n in 2:4){
        
        list.keepX <- list(CyTOF = c(k,l), scRNA = c(m,n))
        predict <- matrix(ncol=12,nrow=0)
        for(rm in 1 : length(unique(sample.AR))){
          sequence <- which(rm == sample.AR)
          trainX <- lapply(A.AR, function(i) i[-sequence, ])
          trainY <- Y.AR[-sequence]
          xtest <- lapply(A.AR, function(i) i[sequence, ])
          
          diablo = block.splsda(X = trainX, Y = trainY, ncomp = 2, keepX = list.keepX, design = design)
          test = predict(diablo, xtest, dist  = "centroids.dist")
          predict <- rbind(predict, c(
            test$AveragedPredict.class$max.dist[,2],
            test$AveragedPredict.class$max.dist[,2],
            test$WeightedVote$centroids.dist[,2],
            test$AveragedPredict.class$max.dist[,1],
            test$AveragedPredict.class$max.dist[,1],
            test$WeightedVote$centroids.dist[,1]
          ))
        }
        accuracy0 <- (sum(predict[,1] == 'AR_t0')+sum(predict[,2] == 'AR_t1'))/length(sample.AR)
        accuracy1 <- (sum(predict[,3] == 'AR_t0')+sum(predict[,4] == 'AR_t1'))/length(sample.AR)
        accuracy2 <- (sum(predict[,5] == 'AR_t0')+sum(predict[,6] == 'AR_t1'))/length(sample.AR)
        accuracy3 <- (sum(predict[,7] == 'AR_t0')+sum(predict[,8] == 'AR_t1'))/length(sample.AR)
        accuracy4 <- (sum(predict[,9] == 'AR_t0')+sum(predict[,10] == 'AR_t1'))/length(sample.AR)
        accuracy5 <- (sum(predict[,11] == 'AR_t0')+sum(predict[,12] == 'AR_t1'))/length(sample.AR)
        
        
        result <- c(k,l,m,n, accuracy0,accuracy1, accuracy2,  accuracy3,accuracy4, accuracy5)
        predict2 <- rbind(predict2, result)
      }
    }
  }
}

predict.AR <- predict2




#make model with best tuning parameters
list.keepAAR = tune.AAR$choice.keepX



list.keepAAR$CyTOF = c(4,3)
list.keepAAR$scRNA = c(4,4)

list.keepAAR$CyTOF = c(2,2)
list.keepAAR$scRNA = c(3,2)


final.diablo.model.A.AR = block.splsda(X = A.AR, Y = Y.AR, ncomp = ncomp, 
                                         keepX = list.keepAAR, design = design)


perf.ar <- perf(final.diablo.model.A.AR, validation = 'loo')
plot(perf.ar)
ar_cons <-perf.ar$features


#make the plots to look at the models
pdf('ARs_arrowplot comp1_2_A.pdf', width = 6, height = 6)
plotArrow(final.diablo.model.A.AR, comp = c(1,2), ind.names = FALSE, legend = TRUE, 
          title = 'DIABLO')
dev.off()


pdf('ARs_loadings_A.pdf', width = 8, height = 5)
plotLoadings(final.diablo.model.A.AR, comp = 1, contrib = 'max', method = 'median')
plotLoadings(final.diablo.model.A.AR, comp = 2, contrib = 'max', method = 'median')
dev.off()



dat.ar <- as.data.frame(Reduce("+", final.diablo.model.A.AR$variates)/length(X)) %>% 
  mutate(time = Y.AR, sample = factor(sample.AR))
ggplot(dat.ar, aes(x = `comp 1`, y = `comp 2`, group = time, color = time)) +
  geom_point(size = 4) +
  stat_ellipse() +
  geom_line(aes(group = sample), color = "gray")+
  xlab("Consensus Component 1") +
  ylab("Consensus Component 2") +
  ggtitle("Allergic Rhinitis") +
  theme(legend.position = c(0.1, 0.9)) + 
  theme_bw() + 
  theme(aspect.ratio = 1) +
  scale_color_manual(values=c("#388ECC", "darkgrey"))
ggsave('AR_indiv.pdf', width = 10, height = 10, units = 'cm')



#LPO results
LPO_results <- melt(data.frame(HD = max(predict.ctrl[,10])*100,
                          AR = max(predict.AR[,10])*100))
ggplot(LPO_results, aes(x=variable, y = value, fill =variable)) + 
  geom_bar(stat='identity') + ylab('Leave person out Accuracy (%)') + 
  theme_bw() + xlab('')  + scale_fill_manual(values = c('darkgreen', 'orange'))
ggsave('accuracry_both.pdf', width = 6, height = 10, units= 'cm')



#get the loadings out per group
loadings_ctrl <- rbind(selectVar(final.diablo.model.A.ctrl, comp = 1)$CyTOF$value,
                   selectVar(final.diablo.model.A.ctrl, comp = 1)$scRNA$value)
loadings_ctrl$Feature <- rownames(loadings_ctrl)
loadings_ctrl$Feature <- factor(loadings_ctrl$Feature, 
                                levels = loadings_ctrl$Feature[order(loadings_ctrl$value.var, decreasing = T)])
loadings_ctrl$Type <- 'CyTOF'
loadings_ctrl$Type[grepl('RNA', loadings_ctrl$Feature)] <- 'scRNA'

ggplot(loadings_ctrl, aes(x=value.var, y=Feature, fill = Type)) + 
  geom_bar(stat = 'identity') + theme_bw() + 
  geom_vline(xintercept = 0, linetype = 'dashed') + 
  ggtitle('Ctrls')
ggsave('loadings_ctrls.pdf', width = 14, height= 10, units= 'cm')


loadings_AR <- rbind(selectVar(final.diablo.model.A.AR, comp = 1)$CyTOF$value,
                       selectVar(final.diablo.model.A.AR, comp = 1)$scRNA$value)
loadings_AR$Feature <- rownames(loadings_AR)
loadings_AR$Feature <- factor(loadings_AR$Feature, 
                                levels = loadings_AR$Feature[order(loadings_AR$value.var, decreasing = T)])
loadings_AR$Type <- 'CyTOF'
loadings_AR$Type[grepl('RNA', loadings_AR$Feature)] <- 'scRNA'

ggplot(loadings_AR, aes(x=value.var, y=Feature, fill = Type)) + 
  geom_bar(stat = 'identity') + theme_bw() + 
  geom_vline(xintercept = 0, linetype = 'dashed') + 
  ggtitle('ARs')
ggsave('loadings_ARs.pdf', width = 14, height= 10, units= 'cm')




#get the selected features for correlation heatmap
final.diablo.model.A.ctrl$explained_variance
features_ctrl <- c(selectVar(final.diablo.model.A.ctrl, comp = 1)$CyTOF$name,
                   selectVar(final.diablo.model.A.ctrl, comp = 1)$scRNA$name)

features_AR <-  c(selectVar(final.diablo.model.A.AR, comp = 1)$CyTOF$name,
                  selectVar(final.diablo.model.A.AR, comp = 1)$scRNA$name)

features_all <- cbind(features_ctrl, features_AR)

selected_data <- cbind(A[[1]][,colnames(A[[1]]) %in% features_all],
                       A[[2]][,colnames(A[[2]]) %in% features_all])


pval <- rcorr(as.matrix(scale(selected_data)), type='pearson')$P
corMat <- cor(scale(selected_data), use = 'pairwise', method = 'pearson')
corMat[pval>0.05] <- 0
color.blocks = brewer.pal(n = 8, name = "Set2")[1:3]

names <- data.frame('Type' = rep('scRNA', ncol(selected_data)))
names[,1] <- as.character(names[,1])
rownames(names) <- rownames(corMat)
names[grepl('Cytof', rownames(names)),1] <- 'CyTOF'
names$Group <- 'None'
names$Group[rownames(names) %in% features_AR] <- 'AR'
names$Group[rownames(names) %in% features_ctrl] <- 'HD'
names$Group[rownames(names) %in% features_AR &
              rownames(names) %in% features_ctrl] <- 'both'


column_dend = hclust(dist(corMat))

breaksList = seq(-1, 1, by = 0.01)
pdf('correlation heatmap all sel.pdf', width = 10,  height = 10)
pheatmap(corMat, annotation_row = names, annotation_col = names,
         annotation_colors =  list(Type = c(CyTOF = color.blocks[1],
                                            scRNA = color.blocks[3])),
         border_color = NA,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),
         breaks = breaksList
)
dev.off()






save.image('sep_datasets_A.RData')

