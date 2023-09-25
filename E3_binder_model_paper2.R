install.packages("tidyverse")
install.packages("rpart")
install.packages("rpart.plot")
install.packages("caret")
install.packages("zoo")
#install.packages("ggbiplot")
install.packages("reshape2")
install.packages("mlr")
install.packages("randomForest")
install.packages("readxl")
install.packages("ggrepel")

library(tidyverse)
library(rpart)
library(rpart.plot)
library(caret)
library(zoo)
library(ggbiplot)
library(reshape2)
library(mlr)
library(randomForest)
library(readxl)
library(ggrepel)


R.version
install.packages("installr")
library(installr)
check.for.updates.R()

install.packages("remotes")
remotes::install_github("vqv/ggbiplot")

setwd("../Desktop/R")
getwd()

#dataset_orig <- read_delim ("E3binder_non_db_v1_ErG.txt", delim = ";")
#dataset_orig <- read_delim ("E3binder_non_db_v1_ErG_nodupl.txt", delim = ";")
#dataset_orig <- read_delim ("E3ligase_v5_unique_ErG.csv", delim = ";")
#dataset_orig <- read_excel("novel_dataset_bind_nonbind_chembl.xlsx", 
#                           sheet = "ErG")

dataset_orig <- read_excel("novel_dataset_bind_nonbind_chembl_nostructure.xlsx")

#names(dataset_orig) <- make.names(names(dataset_orig))

#dataset_orig <- dataset_orig %>% select(ID, Class, everything())

#dataset_pca <- dataset_orig[,-c(2,3,5)]

toremove <- nearZeroVar(dataset_orig[,-c(1:2)])

toremove <- toremove + 2

dataset_pca <- dataset_orig[,-toremove]
dataset_pca <- dataset_orig

# medians on lipinski like descriptors



############################################################
#                                                          #
#                           PCA                            #
#                                                          #
############################################################
dataset_pca
dataset_pca2 <-dataset_pca[ , which(apply(dataset_pca, 2, var) != 0)]

pca_ErG <- prcomp(dataset_pca2[,-c(1:2)], scale. = T, center = T)

pca_ErG_df <- as.data.frame(pca_ErG$x[,c(1,2)], stringsAsFactors = F)
pca_ErG_df$ID <- dataset_orig$ID
pca_ErG_df$Class <- dataset_orig$E3_binder
pca_ErG_df$ID[pca_ErG_df[,1] > 0 | pca_ErG_df[,2] > -3] <- "" 

#candidate_VHL <- pca_ErG_df %>% select(ID) %>% filter(ID != "") %>% filter(!startsWith(ID, "JR")) %>% unlist() %>% unname()
#candidate_VHL <- c(candidate_VHL, "1781","1093","1598","51","1281","250","881", "1763")

p <- ggbiplot(pca_ErG,  groups = pca_ErG_df$Class, circle = F, ellipse = F, var.axes = F)#, labels = pca_ErG_df$ID)
p
p + geom_point(aes(color = as.factor(pca_ErG_df$Class))) +
  theme(legend.position = "top")+
  guides(color=guide_legend(title="Class"))#+
#geom_text_repel(box.padding = 0.5, max.overlaps = 50, label = pca_ErG_df$ID, color = as.numeric(as.factor(pca_ErG_df$Class)) )

pca_ErG_df <- as.data.frame(pca_ErG$x[,c(1,2)], stringsAsFactors = F)
pca_ErG_df$ID <- dataset_orig$ID
pca_ErG_df$Class <- dataset_orig$E3_binder
pca_ErG_df$ID[pca_ErG_df[,1] > 0 | pca_ErG_df[,2] < 0 ] <- "" 




############################################################
#                                                          #
#                           UMAP                           #
#                                                          #
############################################################

library(umap)

custom.config <- umap.defaults
custom.config$n_epochs <- 10000
custom.config$n_neighbors <- 30
custom.config$random_state <- 1206

dataset.umap <- umap(dataset_pca[,-c(1:2)], config=custom.config)

plot.az <- function(x, labels,
                    main="A UMAP visualization of the ErG E3 binder dataset",
                    colors=c("#ff7f00", "#e377c2", "#17becf"),
                    pad=0.1, cex=1, pch=19, add= FALSE, legend.suffix="",
                    cex.main=1, cex.legend=0.85) {
  layout <- x
  if (is(x, "umap")) {
    layout <- x$layout
  }
  
  xylim <- range(layout)
  xylim <- xylim + ((xylim[2]-xylim[1])*pad)*c(-0.5, 0.5)
  if (!add) {
    par(mar=c(0.2,0.7,1.2,0.7), ps=10)
    plot(xylim, xylim, type="n", axes=F, frame=F)
    rect(xylim[1], xylim[1], xylim[2], xylim[2], border="#aaaaaa", lwd=0.25)
  }
  points(layout[,1], layout[,2], col=colors[as.integer(labels)],
         cex=cex, pch=pch)
  mtext(side=3, main, cex=cex.main)
  
  labels.u <- unique(labels)
  legend.pos <- "topleft"
  legend.text <- as.character(labels.u)
  if (add) {
    legend.pos <- "bottomleft"
    legend.text <- paste(as.character(labels.u), legend.suffix)
  }
  
  legend(legend.pos, legend=legend.text, inset=0.03,
         col=colors[as.integer(labels.u)],
         bty="n", pch=pch, cex=cex.legend)
}

#range(dataset.umap$layout)

plot.az(dataset.umap, as.factor(pca_ErG_df$Class))

############################################################
#                                                          #
#                          t-SNE                           #
#                                                          #
############################################################

library(Rtsne)

fortsne <- dataset_pca[!duplicated(dataset_pca[,-c(1:2)]),]
fortsne2 <- distinct(fortsne[,-c(1,2)])

tsne_out <- Rtsne(as.matrix(fortsne2), check_duplicates = FALSE)

# Conversion of matrix to dataframe
tsne_plot <- data.frame(Dim1 = tsne_out$Y[,1],
                        Dim2 = tsne_out$Y[,2],
                        ID = fortsne$ID,
                        Class = fortsne$E3_binder, 
                        stringsAsFactors = F)

#tsne_plot <- tsne_plot %>% mutate(lbl = ifelse(Dim1 > -10 | Dim2 > -10, "", ID))
#tsne_plot <- tsne_plot %>% mutate(lbl = ifelse(Dim1 < 20 | Dim2 < 25, "", ID))
#tsne_plot <- tsne_plot %>% mutate(lbl = ifelse(Dim1 > 22 & Dim2 >0 , "", ID))
tsne_plot <- tsne_plot %>% mutate(Class2 = ifelse(Class == "NO", "Non-binder", "Binder"))
tsne_plot <- tsne_plot %>% mutate(lbl = ifelse(ID == "E3binder_v6_65", ID, ""))
tsne_plot <- tsne_plot %>% mutate(lbl = ifelse(Dim1 > 3 & Dim1 < 6 & Dim2 < 5 & Dim2 >2 , ID, ""))

# Plotting the plot using ggplot() function
ggplot2::ggplot(tsne_plot) + geom_point(aes(x=Dim1,y=Dim2, color=Class2)) +
  ggtitle("t-SNE Plot for the E3binder/nonbinder collection") +
  #geom_text_repel(aes(x=Dim1,y=Dim2),max.overlaps = 50,  label = tsne_plot$lbl)+
  theme(legend.position = "top")+
  guides(color=guide_legend(title="E3 Ligase"))


# look into the ErG more relevant descriptors

relevant <- c("Hf_Ac_d2", "Hf_Ar_d6", "Ac_Ac_d6",
              "Ar_Ac_d7", "Hf_Ar_d3", "Hf_Ar_d5",
              "Hf_Ar_d4", "Hf_Ac_d5", "Ac_Ac_d4")

relevant_IMID <- c("Ac_D_d2","Ac_D_d5","Ac_Ac_d4","Ac_Ac_d5","Ac_Ac_d6","Ac_Ac_d7")

toplot <- dataset_pca %>% dplyr::select(any_of(c("ID","E3_binder", relevant))) %>% melt()
toplot <- dataset_pca %>% dplyr::select(any_of(c("ID","Class", relevant_IMID))) %>% melt()
toplot <- dataset_pca %>% melt()

calc_boxplot_stat <- function(x) {
  coef <- 1.5
  n <- sum(!is.na(x))
  # calculate quantiles
  stats <- quantile(x, probs = c(0.0, 0.25, 0.5, 0.75, 1.0))
  names(stats) <- c("ymin", "lower", "middle", "upper", "ymax")
  iqr <- diff(stats[c(2, 4)])
  # set whiskers
  outliers <- x < (stats[2] - coef * iqr) | x > (stats[4] + coef * iqr)
  if (any(outliers)) {
    stats[c(1, 5)] <- range(c(stats[2:4], x[!outliers]), na.rm = TRUE)
  }
  return(stats)
}

library(ggsignif)

ggplot(toplot, aes(x=E3_binder, y=value, fill=E3_binder)) + 
  stat_summary(fun.data = calc_boxplot_stat, geom="boxplot") + 
  geom_signif(comparisons =list(c("NO", "YES")), test='t.test', map_signif_level = T, vjust = 1.5)+
  facet_wrap(~variable, ncol = 3,scales="free_y")+
  theme(legend.position = "none") +
  ggtitle("E3 Ligase Binder distribution of ErG relevant bits distribution")



