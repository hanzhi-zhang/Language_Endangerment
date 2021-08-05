library(diversitree)
library(ape)
library(nlme)
library(caper)
library(phytools)
library(corHMM)
library(geiger)
library(phangorn)
library(ggplot2)
library(readxl)
library(usedist)

getwd()
setwd("~/Dropbox/HRAF Summer Institute/project")


### Sino-Tibetan languages ###

ST_tree <- read.nexus("./ST_MCC.trees")
rownames(ST_data) <- ST_data$Language
ST_data[ST_tree$tip.label,] -> ST_data #order the data frame by the order of tree tiplabels

#write.csv(ST,"./ST_data.csv")

plot.phylo(ST_tree, show.node.label = F, show.tip.label = T, use.edge.length=T, cex=0.5)
tip_colour<-ST_data$`Level of endangerment`
tip_colour[which(tip_colour==0)]="steelblue1"
tip_colour[which(tip_colour==1)]="lightpink"
tip_colour[which(tip_colour==2)]="palevioletred1"
tip_colour[which(tip_colour==3)]= "orangered"
tip_colour[which(tip_colour==4)]="indianred4"
tip_colour[which(tip_colour==5)]="black"
tip_colour[is.na(tip_colour)]="grey80"
tiplabels(pch=19, cex=0.5, adj=0.5, col=tip_colour)
legend("bottomleft", legend=c("Safe",
                           "Vulnerable",
                           "Definitely endangered",
                           "Severely endangered",
                           "Critically endangered",
                           "Exinct",
                           "Missing data"),fill=c("steelblue1","lightpink","palevioletred1","orangered","indianred4","black","grey80"),bty="n",cex=0.7)

ST_tree_pruned<-drop.tip(ST_tree,ST_tree$tip.label[is.na(ST_data$`Community Marriage Organisation`)])
ST_data_pruned<-subset(ST, !is.na(ST_data$`Community Marriage Organisation`))
rownames(ST_data_pruned) <- ST_data_pruned$Language
ST_data_pruned[ST_tree_pruned$tip.label,] -> ST_data_pruned

plot.phylo(ST_tree_pruned, show.node.label = F, no.margin=T, show.tip.label = T, use.edge.length=T, cex=1)
tip_colour<-ST_data_pruned$`Community Marriage Organisation`
tip_colour[which(tip_colour=='Agamous')]="gold"
tip_colour[which(tip_colour=='Endogamy')]="lightseagreen"
tip_colour[which(tip_colour=='Exogamy')]="plum"

tiplabels(pch=19, cex=1, adj=0.5, col=tip_colour)
legend("bottomleft", legend=c("Agamous",
                              "Endoamous",
                              "Exogamous"),fill=c("gold","lightseagreen", "plum"),bty="n",cex=1)



phylosig(ST_tree, ST_data$`Level of endangerment`, method="lambda", test=T)
ST_data$'Log_speakers' <- log(ST_data$`Number of speakers`)
phylosig(ST_tree, ST_data$`Log_speakers`, method="lambda", test=T)
phylosig(ST_tree, ST_data$`Community Marriage Organisation`, method="lambda", test=T)


### Pama-Nyungan languages ###

PN_tree <- read.nexus("./PN_MCC.trees")
rownames(PN) <- PN$X1
PN[PN_tree$tip.label,] -> PN #order the data frame by the order of tree tiplabels

PN_tree_pruned<-drop.tip(PN_tree,PN_tree$tip.label[is.na(PN$`Level of endangerment`)])

PN_data_pruned<-subset(PN, !is.na(PN$`Level of endangerment`))
rownames(PN_data_pruned) <- PN_data_pruned$X1
PN_data_pruned[PN_tree_pruned$tip.label,] -> PN_data_pruned
#write.csv(ST,"./ST_data.csv")

plot.phylo(PN_tree_pruned, show.node.label = F, no.margin=T, show.tip.label = T, use.edge.length=T, cex=1)
tip_colour<-PN$`Level of endangerment`[!is.na(PN$`Level of endangerment`)]
tip_colour[which(tip_colour==0)]="steelblue1"
tip_colour[which(tip_colour==1)]="lightpink"
tip_colour[which(tip_colour==2)]="palevioletred1"
tip_colour[which(tip_colour==3)]= "orangered"
tip_colour[which(tip_colour==4)]="indianred4"
tip_colour[which(tip_colour==5)]="black"
tip_colour[is.na(tip_colour)]="grey80"
tiplabels(pch=19, cex=1, adj=0.6, col=tip_colour)
legend("bottomleft", legend=c("Safe",
                              "Vulnerable",
                              "Definitely endangered",
                              "Severely endangered",
                              "Critically endangered",
                              "Exinct"),fill=c("steelblue1","lightpink","palevioletred1","orangered","indianred4","black"),bty="n",cex=0.8)

phylosig(PN_tree_pruned, PN$`Level of endangerment`[!is.na(PN$`Level of endangerment`)], method="lambda", test=T)

PN_data_pruned$`Community marriage organisation`<-as.factor(PN_data_pruned$`Community marriage organisation`)
phylosig(PN_tree_pruned, PN_data_pruned$`Number of speakers`, method="lambda", test=T)



PN_tree_pruned_1<-drop.tip(PN_tree,PN_tree$tip.label[is.na(PN$`Community marriage organisation`)])
PN_data_pruned_1<-subset(PN, !is.na(PN$`Community marriage organisation`))
rownames(PN_data_pruned_1) <- PN_data_pruned_1$X1
PN_data_pruned_1[PN_tree_pruned_1$tip.label,] -> PN_data_pruned_1
PN_data_pruned_1$`Community marriage organisation`<-as.character(PN_data_pruned_1$`Community marriage organisation`)

plot.phylo(PN_tree_pruned_1, show.node.label = F, no.margin=T, show.tip.label = T, use.edge.length=T, cex=1)
tip_colour<-PN_data_pruned_1$`Community marriage organisation`
tip_colour[which(tip_colour=='Agamous')]="gold"
tip_colour[which(tip_colour=='Endogamous')]="lightseagreen"
tip_colour[which(tip_colour=='Exogamous')]="plum"

tiplabels(pch=19, cex=1, adj=2, col=tip_colour)
?tiplabels
legend("bottomleft", legend=c("Agamous",
                              "Endoamous",
                              "Exogamous"),fill=c("gold","lightseagreen", "plum"),bty="n",cex=1)

phylosig(PN_tree_pruned_1, PN_data_pruned_1$`Community marriage organisation`, method="lambda", test=T)
phylo.heatmap(PN_tree_pruned_1,as.matrix(PN_data_pruned_1$`Number of speakers`),standardize=TRUE)


PN_tree_pruned_2<-drop.tip(PN_tree,PN_tree$tip.label[is.na(PN$Latitude)])
PN_data_pruned_2<-subset(PN, !is.na(PN$Latitude))
rownames(PN_data_pruned_2) <- PN_data_pruned_2$X1
PN_data_pruned_2[PN_tree_pruned_2$tip.label,] -> PN_data_pruned_2

PN_coord<-as.matrix(PN_data_pruned_2[,c(12,13)])
rownames(PN_coord)<-PN_data_pruned_2$X1
obj<-phylo.to.map(PN_tree_pruned_2,PN_coord,plot=FALSE)
plot(obj,type="direct",ylim=c(-40,-10), xlim=c(110,155))


ST_tree_pruned_2<-drop.tip(ST_tree,ST_tree$tip.label[is.na(ST$Latitude)])
ST_data_pruned_2<-subset(ST, !is.na(ST$Latitude))
rownames(ST_data_pruned_2) <- ST_data_pruned_2$Language
ST_data_pruned_2[ST_tree_pruned_2$tip.label,] -> ST_data_pruned_2

ST_coord<-as.matrix(ST_data_pruned_2[,c(5,6)])
rownames(ST_coord)<-ST_data_pruned_2$X1
obj<-phylo.to.map(ST_tree_pruned_2,ST_coord,plot=FALSE)
plot(obj,type="direct",ylim=c(17,40), xlim=c(80,120))
