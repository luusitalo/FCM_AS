#############################################
# Archipelago Sea mental models
#
#Laura Uusitalo, Patrik Korn, Annaliina Koskinen 2018-2019
#############################################

setwd("D:/FCMs")

# R packages

library(FCMapper) #https://cran.r-project.org/web/packages/FCMapper/FCMapper.pdf
library(MASS)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(data.table)
library(pheatmap)

#https://cran.r-project.org/web/packages/fcm/fcm.pdf

#read in the Fuzzy Cognitive Maps
load("FCMs.RData") #This is the file including the coded maps
nmaps <- length(maps)

#Stakeholder groups of these maps
Group <- c("Recreation", "Science", "Recreation", "Science", 
           "Science", "Policy maker", "eNGO", "Policy maker", 
           "Policy maker", "Recreation", "Science", "Science", 
           "Recreation", "Policy maker", "Policy maker", "eNGO", 
           "Science", "eNGO", "eNGO", "Policy maker", "eNGO", "eNGO", 
           "Recreation", "Recreation", "Fishery", "Fishery", 
           "Fishery", "Fishery", "Fishery", "Fishery")
#################################################################


################################
# Matrix indices; number of concepts, density

mind<-t(matrix.indices(maps[[1]])[c(1:3),2])
colnames(mind) <- c("Connections", "Density", "Variables")

for (i in 2:nmaps){
  mind <- rbind(mind, t(matrix.indices(maps[[i]])[c(1:3),2]))
}


#Fix the connection density: this has been calculated assuming self-loops but they are not used here
#the equation used here is C/N^2
#we want to use C/(N*(N-1))
#therefore, multiply the current value by N^2/(N*(N-1))
for (i in 1:nrow(mind)){
  mind[i,2] <- mind[i,2]*(mind[i,3]/(mind[i,3]-1))
}

#add stakeholder group
mind <- as.data.frame(mind)
mind <- cbind(mind, Group)


#Plots of number of connections and connection density
grouplev<-c("Science", "Policy maker", "eNGO", "Recreation","Fishery")

mind$Group <- factor(mind$Group, levels=grouplev)

p1<-ggplot(mind, aes(x=Group, y=Variables, fill=Group)) +
  geom_boxplot() +
  theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_blank()) +
  labs(y="Number of variables", xlab=" ")

#anova for the differences for no of variables
var.mod<-lm(Variables~Group, data = mind)
summary(var.mod) #no statistical differences in the number of variables between the groups 

colnames(mind) <- c("Connections", "Density", "Variables", "Group")
mind$Group <- factor(mind$Group, levels=grouplev)

p2<-ggplot(mind, aes(x=Group, y=Density, fill=Group)) +
  geom_boxplot() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_blank()) +
  theme(legend.position = "none") +
  labs(y="Connection density", xlab=" ")

#anova for the density
den.mod<-lm(Density~Group, data = mind)
summary(den.mod) #no statistical differences in the connection density between the groups 

png(file="VarNo-Density.png", units="cm", width = 20, height = 10, res=200)
ggarrange(p1, p2, ncol = 2)
dev.off()



## Draw a scatterplot showing each map's number of concepts and density
## This is not in the paper

ggplot(mind, aes(x=Variables, y=Density, shape=Group, color=Group)) +
  geom_point(size=4)

#correlation between #variables and conn.density =~ -0.12
cor(mind$Variables, mind$Density)


###########################
# Centrality analyses

# Harmonize all matrices into the same format for analysis

concept.names <- colnames(maps[[1]])

for (i in 2:nmaps) {
  concept.names <- c(concept.names, colnames(maps[[i]]))
}

all.concept.names <- sort(unique(concept.names))

#empty "all variables" matrix
all_zeros<-matrix(0, nrow=length(all.concept.names), ncol=length(all.concept.names))
colnames(all_zeros)<-all.concept.names
rownames(all_zeros)<-all.concept.names

#combine every map with the all_zeros map
fullmaps <- list(combine.maps(all_zeros, maps[[1]], concept.names1 = all.concept.names, concept.names2 = colnames(maps[[1]])))

for (i in 2:nmaps){
  fullmaps <- c(fullmaps, list(combine.maps(all_zeros, maps[[i]], concept.names1 = all.concept.names, concept.names2 = colnames(maps[[i]])))) 
}

#Concept indices show the importance of each variable in the model: indegree, outdegree, centrality...
CI <- list(concept.indices(fullmaps[[1]], all.concept.names))

for (i in 2:nmaps){
  CI <- c(CI, list(concept.indices(fullmaps[[i]], all.concept.names)))
}


#make a matrix of concept centralities per map
centralities <- CI[[1]][,c(1,4)]
centralities <-cbind(centralities, rep(Group[1],nrow(CI[[1]])))
centralities <-cbind(centralities, rep(1,nrow(CI[[1]])))
colnames(centralities)[3]<-"Group"
colnames(centralities)[4]<-"map"


for (i in 2:nmaps){
  tmp <- CI[[i]][,c(1,4)]
  tmp <- cbind(tmp, rep(Group[i],nrow(CI[[i]])))
  tmp <- cbind(tmp, rep(i,nrow(CI[[i]])))
    colnames(tmp)[3]<-"Group"
    colnames(tmp)[4] <- "map"
  centralities <- rbind(centralities, tmp)
}

cenOrder <- setDT(centralities)[, lapply(.SD, mean),.SDcols=c("Centrality"), by=Concept]
cenOrder <- cenOrder[order(-cenOrder$Centrality),] #38 concepts (=variables) in total

conlev<-cenOrder$Concept
centralities$Concept <- factor(centralities$Concept, levels=conlev)

centralitiesSci <- centralities[centralities$Group=="Science",]
centralitieseNGO <- centralities[centralities$Group=="eNGO",]
centralitiesPol <- centralities[centralities$Group=="Policy maker",]
centralitiesRec <- centralities[centralities$Group=="Recreation",]
centralitiesFish <- centralities[centralities$Group=="Fishery",]


####################

cenAllforPresentation <- ggplot(centralities, aes(x=Concept, y=Centrality, fill=Concept)) + 
  geom_boxplot() +
  ylim(0,13) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size=20), axis.title.x = element_blank(), legend.position = "none") +  
  labs(y="Variable centrality", xlab=" ") + annotate("text", x=1, y=13, label="All maps", hjust=0) + labs(title=element_blank())




cenAll<-ggplot(centralities, aes(x=Concept, y=Centrality, fill=Concept)) + 
  geom_boxplot() +
  ylim(0,13) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text.x = element_blank(), axis.title.x = element_blank(), legend.position = "none") + 
  labs(y="Variable centrality", xlab=" ") + annotate("text", x=1, y=13, label="a) All maps", hjust=0) + labs(title=element_blank())


cenSci<-ggplot(centralitiesSci, aes(x=Concept, y=Centrality, fill=Concept)) + 
  geom_boxplot() +
  ylim(0,13) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text.x = element_blank(), axis.title.x = element_blank(), legend.position = "none") + 
  labs(y=" ", xlab=" ") + annotate("text", x=1, y=13, label="b) Science, n=6", hjust=0) + labs(title=element_blank())

ceneNGO<-ggplot(centralitieseNGO, aes(x=Concept, y=Centrality, fill=Concept)) + 
  geom_boxplot() +
  ylim(0,13) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text.x = element_blank(), axis.title.x = element_blank(), legend.position = "none") + 
  labs(y="", xlab=" ") + annotate("text", x=1, y=13, label="d) eNGO, n=6", hjust=0) + labs(title=element_blank())

cenPol<-ggplot(centralitiesPol, aes(x=Concept, y=Centrality, fill=Concept)) + 
  geom_boxplot() +
  ylim(0,13) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text.x = element_blank(), axis.title.x = element_blank(), legend.position = "none") + 
  labs(y="Variable centrality", xlab=" ")+ annotate("text", x=1, y=13, label="c) Policy makers n=6", hjust=0) + labs(title=element_blank())

cenRec<-ggplot(centralitiesRec, aes(x=Concept, y=Centrality, fill=Concept)) + 
  geom_boxplot() +
  ylim(0,13) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_blank(), legend.position = "none") + 
  labs(y="Variable centrality", xlab=" ") + annotate("text", x=1, y=13, label="e) Recretion, n=6", hjust=0) + labs(title=element_blank())

cenFish<-ggplot(centralitiesFish, aes(x=Concept, y=Centrality, fill=Concept)) +
  geom_boxplot() +
  ylim(0,13) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_blank(), legend.position = "none") + 
  labs(y=" ", xlab=" ") + annotate("text", x=1, y=13, label="f) Fishery, n=6", hjust=0) + labs(title=element_blank())


png(file="centralityplot.png", units="cm", width = 32, height = 30, res=200)
ggarrange(cenAll, cenSci, cenPol,ceneNGO,cenRec,cenFish, ncol = 2, nrow = 3)
dev.off()

###
#cluster the models by variable centrality
###

#format the matrix
centralities2 <- transpose(centralities[1:length(conlev),2]) 
for (i in 2:nmaps){
  centralities2 <- rbind(centralities2, transpose(centralities[centralities$map==i,2]))
}
colnames(centralities2)<-as.character(transpose(centralities[1:length(conlev),1]))

#Multidimensional scaling:
met<- "euclidean" 
d <- dist(centralities2[1:nmaps,], method = met) # distance matrix
malli.mds <- isoMDS(d)
MDS_xy <- data.frame(malli.mds$points)
MDS_xy$Group <-  Group
group <- as.factor(MDS_xy$Group)
malli.mds$stress

#malli.cmd<-cmdscale(d, k=2, eig=TRUE)
#malli.cmd$GOF

ClustCentr <- ggplot(MDS_xy, aes(x=X1, y=X2, shape=Group, color=Group)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.title.x = element_blank(), legend.position = "none") + 
        scale_shape_manual(values = c(15, 16, 17, 18, 20)) +
        scale_color_manual(values = c("coral2", "springgreen2", "cornflowerblue", "orange1", "mediumpurple4")) +
  geom_point(size=4) + labs(title="Clustering by variable centrality", x="", y="") +
  theme(legend.position = "none")
  

#####################
# clustering by the linkages

#turn these matrices into vectors and combine into one matrix for clustering
# one row <- one map; one column is an interaction between two elements

all_models <- as.vector(fullmaps[[1]])
for (i in 2:nmaps){
  all_models <- rbind(all_models, as.vector(fullmaps[[i]]))
}

#Multidimensional scaling:
d2 <- dist(all_models, method = met) # distance matrix

malli2.mds <- isoMDS(d2)
MDS_xy2 <- data.frame(malli2.mds$points)
MDS_xy2$Group <-  Group
malli2.mds$stress

ClustLinkages <- ggplot(MDS_xy2, aes(x=X1, y=X2, shape=Group, color=Group)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.title.x = element_blank(), legend.position = "right", legend.key = element_blank()) + 
        scale_shape_manual(values = c(15, 16, 17, 18, 20)) +
        scale_color_manual(values = c("coral2", "springgreen2", "cornflowerblue", "orange1", "mediumpurple4")) +
  geom_point(size=4) + labs(title="Clustering by linkages", x="", y="")


png(file="MDSplot.png", units="cm", width = 32, height = 15, res=200)
ggarrange(ClustCentr, ClustLinkages, ncol = 2)
dev.off()


###########
# The most important links in the maps
##This is not in the paper
allm<-as.array(fullmaps)
#there's probably a function that does this prettily...
r <- nrow(fullmaps[[1]]) #square matrix; nrow and ncol are equal

#the most important links across all the maps
#initialize a matrix to compute the means
linksAll <- fullmaps[[1]]
linksAll[] <- 0

#compute means
for (i in 1:length(fullmaps)){
  linksAll <- linksAll + fullmaps[[i]]
}
linksAll <- linksAll/length(fullmaps)

colors <- inlmisc::GetColors(9, scheme = "BuRd")
meanLinks <- pheatmap(linksAll, cluster_rows=FALSE, cluster_cols=FALSE, col=colors)

png(file="meanLinks.png", units="cm", width = 20, height = 20, res=200)
meanLinks
dev.off()


######################################
# Simulating the scenarios
######################################

#nochanges
scenarioruns <- nochanges.scenario(maps[1][[1]], colnames(maps[[1]]), 30)  
scenarioruns$mapID <- 1
scenarioruns$Scenario <- "No changes"
for (i in (2:nmaps)) {
  tmp <- nochanges.scenario(maps[i][[1]], colnames(maps[i][[1]]), 30)
  tmp$mapID <- i
  tmp$Scenario <- "No changes"
  scenarioruns<-rbind(scenarioruns, tmp)
}

#increased fishing
for (i in (1:nmaps)) {
  tmp <- changes.scenario(maps[i][[1]], colnames(maps[i][[1]]), 30, set.concepts = "Fishing", 1)
  tmp$mapID <- i
  tmp$Scenario <- "Fishing+"
  scenarioruns<-rbind(scenarioruns, tmp)
}

#decreased fishing
for (i in (1:nmaps)) {
  tmp <- changes.scenario(maps[i][[1]], colnames(maps[i][[1]]), 30, set.concepts = "Fishing", -1)
  tmp$mapID <- i
  tmp$Scenario <- "Fishing-"
  scenarioruns<-rbind(scenarioruns, tmp)
}

#increased eutrophication 

for (i in (1:nmaps)) {
  tmp <- changes.scenario(maps[i][[1]], colnames(maps[i][[1]]), 30, set.concepts = "Eutrophication", 1)
  tmp$mapID <- i
  tmp$Scenario <- "Eutrophication+"
  scenarioruns<-rbind(scenarioruns, tmp)
}

#decreased eutrophication
for (i in (1:nmaps)) {
  tmp <- changes.scenario(maps[i][[1]], colnames(maps[i][[1]]), 30, set.concepts = "Eutrophication", -1)
  tmp$mapID <- i
  tmp$Scenario <- "Eutrophication-"
  scenarioruns<-rbind(scenarioruns, tmp)
}

#increasing salinity
#Salinity variable was defined in the interviews as "the decrease of salinity", hence "salinity 1" means decrease of salinity stengthening
for (i in (1:nmaps)) {
  tmp <- changes.scenario(maps[i][[1]], colnames(maps[i][[1]]), 30, set.concepts = "Salinity", -1)
  tmp$mapID <- i
  tmp$Scenario <- "Salinity+"
  scenarioruns<-rbind(scenarioruns, tmp)
}

#decreasing salinity
#Salinity variable was defined in the interviews as "the decrease of salinity", hence "salinity 1" means decrease of salinity stengthening
for (i in (1:nmaps)) {
  tmp <- changes.scenario(maps[i][[1]], colnames(maps[i][[1]]), 30, set.concepts = "Salinity", 1)
  tmp$mapID <- i
  tmp$Scenario <- "Salinity-"
  scenarioruns<-rbind(scenarioruns, tmp)
}

#increasing temperature
for (i in (1:nmaps)) {
  tmp <- changes.scenario(maps[i][[1]], colnames(maps[i][[1]]), 30, set.concepts = "Temperature", 1)
  tmp$mapID <- i
  tmp$Scenario <- "Temperature+"
  scenarioruns<-rbind(scenarioruns, tmp)
}

#decreasing temperature
for (i in (1:nmaps)) {
  tmp <- changes.scenario(maps[i][[1]], colnames(maps[i][[1]]), 30, set.concepts = "Temperature", -1)
  tmp$mapID <- i
  tmp$Scenario <- "Temperature-"
  scenarioruns<-rbind(scenarioruns, tmp)
}

#add the stakeholder group info
scenarioruns$Group[scenarioruns$mapID %in% c(1, 3, 10, 13)] <- "Recreation"
scenarioruns$Group[scenarioruns$mapID %in% c(2,4,5,11,12,17)] <- "Science"
scenarioruns$Group[scenarioruns$mapID %in% c(6,8,9,14,15, 20)] <- "Policy maker"
scenarioruns$Group[scenarioruns$mapID %in% c(7,16, 18, 19, 21, 22)] <- "eNGO"
scenarioruns$Group[scenarioruns$mapID %in% c(25:30)] <- "Fishery"

write.csv(scenarioruns, "scenario.runs.csv")


# compute & plot the difference from the nochanges scenario for each map & each scenario!
# the scenarioruns table doesn't include information about which map is which - 
# need to modify the scenario simulations so that the map number info is included in the 
# final table

scenariodiff <- scenarioruns
noc.number <- nrow(scenarioruns [scenarioruns$Scenario=="No changes",]) # how many nochange.scenarios rows
noc <- rep(scenarioruns$Equilibrium_value[1:noc.number], 9) #dummy
for (i in 1 : nrow(scenarioruns)){
  #9 scenarios (including the nochanges)
  scenariodiff$difference[i] <- scenarioruns$Equilibrium_value[i] - noc[i]
}

######################
## Plot the 10 most central ecosystem components results per each scenario
## for all maps

tencomp<-c("Cod", "Herring", "Perch",  "Zooplankton", "Benthos", "Zander", "Phytoplankton", "Pike",  "Birds", "Mammals")


#strip the data so that it includes only the relevant ecosystem components
scenariodiff2 <- scenariodiff[scenariodiff$Concept %in% tencomp,]
lev <- tencomp
scenariodiff2$Concept <- factor(scenariodiff2$Concept, levels=lev)


#Fishing +
fishP<-subset(scenariodiff2, Scenario=="Fishing+")
allfPp<-ggplot(fishP, aes(x=Concept, y=difference, fill=Concept)) + 
  geom_boxplot() +
  #geom_violin()+
  geom_hline(yintercept=0) +
  ylim(-1,1) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_blank(), 
        legend.position = "none") + 
  labs(y="difference to base run", xlab=" ", title="High fishing pressure")

#allfPp

#Fishing -
fishM<-subset(scenariodiff2, Scenario=="Fishing-")
allfMp<-ggplot(fishM, aes(x=Concept, y=difference, fill=Concept)) + 
  geom_boxplot() + 
  #geom_violin() +
  geom_hline(yintercept=0) +
  ylim(-1,1) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_blank(), 
        legend.position = "none") +   
  labs(y="difference to base run", xlab=" ", title="Low fishing pressure") 

#allfMp

#Eutrophication +
eutroP<-subset(scenariodiff2, Scenario=="Eutrophication+")
allePp<-ggplot(eutroP, aes(x=Concept, y=difference, fill=Concept)) + 
  geom_boxplot() + 
  geom_hline(yintercept=0) +
  ylim(-1,1) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_blank(), 
        legend.position = "none") + 
  labs(y="difference to base run", xlab=" ", title="High eutrophication") 

#allePp

#Eutroophication -
eutroM<-subset(scenariodiff2, Scenario=="Eutrophication-")
alleMp<-ggplot(eutroM, aes(x=Concept, y=difference, fill=Concept)) + 
  geom_boxplot() + 
  geom_hline(yintercept=0) +
  ylim(-1,1) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_blank(), 
        legend.position = "none") + 
  labs(y="difference to base run", xlab=" ", title="Low eutrophication") 

#alleMp

#Salinity +
salP<-subset(scenariodiff2, Scenario=="Salinity+")
allsPp<-ggplot(salP, aes(x=Concept, y=difference, fill=Concept)) + 
  geom_boxplot() + 
  geom_hline(yintercept=0) +
  ylim(-1,1) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_blank(), 
        legend.position = "none") + 
  labs(y="difference to base run", xlab=" ", title="High salinity") 

#allsPp

#Salinity -
salM<-subset(scenariodiff2, Scenario=="Salinity-")
allsMp<-ggplot(salM, aes(x=Concept, y=difference, fill=Concept)) + 
  geom_boxplot() + 
  geom_hline(yintercept=0) +
  ylim(-1,1) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_blank(), 
        legend.position = "none") + 
  labs(y="difference to base run", xlab=" ", title="Low salinity") 

#allsMp


#Temperature +
tempP<-subset(scenariodiff2, Scenario=="Temperature+")
alltPp<-ggplot(tempP, aes(x=Concept, y=difference, fill=Concept)) + 
  geom_boxplot() + 
  geom_hline(yintercept=0) +
  ylim(-1,1) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_blank(), 
        legend.position = "none") + 
  labs(y="difference to base run", xlab=" ", title="High temperature") 

#alltPp

#Temp-
tempM<-subset(scenariodiff2, Scenario=="Temperature-")
alltMp<-ggplot(tempM, aes(x=Concept, y=difference, fill=Concept)) + 
  geom_boxplot() + 
  geom_hline(yintercept=0) +
  ylim(-1,1) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_blank(), 
        legend.position = "none") + 
    labs(y="difference to base run", xlab=" ", title="Low temperature") 

#alltMp

#the following plot gets super small, have to do something about it
png(file="diffplot.png", units="cm", width = 20, height = 30, res=200)

plot_grid(allfPp + annotate("text", x=1, y=1, label="a) High fishing pressure", hjust=0) + labs(title=element_blank()),
           allfMp+ annotate("text", x=1, y=1, label="b) Low fishing pressure", hjust=0) + labs(title=element_blank()), 
          allePp+ annotate("text", x=1, y=1, label="c) High eutrophication", hjust=0) + labs(title=element_blank()), 
          alleMp+ annotate("text", x=1, y=1, label="d) Low eutrophication", hjust=0) + labs(title=element_blank()), 
          alltPp+ annotate("text", x=1, y=1, label="e) High temperature", hjust=0) + labs(title=element_blank()), 
          alltMp+ annotate("text", x=1, y=1, label="f) Low temperature", hjust=0) + labs(title=element_blank()), 
          allsPp+ annotate("text", x=1, y=1, label="g) High salinity", hjust=0) + labs(title=element_blank()),
          allsMp+ annotate("text", x=1, y=1, label="h) Low salinity", hjust=0) + labs(title=element_blank()),
                   align = 'vh',
                   hjust = -1,
                   nrow = 4)

dev.off()


######################
## Plot the 10 most central ecosystem components resutls per each scenario
## FOR SCIENTISTS

#strip the data so that it includes only the relevant ecosystem components & stakeholder group
scenariodiff3 <- scenariodiff2[scenariodiff2$Group=="Science",]

#Fishing +
fishP<-subset(scenariodiff3, Scenario=="Fishing+")
scifPp<-ggplot(fishP, aes(x=Concept, y=difference, fill=Concept)) + 
  geom_boxplot() + 
  geom_hline(yintercept=0) +
  ylim(-1,1) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_blank(), 
        legend.position = "none") + 
   labs(y="difference to base run", xlab=" ", title="High fishing pressure")

#scifPp

#Fishing -
fishM<-subset(scenariodiff3, Scenario=="Fishing-")
scifMp<-ggplot(fishM, aes(x=Concept, y=difference, fill=Concept)) + 
  geom_boxplot() + 
  geom_hline(yintercept=0) +
  ylim(-1,1) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_blank(), 
        legend.position = "none") + 
    labs(y="difference to base run", xlab=" ", title="Low fishing pressure") 

#scifMp

#Eutrophication +
eutroP<-subset(scenariodiff3, Scenario=="Eutrophication+")
sciePp<-ggplot(eutroP, aes(x=Concept, y=difference, fill=Concept)) + 
  geom_boxplot() + 
  geom_hline(yintercept=0) +
  ylim(-1,1) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_blank(), 
        legend.position = "none") + 
  labs(y="difference to base run", xlab=" ", title="High eutrophication") 

#sciePp

#Eutroophication -
eutroM<-subset(scenariodiff3, Scenario=="Eutrophication-")
scieMp<-ggplot(eutroM, aes(x=Concept, y=difference, fill=Concept)) + 
  geom_boxplot() + 
  geom_hline(yintercept=0) +
  ylim(-1,1) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_blank(), 
        legend.position = "none") + 
  labs(y="difference to base run", xlab=" ", title="Low eutrophication") 

#scieMp

#Salinity +
salP<-subset(scenariodiff3, Scenario=="Salinity+")
scisPp<-ggplot(salP, aes(x=Concept, y=difference, fill=Concept)) + 
  geom_boxplot() + 
  geom_hline(yintercept=0) +
  ylim(-1,1) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_blank(), 
        legend.position = "none") + 
  labs(y="difference to base run", xlab=" ", title="High salinity") 

#scisPp

#Salinity -
salM<-subset(scenariodiff3, Scenario=="Salinity-")
scisMp<-ggplot(salM, aes(x=Concept, y=difference, fill=Concept)) + 
  geom_boxplot() + 
  geom_hline(yintercept=0) +
  ylim(-1,1) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_blank(), 
        legend.position = "none") + 
  labs(y="difference to base run", xlab=" ", title="Low salinity") 

#scisMp


#Temperature +
tempP<-subset(scenariodiff3, Scenario=="Temperature+")
scitPp<-ggplot(tempP, aes(x=Concept, y=difference, fill=Concept)) + 
  geom_boxplot() + 
  geom_hline(yintercept=0) +
  ylim(-1,1) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_blank(), 
        legend.position = "none") + 
  labs(y="difference to base run", xlab=" ", title="High temperature") 

#scitPp

#Temperature -
tempM<-subset(scenariodiff3, Scenario=="Temperature-")
scitMp<-ggplot(tempM, aes(x=Concept, y=difference, fill=Concept)) + 
  geom_boxplot() + 
  geom_hline(yintercept=0) +
  ylim(-1,1) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_blank(), 
        legend.position = "none") + 
  labs(y="difference to base run", xlab=" ", title="Low temperature") 

#scitMp

png(file="diffplotScience.png", units="cm", width = 20, height = 30, res=200)

plot_grid(scifPp + annotate("text", x=1, y=1, label="a) High fishing pressure", hjust=0) + labs(title=element_blank()),
          scifMp+ annotate("text", x=1, y=1, label="b) Low fishing pressure", hjust=0) + labs(title=element_blank()), 
          sciePp+ annotate("text", x=1, y=1, label="c) High eutrophication", hjust=0) + labs(title=element_blank()), 
          scieMp+ annotate("text", x=1, y=1, label="d) Low eutrophication", hjust=0) + labs(title=element_blank()), 
          scitPp+ annotate("text", x=1, y=1, label="e) High temperature", hjust=0) + labs(title=element_blank()), 
          scitMp+ annotate("text", x=1, y=1, label="f) Low temperature", hjust=0) + labs(title=element_blank()), 
          scisPp+ annotate("text", x=1, y=1, label="g) High salinity", hjust=0) + labs(title=element_blank()),
          scisMp+ annotate("text", x=1, y=1, label="h) Low salinity", hjust=0) + labs(title=element_blank()),
          align = 'vh',
          hjust = -1,
          nrow = 4)

dev.off()

######################
## Plot the 10 most central ecosystem components resutls per each scenario
## FOR POLICY MAKERS

#strip the data so that it includes only the relevant ecosystem components & stakeholder group
scenariodiff3 <- scenariodiff2[scenariodiff2$Group=="Policy maker",]

#Fishing +
fishP<-subset(scenariodiff3, Scenario=="Fishing+")
polfPp<-ggplot(fishP, aes(x=Concept, y=difference, fill=Concept)) + 
  geom_boxplot() + 
  geom_hline(yintercept=0) +
  ylim(-1,1) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_blank(), 
        legend.position = "none") + 
    labs(y="difference to base run", xlab=" ", title="High fishing pressure")

#polfPp

#Fishing -
fishM<-subset(scenariodiff3, Scenario=="Fishing-")
polfMp<-ggplot(fishM, aes(x=Concept, y=difference, fill=Concept)) + 
  geom_boxplot() + 
  geom_hline(yintercept=0) +
  ylim(-1,1) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_blank(), 
        legend.position = "none") + 
  labs(y="difference to base run", xlab=" ", title="Low fishing pressure") 

#polfMp

#Eutrophication +
eutroP<-subset(scenariodiff3, Scenario=="Eutrophication+")
polePp<-ggplot(eutroP, aes(x=Concept, y=difference, fill=Concept)) + 
  geom_boxplot() + 
  geom_hline(yintercept=0) +
  ylim(-1,1) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_blank(), 
        legend.position = "none") + 
  labs(y="difference to base run", xlab=" ", title="High eutrophication") 

#polePp

#Eutroophication -
eutroM<-subset(scenariodiff3, Scenario=="Eutrophication-")
poleMp<-ggplot(eutroM, aes(x=Concept, y=difference, fill=Concept)) + 
  geom_boxplot() + 
  geom_hline(yintercept=0) +
  ylim(-1,1) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_blank(), 
        legend.position = "none") + 
  labs(y="difference to base run", xlab=" ", title="Low eutrophication") 

#poleMp

#Salinity +
salP<-subset(scenariodiff3, Scenario=="Salinity+")
polsPp<-ggplot(salP, aes(x=Concept, y=difference, fill=Concept)) + 
  geom_boxplot() + 
  geom_hline(yintercept=0) +
  ylim(-1,1) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_blank(), 
        legend.position = "none") + 
  labs(y="difference to base run", xlab=" ", title="High salinity") 

#polsPp

#Salinity -
salM<-subset(scenariodiff3, Scenario=="Salinity-")
polsMp<-ggplot(salM, aes(x=Concept, y=difference, fill=Concept)) + 
  geom_boxplot() + 
  geom_hline(yintercept=0) +
  ylim(-1,1) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_blank(), 
        legend.position = "none") + 
  labs(y="difference to base run", xlab=" ", title="Low salinity") 

#polsMp


#Temperature +
tempP<-subset(scenariodiff3, Scenario=="Temperature+")
poltPp<-ggplot(tempP, aes(x=Concept, y=difference, fill=Concept)) + 
  geom_boxplot() + 
  geom_hline(yintercept=0) +
  ylim(-1,1) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_blank(), 
        legend.position = "none") + 
  labs(y="difference to base run", xlab=" ", title="High temperature") 

#poltPp

#Temperature -
tempM<-subset(scenariodiff3, Scenario=="Temperature-")
poltMp<-ggplot(tempM, aes(x=Concept, y=difference, fill=Concept)) + 
  geom_boxplot() + 
  geom_hline(yintercept=0) +
  ylim(-1,1) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_blank(), 
        legend.position = "none") + 
  labs(y="difference to base run", xlab=" ", title="Low temperature") 

#poltMp

png(file="diffplotPolicy.png", units="cm", width = 20, height = 30, res=200)

plot_grid(polfPp + annotate("text", x=1, y=1, label="a) High fishing pressure", hjust=0) + labs(title=element_blank()),
          polfMp+ annotate("text", x=1, y=1, label="b) Low fishing pressure", hjust=0) + labs(title=element_blank()), 
          polePp+ annotate("text", x=1, y=1, label="c) High eutrophication", hjust=0) + labs(title=element_blank()), 
          poleMp+ annotate("text", x=1, y=1, label="d) Low eutrophication", hjust=0) + labs(title=element_blank()), 
          poltPp+ annotate("text", x=1, y=1, label="e) High temperature", hjust=0) + labs(title=element_blank()), 
          poltMp+ annotate("text", x=1, y=1, label="f) Low temperature", hjust=0) + labs(title=element_blank()), 
          polsPp+ annotate("text", x=1, y=1, label="g) High salinity", hjust=0) + labs(title=element_blank()),
          polsMp+ annotate("text", x=1, y=1, label="h) Low salinity", hjust=0) + labs(title=element_blank()),
          align = 'vh',
          hjust = -1,
          nrow = 4)

dev.off()


######################
## Plot the 10 most central ecosystem components resutls per each scenario
## FOR eNGOs

#strip the data so that it includes only the relevant ecosystem components & stakeholder group
scenariodiff3 <- scenariodiff2[scenariodiff2$Group=="eNGO",]

#Fishing +
fishP<-subset(scenariodiff3, Scenario=="Fishing+")
engofPp<-ggplot(fishP, aes(x=Concept, y=difference, fill=Concept)) + 
  geom_boxplot() + 
  geom_hline(yintercept=0) +
  ylim(-1,1) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_blank(), 
        legend.position = "none") + 
  labs(y="difference to base run", xlab=" ", title="High fishing pressure")

#engofPp

#Fishing -
fishM<-subset(scenariodiff3, Scenario=="Fishing-")
engofMp<-ggplot(fishM, aes(x=Concept, y=difference, fill=Concept)) + 
  geom_boxplot() + 
  geom_hline(yintercept=0) +
  ylim(-1,1) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_blank(), 
        legend.position = "none") + 
    labs(y="difference to base run", xlab=" ", title="Low fishing pressure") 

#engofMp

#Eutrophication +
eutroP<-subset(scenariodiff3, Scenario=="Eutrophication+")
engoePp<-ggplot(eutroP, aes(x=Concept, y=difference, fill=Concept)) + 
  geom_boxplot() + 
  geom_hline(yintercept=0) +
  ylim(-1,1) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_blank(), 
        legend.position = "none") + 
    labs(y="difference to base run", xlab=" ", title="High eutrophication") 

#engoePp

#Eutrophication -
eutroM<-subset(scenariodiff3, Scenario=="Eutrophication-")
engoeMp<-ggplot(eutroM, aes(x=Concept, y=difference, fill=Concept)) + 
  geom_boxplot() + 
  geom_hline(yintercept=0) +
  ylim(-1,1) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_blank(), 
        legend.position = "none") + 
    labs(y="difference to base run", xlab=" ", title="Low eutrophication") 

#engoeMp

#Salinity +
salP<-subset(scenariodiff3, Scenario=="Salinity+")
engosPp<-ggplot(salP, aes(x=Concept, y=difference, fill=Concept)) + 
  geom_boxplot() + 
  geom_hline(yintercept=0) +
  ylim(-1,1) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_blank(), 
        legend.position = "none") + 
    labs(y="difference to base run", xlab=" ", title="High salinity") 

#engosPp

#Salinity -
salM<-subset(scenariodiff3, Scenario=="Salinity-")
engosMp<-ggplot(salM, aes(x=Concept, y=difference, fill=Concept)) + 
  geom_boxplot() + 
  geom_hline(yintercept=0) +
  ylim(-1,1) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_blank(), 
        legend.position = "none") + 
    labs(y="difference to base run", xlab=" ", title="Low salinity") 

#engosMp


#Temperature +
tempP<-subset(scenariodiff3, Scenario=="Temperature+")
engotPp<-ggplot(tempP, aes(x=Concept, y=difference, fill=Concept)) + 
  geom_boxplot() + 
  geom_hline(yintercept=0) +
  ylim(-1,1) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_blank(), 
        legend.position = "none") + 
    labs(y="difference to base run", xlab=" ", title="High temperature") 

#engotPp

#Temperature -
tempM<-subset(scenariodiff3, Scenario=="Temperature-")
engotMp<-ggplot(tempM, aes(x=Concept, y=difference, fill=Concept)) + 
  geom_boxplot() + 
  geom_hline(yintercept=0) +
  ylim(-1,1) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_blank(), 
        legend.position = "none") + 
    labs(y="difference to base run", xlab=" ", title="Low temperature") 

#engotMp

png(file="diffploteNGO.png", units="cm", width = 20, height = 30, res=200)

plot_grid(engofPp + annotate("text", x=1, y=1, label="a) High fishing pressure", hjust=0) + labs(title=element_blank()),
          engofMp+ annotate("text", x=1, y=1, label="b) Low fishing pressure", hjust=0) + labs(title=element_blank()), 
          engoePp+ annotate("text", x=1, y=1, label="c) High eutrophication", hjust=0) + labs(title=element_blank()), 
          engoeMp+ annotate("text", x=1, y=1, label="d) Low eutrophication", hjust=0) + labs(title=element_blank()), 
          engotPp+ annotate("text", x=1, y=1, label="e) High temperature", hjust=0) + labs(title=element_blank()), 
          engotMp+ annotate("text", x=1, y=1, label="f) Low temperature", hjust=0) + labs(title=element_blank()), 
          engosPp+ annotate("text", x=1, y=1, label="g) High salinity", hjust=0) + labs(title=element_blank()),
          engosMp+ annotate("text", x=1, y=1, label="h) Low salinity", hjust=0) + labs(title=element_blank()),
          align = 'vh',
          hjust = -1,
          nrow = 4)

dev.off()


######################
## Plot the 10 most central ecosystem components resutls per each scenario
## FOR recreational users

#strip the data so that it includes only the relevant ecosystem components & stakeholder group
scenariodiff3 <- scenariodiff2[scenariodiff2$Group=="Recreation",]


#Fishing +
fishP<-subset(scenariodiff3, Scenario=="Fishing+")
recfPp<-ggplot(fishP, aes(x=Concept, y=difference, fill=Concept)) + 
  geom_boxplot() + 
  geom_hline(yintercept=0) +
  ylim(-1,1) + 
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),       
      panel.background = element_blank(), axis.line = element_line(colour = "black"),        
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_blank(),        
      legend.position = "none") + 
  labs(y="difference to base run", xlab=" ", title="High fishing pressure")



#recfPp

#Fishing -
fishM<-subset(scenariodiff3, Scenario=="Fishing-")
recfMp<-ggplot(fishM, aes(x=Concept, y=difference, fill=Concept)) + 
  geom_boxplot() + 
  geom_hline(yintercept=0) +
  ylim(-1,1) + 
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),       
      panel.background = element_blank(), axis.line = element_line(colour = "black"),        
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_blank(),        
      legend.position = "none") + 
  labs(y="difference to base run", xlab=" ", title="Low fishing pressure") 

recfMp

#Eutrophication +
eutroP<-subset(scenariodiff3, Scenario=="Eutrophication+")
recePp<-ggplot(eutroP, aes(x=Concept, y=difference, fill=Concept)) + 
  geom_boxplot() + 
  geom_hline(yintercept=0) +
  ylim(-1,1) + 
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),       
      panel.background = element_blank(), axis.line = element_line(colour = "black"),        
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_blank(),        
      legend.position = "none") + 
  labs(y="difference to base run", xlab=" ", title="High eutrophication") 

recePp

#Eutroophication -
eutroM<-subset(scenariodiff3, Scenario=="Eutrophication-")
receMp<-ggplot(eutroM, aes(x=Concept, y=difference, fill=Concept)) + 
  geom_boxplot() + 
  geom_hline(yintercept=0) +
  ylim(-1,1) + 
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),       
      panel.background = element_blank(), axis.line = element_line(colour = "black"),        
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_blank(),        
      legend.position = "none") + 
  labs(y="difference to base run", xlab=" ", title="Low eutrophication") 

receMp

#Salinity +
salP<-subset(scenariodiff3, Scenario=="Salinity+")
recsPp<-ggplot(salP, aes(x=Concept, y=difference, fill=Concept)) + 
  geom_boxplot() + 
  geom_hline(yintercept=0) +
  ylim(-1,1) + 
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),       
      panel.background = element_blank(), axis.line = element_line(colour = "black"),        
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_blank(),        
      legend.position = "none") + 
  labs(y="difference to base run", xlab=" ", title="High salinity") 

recsPp

#Salinity -
salM<-subset(scenariodiff3, Scenario=="Salinity-")
recsMp<-ggplot(salM, aes(x=Concept, y=difference, fill=Concept)) + 
  geom_boxplot() + 
  geom_hline(yintercept=0) +
  ylim(-1,1) + 
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),       
      panel.background = element_blank(), axis.line = element_line(colour = "black"),        
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_blank(),        
      legend.position = "none") + 
  labs(y="difference to base run", xlab=" ", title="Low salinity") 

recsMp


#Temperature +
tempP<-subset(scenariodiff3, Scenario=="Temperature+")
rectPp<-ggplot(tempP, aes(x=Concept, y=difference, fill=Concept)) + 
  geom_boxplot() + 
  geom_hline(yintercept=0) +
  ylim(-1,1) + 
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),       
      panel.background = element_blank(), axis.line = element_line(colour = "black"),        
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_blank(),        
      legend.position = "none") + 
  labs(y="difference to base run", xlab=" ", title="High temperature") 

rectPp

#Temperature -
tempM<-subset(scenariodiff3, Scenario=="Temperature-")
rectMp<-ggplot(tempM, aes(x=Concept, y=difference, fill=Concept)) + 
  geom_boxplot() + 
  geom_hline(yintercept=0) +
  ylim(-1,1) + 
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),       
      panel.background = element_blank(), axis.line = element_line(colour = "black"),        
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_blank(),        
      legend.position = "none") + 
  labs(y="difference to base run", xlab=" ", title="Low temperature") 

rectMp

png(file="diffplotRecreation.png", units="cm", width = 20, height = 30, res=200)

plot_grid(recfPp + annotate("text", x=1, y=1, label="a) High fishing pressure", hjust=0) + labs(title=element_blank()),
          recfMp+ annotate("text", x=1, y=1, label="b) Low fishing pressure", hjust=0) + labs(title=element_blank()), 
          recePp+ annotate("text", x=1, y=1, label="c) High eutrophication", hjust=0) + labs(title=element_blank()), 
          receMp+ annotate("text", x=1, y=1, label="d) Low eutrophication", hjust=0) + labs(title=element_blank()), 
          rectPp+ annotate("text", x=1, y=1, label="e) High temperature", hjust=0) + labs(title=element_blank()), 
          rectMp+ annotate("text", x=1, y=1, label="f) Low temperature", hjust=0) + labs(title=element_blank()), 
          recsPp+ annotate("text", x=1, y=1, label="g) High salinity", hjust=0) + labs(title=element_blank()),
          recsMp+ annotate("text", x=1, y=1, label="h) Low salinity", hjust=0) + labs(title=element_blank()),
          align = 'vh',
          hjust = -1,
          nrow = 4)

dev.off()


######################
## Plot the 10 most central ecosystem components resutls per each scenario
## FOR fishery

#strip the data so that it includes only the relevant ecosystem components & stakeholder group
scenariodiff3 <- scenariodiff2[scenariodiff2$Group=="Fishery",]


#Fishing +
fishP<-subset(scenariodiff3, Scenario=="Fishing+")
fisfPp<-ggplot(fishP, aes(x=Concept, y=difference, fill=Concept)) + 
  geom_boxplot() + 
  geom_hline(yintercept=0) +
  ylim(-1,1) + 
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),       
      panel.background = element_blank(), axis.line = element_line(colour = "black"),        
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_blank(),        
      legend.position = "none") + 
  labs(y="difference to base run", xlab=" ", title="High fishing pressure")

fisfPp

#Fishing -
fishM<-subset(scenariodiff3, Scenario=="Fishing-")
fisfMp<-ggplot(fishM, aes(x=Concept, y=difference, fill=Concept)) + 
  geom_boxplot() + 
  geom_hline(yintercept=0) +
  ylim(-1,1) + 
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),       
      panel.background = element_blank(), axis.line = element_line(colour = "black"),        
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_blank(),        
      legend.position = "none") + 
  labs(y="difference to base run", xlab=" ", title="Low fishing pressure") 

fisfMp


#Eutrophication +
eutroP<-subset(scenariodiff3, Scenario=="Eutrophication+")
fisePp<-ggplot(eutroP, aes(x=Concept, y=difference, fill=Concept)) + 
  geom_boxplot() + 
  geom_hline(yintercept=0) +
  ylim(-1,1) + 
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),       
      panel.background = element_blank(), axis.line = element_line(colour = "black"),        
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_blank(),        
      legend.position = "none") + 
  labs(y="difference to base run", xlab=" ", title="High eutrophication") 

fisePp

#Eutrophication -
eutroM<-subset(scenariodiff3, Scenario=="Eutrophication-")
fiseMp<-ggplot(eutroM, aes(x=Concept, y=difference, fill=Concept)) + 
  geom_boxplot() + 
  geom_hline(yintercept=0) +
  ylim(-1,1) + 
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),       
      panel.background = element_blank(), axis.line = element_line(colour = "black"),        
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_blank(),        
      legend.position = "none") + 
  labs(y="difference to base run", xlab=" ", title="Low eutrophication") 

fiseMp

#Salinity +
salP<-subset(scenariodiff3, Scenario=="Salinity+")
fissPp<-ggplot(salP, aes(x=Concept, y=difference, fill=Concept)) + 
  geom_boxplot() + 
  geom_hline(yintercept=0) +
  ylim(-1,1) + 
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),       
      panel.background = element_blank(), axis.line = element_line(colour = "black"),       
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_blank(),        
      legend.position = "none") + 
  labs(y="difference to base run", xlab=" ", title="High salinity") 

fissPp

#Salinity -
salM<-subset(scenariodiff3, Scenario=="Salinity-")
fissMp<-ggplot(salM, aes(x=Concept, y=difference, fill=Concept)) + 
  geom_boxplot() + 
  geom_hline(yintercept=0) +
  ylim(-1,1) + 
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),       
      panel.background = element_blank(), axis.line = element_line(colour = "black"),        
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_blank(),        
      legend.position = "none") + 
  labs(y="difference to base run", xlab=" ", title="Low salinity") 

fissMp


#Temperature +
tempP<-subset(scenariodiff3, Scenario=="Temperature+")
fistPp<-ggplot(tempP, aes(x=Concept, y=difference, fill=Concept)) + 
  geom_boxplot() + 
  geom_hline(yintercept=0) +
  ylim(-1,1) + 
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),       
      panel.background = element_blank(), axis.line = element_line(colour = "black"),        
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_blank(),        
      legend.position = "none") + 
  labs(y="difference to base run", xlab=" ", title="High temperature") 

fistPp

#temp -
tempM<-subset(scenariodiff3, Scenario=="Temperature-")
fistMp<-ggplot(tempM, aes(x=Concept, y=difference, fill=Concept)) + 
  geom_boxplot() + 
  geom_hline(yintercept=0) +
  ylim(-1,1) + 
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),       
      panel.background = element_blank(), axis.line = element_line(colour = "black"),        
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_blank(),        
      legend.position = "none") + 
  labs(y="difference to base run", xlab=" ", title="Low temperature") 

fistMp

png(file="diffplotFishery.png", units="cm", width = 20, height = 30, res=200)

plot_grid(fisfPp+ annotate("text", x=1, y=1, label="a) High fishing pressure", hjust=0) + labs(title=element_blank()),
          fisfMp+ annotate("text", x=1, y=1, label="b) Low fishing pressure", hjust=0) + labs(title=element_blank()), 
          fisePp+ annotate("text", x=1, y=1, label="c) High eutrophication", hjust=0) + labs(title=element_blank()), 
          fiseMp+ annotate("text", x=1, y=1, label="d) Low eutrophication", hjust=0) + labs(title=element_blank()), 
          fistPp+ annotate("text", x=1, y=1, label="e) High temperature", hjust=0) + labs(title=element_blank()), 
          fistMp+ annotate("text", x=1, y=1, label="f) Low temperature", hjust=0) + labs(title=element_blank()), 
          fissPp+ annotate("text", x=1, y=1, label="g) High salinity", hjust=0) + labs(title=element_blank()),
          fissMp+ annotate("text", x=1, y=1, label="h) Low salinity", hjust=0) + labs(title=element_blank()),
          align = 'vh',
          hjust = -1,
          nrow = 4)

dev.off()

###############################################
# Plots for presentation

#all respondents, different scenarios

allfPp2 <- allfPp + annotate("text", x=1, y=1, label="High fishing pressure", hjust=0) + labs(title=element_blank()) + theme(axis.title.y = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank())
allfMp2 <- allfMp + annotate("text", x=1, y=1, label="Low fishing pressure", hjust=0) + labs(title=element_blank()) + theme(axis.title.y = element_blank())
allePp2 <- allePp + annotate("text", x=1, y=1, label="High eutrophication", hjust=0) + labs(title=element_blank()) + theme(axis.title.y = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank()) 
alleMp2 <- alleMp + annotate("text", x=1, y=1, label="Low eutrophication", hjust=0) + labs(title=element_blank()) + theme(axis.title.y = element_blank()) 
alltPp2 <- alltPp + annotate("text", x=1, y=1, label="High temperature", hjust=0) + labs(title=element_blank()) + theme(axis.title.y = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank()) 
alltMp2 <- alltMp + annotate("text", x=1, y=1, label="Low temperature", hjust=0) + labs(title=element_blank()) + theme(axis.title.y = element_blank())
allsPp2 <- allsPp + annotate("text", x=1, y=1, label="High salinity", hjust=0) + labs(title=element_blank())+ theme(axis.title.y = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank())
allsMp2 <- allsMp + annotate("text", x=1, y=1, label="Low salinity", hjust=0) + labs(title=element_blank()) + theme(axis.title.y = element_blank())

png(file="diffplot_h.png", units="cm", width = 30, height = 20, res=200)
ggarrange(allfPp2, allePp2, alltPp2, allsPp2, 
          allfMp2, alleMp2, alltMp2, allsMp2, ncol = 4, nrow = 2)
dev.off()

# decreased fishing scenario, different groups
allfMp2 <- allfMp + annotate("text", x=1, y=1, label="All respondents", hjust=0) + labs(title=element_blank()) + theme(axis.title.y = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank())
fisfMp2 <- fisfMp + annotate("text", x=1, y=1, label="Fishery", hjust=0) + labs(title=element_blank()) + theme(axis.title.y = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank()) 
engofMp2 <- engofMp + annotate("text", x=1, y=1, label="eNGO", hjust=0) + labs(title=element_blank()) + theme(axis.title.y = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank()) 
recfMp2 <- recfMp + annotate("text", x=1, y=1, label="Recreation", hjust=0) + labs(title=element_blank()) + theme(axis.title.y = element_blank()) 
polfMp2 <- polfMp + annotate("text", x=1, y=1, label="Policy makers", hjust=0) + labs(title=element_blank()) + theme(axis.title.y = element_blank())
scifMp2 <- scifMp + annotate("text", x=1, y=1, label="Science", hjust=0) + labs(title=element_blank())+ theme(axis.title.y = element_blank())

png(file="fishingMinusScen_h.png", units="cm", width = 30, height = 20, res=200)
ggarrange(allfMp2, fisfMp2, engofMp2,
          recfMp2, polfMp2, scifMp2, ncol = 3, nrow = 2)
dev.off()


# decreased eutrophication scenario, different groups
alleMp2 <- alleMp + annotate("text", x=1, y=1, label="All respondents", hjust=0) + labs(title=element_blank()) + theme(axis.title.y = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank())
fiseMp2 <- fiseMp + annotate("text", x=1, y=1, label="Fishery", hjust=0) + labs(title=element_blank()) + theme(axis.title.y = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank()) 
engoeMp2 <- engoeMp + annotate("text", x=1, y=1, label="eNGO", hjust=0) + labs(title=element_blank()) + theme(axis.title.y = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank()) 
receMp2 <- receMp + annotate("text", x=1, y=1, label="Recreation", hjust=0) + labs(title=element_blank()) + theme(axis.title.y = element_blank()) 
poleMp2 <- poleMp + annotate("text", x=1, y=1, label="Policy makers", hjust=0) + labs(title=element_blank()) + theme(axis.title.y = element_blank())
scieMp2 <- scieMp + annotate("text", x=1, y=1, label="Science", hjust=0) + labs(title=element_blank())+ theme(axis.title.y = element_blank())

png(file="eutroMinusScen_h.png", units="cm", width = 30, height = 20, res=200)
ggarrange(alleMp2, fiseMp2, engoeMp2,
          receMp2, poleMp2, scieMp2, ncol = 3, nrow = 2)
dev.off()



######################
## Plot the equilibrium values of the drivers in the base scenario
## for all maps

drivers<-scenarioruns[scenarioruns$Scenario=="No changes",]

drivers <- drivers[drivers$Concept %in% c("Fishing", "Eutrophication", "Temperature", "Salinity"),]
levD<- c("Fishing", "Eutrophication", "Temperature", "Salinity")
drivers$Concept <- factor(drivers$Concept, levels=levD)

# Salinity variable was defined in the interviews as "the decrease of salinity",
# therefore salinity needs to be inverted. Doing that here.
drivers$Equilibrium_value[drivers$Concept=="Salinity"] <- -1 * drivers$Equilibrium_value[drivers$Concept=="Salinity"]  

driverplot<-ggplot(drivers, aes(x=Concept, y=Equilibrium_value, fill=Concept)) + 
  geom_boxplot() + 
  geom_hline(yintercept=0) +
  ylim(-1,1) + 
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),       panel.background = element_blank(), axis.line = element_line(colour = "black"),        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_blank(),        legend.position = "none") + 
  labs(y="Equilibrium value", xlab=" ")

png(file="equilibrium.png", units="cm", width = 20, height = 20, res=200)
driverplot
dev.off()

