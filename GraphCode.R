# plot(clusterVCtr,tp,type = "b",main = "True positive vs Diff cluster numbers")
# plot(clusterVCtr,fp,type = "b",main = "False positive vs Diff cluster numbers")
# plot(clusterVCtr,tn,type = "b",main = "True Negative vs Diff cluster numbers")
# plot(clusterVCtr,fn,type = "b",main = "False Negative vs Diff cluster numbers")
# plot(clusterVCtr,timeTaken,type = "b",main = "Amount of time vs Diff cluster numbers")
# 
# plot(clusterVCtr,precision,type = "b",main = "Precision vs Diff cluster numbers")
# plot(clusterVCtr,recall,type = "b",main = "Recall vs Diff cluster numbers")
# 
# plot(recall,precision,type = "b",main = "Precision vs Recall",xlim=c(0,1), lwd=2)
# # plot(recall,precision,type = "b",main = "Precision vs Recall",xlim=c(0,1), ylim=c(0,1), lwd=2)




df <- data.frame(TP=tp,FP=fp,TN=tn,FN=fn,PRECISION=precision,RECALL=recall,TIME=timeTaken)

#df <- data.frame(TP=tp,FP=fp,TN=tn,FN=fn,PRECISION=precision,RECALL=recall)

write.table(df, file = "C:\\Users\\USER\\Documents\\R\\mergednet1WT0dot1ALLk3to10TimeHrs.tsv", sep = "\t", row.names = FALSE)

plot(thresholdVctr,tp,type = "b",main = "True positive vs Diff Threshold")


plot(clusterVCtr,tp,type = "b",main = "True positive vs Diff cluster numbers")
plot(clusterVCtr,fp,type = "b",main = "False positive vs Diff cluster numbers")
plot(clusterVCtr,tn,type = "b",main = "True Negative vs Diff cluster numbers")
plot(clusterVCtr,fn,type = "b",main = "False Negative vs Diff cluster numbers")
plot(clusterVCtr,timeTaken,type = "b",main = "Amount of time vs Diff cluster numbers")

plot(clusterVCtr,precision,type = "b",main = "Precision vs Diff cluster numbers")
plot(clusterVCtr,recall,type = "b",main = "Recall vs Diff cluster numbers")

plot(recall,precision,type = "b",main = "Precision vs Recall",xlim=c(0,1), lwd=2)
# plot(recall,precision,type = "b",main = "Precision vs Recall",xlim=c(0,1), ylim=c(0,1), lwd=2)


#================================================

# DOING FOR UNMERGED

# COLUMNS "THRESHOLD"	"TP"	"FP"	"TN"	"FN"	"PRECISION"	"RECALL"
#             1        2     3      4     5       6           7



clusterVCtr <- c(3,4,5,6,7,8,9)

thresholdVctr <- c(0,0.01,0.03,0.05,0.1,0.2)

imgfolder <- "C:\\Users\\USER\\Desktop\\r images\\"
rStoreDataPath <- "C:\\Users\\USER\\Documents\\R\\"


ncluster <- 9
tpBeforeDPI <- c()
tpAfterDPI <- c()
# before DPI
fname <- paste(rStoreDataPath,"merged result\\","Net1MergedCluster",ncluster,"TPFPsForAllThreshold.tsv",sep="")

tmp <- read.table(fname, header = T, sep = "\t")
tmp <- as.matrix(tmp)
tpBeforeDPI <- tmp[,3]


# after DPI
fname <- paste(rStoreDataPath,"merged result DPI\\","Net1MergedCluster",ncluster,"TPFPsForAllThresholdDPI.tsv",sep="")
tmp <- read.table(fname, header = T, sep = "\t")
tmp <- as.matrix(tmp)
tpAfterDPI <- tmp[,3]

m <- paste("False Positive vs Threshold (Before and After DPI for Merged,Cluster = ",ncluster,")", sep="")
jpeg(file = paste(imgfolder,"imgMK",ncluster," FPvsTH.jpeg", sep=""),width = 600,height = 600)
plot(thresholdVctr,tpBeforeDPI, type = "b",col = "red", xlab = "Threshold Values", ylab = "False Positive", main=m,cex=1.5)
lines(thresholdVctr,tpAfterDPI, type = "b",col = "blue",cex=1.5)
legend("topright", c("Before DPI","After DPI"), cex=1.5, bty="n", fill=c("red","blue"))
dev.off()

#===========================================


a <- as.integer(abs(runif(15))*100)
a <- matrix(a, nrow = 3)
barplot(a, beside = T, names.arg = c("A","B","C","D","E"))




dat <- read.table(text = "Cluster1   Cluster2   Cluster3   Cluster4   Cluster5   Cluster6    Cluster7
                  1 480 780 431 295 670 360  190
                  2 720 350 377 255 340 615  345
                  3 460 480 179 560  60 735 1260
                  4 220 240 876 789 820 100   75", header = TRUE)


barplot(as.matrix(dat), beside = T)



#================================================               ====================================
#================================================               ====================================
#================================================ TRUE POSITIVE ====================================
#================================================               ====================================

# COLUMNS "THRESHOLD"	"TP"	"FP"	"TN"	"FN"	"PRECISION"	"RECALL"
#             1        2     3      4     5       6           7



## tp(Y) vs cluster number(X) FOR different threshold in same graph
# 
#                                   BEFORE DPI


tpBeforeDPI <- matrix(0,nrow = length(thresholdVctr), ncol = length(clusterVCtr))

for(i in 1:length(clusterVCtr)){
  ncluster <- clusterVCtr[i]
  
  # before DPI
  fname <- paste(rStoreDataPath,"resultSelectedMergedFromUnMergedCluster=",ncluster,"TPFPsForAllThreshold.tsv",sep="")
  tmp <- read.table(fname, header = T, sep = "\t")
  tmp <- as.matrix(tmp)
  tpBeforeDPI[,i] <- tmp[,2]
  
}

m <- "True Positive VS Cluster Number (Before DPI)"
subT <- paste("Number of Genes = ",ngenes)

jpeg(file = paste(imageDataPath,subDirPath," ",m,".jpeg", sep=""),width = wh,height = wh)
# first plot is drawn outside the loop
par(mar=marVector, xpd=TRUE)
plot(clusterVCtr, tpBeforeDPI[1,], main = m, sub = subT, col = colVctr[1], ylim = c(0,max(tpBeforeDPI)),xlab = "Cluster Number", ylab = "True Positive", type = "b", cex=1.5,pch=pchVctr[1])

for(i in 2:length(thresholdVctr)){
  lines(clusterVCtr,tpBeforeDPI[i,], type = "b", col = colVctr[i], cex=1.5,pch=pchVctr[i])
  
}

legend("topright", inset=insetVector, legend=thresholdCharVctr, cex=1.5, bty="n",pch=pchVctr,col = colVctr, title="Threshold")

dev.off()


#=============================================
rm(list=ls())

#=============================================

## tp(Y) vs cluster number(X) FOR different threshold in same graph
#                                   AFTER DPI

# +1 in the matrix is for no clustr col at the beginning
tpAfterDPI <- matrix(0,nrow = length(thresholdVctr), ncol = length(clusterVCtr))

for(i in 1:length(clusterVCtr)){
  ncluster <- clusterVCtr[i]
  
  # after DPI
  fname <- paste(rStoreDataPath,"resultSelectedMergedFromUnMergedCluster=",ncluster,"TPFPsForAllThresholdDPI.tsv",sep="")
  tmp <- read.table(fname, header = T, sep = "\t")
  tmp <- as.matrix(tmp)
  tpAfterDPI[,i] <- tmp[,2]
  
}

m <- "True Positive VS Cluster Number (After DPI)"
subT <- paste("Number of Genes = ",ngenes)

jpeg(file = paste(imageDataPath,subDirPath," ",m,".jpeg", sep=""),width = wh,height = wh)
# first plot is drawn outside the loop
par(mar=marVector, xpd=TRUE)
plot(clusterVCtr, tpAfterDPI[1,], main = m, sub = subT, col = colVctr[1], ylim = c(0,max(tpBeforeDPI)),xlab = "Cluster Number", ylab = "True Positive", type = "b", cex=1.5,pch=pchVctr[1])

for(i in 2:length(thresholdVctr)){
  lines(clusterVCtr,tpAfterDPI[i,], type = "b", col = colVctr[i], cex=1.5,pch=pchVctr[i])
  
}

legend("topright", inset=insetVector, legend=thresholdCharVctr, cex=1.5, bty="n" ,pch=pchVctr,col = colVctr, title="Threshold")

dev.off()


#=============================================
#=============================================

# COLUMNS "THRESHOLD"	"TP"	"FP"	"TN"	"FN"	"PRECISION"	"RECALL"
#             1        2     3      4     5       6           7


## tp(Y) vs threshold(X, including no cluster) FOR different cluster number in same graph
#                                   BEFORE DPI

tpBeforeDPI <- matrix(0,nrow = length(thresholdVctr), ncol = (length(clusterVCtr)+1))

# no cluster
fname <- paste(rStoreDataPath,"resultNoClusterTPFPsForAllThreshold.tsv",sep="")
tmp <- read.table(fname, header = T, sep = "\t")
tmp <- as.matrix(tmp)
tpBeforeDPI[,1] <- tmp[,2]



for(i in 1:length(clusterVCtr)){
  ncluster <- clusterVCtr[i]
  
  # before DPI
  fname <- paste(rStoreDataPath,"resultSelectedMergedFromUnMergedCluster=",ncluster,"TPFPsForAllThreshold.tsv",sep="")
  tmp <- read.table(fname, header = T, sep = "\t")
  tmp <- as.matrix(tmp)
  tpBeforeDPI[,(i+1)] <- tmp[,2]
  
}

# thresholdVctr
# clusterVCtr
m <- "True Positive VS Threshold Value (Before DPI)"
subT <- paste("Number of Genes = ",ngenes)

jpeg(file = paste(imageDataPath,subDirPath," ",m,".jpeg", sep=""),width = wh,height = wh)
# first plot is drawn outside the loop
par(mar=marVector, xpd=TRUE)
plot(thresholdVctr, tpBeforeDPI[,1], main = m, sub = subT, col = colVctr[1], ylim = c(0,max(tpBeforeDPI)),xlab = "Threshold Value", ylab = "True Positive", type = "b", cex=1.5,pch=pchVctr[1])

for(i in 1:length(clusterVCtr)){
  lines(thresholdVctr,tpBeforeDPI[,(i+1)], type = "b", col = colVctr[(i+1)], cex=1.5,pch=pchVctr[(i+1)])
  
}

legend("topright", inset=insetVector, legend=clusterCharVCtr, cex=1.5, bty="n", title="Cluster",pch=pchVctr,col = colVctr)

dev.off()


#=============================================

## tp(Y) vs threshold(X, including no cluster) FOR different cluster number in same graph
#                                   AFTER DPI

tpAftereDPI <- matrix(0,nrow = length(thresholdVctr), ncol = (length(clusterVCtr)+1))

# no cluster
fname <- paste(rStoreDataPath,"resultNoClusterTPFPsForAllThresholdDPI.tsv",sep="")
tmp <- read.table(fname, header = T, sep = "\t")
tmp <- as.matrix(tmp)
tpAftereDPI[,1] <- tmp[,2]



for(i in 1:length(clusterVCtr)){
  ncluster <- clusterVCtr[i]
  
  # after DPI
  fname <- paste(rStoreDataPath,"resultSelectedMergedFromUnMergedCluster=",ncluster,"TPFPsForAllThresholdDPI.tsv",sep="")
  tmp <- read.table(fname, header = T, sep = "\t")
  tmp <- as.matrix(tmp)
  tpAftereDPI[,(i+1)] <- tmp[,2]
  
}

# thresholdVctr
# clusterVCtr
m <- "True Positive VS Threshold Value (After DPI)"
subT <- paste("Number of Genes = ",ngenes)


jpeg(file = paste(imageDataPath,subDirPath," ",m,".jpeg", sep=""),width = wh,height = wh)
# first plot is drawn outside the loop
par(mar=marVector, xpd=TRUE)
plot(thresholdVctr, tpAftereDPI[,1], main = m, sub = subT, col = colVctr[1], ylim = c(0,max(tpBeforeDPI)),xlab = "Threshold Value", ylab = "True Positive", type = "b", cex=1.5,pch=pchVctr[1])

for(i in 1:length(clusterVCtr)){
  lines(thresholdVctr,tpAftereDPI[,(i+1)], type = "b", col = colVctr[(i+1)], cex=1.5,pch=pchVctr[(i+1)])
  
}

legend("topright", inset=insetVector, legend=clusterCharVCtr, cex=1.5, bty="n", title="Cluster",pch=pchVctr,col = colVctr)

dev.off()


#=============================================  BAR PLOT


# barplot()
# number of groups in x axis   = number of cols in a matrix
# number of bars in each group = number of rows


#==================     TP of diff cluster in x-axis as groups      BEFORE DPI

# make tp/fp bar plot for diff cluster number as group

# COLUMNS "THRESHOLD"	"TP"	"FP"	"TN"	"FN"	"PRECISION"	"RECALL"
#             1        2     3      4     5       6           7


# bar labels
barplotClusterName <- c()
barplotClusterName[1] <- "ARC"
barplotClusterName[2] <- "NC"

for(a in 1:length(clusterVCtr)){
  numc <- clusterVCtr[a]
  barplotClusterName[a+2] <- paste("",numc)
  
}


# getting the maxY for ylimit from no cluster
fname <- paste(rStoreDataPath,"resultNoClusterTPFPsForAllThreshold.tsv",sep="")
tmp <- read.table(fname, header = T, sep = "\t")
tmp <- as.matrix(tmp)
rownames(tmp) <-NULL
colnames(tmp) <-NULL
maxY <- max(tmp[,2])  # maximum true positive value for ylimit


tempthresholdVctr <- c(0.0,0.01,0.03,0.05,0.08)

for(t in 1:length(tempthresholdVctr)){
  
  tpMatrix <- matrix(0, nrow = 2, ncol = (length(clusterVCtr)+2))   #+2 for aracne & no cluster
  
  # before DPI
  
  #aracne
  fname <- paste(rStoreDataPath,"aracne.tsv",sep="")
  tmp <- read.table(fname, header = T, sep = "\t")
  tmp <- as.matrix(tmp)
  rownames(tmp) <-NULL
  colnames(tmp) <-NULL
  # only considering true positives
  tpos <- tmp[,1]  # tp in 1st col for aracne
  tpMatrix[2,1] <- tpos   # 1 for 1st col and 2nd row to position it in the middle
  
  
  # no cluster
  fname <- paste(rStoreDataPath,"resultNoClusterTPFPsForAllThreshold.tsv",sep="")
  tmp <- read.table(fname, header = T, sep = "\t")
  tmp <- as.matrix(tmp)
  rownames(tmp) <-NULL
  colnames(tmp) <-NULL
  # only considering true positives
  tpos <- tmp[t,2]  # 1 = 1st threshold, 2 = 2nd col for TP
  tpMatrix[2,2] <- tpos    # col 2 for 2nd col and 2nd row to position it in the middle
  
  
  
  for(m in 1:length(clusterVCtr)){
    
    ncluster <- clusterVCtr[m]
    
    # before DPI
    
    # unmerged
    fname <- paste(rStoreDataPath,"resultUnMergedCluster=",ncluster,"TPFPsForAllThreshold.tsv",sep="")
    tmp <- read.table(fname, header = T, sep = "\t")
    tmp <- as.matrix(tmp)
    rownames(tmp) <-NULL
    colnames(tmp) <-NULL
    # only considering true positives
    tpos <- tmp[t,2]  
    tpMatrix[1,(m+2)] <- tpos    
    
    
    # selected merged
    fname <- paste(rStoreDataPath,"resultSelectedMergedFromUnMergedCluster=",ncluster,"TPFPsForAllThreshold.tsv",sep="")
    tmp <- read.table(fname, header = T, sep = "\t")
    tmp <- as.matrix(tmp)
    rownames(tmp) <-NULL
    colnames(tmp) <-NULL
    # only considering true positives
    tpos <- tmp[t,2]  
    tpMatrix[2,(m+2)] <- tpos    
    
  }
  
  subT <- paste("Number of Genes = ",ngenes)
  
  m <- paste("True Positive VS Cluster (Before DPI),Threshold=",tempthresholdVctr[t],sep="")
  jpeg(file = paste(imageDataPath,subDirPath," ",m,".jpeg", sep=""),width = wh,height = wh)
  
  barplot(tpMatrix, beside = T, names.arg = barplotClusterName, main = m, sub = subT,ylim = c(0,maxY) ,xlab = "Cluster Number", ylab = "True Positive", col=c("darkgreen","lightgreen"), legend = c("Unmerged","Selected Merged"))
  
  dev.off()
  
  
}




rm(list=ls())





#=============================================

#==================     TP of diff cluster in x-axis as groups      AFTER DPI

# make tp/fp bar plot for diff cluster number as group

# COLUMNS "THRESHOLD"	"TP"	"FP"	"TN"	"FN"	"PRECISION"	"RECALL"
#             1        2     3      4     5       6           7


# bar labels
barplotClusterName <- c()
barplotClusterName[1] <- "ARC"
barplotClusterName[2] <- "NC"

for(a in 1:length(clusterVCtr)){
  numc <- clusterVCtr[a]
  barplotClusterName[a+2] <- paste("",numc)
  
}


# getting the maxY for ylimit from no cluster
fname <- paste(rStoreDataPath,"resultNoClusterTPFPsForAllThreshold.tsv",sep="")
tmp <- read.table(fname, header = T, sep = "\t")
tmp <- as.matrix(tmp)
rownames(tmp) <-NULL
colnames(tmp) <-NULL
maxY <- max(tmp[,2])  # maximum true positive value for ylimit


tempthresholdVctr <- c(0.0,0.01,0.03,0.05,0.08)

for(t in 1:length(tempthresholdVctr)){
  
  tpMatrix <- matrix(0, nrow = 2, ncol = (length(clusterVCtr)+2))   #+2 for aracne & no cluster
  
  # after DPI
  
  #aracne
  fname <- paste(rStoreDataPath,"aracne.tsv",sep="")
  tmp <- read.table(fname, header = T, sep = "\t")
  tmp <- as.matrix(tmp)
  rownames(tmp) <-NULL
  colnames(tmp) <-NULL
  # only considering true positives
  tpos <- tmp[,1]  # tp in 1st col for aracne
  tpMatrix[2,1] <- tpos   # 1 for 1st col and 2nd row to position it in the middle
  
  
  # no cluster
  fname <- paste(rStoreDataPath,"resultNoClusterTPFPsForAllThresholdDPI.tsv",sep="")
  tmp <- read.table(fname, header = T, sep = "\t")
  tmp <- as.matrix(tmp)
  rownames(tmp) <-NULL
  colnames(tmp) <-NULL
  # only considering true positives
  tpos <- tmp[t,2]  # 1 = 1st threshold, 2 = 2nd col for TP
  tpMatrix[2,2] <- tpos    # col 2 for 2nd col and 2nd row to position it in the middle
  
  
  for(m in 1:length(clusterVCtr)){
    
    ncluster <- clusterVCtr[m]
    
    # before DPI
    
    # unmerged
    fname <- paste(rStoreDataPath,"resultUnMergedCluster=",ncluster,"TPFPsForAllThresholdDPI.tsv",sep="")
    tmp <- read.table(fname, header = T, sep = "\t")
    tmp <- as.matrix(tmp)
    rownames(tmp) <-NULL
    colnames(tmp) <-NULL
    # only considering true positives
    tpos <- tmp[t,2]  
    tpMatrix[1,(m+2)] <- tpos    
    
    
    # selected merged
    fname <- paste(rStoreDataPath,"resultSelectedMergedFromUnMergedCluster=",ncluster,"TPFPsForAllThresholdDPI.tsv",sep="")
    tmp <- read.table(fname, header = T, sep = "\t")
    tmp <- as.matrix(tmp)
    rownames(tmp) <-NULL
    colnames(tmp) <-NULL
    # only considering true positives
    tpos <- tmp[t,2]  
    tpMatrix[2,(m+2)] <- tpos    
    
  }
  
  subT <- paste("Number of Genes = ",ngenes)
  m <- paste("True Positive VS Cluster (After DPI),Threshold=",tempthresholdVctr[t],sep="")
  jpeg(file = paste(imageDataPath,subDirPath," ",m,".jpeg", sep=""),width = wh,height = wh)
  
  barplot(tpMatrix, beside = T, names.arg = barplotClusterName, main = m,sub = subT,ylim = c(0,maxY) ,xlab = "Cluster Number", ylab = "True Positive", col=c("darkgreen","lightgreen"), legend = c("Unmerged","Selected Merged"))
  
  dev.off()
  
  
}




rm(list=ls())



#==========================================================================================
#==========================================================================================
#==========================================================================================

#=============================================
#=============================================
#=============================================  FALSE POSITIVE
#=============================================
#==========================================================================================

#==========================================================================================


# COLUMNS "THRESHOLD"	"TP"	"FP"	"TN"	"FN"	"PRECISION"	"RECALL"
#             1        2     3      4     5       6           7



## fp(Y) vs cluster number(X) FOR different threshold in same graph
# 
#                                   BEFORE DPI


tpBeforeDPI <- matrix(0,nrow = length(thresholdVctr), ncol = length(clusterVCtr))

for(i in 1:length(clusterVCtr)){
  ncluster <- clusterVCtr[i]
  
  # before DPI
  fname <- paste(rStoreDataPath,"resultSelectedMergedFromUnMergedCluster=",ncluster,"TPFPsForAllThreshold.tsv",sep="")
  tmp <- read.table(fname, header = T, sep = "\t")
  tmp <- as.matrix(tmp)
  tpBeforeDPI[,i] <- tmp[,3]
  
}

m <- "False Positive VS Cluster Number (Before DPI)"
subT <- paste("Number of Genes = ",ngenes)

jpeg(file = paste(imageDataPath,subDirPath," ",m,".jpeg", sep=""),width = wh,height = wh)
# first plot is drawn outside the loop
par(mar=marVector, xpd=TRUE)
plot(clusterVCtr, tpBeforeDPI[1,], main = m, sub = subT, col = colVctr[1], ylim = c(0,max(tpBeforeDPI)),xlab = "Cluster Number", ylab = "False Positive", type = "b", cex=1.5,pch=pchVctr[1])

for(i in 2:length(thresholdVctr)){
  lines(clusterVCtr,tpBeforeDPI[i,], type = "b", col = colVctr[i], cex=1.5,pch=pchVctr[i])
  
}

legend("topright", inset=insetVector, legend=thresholdCharVctr, cex=1.5, bty="n",pch=pchVctr,col = colVctr, title="Threshold")

dev.off()


#=============================================
rm(list=ls())

#=============================================

## fp(Y) vs cluster number(X) FOR different threshold in same graph
#                                   AFTER DPI

# +1 in the matrix is for no clustr col at the beginning
tpAfterDPI <- matrix(0,nrow = length(thresholdVctr), ncol = length(clusterVCtr))

for(i in 1:length(clusterVCtr)){
  ncluster <- clusterVCtr[i]
  
  # after DPI
  fname <- paste(rStoreDataPath,"resultSelectedMergedFromUnMergedCluster=",ncluster,"TPFPsForAllThresholdDPI.tsv",sep="")
  tmp <- read.table(fname, header = T, sep = "\t")
  tmp <- as.matrix(tmp)
  tpAfterDPI[,i] <- tmp[,3]
  
}

m <- "False Positive VS Cluster Number (After DPI)"
subT <- paste("Number of Genes = ",ngenes)

jpeg(file = paste(imageDataPath,subDirPath," ",m,".jpeg", sep=""),width = wh,height = wh)
# first plot is drawn outside the loop
par(mar=marVector, xpd=TRUE)
plot(clusterVCtr, tpAfterDPI[1,], main = m, sub = subT, col = colVctr[1], ylim = c(0,max(tpBeforeDPI)),xlab = "Cluster Number", ylab = "False Positive", type = "b", cex=1.5,pch=pchVctr[1])

for(i in 2:length(thresholdVctr)){
  lines(clusterVCtr,tpAfterDPI[i,], type = "b", col = colVctr[i], cex=1.5,pch=pchVctr[i])
  
}

legend("topright", inset=insetVector, legend=thresholdCharVctr, cex=1.5, bty="n" ,pch=pchVctr,col = colVctr, title="Threshold")

dev.off()


#=============================================
#=============================================

# COLUMNS "THRESHOLD"	"TP"	"FP"	"TN"	"FN"	"PRECISION"	"RECALL"
#             1        2     3      4     5       6           7


## fp(Y) vs threshold(X, including no cluster) FOR different cluster number in same graph
#                                   BEFORE DPI

tpBeforeDPI <- matrix(0,nrow = length(thresholdVctr), ncol = (length(clusterVCtr)+1))

# no cluster
fname <- paste(rStoreDataPath,"resultNoClusterTPFPsForAllThreshold.tsv",sep="")
tmp <- read.table(fname, header = T, sep = "\t")
tmp <- as.matrix(tmp)
tpBeforeDPI[,1] <- tmp[,3]



for(i in 1:length(clusterVCtr)){
  ncluster <- clusterVCtr[i]
  
  # before DPI
  fname <- paste(rStoreDataPath,"resultSelectedMergedFromUnMergedCluster=",ncluster,"TPFPsForAllThreshold.tsv",sep="")
  tmp <- read.table(fname, header = T, sep = "\t")
  tmp <- as.matrix(tmp)
  tpBeforeDPI[,(i+1)] <- tmp[,3]
  
}

# thresholdVctr
# clusterVCtr
m <- "Flase Positive VS Threshold Value (Before DPI)"
subT <- paste("Number of Genes = ",ngenes)

jpeg(file = paste(imageDataPath,subDirPath," ",m,".jpeg", sep=""),width = wh,height = wh)
# first plot is drawn outside the loop
par(mar=marVector, xpd=TRUE)
plot(thresholdVctr, tpBeforeDPI[,1], main = m, sub = subT, col = colVctr[1], ylim = c(0,max(tpBeforeDPI)),xlab = "Threshold Value", ylab = "False Positive", type = "b", cex=1.5,pch=pchVctr[1])

for(i in 1:length(clusterVCtr)){
  lines(thresholdVctr,tpBeforeDPI[,(i+1)], type = "b", col = colVctr[(i+1)], cex=1.5,pch=pchVctr[(i+1)])
  
}

legend("topright", inset=insetVector, legend=clusterCharVCtr, cex=1.5, bty="n", title="Cluster",pch=pchVctr,col = colVctr)

dev.off()


#=============================================

## fp(Y) vs threshold(X, including no cluster) FOR different cluster number in same graph
#                                   AFTER DPI

tpAftereDPI <- matrix(0,nrow = length(thresholdVctr), ncol = (length(clusterVCtr)+1))

# no cluster
fname <- paste(rStoreDataPath,"resultNoClusterTPFPsForAllThresholdDPI.tsv",sep="")
tmp <- read.table(fname, header = T, sep = "\t")
tmp <- as.matrix(tmp)
tpAftereDPI[,1] <- tmp[,3]



for(i in 1:length(clusterVCtr)){
  ncluster <- clusterVCtr[i]
  
  # after DPI
  fname <- paste(rStoreDataPath,"resultSelectedMergedFromUnMergedCluster=",ncluster,"TPFPsForAllThresholdDPI.tsv",sep="")
  tmp <- read.table(fname, header = T, sep = "\t")
  tmp <- as.matrix(tmp)
  tpAftereDPI[,(i+1)] <- tmp[,3]
  
}

# thresholdVctr
# clusterVCtr
m <- "False Positive VS Threshold Value (After DPI)"
subT <- paste("Number of Genes = ",ngenes)


jpeg(file = paste(imageDataPath,subDirPath," ",m,".jpeg", sep=""),width = wh,height = wh)
# first plot is drawn outside the loop
par(mar=marVector, xpd=TRUE)
plot(thresholdVctr, tpAftereDPI[,1], main = m, sub = subT, col = colVctr[1], ylim = c(0,max(tpBeforeDPI)),xlab = "Threshold Value", ylab = "False Positive", type = "b", cex=1.5,pch=pchVctr[1])

for(i in 1:length(clusterVCtr)){
  lines(thresholdVctr,tpAftereDPI[,(i+1)], type = "b", col = colVctr[(i+1)], cex=1.5,pch=pchVctr[(i+1)])
  
}

legend("topright", inset=insetVector, legend=clusterCharVCtr, cex=1.5, bty="n", title="Cluster",pch=pchVctr,col = colVctr)

dev.off()


#=============================================  BAR PLOT


# barplot()
# number of groups in x axis   = number of cols in a matrix
# number of bars in each group = number of rows


#==================     FP of diff cluster in x-axis as groups      BEFORE DPI

# COLUMNS "THRESHOLD"	"TP"	"FP"	"TN"	"FN"	"PRECISION"	"RECALL"
#             1        2     3      4     5       6           7


# bar labels
barplotClusterName <- c()
barplotClusterName[1] <- "ARC"
barplotClusterName[2] <- "NC"

for(a in 1:length(clusterVCtr)){
  numc <- clusterVCtr[a]
  barplotClusterName[a+2] <- paste("",numc)
  
}


# getting the maxY for ylimit from no cluster
fname <- paste(rStoreDataPath,"resultNoClusterTPFPsForAllThreshold.tsv",sep="")
tmp <- read.table(fname, header = T, sep = "\t")
tmp <- as.matrix(tmp)
rownames(tmp) <-NULL
colnames(tmp) <-NULL
maxY <- max(tmp[,3])  # maximum false positive value for ylimit


tempthresholdVctr <- c(0.0,0.01,0.03,0.05,0.08)

for(t in 1:length(tempthresholdVctr)){
  
  tpMatrix <- matrix(0, nrow = 2, ncol = (length(clusterVCtr)+2))   #+2 for aracne & no cluster
  
  # before DPI
  
  #aracne
  fname <- paste(rStoreDataPath,"aracne.tsv",sep="")
  tmp <- read.table(fname, header = T, sep = "\t")
  tmp <- as.matrix(tmp)
  rownames(tmp) <-NULL
  colnames(tmp) <-NULL
  # only considering false positives
  tpos <- tmp[,2]  # fp in 2nd col for aracne
  tpMatrix[2,1] <- tpos   # 1 for 1st col and 2nd row to position it in the middle
  
  
  # no cluster
  fname <- paste(rStoreDataPath,"resultNoClusterTPFPsForAllThreshold.tsv",sep="")
  tmp <- read.table(fname, header = T, sep = "\t")
  tmp <- as.matrix(tmp)
  rownames(tmp) <-NULL
  colnames(tmp) <-NULL
  # only considering true positives
  tpos <- tmp[t,3]  # 1 = 1st threshold, 3 = 3nd col for FP
  tpMatrix[2,2] <- tpos    # col 2 for 2nd col and 2nd row to position it in the middle
  
  
  
  for(m in 1:length(clusterVCtr)){
    
    ncluster <- clusterVCtr[m]
    
    # before DPI
    
    # unmerged
    fname <- paste(rStoreDataPath,"resultUnMergedCluster=",ncluster,"TPFPsForAllThreshold.tsv",sep="")
    tmp <- read.table(fname, header = T, sep = "\t")
    tmp <- as.matrix(tmp)
    rownames(tmp) <-NULL
    colnames(tmp) <-NULL
    # only considering false positives
    tpos <- tmp[t,3]  
    tpMatrix[1,(m+2)] <- tpos    
    
    
    # selected merged
    fname <- paste(rStoreDataPath,"resultSelectedMergedFromUnMergedCluster=",ncluster,"TPFPsForAllThreshold.tsv",sep="")
    tmp <- read.table(fname, header = T, sep = "\t")
    tmp <- as.matrix(tmp)
    rownames(tmp) <-NULL
    colnames(tmp) <-NULL
    # only considering false positives
    tpos <- tmp[t,3]  
    tpMatrix[2,(m+2)] <- tpos    
    
  }
  
  subT <- paste("Number of Genes = ",ngenes)
  
  m <- paste("False Positive VS Cluster (Before DPI),Threshold=",tempthresholdVctr[t],sep="")
  jpeg(file = paste(imageDataPath,subDirPath," ",m,".jpeg", sep=""),width = wh,height = wh)
  
  barplot(tpMatrix, beside = T, names.arg = barplotClusterName, main = m, sub = subT,ylim = c(0,maxY) ,xlab = "Cluster Number", ylab = "False Positive", col=c("darkgreen","lightgreen"), legend = c("Unmerged","Selected Merged"))
  
  dev.off()
  
  
}




rm(list=ls())





#=============================================

#==================     TP of diff cluster in x-axis as groups      AFTER DPI

# make tp/fp bar plot for diff cluster number as group

# COLUMNS "THRESHOLD"	"TP"	"FP"	"TN"	"FN"	"PRECISION"	"RECALL"
#             1        2     3      4     5       6           7


# bar labels
barplotClusterName <- c()
barplotClusterName[1] <- "ARC"
barplotClusterName[2] <- "NC"

for(a in 1:length(clusterVCtr)){
  numc <- clusterVCtr[a]
  barplotClusterName[a+2] <- paste("",numc)
  
}


# getting the maxY for ylimit from no cluster
fname <- paste(rStoreDataPath,"resultNoClusterTPFPsForAllThreshold.tsv",sep="")
tmp <- read.table(fname, header = T, sep = "\t")
tmp <- as.matrix(tmp)
rownames(tmp) <-NULL
colnames(tmp) <-NULL
maxY <- max(tmp[,3])  # maximum false positive value for ylimit


tempthresholdVctr <- c(0.0,0.01,0.03,0.05,0.08)

for(t in 1:length(tempthresholdVctr)){
  
  tpMatrix <- matrix(0, nrow = 2, ncol = (length(clusterVCtr)+2))   #+2 for aracne & no cluster
  
  # after DPI
  
  #aracne
  fname <- paste(rStoreDataPath,"aracne.tsv",sep="")
  tmp <- read.table(fname, header = T, sep = "\t")
  tmp <- as.matrix(tmp)
  rownames(tmp) <-NULL
  colnames(tmp) <-NULL
  # only considering true positives
  tpos <- tmp[,2]  # fp in 2nd col for aracne
  tpMatrix[2,1] <- tpos   # 1 for 1st col and 2nd row to position it in the middle
  
  
  # no cluster
  fname <- paste(rStoreDataPath,"resultNoClusterTPFPsForAllThresholdDPI.tsv",sep="")
  tmp <- read.table(fname, header = T, sep = "\t")
  tmp <- as.matrix(tmp)
  rownames(tmp) <-NULL
  colnames(tmp) <-NULL
  
  tpos <- tmp[t,3]  # 1 = 1st threshold, 3 = 3rd col for FP
  tpMatrix[2,2] <- tpos    # col 2 for 2nd col and 2nd row to position it in the middle
  
  
  for(m in 1:length(clusterVCtr)){
    
    ncluster <- clusterVCtr[m]
    
    # before DPI
    
    # unmerged
    fname <- paste(rStoreDataPath,"resultUnMergedCluster=",ncluster,"TPFPsForAllThresholdDPI.tsv",sep="")
    tmp <- read.table(fname, header = T, sep = "\t")
    tmp <- as.matrix(tmp)
    rownames(tmp) <-NULL
    colnames(tmp) <-NULL
    # only considering false positives
    tpos <- tmp[t,3]  
    tpMatrix[1,(m+2)] <- tpos    
    
    
    # selected merged
    fname <- paste(rStoreDataPath,"resultSelectedMergedFromUnMergedCluster=",ncluster,"TPFPsForAllThresholdDPI.tsv",sep="")
    tmp <- read.table(fname, header = T, sep = "\t")
    tmp <- as.matrix(tmp)
    rownames(tmp) <-NULL
    colnames(tmp) <-NULL
    # only considering false positives
    tpos <- tmp[t,3]  
    tpMatrix[2,(m+2)] <- tpos    
    
  }
  
  subT <- paste("Number of Genes = ",ngenes)
  m <- paste("False Positive VS Cluster (After DPI),Threshold=",tempthresholdVctr[t],sep="")
  jpeg(file = paste(imageDataPath,subDirPath," ",m,".jpeg", sep=""),width = wh,height = wh)
  
  barplot(tpMatrix, beside = T, names.arg = barplotClusterName, main = m,sub = subT,ylim = c(0,maxY) ,xlab = "Cluster Number", ylab = "False Positive", col=c("darkgreen","lightgreen"), legend = c("Unmerged","Selected Merged"))
  
  dev.off()
  
  
}




rm(list=ls())





#=========================================                      ==================================
#=========================================                      ==================================
#=========================================  PRECISION VS RECALL ==================================
#=========================================                      ==================================
#=========================================                      ==================================

# COLUMNS "THRESHOLD"	"TP"	"FP"	"TN"	"FN"	"PRECISION"	"RECALL"
#             1        2     3      4     5       6           7



## pr(Y) vs rc(X, including no cluster) FOR different cluster number in same graph
#                                   BEFORE DPI

prBeforeDPI <- matrix(0,nrow = length(thresholdVctr), ncol = (length(clusterVCtr)+1))

rcBeforeDPI <- matrix(0,nrow = length(thresholdVctr), ncol = (length(clusterVCtr)+1))


# no cluster
fname <- paste(rStoreDataPath,"resultNoClusterTPFPsForAllThreshold.tsv",sep="")
tmp <- read.table(fname, header = T, sep = "\t")
tmp <- as.matrix(tmp)
prBeforeDPI[,1] <- tmp[,6]
rcBeforeDPI[,1] <- tmp[,7]


for(i in 1:length(clusterVCtr)){
  ncluster <- clusterVCtr[i]
  
  # before DPI
  fname <- paste(rStoreDataPath,"resultSelectedMergedFromUnMergedCluster=",ncluster,"TPFPsForAllThreshold.tsv",sep="")
  tmp <- read.table(fname, header = T, sep = "\t")
  tmp <- as.matrix(tmp)
  prBeforeDPI[,(i+1)] <- tmp[,6]
  rcBeforeDPI[,(i+1)] <- tmp[,7]
  
}


m <- "Precision VS Recall (Before DPI)"
subT <- paste("Number of Genes = ",ngenes)

jpeg(file = paste(imageDataPath,subDirPath," ",m,".jpeg", sep=""),width = wh,height = wh)
# first plot is drawn outside the loop
par(mar=marVector, xpd=TRUE)
plot(rcBeforeDPI[,1], prBeforeDPI[,1], main = m, sub = subT, col = colVctr[1], xlim = c(0,1), ylim = c(0,1),xlab = "Recall", ylab = "Precision", type = "b", cex=1.5,pch=pchVctr[1])

for(i in 1:length(clusterVCtr)){
  lines(rcBeforeDPI[,(i+1)],prBeforeDPI[,(i+1)], type = "b", col = colVctr[(i+1)], cex=1.5,pch=pchVctr[(i+1)])
  
}

legend("topright", inset=insetVector, legend=clusterCharVCtr, cex=1.5, bty="n", title="Cluster",pch=pchVctr,col = colVctr)

dev.off()


#================================================================================================


## pr(Y) vs rc(X, including no cluster) FOR different cluster number in same graph
#                                   AFTER DPI

prBeforeDPI <- matrix(0,nrow = length(thresholdVctr), ncol = (length(clusterVCtr)+1))

rcBeforeDPI <- matrix(0,nrow = length(thresholdVctr), ncol = (length(clusterVCtr)+1))


# no cluster
fname <- paste(rStoreDataPath,"resultNoClusterTPFPsForAllThresholdDPI.tsv",sep="")
tmp <- read.table(fname, header = T, sep = "\t")
tmp <- as.matrix(tmp)
prBeforeDPI[,1] <- tmp[,6]
rcBeforeDPI[,1] <- tmp[,7]


for(i in 1:length(clusterVCtr)){
  ncluster <- clusterVCtr[i]
  
  # before DPI
  fname <- paste(rStoreDataPath,"resultSelectedMergedFromUnMergedCluster=",ncluster,"TPFPsForAllThresholdDPI.tsv",sep="")
  tmp <- read.table(fname, header = T, sep = "\t")
  tmp <- as.matrix(tmp)
  prBeforeDPI[,(i+1)] <- tmp[,6]
  rcBeforeDPI[,(i+1)] <- tmp[,7]
  
}


m <- "Precision VS Recall (After DPI)"
subT <- paste("Number of Genes = ",ngenes)

jpeg(file = paste(imageDataPath,subDirPath," ",m,".jpeg", sep=""),width = wh,height = wh)
# first plot is drawn outside the loop
par(mar=marVector, xpd=TRUE)
plot(rcBeforeDPI[,1], prBeforeDPI[,1], main = m, sub = subT, col = colVctr[1], xlim = c(0,1), ylim = c(0,1),xlab = "Recall", ylab = "Precision", type = "b", cex=1.5,pch=pchVctr[1])

for(i in 1:length(clusterVCtr)){
  lines(rcBeforeDPI[,(i+1)],prBeforeDPI[,(i+1)], type = "b", col = colVctr[(i+1)], cex=1.5,pch=pchVctr[(i+1)])
  
}

legend("topright", inset=insetVector, legend=clusterCharVCtr, cex=1.5, bty="n", title="Cluster",pch=pchVctr,col = colVctr)

dev.off()



#==============================================================================



for(i in 1: dim(rcBeforeDPI)[2]){
  x <- rcBeforeDPI[,i]
  y <- prBeforeDPI[,i]
  area <- trapz(y,x)
  print(paste("area of ", clusterCharVCtr[i]," = ", area))
}

plot(x,y)
plot(y,x)



#====================================               pr vs rc barplot         =================

# ==================================               PR RC BEFORE DPI

# bar labels
barplotClusterName <- c()
barplotClusterName[1] <- "ARC"
barplotClusterName[2] <- "NC"

for(a in 1:length(clusterVCtr)){
  numc <- clusterVCtr[a]
  barplotClusterName[a+2] <- paste("",numc)
  
}

# thresholdVctr

for(t in 1:length(thresholdVctr)){
  
  prMatrix <- matrix(0, nrow = 2, ncol = (length(clusterVCtr)+2))   #+2 for aracne & no cluster
  
  # row 1 = precision
  # row 2 = recall
  
  #aracne
  fname <- paste(rStoreDataPath,"aracne.tsv",sep="")
  tmp <- read.table(fname, header = T, sep = "\t")
  tmp <- as.matrix(tmp)
  rownames(tmp) <-NULL
  colnames(tmp) <-NULL
  # pr
  prMatrix[1,1] <- tmp[,5]
  # rc
  prMatrix[2,1] <- tmp[,6]
  
  
  # no cluster
  
  # before DPI
  fname <- paste(rStoreDataPath,"resultNoClusterTPFPsForAllThreshold.tsv",sep="")
  tmp <- read.table(fname, header = T, sep = "\t")
  tmp <- as.matrix(tmp)
  rownames(tmp) <-NULL
  colnames(tmp) <-NULL
  # pr
  prMatrix[1,2] <- tmp[t,6]
  # rc
  prMatrix[2,2] <- tmp[t,7]
  
  
  
  for(m in 1:length(clusterVCtr)){
    
    ncluster <- clusterVCtr[m]
    
    # before DPI
    fname <- paste(rStoreDataPath,"resultSelectedMergedFromUnMergedCluster=",ncluster,"TPFPsForAllThreshold.tsv",sep="")
    tmp <- read.table(fname, header = T, sep = "\t")
    tmp <- as.matrix(tmp)
    rownames(tmp) <-NULL
    colnames(tmp) <-NULL
    # pr
    prMatrix[1,(m+2)] <- tmp[t,6]
    # rc
    prMatrix[2,(m+2)] <- tmp[t,7]
    
    
  }
  
  subT <- paste("Number of Genes = ",ngenes)
  
  m <- paste("Precisoin and Recall (Before DPI) for Threshold=",thresholdVctr[t],sep="")
  jpeg(file = paste(imageDataPath,subDirPath," ",m,".jpeg", sep=""),width = wh,height = wh)
  
  barplot(prMatrix, beside = T, names.arg = barplotClusterName, main = m, sub = subT ,ylim = c(0,1) ,xlab = "Cluster Number", ylab = "Precision & Recall", col=c("darkgreen","lightgreen"), legend = c("Precision","Recall"))
  
  dev.off()
  
  
}




rm(list=ls())



# ==================================               PR RC AFTER DPI

# bar labels
barplotClusterName <- c()
barplotClusterName[1] <- "ARC"
barplotClusterName[2] <- "NC"

for(a in 1:length(clusterVCtr)){
  numc <- clusterVCtr[a]
  barplotClusterName[a+2] <- paste("",numc)
  
}

# thresholdVctr

for(t in 1:length(thresholdVctr)){
  
  prMatrix <- matrix(0, nrow = 2, ncol = (length(clusterVCtr)+2))   #+2 for aracne & no cluster
  
  # row 1 = precision
  # row 2 = recall
  
  #aracne
  fname <- paste(rStoreDataPath,"aracne.tsv",sep="")
  tmp <- read.table(fname, header = T, sep = "\t")
  tmp <- as.matrix(tmp)
  rownames(tmp) <-NULL
  colnames(tmp) <-NULL
  # pr
  prMatrix[1,1] <- tmp[,5]
  # rc
  prMatrix[2,1] <- tmp[,6]
  
  
  # no cluster
  
  # after DPI
  fname <- paste(rStoreDataPath,"resultNoClusterTPFPsForAllThresholdDPI.tsv",sep="")
  tmp <- read.table(fname, header = T, sep = "\t")
  tmp <- as.matrix(tmp)
  rownames(tmp) <-NULL
  colnames(tmp) <-NULL
  # pr
  prMatrix[1,2] <- tmp[t,6]
  # rc
  prMatrix[2,2] <- tmp[t,7]
  
  
  
  for(m in 1:length(clusterVCtr)){
    
    ncluster <- clusterVCtr[m]
    
    # after DPI
    fname <- paste(rStoreDataPath,"resultSelectedMergedFromUnMergedCluster=",ncluster,"TPFPsForAllThresholdDPI.tsv",sep="")
    tmp <- read.table(fname, header = T, sep = "\t")
    tmp <- as.matrix(tmp)
    rownames(tmp) <-NULL
    colnames(tmp) <-NULL
    # pr
    prMatrix[1,(m+2)] <- tmp[t,6]
    # rc
    prMatrix[2,(m+2)] <- tmp[t,7]
    
    
  }
  
  subT <- paste("Number of Genes = ",ngenes)
  
  m <- paste("Precisoin and Recall (After DPI) for Threshold=",thresholdVctr[t],sep="")
  jpeg(file = paste(imageDataPath,subDirPath," ",m,".jpeg", sep=""),width = wh,height = wh)
  
  barplot(prMatrix, beside = T, names.arg = barplotClusterName, main = m, sub = subT ,ylim = c(0,1) ,xlab = "Cluster Number", ylab = "Precision & Recall", col=c("darkgreen","lightgreen"), legend = c("Precision","Recall"))
  
  dev.off()
  
  
}




rm(list=ls())


