require("infotheo")
require("minet")


findConnectionSingleCluster <- function(cid){
  
  clusterId <- cid      #cluster no.
  ngenesInCluster <-  geneNumberClusterList[[clusterId]]     #no of genes in this cluster
  
  cgnames <- geneNameClusterList[[clusterId]]
  
  tempMI <- tempMIClusterList[[clusterId]]
  
  
  dtempMI <- discretize(tempMI)
  
  #=============================== mutual info. AND distance matrix 
  
  MIMatrix <- matrix(0,nrow = ngenesInCluster, ncol = ngenesInCluster)
  NMIMatrix <- matrix(0,nrow = ngenesInCluster, ncol = ngenesInCluster)
  
  DMatrix <- matrix(0,nrow = ngenesInCluster, ncol = ngenesInCluster)
  
  
  for(i in 1:ngenesInCluster){
    for(j in 1:ngenesInCluster){
      MIMatrix[i,j] <- mutinformation(dtempMI[,i], dtempMI[,j])
      NMIMatrix[i,j] <- sqrt(1-exp(-2*MIMatrix[i,j])) #linfoot
      DMatrix[i,j] <- exp(-NMIMatrix[i,j])
    }
    
  }
  
  #============================== single entropy
  
  E1Matrix <- E1MatrixClusterList[[clusterId]]
  
  #============================  conditional entropy H(X|Y)  
  
  CEMatrix <- matrix(0,nrow = ngenesInCluster, ncol = ngenesInCluster)
  
  for(i in 1:ngenesInCluster){
    for(j in 1:ngenesInCluster){
      CEMatrix[i,j] <- condentropy(dtempMI[,i],dtempMI[,j])
    }
    
  }
  
  
  # TO DO->  rel_entr_reduc_2(m,n) = MI(m,n,1)/H1(m);
  ERCE2Matrix <- matrix(0,nrow = ngenesInCluster, ncol = ngenesInCluster)
  
  for(i in 1:ngenesInCluster){
    for(j in 1:ngenesInCluster){
      ERCE2Matrix[i,j] <- (MIMatrix[i,j]) / E1Matrix[i]
    }
    
  }
  
  
  #============================ joint entropy H(X,Y)
  
  # JE2Matrix <- matrix(0,nrow = ngenesInCluster, ncol = ngenesInCluster)
  # 
  # for(i in 1:ngenesInCluster){
  #   for(j in 1:ngenesInCluster){
  #     JE2Matrix[i,j] <- E1Matrix[i] + CEMatrix[j,i]
  #   }
  #   
  # }
  
  
  #===================================   entropy reduction ERT 2
  
  
  # #this only considers the minimum for a given gene
  # diagMatrix <- diag(ngenesInCluster)
  # CEMatrix <- CEMatrix + diagMatrix*(10^3)
  # 
  # ERTMatrix <- matrix(0,nrow = ngenesInCluster, ncol = 1)
  # for(i in 1:ngenesInCluster){
  #   minval <- min(CEMatrix[i,])
  #   pos <- which(CEMatrix[i,] == minval)
  #   ERTMatrix[i,1] <- pos
  #   
  # }
  
  #considers all the genes where cond entropy is less than individual entropy
  diagMatrix <- diag(ngenesInCluster)
  CEMatrix <- CEMatrix + diagMatrix*(10^3)
  
  ERTMatrix <- matrix(0,nrow = ngenesInCluster, ncol = ngenesInCluster)
  for(i in 1:ngenesInCluster){
    for(j in 1:ngenesInCluster){
      if(CEMatrix[i,j] < E1Matrix[i,1]){
        # dependent
        ERTMatrix[i,j] <- 1
      }
      
    }
    
  }
  
  
  
  #===================================   connection matrix
  
  connection <- matrix(0,nrow = ngenesInCluster, ncol = ngenesInCluster)
  
  # # entropy reduction 2 var with only ONE min
  # for(i in 1:ngenesInCluster){
  #   connection[i,ERTMatrix[i]] <- ERCE2Matrix[i,ERTMatrix[i]]
  #   
  # }
  
  
  # entropy reduction 2 var with all mins
  for(i in 1:ngenesInCluster){
    
    for(j in 1:ngenesInCluster){
      if(ERTMatrix[i,j] == 1){
        connection[i,j] <- ERCE2Matrix[i,j]
      }
      
    }
    
  }
  
  
  # entropy reduction 3 var
  # for(i in 1:ngenesInCluster){
  #   connection[i,ERT3Matrix[i,2]] <- ERCE3Matrix[i,ERT3Matrix[i,1],ERT3Matrix[i,2]]
  #   
  # }
  
  
  
  #===================================   redefine connection using threshold
  
  
  # for(i in 1:ngenesInCluster){
  #   for(j in 1:ngenesInCluster){
  #     if(abs(connection[i,j]) < threshold){
  #       connection[i,j] <- 0
  #       
  #     }
  #     
  #   }
  #   
  # }
  
  rownames(connection) <- cgnames
  colnames(connection) <- cgnames
  
  return(connection)
  
}
#======================   FUNCTION ENDS

findConnectionTwoCluster <- function(cid1, cid2){
  cluster1 <- cid1
  cluster2 <- cid2
  
  ngenesInCluster1 <-  geneNumberClusterList[[cluster1]]     #no of genes in cluster 1
  ngenesInCluster2 <-  geneNumberClusterList[[cluster2]]    #no of genes in cluster 2
  
  
  # cluster 1
  tempMI1 <- tempMIClusterList[[cluster1]]
  c1gnames <- geneNameClusterList[[cluster1]] #gene names in a cluster 1
  
  # cluster 2
  tempMI2 <- tempMIClusterList[[cluster2]]
  c2gnames <- geneNameClusterList[[cluster2]] #gene names in a cluster 2
  
  
  dtempMI1 <- discretize(tempMI1)
  dtempMI2 <- discretize(tempMI2)
  
  #=============================== mutual info. AND distance matrix 
  
  #cluster 1 I(x,y)
  cluster1MIMatrix <- matrix(0,nrow = ngenesInCluster1, ncol = ngenesInCluster2)
  cluster1NMIMatrix <- matrix(0,nrow = ngenesInCluster1, ncol = ngenesInCluster2)
  cluster1DMatrix <- matrix(0,nrow = ngenesInCluster1, ncol = ngenesInCluster2)
  
  for(i in 1:ngenesInCluster1){
    for(j in 1:ngenesInCluster2){
      cluster1MIMatrix[i,j] <- mutinformation(dtempMI1[,i], dtempMI2[,j])
      cluster1NMIMatrix[i,j] <- sqrt(1-exp(-2*cluster1MIMatrix[i,j])) #linfoot
      cluster1DMatrix[i,j] <- exp(-cluster1NMIMatrix[i,j])
    }
    
  }
  
  
  #cluster 2 I(y,x)
  cluster2MIMatrix <- matrix(0,nrow = ngenesInCluster2, ncol = ngenesInCluster1)
  cluster2NMIMatrix <- matrix(0,nrow = ngenesInCluster2, ncol = ngenesInCluster1)
  cluster2DMatrix <- matrix(0,nrow = ngenesInCluster2, ncol = ngenesInCluster1)
  
  for(i in 1:ngenesInCluster2){
    for(j in 1:ngenesInCluster1){
      cluster2MIMatrix[i,j] <- mutinformation(dtempMI2[,i], dtempMI1[,j])
      cluster2NMIMatrix[i,j] <- sqrt(1-exp(-2*cluster2MIMatrix[i,j])) #linfoot
      cluster2DMatrix[i,j] <- exp(-cluster2NMIMatrix[i,j])
    }
    
  }
  
  
  
  #============================== single entropy
  
  # cluster 1
  cluster1E1Matrix <-  E1MatrixClusterList[[cluster1]]
  
  
  # cluster 2
  cluster2E1Matrix <- E1MatrixClusterList[[cluster2]]
  
  
  #============================  conditional entropy H(X|Y)  
  
  
  # cluster 1  H(X|Y) 
  
  cluster1CEMatrix <- matrix(0,nrow = ngenesInCluster1, ncol = ngenesInCluster2)
  
  for(i in 1:ngenesInCluster1){
    for(j in 1:ngenesInCluster2){
      cluster1CEMatrix[i,j] <- condentropy(dtempMI1[,i],dtempMI2[,j])
    }
    
  }
  
  
  # TO DO->  rel_entr_reduc_2(m,n) = MI(m,n,1)/H1(m);
  cluster1ERCE2Matrix <- matrix(0,nrow = ngenesInCluster1, ncol = ngenesInCluster2)
  
  for(i in 1:ngenesInCluster1){
    for(j in 1:ngenesInCluster2){
      cluster1ERCE2Matrix[i,j] <- (cluster1MIMatrix[i,j]) / cluster1E1Matrix[i]
    }
    
  }
  
  
  
  
  # cluster 2  H(Y|X) 
  
  cluster2CEMatrix <- matrix(0,nrow = ngenesInCluster2, ncol = ngenesInCluster1)
  
  for(i in 1:ngenesInCluster2){
    for(j in 1:ngenesInCluster1){
      cluster2CEMatrix[i,j] <- condentropy(dtempMI2[,i],dtempMI1[,j])
    }
    
  }
  
  
  # TO DO->  rel_entr_reduc_2(m,n) = MI(m,n,1)/H1(m);
  cluster2ERCE2Matrix <- matrix(0,nrow = ngenesInCluster2, ncol = ngenesInCluster1)
  
  for(i in 1:ngenesInCluster2){
    for(j in 1:ngenesInCluster1){
      cluster2ERCE2Matrix[i,j] <- (cluster2MIMatrix[i,j]) / cluster2E1Matrix[i]
    }
    
  }
  
  
  #===================================   entropy reduction ERT 2
  
  
  # cluster 1
  #considers all the genes where cond entropy is less than individual entropy
  # cluster1diagMatrix <- diag(1, nrow = ngenesInCluster1, ncol = ngenesInCluster2)
  # cluster1CEMatrix <- cluster1CEMatrix + cluster1diagMatrix*(10^3)
  
  cluster1ERTMatrix <- matrix(0,nrow = ngenesInCluster1, ncol = ngenesInCluster2)
  for(i in 1:ngenesInCluster1){
    for(j in 1:ngenesInCluster2){
      
      if(cluster1CEMatrix[i,j] < cluster1E1Matrix[i,1]){
        # dependent
        cluster1ERTMatrix[i,j] <- 1
      }
      
    }
    
  }
  
  # cluster 2
  #considers all the genes where cond entropy is less than individual entropy
  # cluster2diagMatrix <- diag(1, nrow = ngenesInCluster2, ncol = ngenesInCluster1)
  # cluster2CEMatrix <- cluster2CEMatrix + cluster2diagMatrix*(10^3)
  
  cluster2ERTMatrix <- matrix(0,nrow = ngenesInCluster2, ncol = ngenesInCluster1)
  for(i in 1:ngenesInCluster2){
    for(j in 1:ngenesInCluster1){
      
      if(cluster2CEMatrix[i,j] < cluster2E1Matrix[i,1]){
        # dependent
        cluster2ERTMatrix[i,j] <- 1
      }
      
    }
    
  }
  
  
  #===================================   connection matrix
  
  # cluster 1
  cluster1connection <- matrix(0,nrow = ngenesInCluster1, ncol = ngenesInCluster2)
  
  # entropy reduction 2 var with all mins
  for(i in 1:ngenesInCluster1){
    
    for(j in 1:ngenesInCluster2){
      if(cluster1ERTMatrix[i,j] == 1){
        cluster1connection[i,j] <- cluster1ERCE2Matrix[i,j]
      }
      
    }
    
  }
  
  
  
  # cluster 2
  cluster2connection <- matrix(0,nrow = ngenesInCluster2, ncol = ngenesInCluster1)
  
  # entropy reduction 2 var with all mins
  for(i in 1:ngenesInCluster2){
    
    for(j in 1:ngenesInCluster1){
      if(cluster2ERTMatrix[i,j] == 1){
        cluster2connection[i,j] <- cluster2ERCE2Matrix[i,j]
      }
      
    }
    
  }
  
  
  
  connection <- matrix(0,nrow = ngenesInCluster1, ncol = ngenesInCluster2)
  
  for(i in 1:ngenesInCluster1){
    
    for(j in 1:ngenesInCluster2){
      cluster1val <- cluster1connection[i,j]
      cluster2val <- cluster2connection[j,i]
      
      if(cluster1val > cluster2val){
        connection[i,j] <- cluster1val
        
      }else{
        connection[i,j] <- cluster2val
        
      }
      
      
    }
    
  }
  
  
  
  #===================================   redefine connection using threshold
  # threshold is being checked after final matrix created
  
  # for(i in 1:ngenesInCluster1){
  #   for(j in 1:ngenesInCluster2){
  #     if(abs(connection[i,j]) < threshold){
  #       connection[i,j] <- 0
  # 
  #     }
  # 
  #   }
  # 
  # }
  
  rownames(connection) <- c1gnames
  colnames(connection) <- c2gnames
  
  return(connection)
  
  
}
#======================   FUNCTION ENDS

#===================================   path to folders

# sub dir path after root data path

# subDirPath <- "net 1"
# subDirPath <- "net 3"
subDirPath <- "dream 4"
# subDirPath <- "gnw"

# subDirPath <- ""

netDataPath <- paste("C:\\Users\\USER\\Desktop\\r data\\",subDirPath,"\\",sep="")
#folder where to save genereated tsvs graphs
rStoreDataPath <- paste("C:\\Users\\USER\\Desktop\\data\\",subDirPath,"\\",sep="")

imageDataPath <- paste("C:\\Users\\USER\\Desktop\\image\\",subDirPath,"\\",sep="")

dataNumber <- ""

# data file name = data.tsv
# gold standard file name = gs.tsv

dataFileName <- paste("data",dataNumber,".tsv",sep="")
goldStandardFileName <- paste("gs",dataNumber,".tsv",sep="")

wh <- 400
marVector <- c(5.1, 4.1, 4.1, 7.1)
insetVector <- c(-0.4,0)

#===================      Input Data

fname <- paste(netDataPath,dataFileName,sep="")
matA <- read.table(fname, header = T, sep = "\t")

geneNames <- colnames(matA)

tempdata <- t(matA)         #rows = genes & cols = experiments
mydata <- tempdata

tempdata <- scale(tempdata)

ngenes <- nrow(tempdata)
nexp <- ncol(tempdata)
#=================== ENDS

#===================      Gold Standard

fname <- paste(netDataPath,goldStandardFileName,sep="")
gsData <- read.table(fname, header = F, sep = "\t")


gsMat <- as.matrix(gsData)

gsMat <- gsMat[(which(gsMat[,3] == "1")),] #get only pairs between which edge exists


rownames(gsMat) <- NULL
colnames(gsMat) <- NULL

gsConnectionMat <- matrix(0, nrow = ngenes, ncol = ngenes)
rownames(gsConnectionMat) <- geneNames
colnames(gsConnectionMat) <- geneNames

tempRowNames <- rownames(gsConnectionMat)

for(i in 1:dim(gsMat)[1]){
  fgene <- gsMat[i,1] #first gene
  sgene <- gsMat[i,2] #second gene
  g1 <- which(tempRowNames == fgene)
  g2 <- which(tempRowNames == sgene)
  
  gsConnectionMat[g1,g2] <- 1
  gsConnectionMat[g2,g1] <- 1
  
}

#=================== ENDS

clusterVCtr <- c()

clusterCharVCtr <- c()

if(subDirPath == "net 1"){
  clusterVCtr <- c(3,5,10,20,50,100) #for net 1
  clusterCharVCtr <- c("NC","3","5","10","20","50","100")
  
}else if(subDirPath == "net 3"){
  clusterVCtr <- c(10,20,50,100) #for net 3 
  clusterCharVCtr <- c("NC","10","20","50","100")
  
}else if(subDirPath == "dream 4"){
  clusterVCtr <- c(2,5,8,10,15,20) #for dream 4
  clusterCharVCtr <- c("NC","2","5","8","10","15","20")
  
}else if(subDirPath == "gnw"){
  clusterVCtr <- c(3,5,10,20,50,100)   #for gnw
  clusterCharVCtr <- c("NC","3","5","10","20","50","100")
  
}



thresholdVctr <- c(0,0.01,0.03,0.05,0.08,0.1,0.2)

thresholdCharVctr <- c("0","0.01","0.03","0.05","0.08","0.1","0.2")

pchVctr <- c(15,16,17,18,19,20,8)

colVctr <- c("red","green","dark green","orange","purple","black","blue")

clusterMaxIter <- 300
clusterAlgo <- "Lloyd"


# resultSelectedMergedFromUnMergedCluster=2TPFPsForAllThreshold.tsv   ==> selected merged from unmerged
# resultSelectedMergedFromUnMergedCluster=TimeIn.tsv       ==> selected merged from unmerged time in sec
# resultSelectedMergedFromUnMergedCluster=2ConnectoinMatrix.tsv


# resultUnMergedCluster=2TPFPsForAllThreshold.tsv   ==> unmerged tp/fp
# resultUnMergedCluster=TimeIn.tsv                 ==> unmerged time in sec
# resultUnMergedCluster=2ConnectoinMatrix.tsv


# resultNoClusterTPFPsForAllThreshold.tsv           ==> no cluster tp/fp
# resultNoClusterTimeIn.tsv                         ==> no cluster time in sec
# resultNoClusterConnectoinMatrix.tsv





