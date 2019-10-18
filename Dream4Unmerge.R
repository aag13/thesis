
# UNMERGED

geneConnectionMatFileName <- paste(rStoreDataPath,"result",dataNumber,"UnMergedCluster=",sep="")


timeTaken <- c()

for(m in 1:length(clusterVCtr)){
  
  start.time <- as.numeric(Sys.time())
  
  tp <- c()
  fp <- c()
  tn <- c()
  fn <- c()
  precision <- c()
  recall <- c()
  
  ncluster <- clusterVCtr[m]
  
  # seed is used so that for k=a value will always cluster the genes in same way for all
  # merge, selected merge, unmerge
  set.seed(ncluster)
  
  kCluster <- kmeans(tempdata, centers = ncluster,iter.max = clusterMaxIter, nstart = 20, algorithm = clusterAlgo)
  print(paste("calculated the cluster with k = ",ncluster))
  
  #==================== group genes of same cluster into an accessible list
  
  clusterList <- list()
  
  for(i in 1:ncluster){
    clusterList[[i]] <- which(kCluster$cluster == i)
    
  }
  
  # number of genes of a cluster
  geneNumberClusterList <- list()
  
  # gene names of a cluster
  geneNameClusterList <- list()
  
  #data matrix of genes of a cluster
  tempMIClusterList <- list()
  
  #single entropy of genes of a cluster
  E1MatrixClusterList <- list()
  
  for(i in 1:ncluster){
    ngenesInCluster <- length(clusterList[[i]])
    geneNumberClusterList[[i]] <- ngenesInCluster #assigned to a list
    
    tempMI <- matrix(0, nrow = nexp, ncol = ngenesInCluster)
    cgnames <- rep("",ngenesInCluster) #gene names in a cluster
    
    #create the matrix of genes in this cluster
    for(j in 1:ngenesInCluster){
      gnum <- clusterList[[i]][[j]]
      tempMI[,j] <- matA[,gnum]
      cgnames[j] <- geneNames[gnum]
      
    }
    
    geneNameClusterList[[i]] <- cgnames
    tempMIClusterList[[i]] <- tempMI
    
    dtempMI <- discretize(tempMI)
    
    E1Matrix <- matrix(0,nrow = ngenesInCluster, ncol = 1)
    for(e in 1:ngenesInCluster){
      E1Matrix[e,1] <- entropy(dtempMI[,e])
      
    }
    E1MatrixClusterList[[i]] <- E1Matrix
    
    
  }
  
  
  # contains genes and connectoins in matrix format of individual cluster
  connectionList <- list()
  
  for(nofc in 1: ncluster){
    vc <- findConnectionSingleCluster(nofc)
    connectionList[[nofc]] <- vc
    
  }
  
  #======  Function to create single gene interaction matrix by merging genes from cluster
  
  geneConnectionMat <- matrix(0,nrow = ngenes, ncol = ngenes)
  
  rownames(geneConnectionMat) <- geneNames
  colnames(geneConnectionMat) <- geneNames
  
  for(i in 1:ncluster){
    #for each cluster
    clstr <- connectionList[[i]] # connection matrix of this cluster
    clstrGeneNames <- rownames(clstr) # gene names of this cluster
    nofg <- length(clstrGeneNames)
    for(j in 1:nofg){
      #for each row
      rgene <- clstrGeneNames[j]    #name of gene of cluster matrix of row
      matGeneNoR <- which(geneNames == rgene)    #row gene number in the main matrix
      for(k in 1:nofg){
        cgene <- clstrGeneNames[k]    #name of gene of cluster matrix of col
        matGeneNoC <- which(geneNames == cgene)    #col gene number in the main matrix
        
        # put 1 in th row, th col
        #geneConnectionMat[matGeneNoR,matGeneNoC] <- 1
        # instead of 1 just keep the value and use threshold later
        
        geneConnectionMat[matGeneNoR,matGeneNoC] <- clstr[j,k]
        #so the main connection matrix now has entropy reduction values
        
        
      }
      
    }
    
  }
  
  # done with individual cluster connection
  
  # save the connection matrix of this cluster
  df <- data.frame(geneConnectionMat)
  fname <- paste(geneConnectionMatFileName,ncluster,"ConnectoinMatrix.tsv",sep="")
  write.table(df, file = fname, sep = "\t")
  
  
  
  # for each threshold find TPs,FPs......
  
  for(n in 1:length(thresholdVctr)){
    
    threshold <- thresholdVctr[n]
    
    tempGeneConnMat <- geneConnectionMat
    
    for(i in 1:ngenes){
      for(j in 1:ngenes){
        if(abs(tempGeneConnMat[i,j]) <= threshold){
          if(abs(tempGeneConnMat[j,i]) <= threshold){
            
            tempGeneConnMat[i,j] <- 0
            tempGeneConnMat[j,i] <- 0
            
          }else{
            tempGeneConnMat[i,j] <- 1
            tempGeneConnMat[j,i] <- 1
          }
        }else{
          tempGeneConnMat[i,j] <- 1
          tempGeneConnMat[j,i] <- 1
          
        }
        
      }
    }
    
    #=============================== TP, FP ,TN, FN
    
    ctp <- 0
    cfp <- 0
    cfn <- 0
    ctn <- 0
    
    for(e in 1:ngenes){
      #for each row
      for(f in 1:ngenes){
        #for each col
        if((tempGeneConnMat[e,f]==1) & (gsConnectionMat[e,f]==1)){
          #TP
          ctp <- ctp + 1
        }else if((tempGeneConnMat[e,f]==1) & (gsConnectionMat[e,f]==0)){
          #FP
          cfp <- cfp + 1
          
        }else if((tempGeneConnMat[e,f]==0) & (gsConnectionMat[e,f]==1)){
          #FN
          cfn <- cfn + 1
          
        }else{
          #TN
          ctn <- ctn + 1
          
        }
        
      }
    }
    
    cpr <- ctp/(ctp+cfp)
    crc <- ctp/(ctp+cfn)
    
    tp[n] <- ctp
    fp[n] <- cfp
    tn[n] <- ctn
    fn[n] <- cfn
    
    precision[n] <- cpr
    recall[n] <- crc
    
    print(paste("done for cluster = ",ncluster," with threshold = ", threshold))
    
    
  } # for each threshold ends
  
  end.time <- as.numeric(Sys.time())
  time.taken <- end.time - start.time
  
  timeTaken[m] <- time.taken
  print(paste("UNMERGED VERSION...........time taken for cluster = ",ncluster, " is ",time.taken," sec"))
  
  df <- data.frame(THRESHOLD=thresholdVctr,TP=tp,FP=fp,TN=tn,FN=fn,PRECISION=precision,
                   RECALL=recall)
  
  fname <- paste(geneConnectionMatFileName,ncluster,"TPFPsForAllThreshold.tsv",sep="")
  write.table(df, file = fname, sep = "\t", row.names = FALSE)
  
  
} # for loops ends for diff cluster

#==============================================
# store the time needed for each cluster
df <- data.frame(CLUSTER=clusterVCtr ,TIME=timeTaken)
fname <- paste(geneConnectionMatFileName,"TimeInSec.tsv",sep="")
write.table(df, file = fname, sep = "\t", row.names = FALSE)


#==============================================

