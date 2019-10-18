# #======================   calculating MI matrix of a dataset

# dataMat <- as.matrix(matA)
# dataMat <- scale(dataMat)
# dataMat <- discretize(dataMat)
# minetMIM1 <- build.mim(dataMat)
# 
# minetMIM2 <- build.mim(matA)
# 
# all(minetMIM2 == minetMIM3)
# 
# View(minetMIM2)
# View(minetMIM3)

# MIMat <- matrix(0,nrow = ngenes, ncol = ngenes)
# 
# rownames(MIMat) <- geneNames
# colnames(MIMat) <- geneNames
# 
# #start.time <- Sys.time()
# for(a in 1:ngenes){
#   for(b in 1:ngenes){
#     MIMat[a,b] <- mutinformation(dataMat[,a],dataMat[,b])
# 
#   }
# 
# }


# mim <- round(build.mim(matA), 5)
# df <- data.frame(mim)
# fname <- paste(rStoreDataPath,"mim.tsv",sep="")
# write.table(df, file = fname, sep = "\t")
# 
# rm(list=ls())

# end.time <- Sys.time()
# time.taken <- end.time - start.time
# time.taken
# 
# df <- data.frame(MIMat)
# write.table(df, file = "C:\\Users\\USER\\Desktop\\dream4 data\\Dream4-5MIMatrix.tsv", sep = "\t",row.names = F)
# 
# # Dream4-1MIMatrix.tsv
# 
# rm(list=ls())

#=================================================


# getting net 1 MI matrix

fname <- paste(rStoreDataPath,"mim.tsv",sep="")
MIMat <- read.table(fname, header = T, sep = "\t")
MIMat <- as.matrix(MIMat)

# No CLuster DPI

start.time <- as.numeric(Sys.time())

tp <- c()
fp <- c()
tn <- c()
fn <- c()
precision <- c()
recall <- c()

# resultNoClusterConnectoinMatrix.tsv

fname <- paste(rStoreDataPath,"resultNoClusterConnectoinMatrix.tsv",sep="")
tmpgconnMat <- read.table(fname, header = T, sep = "\t")

geneConnectionMat <- as.matrix(tmpgconnMat)


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
  
  dpiTempGeneConnMat <- tempGeneConnMat
  
  # now run DPI
  for(a in 3:ngenes){
    for(b in 2:(a-1)){
      for(c in 1:(b-1)){
        # do DPI only when a,b / b,c / a,c is 1 in geneConnectionMat
        if(dpiTempGeneConnMat[a,b] == 1){
          if(dpiTempGeneConnMat[b,c] == 1){
            if(dpiTempGeneConnMat[a,c] == 1){
              # all three genes have pair connection
              ab <- MIMat[a,b]
              bc <- MIMat[b,c]
              ac <- MIMat[a,c]
              
              if(ab <= min(bc,ac)){
                # ab is indirect, remove it
                tempGeneConnMat[a,b] = 0
                tempGeneConnMat[b,a] = 0
                
              }else if(bc <= min(ab,ac)){
                # bc is indirect, remove it
                tempGeneConnMat[b,c] = 0
                tempGeneConnMat[c,b] = 0
                
              }else if(ac <= min(ab,bc)){
                # ac is indirect, remove it
                tempGeneConnMat[a,c] = 0
                tempGeneConnMat[c,a] = 0
                
              }
              
            }
          }
        }
        
        
        
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
  
  print(paste("done with threshold = ", threshold))
  
  
} # for each threshold ends

end.time <- as.numeric(Sys.time())
time.taken <- end.time - start.time


print(paste("DPI NO CLUSTER VERSION......time taken  = ",time.taken," sec"))

df <- data.frame(THRESHOLD=thresholdVctr,TP=tp,FP=fp,TN=tn,FN=fn,PRECISION=precision,
                 RECALL=recall)
fname <- paste(rStoreDataPath ,"resultNoClusterTPFPsForAllThresholdDPI.tsv",sep="")
write.table(df, file = fname, sep = "\t", row.names = FALSE)


# store the time needed for each cluster
df <- data.frame(TIME=time.taken)
fname <- paste(rStoreDataPath ,"resultNoClusterTimeInSecDPI.tsv",sep="")
write.table(df, file = fname, sep = "\t", row.names = FALSE)



#============================================

# SELECTED MERGED FROM UNMERGED DPI

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
  
  # resultSelectedMergedFromUnMergedCluster=2ConnectoinMatrix.tsv
  
  fname <- paste(rStoreDataPath,"resultSelectedMergedFromUnMergedCluster=", ncluster,"ConnectoinMatrix.tsv",sep="")
  tmpgconnMat <- read.table(fname, header = T, sep = "\t")
  
  geneConnectionMat <- as.matrix(tmpgconnMat)
  
  
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
    
    dpiTempGeneConnMat <- tempGeneConnMat
    
    # now run DPI
    for(a in 3:ngenes){
      for(b in 2:(a-1)){
        for(c in 1:(b-1)){
          # do DPI only when a,b / b,c / a,c is 1 in geneConnectionMat
          if(dpiTempGeneConnMat[a,b] == 1){
            if(dpiTempGeneConnMat[b,c] == 1){
              if(dpiTempGeneConnMat[a,c] == 1){
                # all three genes have pair connection
                ab <- MIMat[a,b]
                bc <- MIMat[b,c]
                ac <- MIMat[a,c]
                
                if(ab <= min(bc,ac)){
                  # ab is indirect, remove it
                  tempGeneConnMat[a,b] = 0
                  tempGeneConnMat[b,a] = 0
                  
                }else if(bc <= min(ab,ac)){
                  # bc is indirect, remove it
                  tempGeneConnMat[b,c] = 0
                  tempGeneConnMat[c,b] = 0
                  
                }else if(ac <= min(ab,bc)){
                  # ac is indirect, remove it
                  tempGeneConnMat[a,c] = 0
                  tempGeneConnMat[c,a] = 0
                  
                }
                
              }
            }
          }
          
          
          
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
  print(paste("DPI SELECTED MERGED VERSION............time taken for cluster = ",ncluster, " is ",time.taken," sec"))
  
  df <- data.frame(THRESHOLD=thresholdVctr,TP=tp,FP=fp,TN=tn,FN=fn,PRECISION=precision,
                   RECALL=recall)
  
  # resultSelectedMergedFromUnMergedCluster=2TPFPsForAllThreshold.tsv
  
  fname <- paste(rStoreDataPath ,"resultSelectedMergedFromUnMergedCluster=",ncluster,"TPFPsForAllThresholdDPI.tsv",sep="")
  
  write.table(df, file = fname, sep = "\t", row.names = FALSE)
  
  
  
  
} # for loops ends for diff cluster

# store the time needed for each cluster
df <- data.frame(CLUSTER=clusterVCtr ,TIME=timeTaken)
fname <- paste(rStoreDataPath ,"resultSelectedMergedFromUnMergedCluster=TimeInSecDPI.tsv",sep="")
write.table(df, file = fname, sep = "\t", row.names = FALSE)


#=======================================

# UNMERGED DPI

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
  
  
  fname <- paste(rStoreDataPath,"resultUnMergedCluster=", ncluster,"ConnectoinMatrix.tsv",sep="")
  tmpgconnMat <- read.table(fname, header = T, sep = "\t")
  
  geneConnectionMat <- as.matrix(tmpgconnMat)
  
  
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
    
    dpiTempGeneConnMat <- tempGeneConnMat
    
    # now run DPI
    for(a in 3:ngenes){
      for(b in 2:(a-1)){
        for(c in 1:(b-1)){
          # do DPI only when a,b / b,c / a,c is 1 in geneConnectionMat
          if(dpiTempGeneConnMat[a,b] == 1){
            if(dpiTempGeneConnMat[b,c] == 1){
              if(dpiTempGeneConnMat[a,c] == 1){
                # all three genes have pair connection
                ab <- MIMat[a,b]
                bc <- MIMat[b,c]
                ac <- MIMat[a,c]
                
                if(ab <= min(bc,ac)){
                  # ab is indirect, remove it
                  tempGeneConnMat[a,b] = 0
                  tempGeneConnMat[b,a] = 0
                  
                }else if(bc <= min(ab,ac)){
                  # bc is indirect, remove it
                  tempGeneConnMat[b,c] = 0
                  tempGeneConnMat[c,b] = 0
                  
                }else if(ac <= min(ab,bc)){
                  # ac is indirect, remove it
                  tempGeneConnMat[a,c] = 0
                  tempGeneConnMat[c,a] = 0
                  
                }
                
              }
            }
          }
          
          
          
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
    
    print(paste("done for cluster = ",ncluster," after DPI with threshold = ", threshold))
    
    
  } # for each threshold ends
  
  end.time <- as.numeric(Sys.time())
  time.taken <- end.time - start.time
  
  
  timeTaken[m] <- time.taken
  print(paste("UNMERGED VERSION............time taken for cluster = ",ncluster, " is ",time.taken," sec"))
  
  df <- data.frame(THRESHOLD=thresholdVctr,TP=tp,FP=fp,TN=tn,FN=fn,PRECISION=precision,
                   RECALL=recall)
  
  fname <- paste(rStoreDataPath ,"resultUnMergedCluster=",ncluster,"TPFPsForAllThresholdDPI.tsv",sep="")
  
  write.table(df, file = fname, sep = "\t", row.names = FALSE)
  
  
  
  
} # for loops ends for diff cluster

# store the time needed for each cluster
df <- data.frame(CLUSTER=clusterVCtr ,TIME=timeTaken)
fname <- paste(rStoreDataPath ,"resultUnMergedCluster=TimeInSecDPI.tsv",sep="")
write.table(df, file = fname, sep = "\t", row.names = FALSE)

#==================================DPI END





