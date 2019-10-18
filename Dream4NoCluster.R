# NO CLUSTER

geneConnectionMatFileName <- paste(rStoreDataPath,"result",dataNumber,"NoCluster",sep="")

start.time <- as.numeric(Sys.time())

tp <- c()
fp <- c()
tn <- c()
fn <- c()
precision <- c()
recall <- c()

# NO CLUSTER

tempd <- matA
dtemp <- discretize(tempd)

#=============================== mutual info. AND distance matrix 

MIMatrix <- matrix(0,nrow = ngenes, ncol = ngenes)
NMIMatrix <- matrix(0,nrow = ngenes, ncol = ngenes)

DMatrix <- matrix(0,nrow = ngenes, ncol = ngenes)


for(i in 1:ngenes){
  for(j in 1:ngenes){
    MIMatrix[i,j] <- mutinformation(dtemp[,i], dtemp[,j])
    NMIMatrix[i,j] <- sqrt(1-exp(-2*MIMatrix[i,j])) #linfoot
    DMatrix[i,j] <- exp(-NMIMatrix[i,j])
  }
  
}

#============================== single entropy

E1Matrix <- matrix(0,nrow = ngenes, ncol = 1)
for(e in 1:ngenes){
  E1Matrix[e,1] <- entropy(dtemp[,e])
  
}


#============================  conditional entropy H(X|Y)  

CEMatrix <- matrix(0,nrow = ngenes, ncol = ngenes)

for(i in 1:ngenes){
  for(j in 1:ngenes){
    CEMatrix[i,j] <- condentropy(dtemp[,i],dtemp[,j])
  }
  
}


# TO DO->  rel_entr_reduc_2(m,n) = MI(m,n,1)/H1(m);
ERCE2Matrix <- matrix(0,nrow = ngenes, ncol = ngenes)

for(i in 1:ngenes){
  for(j in 1:ngenes){
    ERCE2Matrix[i,j] <- (MIMatrix[i,j]) / E1Matrix[i]
  }
  
}


#============================ joint entropy H(X,Y)

# JE2Matrix <- matrix(0,nrow = ngenes, ncol = ngenes)
# 
# for(i in 1:ngenes){
#   for(j in 1:ngenes){
#     JE2Matrix[i,j] <- E1Matrix[i] + CEMatrix[j,i]
#   }
#   
# }


#===================================   entropy reduction ERT 2


# #this only considers the minimum for a given gene
# diagMatrix <- diag(ngenes)
# CEMatrix <- CEMatrix + diagMatrix*(10^3)
# 
# ERTMatrix <- matrix(0,nrow = ngenes, ncol = 1)
# for(i in 1:ngenes){
#   minval <- min(CEMatrix[i,])
#   pos <- which(CEMatrix[i,] == minval)
#   ERTMatrix[i,1] <- pos
#   
# }

#considers all the genes where cond entropy is less than individual entropy
diagMatrix <- diag(ngenes)
CEMatrix <- CEMatrix + diagMatrix*(10^3)

ERTMatrix <- matrix(0,nrow = ngenes, ncol = ngenes)
for(i in 1:ngenes){
  for(j in 1:ngenes){
    if(CEMatrix[i,j] < E1Matrix[i,1]){
      # dependent
      ERTMatrix[i,j] <- 1
    }
    
  }
  
}



#===================================   connection matrix

connection <- matrix(0,nrow = ngenes, ncol = ngenes)

# # entropy reduction 2 var with only ONE min
# for(i in 1:ngenes){
#   connection[i,ERTMatrix[i]] <- ERCE2Matrix[i,ERTMatrix[i]]
#   
# }


# entropy reduction 2 var with all mins
for(i in 1:ngenes){
  
  for(j in 1:ngenes){
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

rownames(connection) <- geneNames
colnames(connection) <- geneNames

geneConnectionMat <- connection


# save the connection matrix
df <- data.frame(geneConnectionMat)
fname <- paste(geneConnectionMatFileName,"ConnectoinMatrix.tsv",sep="")

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
  
  print(paste("done for threshold = ", threshold))
  
  
} # for each threshold ends

end.time <- as.numeric(Sys.time())
time.taken <- end.time - start.time

print(paste("NO CLUSTER VERSION...........time taken = ",time.taken," sec"))

df <- data.frame(THRESHOLD=thresholdVctr,TP=tp,FP=fp,TN=tn,FN=fn,PRECISION=precision,
                 RECALL=recall)

fname <- paste(geneConnectionMatFileName,"TPFPsForAllThreshold.tsv",sep="")

write.table(df, file = fname, sep = "\t", row.names = FALSE)

#==============================================



df <- data.frame(TIME=time.taken)

fname <- paste(geneConnectionMatFileName,"TimeInSec.tsv",sep="")

write.table(df, file = fname, sep = "\t", row.names = FALSE)







