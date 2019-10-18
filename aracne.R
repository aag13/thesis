# running aracne

#========================================== minet aracne

# connection matrix

start.time <- as.numeric(Sys.time())

arc <- aracne(round(build.mim(matA), 5))

# for each threshold find TPs,FPs......
tp <- c()
fp <- c()
tn <- c()
fn <- c()

precision <- c()
recall <- c()

tempGeneConnMat <- matrix(0,nrow = ngenes, ncol = ngenes)

for(i in 1:ngenes){
  for(j in 1:ngenes){
    if(abs(arc[i,j]) <= 0){
      if(abs(arc[j,i]) <= 0){
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

tp <- ctp
fp <- cfp
tn <- ctn
fn <- cfn

precision <- cpr
recall <- crc


df <- data.frame(TP=tp,FP=fp,TN=tn,FN=fn,PRECISION=precision,RECALL=recall)
fname <- paste(rStoreDataPath ,"aracne.tsv",sep="")
write.table(df, file = fname, sep = "\t", row.names = FALSE)

end.time <- as.numeric(Sys.time())
time.taken <- end.time - start.time

print(paste("ARACNE...........time taken = ",time.taken," sec"))


df <- data.frame(TIME=time.taken)

fname <- paste(rStoreDataPath,"aracneTimeInSec.tsv",sep="")

write.table(df, file = fname, sep = "\t", row.names = FALSE)






