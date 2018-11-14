# May 22, 2017
# Dawit A. Yohannes
# dawit.yohannes@helsinki.fi
#
#
### ***************************************************************
### ******** Repseq Da Analysis: differential analysis of TCR repertoires using immune modules/structures 
### ***************************************************************


# General utility functions ....................................................................

addItemToObject <- function(object,item,itemName){
  
  # replace item or add it as new item
  if(itemName %in% names(object)){
    object[[itemName]] <- item
    
    
  }else{
    object[[length(object)+1]] <- item
    names(object)[length(object)] <- itemName 
    
  }
  
  return(object)
  
}


item.exists <- function(object,item){ 
  
  return(as.character(item) %in% names(object)) 
  
}


# Setup immune repertoire samples for differential abundance analysis....................................................................

#' Setup repseq object for differential abudance analysis, CDR3 counts are changed to counts per million for all samples
#' 
#' @param sampleNames is a vector of sample names (that have been already read into R )
#' @param samGroup a vector of group labels, e.g c(0,0,1,1), 0 for samples in first condition and 1 for samples in second condition.
#' @param normToCPM normalize all repseq samples to a total repertoire size of 1e6 (1 million) by default. 
#' @return a repSeq object which is a list containing the raw data (normalized to 1e6 in all samples if normToCPM is set to true),the sample names, and their group status
#' @export

setUp <- function(sampleNames=NULL,samGroup=NULL,normToCPM=T){
  if(is.null(sampleNames))
    stop("Sample names not given in samNames.")
  if(is.null(samGroup))
    stop("Sample names not given in samNames.")
  if(length(sampleNames) != length(samGroup))
    stop("The number of samples must be equal to the number of group labels.")
  
  sd = lapply(sampleNames,function(x) get(x))
  
  repSeqObj <- list(sampleData=sd,samNames=sampleNames,group=samGroup)
  names(repSeqObj$sampleData) <- repSeqObj$samNames
  
  for(i in 1:length(names(repSeqObj$sampleData))){
    
    # remove clonotypes with negative counts..
    repSeqObj$sampleData[[i]] <- repSeqObj$sampleData[[i]][!repSeqObj$sampleData[[i]]$COUNT < 0,]
    
    # remove clonotypes with a length of less than 3
    tooShortAA <- sapply(repSeqObj$sampleData[[i]]$AMINOACID,nchar) # some samples have too short AA clones which are only 2 AA long, shorter than the kmer length 3 we are working with, we just remove them, they are most likely errors.
    repSeqObj$sampleData[[i]] <- repSeqObj$sampleData[[i]][!(tooShortAA < 3),]
    
    
    repSeqObj$sampleData[[i]]$FREQUENCYCOUNT <- repSeqObj$sampleData[[i]]$COUNT / sum(repSeqObj$sampleData[[i]]$COUNT)
    
    if(normToCPM==T){
    # normalizing total repertoire size to counts per million
    tototalReads = 1e6
    repSeqObj$sampleData[[i]]$COUNT <- repSeqObj$sampleData[[i]]$FREQUENCYCOUNT * tototalReads; 
    }
    
  }
  
  return(repSeqObj)
  
}



# Repseq sample scaling and normalization ....................................................................

#' scale all samples to equal number of sequence reads
#' 
#' @param repSeqObj is a repseq data object containing all repertoire sample data
#' @param totalReads the total size of reads per sample to scale our data to. All repertoires are scaled to 100k reads by default.
#' @param resample a boolean argument to decide whether to resample repertoires, needed when data is too big, default is True
#' @param resampleSize Maximum possible number of clonotype reads repertoires will be resampled to, default is 1000.
#' @param useProb a boolean to indicate whether the resampling should be done probablistically or not.
#' @return repSeq data Object containing normalized repertoire data after resampling
#' @export

scaleSamples <- function(repSeqObj,totalReads=1e5,resample=T,resampleSize=3000,useProb=T){
  
  scaledSampleData = list()
  
  for(i in 1:length(names(repSeqObj$sampleData))){
    
    if(resample == T){
      
      # resampling done on actual relative frequency of clonotypes
      
      temp = repSeqObj$sampleData[[i]];
      if(useProb==T){
        resampledTemp = sample(temp$NUCLEOTIDE,resampleSize,replace=T,prob=temp$COUNT/sum(temp$COUNT))
      }else{
        resampledTemp = sample(temp$NUCLEOTIDE,resampleSize,replace=T)
      }
      resampledTempCount = sort(table(resampledTemp),decreasing=T) 
      
      resampledTemp = temp[match(names(resampledTempCount),temp$NUCLEOTIDE),]
      resampledTemp$COUNT = resampledTempCount
      resampledTemp$FREQUENCYCOUNT = resampledTempCount/sum(resampledTempCount)
      
      temp = resampledTemp
      
    }else{
      
      temp = repSeqObj$sampleData[[i]];  
      temp$FREQUENCYCOUNT <- temp$COUNT/sum(temp$COUNT)
    }
    
    # normalizing total repertoire size
    
    temp$COUNT <- temp$FREQUENCYCOUNT * totalReads; 
    scaledSampleData[[length(scaledSampleData)+1]] <- temp;
    

  }
  
  names(scaledSampleData) <- repSeqObj$samNames 
  
  
  
  repSeqObj <- addItemToObject(repSeqObj,scaledSampleData,"scaledSampleData")
  

  
  return(repSeqObj)
  
}



# Repseq sample scaling and normalization for background repertoires....................................................................

#' Pool and resample samples to create reference/background repertoires, normalize samples to equal number of sequence reads
#' 
#' @param repSeqObj is an object containing all repertoire sample data
#' @param totalReads the total size of reads per sample to scale our data to. All repertoires are scaled to 100k reads by default.
#' @param resample a boolean argument to decide whether to resample repertoires, needed when data is too big, default is True. If False, resamples as big as the original datasets are drawn from the pooled repertoire dataset
#' @param resampleSize Maximum possible number of clonotype reads repertoires will be resampled to, default is 1000.
#' @param useProb a boolean to indicate whether the resampling should be done probablistically or not.
#' @return repSeq data object created by downsampling after pooling all datasets  
#' @export

resampleReferenceSamples <- function(repSeqObj,totalReads=1e5,resample=T,resampleSize=1000,useProb=T){
  
  scaledSampleData = list()
  
  # first normalize and pool the data. 
  # Normalization is critical since relative frequency in the pool will be done from the clone sizes in respective samples (which depends on library size)
  pooledRep <- NULL
  
  for(i in 1:length(names(repSeqObj$sampleData))){
    tempForPooling <- repSeqObj$sampleData[[i]]
    tempForPoolingRelFreq <- tempForPooling$COUNT / sum(tempForPooling$COUNT)
    tempForPooling$COUNT <- tempForPoolingRelFreq * totalReads; 
    pooledRep <- rbind(pooledRep,tempForPooling)
  }
  
  pooledProp <- pooledRep$COUNT/sum(pooledRep$COUNT)
 
  for(i in 1:length(names(repSeqObj$sampleData))){
    
    if(resample == T){
      
      # resampling done on actual relative frequency of clonotypes
      
      temp = pooledRep;
      if(useProb==T){
        resampledTemp = sample(temp$NUCLEOTIDE,resampleSize,replace=T,prob=pooledProp)
      }else{
        resampledTemp = sample(temp$NUCLEOTIDE,resampleSize,replace=T)
      }
      resampledTempCount = sort(table(resampledTemp),decreasing=T) 
      
      resampledTemp = temp[match(names(resampledTempCount),temp$NUCLEOTIDE),]
      resampledTemp$COUNT = resampledTempCount
      resampledTemp$FREQUENCYCOUNT = resampledTempCount/sum(resampledTempCount)
      
    }else{
      
      temp = pooledRep;
      repSisze = nrow(repSeqObj$sampleData[[i]])
     
      if(useProb==T){
        resampledTemp = sample(temp$NUCLEOTIDE,repSisze,replace=T,prob=pooledProp)
      }else{
        resampledTemp = sample(temp$NUCLEOTIDE,repSisze,replace=T)
      }
      resampledTempCount = sort(table(resampledTemp),decreasing=T) 
      
      resampledTemp = temp[match(names(resampledTempCount),temp$NUCLEOTIDE),]
      resampledTemp$COUNT = resampledTempCount
      resampledTemp$FREQUENCYCOUNT = resampledTempCount/sum(resampledTempCount)
    }
    
    # normalizing total repertoire size
    
    resampledTemp$COUNT <- resampledTemp$FREQUENCYCOUNT * totalReads; 
    scaledSampleData[[length(scaledSampleData)+1]] <- resampledTemp;
    
  }
  
  names(scaledSampleData) <- repSeqObj$samNames 
  
  
  
  repSeqObj <- addItemToObject(repSeqObj,scaledSampleData,"scaledSampleData")
  

  
  return(repSeqObj)
  
}


#' Pool and resample samples to create reference/background repertoires and return only CDR3s, normalizes samples to equal number of sequence reads
#' 
#' @param repSeqObj is an object containing all repertoire sample data
#' @param nClones number of CDR3 sequences to randomly draw
#' @param cType the type of CDR3 sequence to return for the drawn clones, either AA or NT
#' @param useProb a boolean to indicate whether the resampling should be done probablistically or not.
#' @return returns a vector of CDR3 sequences with size nClones that is randomly sampled from the pooled datasets of all repseq samples in repSeqObj 
#' @export

sampleFromPooled <- function(repSeqObj,nClones,cType="AA",useProb=T){
  
  scaledSampleData = list()
  
  pooledRep <- NULL
  
  for(i in 1:length(names(repSeqObj$sampleData))){
    tempForPooling <- repSeqObj$sampleData[[i]]
    tempForPoolingRelFreq <- tempForPooling$COUNT / sum(tempForPooling$COUNT)
    pooledRep <- rbind(pooledRep,tempForPooling)
  }
  
  pooledProp <- pooledRep$COUNT/sum(pooledRep$COUNT)
  

  temp = pooledRep;
  
  if(cType == "AA"){
    cpool = temp$AMINOACID
  }else{
    cpool = temp$NUCLEOTIDE
  }
  
  if(useProb==T){
    randomDrawnClones = sample(cpool,nClones,replace=T,prob=pooledProp)
  }else{
    randomDrawnClones = sample(cpool,nClones,replace=T)
  }
  
  return(randomDrawnClones)
  
}

#' Pool and resample samples to create reference/background repertoires, normalizes samples to equal number of sequence reads
#' 
#' @param repSeqObj is an object containing all repertoire sample data
#' @param nClones number of CDR3 sequences to randomly draw
#' @param useProb a boolean to indicate whether the resampling should be done probablistically or not.
#' @return returns a single repseq data of randomly selected CDR3 sequences and their associated features, that is randomly sampled from the pooled datasets of all repseq samples in repSeqObj 
#' @export

sampleWithFeatuersFromPooled <- function(repSeqObj,nClones,useProb=T){
  
  scaledSampleData = list()
  
  pooledRep <- NULL
  
  for(i in 1:length(names(repSeqObj$sampleData))){
    tempForPooling <- repSeqObj$sampleData[[i]]
    tempForPoolingRelFreq <- tempForPooling$COUNT / sum(tempForPooling$COUNT)
    pooledRep <- rbind(pooledRep,tempForPooling)
  }
  
  pooledProp <- pooledRep$COUNT/sum(pooledRep$COUNT)
  
  
  temp = pooledRep;
  
  if(useProb==T){
    randomDrawnClonesidx = sample(1:nrow(temp),nClones,replace=T,prob=pooledProp)
    randomDrawnClones <- temp[randomDrawnClonesidx,]
    
  }else{
    randomDrawnClonesidx = sample(1:nrow(temp),nClones,replace=T)
    randomDrawnClones <- temp[randomDrawnClonesidx,]
    
  }
  
  return(randomDrawnClones)
  
}


#' Pool all repertoires
#' 
#' @param repSeqObj is an object containing all repertoire sample data
#' @param g a vector indicating the group to sample from, default is NULL and pools all samples from all groups, to pool samples only in group 1, g should be c(1) 
#' @return a single resepse dataset that has all CDR3 sequences from all samples, the frequency of all CDR3s is updated depending o the original counts in the original data.
#' @export

getPooledSamples <- function(repSeqObj,g=NULL){
  
  scaledSampleData = list()
  
  if(is.null(g)){
    g=repObj$group
  }
  pooledRep <- NULL
  
  for(i in names(repSeqObj$sampleData)[repSeqObj$group %in% g]){
    tempForPooling <- repSeqObj$sampleData[[i]]
    tempForPoolingRelFreq <- tempForPooling$COUNT / sum(tempForPooling$COUNT)
    pooledRep <- rbind(pooledRep,tempForPooling)
  }
  
  pooledRep$FREQUENCYCOUNT <- pooledRep$COUNT/sum(pooledRep$COUNT)
  
  return(pooledRep)
  
}


# Clonotype clustering and related ....................................................................


#' This function performs pairwise 4-mer nt usage clustering of clonotypes in repSeq data
#' 
#' @param sam a repertoire data (of only productive clonotypes,immunoseq format) for whose clonotypes are going to be clustered based on 4-mer usage frequencies.
#' @return a list containing the clonotypes,clonotypes with 4-mer usage frequencies,distance matrix between clonotypes,cluster labels for each clonotype,and cluster centroids.
#' @export

getClusterLables <- function(sam,k=10,clusterby="NT",kmerWidth=4,posWt=F,distMethod="euclidean",useDynamicTreeCut=T){
  
 
  if(clusterby=="NT"){
    seqs <- unique(sam$CDR3NT) # use extracted CDR3 NT sequence 
    seq_mers <- getKmerFrequency(seqs,type="NT",k=kmerWidth)
     
  }else{
    seqs <- unique(sam$AMINOACID) # use AA
    seq_mers <- getKmerFrequency(seqs,type="AA",k=kmerWidth)
  }
  
  kmers <- colnames(seq_mers)
  
  if(posWt==T){
    
    cat("\t...calculating positional weights for kmer frequencies \n")
    
    # give position based weights to kmers, and normalize kmer frequence matrix
    weighted.seq.mers = NULL
    
    for(seq in seqs){
      
      kmerWeights = sapply(kmers,function(x) determineWeight(x,as.character(seq))) 
      weighted.seq.mers = rbind(weighted.seq.mers,kmerWeights)  
    }
    
    rownames(weighted.seq.mers) <- seqs
    
    #dim(weighted.seq.mers)
    
    seq_mers <- weighted.seq.mers
    
    # normalize the kmer counts for sequence length
    #seq_mers <- t(apply(seq_mers,1,function(x) x/sum(x)))
    seq_mers <- as(seq_mers, "sparseMatrix")
    
  }
  
 
  #print(head(seq_mers))
  
  
  #dcalculated <- dist(seq_mers)
  dcalculated <- dist.matrix(seq_mers, method=distMethod, convert=T, as.dist=TRUE)
  #dcalculated <- dist.matrix(seq_mers, method="cosine", convert=T, as.dist=TRUE)
  

  #print(head(dcalculated))
  hc<-hclust(dcalculated,method="complete")
  

  #clusters=cutreeDynamicTree(hc, maxTreeHeight = max(dcalculated), minModuleSize = 50)
  #cls= cutreeHybrid(hc,distM=as.matrix(dcalculated))
  #cls <- cutree(hc,k)
  if(useDynamicTreeCut==T){
    cls <-cutreeDynamic(hc,distM=as.matrix(dcalculated),minClusterSize=20,verbose=0) 
  }else{
    cls <- cutree(hc,k)
  }
  

 
  clsCenters <- getCenters(seq_mers,cls)
  #distanceAndClusters <- list(seqs=as.character(seqs),seqmers=seq_mers,distmatrix=as.matrix(dcalculated),clusters=cls,clusterCenters=clsCenters)
  
  # without returning the distance matrix
  distanceAndClusters <- list(seqs=as.character(seqs),seqmers=seq_mers,clusters=cls,clusterCenters=clsCenters)
  
  
 
  return(distanceAndClusters)
  
}


#' Compute cluster centroids
#' 
#' This function calculates cluster centroids given sequences (rows) and k-mer usage (columns), and cluster assignment labels for each sequence. Centroid vector for each cluster is the average of the k-mer usage frequencies of clonotypes in the cluster.
#' 
#' @param seqmers is a dataframe/matrix with observations (sequences) on the rows and features (k-mer usage frequences) on the columns
#' @param clslabels is a vector containing cluster assignment labels for every sequence in seqmers
#' @return a matrix with rows containing the cluster labels and columns containing an average usage frequency for each possible k-mer representing the cluster centroid.
#' 

getCenters <- function(seqmers,clslabels){
  clusterCenters <- NULL # holds the cluster centroids
  clusterCenterNames <- c()
  
  clusters <- table(clslabels)
  
  for(i in names(clusters)){
    # drop clusters labeled zero which contain clearly unassigned members when using dynamic tree cut
    if(i > 0){
      selected <- seqmers[which(clslabels==i),,drop=F]
      centroid <- Matrix::colMeans(selected)
      clusterCenterNames <- c(clusterCenterNames,paste(i,sep=""))
      clusterCenters <- rbind(clusterCenters,centroid)
    }
  }
  
  rownames(clusterCenters) <- clusterCenterNames
  return(clusterCenters)
  
}

#' Finding optimal number of clusters
#' 
#' Finds optimal k as the average optimal k detected from the within sample clustering of selected repertoire samples.
#' 
#' @param repSeqObj is an object containing all repertoire sample data
#' @param distMethod the distance method used for determining distance between CDR3 feature vectors, default "euclidean"
#' @param posWt boolean to give weights to kmer frequencies depending on their position in the CDR3
#' @return returns an optimal k for dividing unsupervised clustering results into k compact clusters.
#' 

findOptimalK <- function(repSeqObj,nSamEval=2,clusterby,minCSizePerc = 0.1,minNClonesPerCluster=20,kmerWidth=4,posWt=T,distMethod="euclidean"){
  
  # select samples on which to evaluate possible ks
  
  #cat("Searching for optimal k...\n")
  
  
  numberOfSamples = length(repSeqObj$samNames)
  numberOfPairs = numberOfSamples/2
  # 
  # if(nSamEval > numberOfPairs){
  #   sizeOfRandomSamples = numberOfPairs/2
  # }else{
  #   sizeOfRandomSamples = nSamEval
  # }
  # 
  # grpLevels = levels(factor(repSeqObj$group))
  # 
  # selectedSamplesIndex = sample(1:numberOfPairs,sizeOfRandomSamples)
  # selectedSampleNames = repSeqObj$samNames[which(repSeqObj$group == grpLevels[1])][selectedSamplesIndex]
  # 
  selectedSamplesIndex = sample(1:numberOfSamples,nSamEval)
  selectedSampleNames = repSeqObj$samNames[selectedSamplesIndex]
  
  repSeqObj  = addItemToObject(repSeqObj,selectedSampleNames,"selectedSampleNames")
  repSeqObj  = addItemToObject(repSeqObj,numberOfSamples,"numberOfSamples")
  
  optimalK = c()
  
  # find optimal k for each randomly selected sample
  
  
  for(sam in repSeqObj$selectedSampleNames){
    
    sam1 =  repSeqObj$scaledSampleData[[sam]]
    
    
    # decide to use either AA or CDR3nt
    if(clusterby=="NT"){
      seqs <- unique(sam1$CDR3NT) # use extracted CDR3 NT sequence 
      seq_mers <- getKmerFrequency(seqs,type="NT",k=kmerWidth)
      
    }else{
      seqs <- unique(sam1$AMINOACID) # use AA
      seq_mers <- getKmerFrequency(seqs,type="AA",k=kmerWidth)
      
    }
    
   
    kmers <- colnames(seq_mers)
    
    
    if(posWt==T){
      
      cat("\t...calculating positional weights for kmers frequencies... \n")
      
      # give position based weights to kmers, and normalize kmer frequence matrix
      weighted.seq.mers = NULL
      
      for(seq in seqs){
        
        kmerWeights = sapply(kmers,function(x) determineWeight(x,as.character(seq))) 
        weighted.seq.mers = rbind(weighted.seq.mers,kmerWeights)  
      }
      
      rownames(weighted.seq.mers) <- seqs
      
      dim(weighted.seq.mers)
      
      seq_mers <- weighted.seq.mers
      
      # normalize the kmer counts for sequence length
      #seq_mers <- t(apply(seq_mers,1,function(x) x/sum(x)))
      seq_mers <- as(seq_mers, "sparseMatrix")
    }
    
    #look for optimal k, from k that allows clusters that contain 2% of the repertoire to k that allows clusters to contain 20 clones per cluster
    minNClust = round(nrow(seq_mers)/(minCSizePerc * nrow(seq_mers)))
    #maxK = round(nrow(seq_mers)/4)
    maxK = round(nrow(seq_mers)/minNClonesPerCluster)
    
    # normalize the kmer counts for sequence length
    #seq_mers <- t(apply(seq_mers,1,function(x) x/sum(x)))
    
    
    doParallel::registerDoParallel(cores=detectCores())  
    
    
    # evaluate ks from 5-maxK based on silhouette and data compression 
    
    ks_sils = foreach(i=1:10,.export=c('silhouette','summary','cosineDist'),.packages='cluster',.combine=rbind) %dopar% {
      
      #seqmersResampled <- seq_mers[sample(1:nrow(seq_mers),size=round(nrow(seq_mers)/2),replace=F),]
      #dcalculated <- cosineDist(seqmersResampled)
      #dcalculated <- dist(seqmersResampled)
      
     
       
      seqmersResampled <- seq_mers
      dcalculated <- dist.matrix(seqmersResampled, method=distMethod, convert=T, as.dist=TRUE)
        
      hc<-hclust(dcalculated,method="complete")  # it was average linkage, but it did not work well
      
     
      sils <- c()
      for(j in minNClust:maxK){
        
       
        if(useCmeans==T){
          fuzzyc <- e1071::cmeans(seqmersResampled, centers=50, iter.max = 100, verbose = FALSE,dist = "euclidean", method = "cmeans")
          silmuke <- fclust::SIL(as.matrix(seqmersResampled), fuzzyc$membership)
        }else{
          cls <- cutree(hc,j)
          #cls <-cutreeDynamic(hc,distM=as.matrix(dcalculated),minClusterSize=j,verbose=0)
          siSeqMers <- cluster::silhouette(cls,dcalculated)
        }
        avgSilWidth <- summary(siSeqMers)$avg.width
        
        nClusSilWidthAboveAverage <- sum(summary(siSeqMers)$clus.avg.widths > avgSilWidth)/j # between 0-1
        #compressionGain <- 1-j/nrow(seqmersResampled) # between 0-1
        compressionGain <-j/nrow(seqmersResampled)
        avgSil <- (avgSilWidth + 1)/2 # add 1 to avg silhouette and divide by 2 to make it between 0 and 1
        
         
        kscore <- avgSil
        sils <- c(sils,kscore)
        
      }  
      
      sils 
      
    }
    
    maxScoreK = minNClust + which.max(apply(ks_sils,2,mean)) - 1
    
  
    plot(minNClust:maxK,apply(ks_sils,2,mean),xlab="K",ylab="mean K evaulation score",type="l",main=paste("K evaluation scores in sample:",sam))
    abline(v=maxScoreK,col="red")
    axis(1, at=maxScoreK,labels=maxScoreK,col="red")
  
    

  
    optimalK = c(optimalK,maxScoreK)
    
  }
  
  
  optk = mean(optimalK)
  return(ceiling(optk))
  
}



CountkmerFrequency <- function(seqs,type="NT",k=4){
  
  DNAletters=c("A","C","G","T")
  AAletters=c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
  
  if(type=="NT"){
    allCombinations = do.call(expand.grid, rep(list(DNAletters), k))
    kmers = apply(allCombinations,1,paste,collapse="")
    
  }else{
    allCombinations = do.call(expand.grid, rep(list(AAletters), k))
    kmers = apply(allCombinations,1,paste,collapse="")
    
  }
  
  kmerCount <- sapply(seqs,function(x){
    
    sapply(kmers,function(y){sum(gregexpr(paste("(?=",y,")",sep=""),x,perl=T)[[1]] > 0)}) #with look-ahead assertion to allow overlapping hits
    
  })
  
  kmerCount <- t(kmerCount)
  
  kmerCount <- kmerCount[,order(colnames(kmerCount))]
  
  return(kmerCount)
  
}


#' Count kmer frequencies in a CDR3
#' 
#' @param seqs a vector of all CDR3 sequences
#' @param type the type of kmers, NT or AA
#' @param k the size of k, default is 4
#' @param normForLength normalize the kmer frequencies by the length of the CDR3, default is False
#' @return a sparse matrix of sequences versus kmer counts
getKmerFrequency <- function(seqs,type="NT",k=4,normForLength=F){
  
  if(type=="NT"){
    seqSet <- Biostrings::DNAStringSet(seqs) 
    seq_mers <- Biostrings::oligonucleotideFrequency(seqSet,width=k,as.prob=F)
    rownames(seq_mers) <- seqs   
    
  }else{
    #cat("\t...Determining AA k-mers and their frequency in clonotypes may take a long time for k >= 4... \n")
    
    seqList <- as.list(seqs)
    seqAAbin <- as.AAbin(seqList)
    seq_mers <- kmer::kcount(seqAAbin,k=k)
    rownames(seq_mers) <- seqs
    
    ## selecting the most variable kmers only since AA kmers are numorous
    #vars = apply(seq_mers,2,var)
    #seq_mers <- seq_mers[,vars > 0]
    
  }
  
  # normalize the kmer counts for sequence length
  if(normForLength == T) 
    seq_mers <- t(apply(seq_mers,1,function(x) x/sum(x)))
  
  seq_mers <- as(seq_mers, "sparseMatrix")
  return(seq_mers)
  
}



determineWeight <- function(kmer,seq){
  
  kmerPositions = as.numeric(gregexpr(kmer,seq)[[1]])
  seqLength <- nchar(seq)
  weightGroups <- seq(seqLength/3,seqLength,seqLength/3)
  
  kmerWts =c()
  for(i in kmerPositions){
    if(i == -1){
      wt=0
    }else{
      # The kmer is observed in the CDR3. We increase the weight for its frequency by increasing it depending on where it's found in the sequence
      wps=sum(weightGroups > i)
      if(wps == 3){wt=5} #v-area weight
      if(wps == 2){wt=10} #middle area weight
      if(wps == 1){wt=1} #j-area weight
    }
    
    kmerWts = c(kmerWts,wt)
  }
  
  return(sum(kmerWts))
  
}


 
cosineDist <- function(x){
  as.dist(1 - x%*%t(x)/(sqrt(rowSums(x^2) %*% t(rowSums(x^2))))) 
}


shannonEntropy <- function(freqs){
  return(-sum(freqs * (log2(freqs)/log2(length(freqs)))))
}


 
getClusterFoldChanges <- function(sam1,sam2,consensusT,s1cls,s2cls){
  
  clusterFc=c()
  for (ck in 1:nrow(consensusT)){
    
    consensusK <- consensusT[ck,]
    
    sam1selected <- sam1[s1cls$clusters==consensusK[1],]
    sam2selected <- sam2[s2cls$clusters==consensusK[2],]
    

    sam1CountPerClone <- sum(sam1$COUNT) / nrow(sam1)
    sam2CountPerClone <- sum(sam2$COUNT) / nrow(sam2)
    
    sam1SelectedCountPerClone <- sum(sam1selected$COUNT) / nrow(sam1selected)  
    sam2SelectedCountPerClone <- sum(sam2selected$COUNT) / nrow(sam2selected)
    
    sam1SelectedVsTotal <- sam1SelectedCountPerClone/sam1CountPerClone
    sam2SelectedVsTotal <- sam2SelectedCountPerClone/sam2CountPerClone
    
    sam1vssam2Selected <- round(sam2SelectedVsTotal/sam1SelectedVsTotal,digits = 2)
    
    clusterFc <- c(clusterFc,sam1vssam2Selected)
    
  }
  
  return(clusterFc)
  
}


#' finds optimal clusters within samples and matches them across samples
#' @keywords internal
#' 
findOptimalClusters <- function(repSeqObj,k,clusterby="NT",kmerWidth=4,posWt=F,distMethod="euclidean",useDynamicTreeCut=T,matchingMethod=c("hc","km","og")){
  
  # First do clustering for all samples
  withinSampleClusters = list()
  
  #cat("Performing within sample clustering of clonotypes :\n")
  
  
  # start progress bar
  #pb <- txtProgressBar(min = 0, max = length(repSeqObj$samNames), style = 3)
  
  for(i in 1:length(repSeqObj$samNames)){
    
    cResult = getClusterLables(repSeqObj$scaledSampleData[[i]],k,clusterby,kmerWidth,posWt,distMethod,useDynamicTreeCut=useDynamicTreeCut)
    
    withinSampleClusters[[length(withinSampleClusters)+1]] <- cResult
    
    #setTxtProgressBar(pb, i)
    
  }
  
  #close(pb)
  
  names(withinSampleClusters) <- repSeqObj$samNames 
  
  
  repSeqObj <- addItemToObject(repSeqObj,withinSampleClusters,"withinSampleClusters")
  
  
  # Next across sample clustering of centroid - to find matching clusters across samples
  
  # Do matching using kmeans, and allow for some unmatched clusters 
  
  
      # #first collect all centroids into one matrix
      # 
      # combinedCentroids = NULL
      # for(sam in repSeqObj$samNames){
      # 
      #   samClusterCentroids <- repSeqObj$withinSampleClusters[[sam]]$clusterCenters
      #   ntimes <- nrow(samClusterCentroids)
      #   #print(nrow(samClusterCentroids))
      # 
      #   rownames(samClusterCentroids) <- paste(rep(sam,ntimes),"_",rownames(samClusterCentroids),sep="")
      #   combinedCentroids <- rbind(combinedCentroids,samClusterCentroids)
      # 
      # }
      # 
      # maxKacrossSamples = max(sapply(repSeqObj$samNames,function(x) length(table(repSeqObj$withinSampleClusters[[x]]$clusters))))
      # minKacrossSamples = min(sapply(repSeqObj$samNames,function(x) length(table(repSeqObj$withinSampleClusters[[x]]$clusters))))
      # 
      # # find good k across samples
      # sils <- c()
      # combinedCentroidsdist <- dist(combinedCentroids)
      # 
      # for(j in minKacrossSamples:maxKacrossSamples){
      #   
      #   centroidKmcls <- kmeans(combinedCentroids,maxKacrossSamples,iter.max = 50,nstart = 50)
      #   
      #   sicentroid <- silhouette(centroidKmcls$cluster,combinedCentroidsdist)
      #   avgSilWidth <- summary(sicentroid)$avg.width
      #   
      #   nClusSilWidthAboveAverage <- sum(summary(sicentroid)$clus.avg.widths > avgSilWidth)/j # between 0-1
      #   #compressionGain <- 1-j/nrow(seqmersResampled) # between 0-1
      #   compressionGain <-j/nrow(combinedCentroids)
      #   avgSil <- (avgSilWidth + 1)/2 # add 1 to avg silhouette and divide by 2 to make it between 0 and 1
      #   
      #   #kscore <- (nClusSilWidthAboveAverage + compressionGain + avgSil)/3 # between 0 and 1
      #   
      #   kscore <- (nClusSilWidthAboveAverage + avgSil)/2
      #   sils <- c(sils,avgSilWidth)
      #   
      # }  
      # 
      #  
      # 
      # goodK = minKacrossSamples + which.max(sils) - 1
      # 
      # 
      # centroidKmcls <- kmeans(combinedCentroids,goodK,iter.max = 50,nstart = 50)
      # matchTable = matrix(data = NA, nrow = goodK, ncol = length(repSeqObj$samNames),dimnames = list(1:goodK,repSeqObj$samNames))
      # 
      # for(i in 1:length(table(centroidKmcls$cluster))){
      #   tempCentroidClustMem <- centroidKmcls$cluster[centroidKmcls$cluster==i]
      #   availSamples <- sapply(names(tempCentroidClustMem),function(x) strsplit(x,"_")[[1]][1])
      #   matchingClustersInSamples <- sapply(names(tempCentroidClustMem),function(x) strsplit(x,"_")[[1]][2])
      #   # merge within sample clusters if they are matching the same other cluster in another sample
      #   # so merging happens if matching clusters across samples involve more than one cluster per sample at least in one case
      #   removedLables <- c()
      #   if(sum(table(availSamples) > 1) > 0){
      #     for(name in names(which(table(availSamples)>1))){
      #       clustersToBeMerged <- matchingClustersInSamples[availSamples == name]
      #       idxTobeRemoved <- which(availSamples == name)[-1]
      # 
      #       repSeqObj$withinSampleClusters[[name]]$clusters[repSeqObj$withinSampleClusters[[name]]$clusters %in% clustersToBeMerged] <- as.character(clustersToBeMerged[1])
      # 
      #       updatedCentroid <- getCenters(repSeqObj$withinSampleClusters[[name]]$seqmers[repSeqObj$withinSampleClusters[[name]]$clusters %in% clustersToBeMerged[1],],
      #                  repSeqObj$withinSampleClusters[[name]]$clusters[repSeqObj$withinSampleClusters[[name]]$clusters %in% clustersToBeMerged[1]])
      # 
      #       repSeqObj$withinSampleClusters[[name]]$clusterCenters[as.character(clustersToBeMerged[1]),] <- updatedCentroid
      #       removedLables <- c(removedLables,idxTobeRemoved)
      #     }
      #   }
      #   #print(i)
      #   #print(availSamples)
      #   #print(matchingClustersInSamples)
      # 
      #   if(length(removedLables) > 0)
      #     matchingClustersInSamples <- matchingClustersInSamples[- removedLables]
      #   availSamples <- unique(availSamples)
      # 
      #   matchTable[i, match(availSamples, repSeqObj$samNames)] <- as.numeric(matchingClustersInSamples)
      # }

  

  if(matchingMethod == "og"){    
    matchTable <- getClusterMatches(repSeqObj,matchingMethod="og") # og is the own greedy method of matching
  }else if(matchingMethod == "km"){
    matchTable <- getClusterMatches(repSeqObj,matchingMethod="km")
  }else if(matchingMethod == "hc"){
    matchTable <- getClusterMatches(repSeqObj,matchingMethod="hc",distMethod=distMethod)
  }else{
    cat("Matching method can only be one of hc, km or og.hc method will be used.\n") 
    matchTable <- getClusterMatches(repSeqObj,matchingMethod="hc",distMethod=distMethod)
    
  }
  
  
  # cluster matching using sample 1 only
  # This is not the optimal way to match the cluster centroids across samples. A balanced k-means (as commented above)
  # is probably the best way to go about this. In this bit, identification of matching cluster centroids in other samples is done for each 
  # centroid of sample 1. If a cluster centroid (of for example sample 2) has already been matched to
  # one of the centroids in sample 1, it's removed from further consideration and the next best match is selected (until)
  # none of the centroids have been matched previously to any other sample 1 centroid.  
  
  
    # matchTable = NULL
    # 
    # clMatches = c()
    # samClusterCentroids = repSeqObj$withinSampleClusters[[1]]$clusterCenters
    # subRepertoireNames = 1:nrow(samClusterCentroids)
    # 
    # cat("Finding matching clusters across samples: \n")
    # 
    # 
    # for (i in 1:nrow(samClusterCentroids)){
    #   
    #   for(s2 in 2:length(repSeqObj$samNames)){
    #     
    #     combinedsams1 <- rbind(samClusterCentroids[i,],repSeqObj$withinSampleClusters[[s2]]$clusterCenters)
    #     
    #     imin <- getMatchingCluster(combinedsams1)
    #     
    #     if(!is.null(matchTable) & (imin %in% matchTable[,s2])){
    #       
    #       cimins <- c()
    #       while(imin %in% matchTable[,s2]){
    #         
    #         cimins <- c(cimins,imin)
    #         #print(cimins)
    #         combinedsams1 <- rbind(samClusterCentroids[i,],repSeqObj$withinSampleClusters[[s2]]$clusterCenters[!rownames(repSeqObj$withinSampleClusters[[s2]]$clusterCenters) %in% as.numeric(cimins),])
    #         
    #         if((length(subRepertoireNames)- length(cimins)) > 1){
    #           imin <- getMatchingCluster(combinedsams1)
    #         }else{
    #           # cluster is more similar to all other clusters but this one. But all the other clusters have found matches before
    #           imin <- subRepertoireNames[!(subRepertoireNames %in% cimins)]
    #         }
    #         
    #       }
    #       
    #     }
    #     
    #     
    #     clMatches <- c(clMatches,as.numeric(imin))
    #     
    #   }
    #   
    #   #print(clMatches)
    #   matchTable <- rbind(matchTable,c(i,clMatches))
    #   clMatches = c()
    # }
    # 
    # colnames(matchTable) <- repSeqObj$samNames
  
  repSeqObj <- addItemToObject(repSeqObj,matchTable,"clusterMatchTable")
  
  
  return(repSeqObj)
  
  
  
}



getClusterMatches<- function(repSeqObj,matchingMethod=c("hc","km","og"),distMethod="euclidean"){
  

  if(matchingMethod == "og"){
    # cluster matching using sample 1 only
    # This is not the optimal way to match the cluster centroids across samples. A balanced k-means (as commented above)
    # is probably the best way to go about this. In this bit, identification of matching cluster centroids in other samples is done for each 
    # centroid of sample 1. If a cluster centroid (of for example sample 2) has already been matched to
    # one of the centroids in sample 1, it's removed from further consideration and the next best match is selected (until)
    # none of the centroids have been matched previously to any other sample 1 centroid.  
  
    matchTable = NULL
    
    clMatches = c()
    samClusterCentroids = repSeqObj$withinSampleClusters[[1]]$clusterCenters
    subRepertoireNames = 1:nrow(samClusterCentroids)
    
    #cat("Finding matching clusters across samples: \n")
    
    
    for (i in 1:nrow(samClusterCentroids)){
      
      for(s2 in 2:length(repSeqObj$samNames)){
        
        combinedsams1 <- rbind(samClusterCentroids[i,],repSeqObj$withinSampleClusters[[s2]]$clusterCenters)
        
        imin <- getMatchingCluster(combinedsams1)
        
        if(!is.null(matchTable) & (imin %in% matchTable[,s2])){
          
          cimins <- c()
          while(imin %in% matchTable[,s2]){
            
            cimins <- c(cimins,imin)
            #print(cimins)
            combinedsams1 <- rbind(samClusterCentroids[i,],repSeqObj$withinSampleClusters[[s2]]$clusterCenters[!rownames(repSeqObj$withinSampleClusters[[s2]]$clusterCenters) %in% as.numeric(cimins),])
            
            if((length(subRepertoireNames)- length(cimins)) > 1){
              imin <- getMatchingCluster(combinedsams1)
            }else{
              # cluster is more similar to all other clusters but this one. But all the other clusters have found matches before
              imin <- subRepertoireNames[!(subRepertoireNames %in% cimins)]
            }
            
          }
          
        }
        
        
        clMatches <- c(clMatches,as.numeric(imin))
        
      }
      
      #print(clMatches)
      matchTable <- rbind(matchTable,c(i,clMatches))
      clMatches = c()
    }
    
    colnames(matchTable) <- repSeqObj$samNames
    
  }else{
    
    # If we expect different number of clusters in the  
    # #first collect all centroids into one matrix
    #cat("Finding matching clusters across samples: \n")
    
    combinedCentroids = NULL
    for(sam in repSeqObj$samNames){
      
      samClusterCentroids <- repSeqObj$withinSampleClusters[[sam]]$clusterCenters
      ntimes <- nrow(samClusterCentroids)
     
      
      rownames(samClusterCentroids) <- paste(rep(sam,ntimes),"_",rownames(samClusterCentroids),sep="")
      combinedCentroids <- rbind(combinedCentroids,samClusterCentroids)
      
    }
    
    if(matchingMethod == "km"){
    maxKacrossSamples = max(sapply(repSeqObj$samNames,function(x) length(table(repSeqObj$withinSampleClusters[[x]]$clusters))))
    minKacrossSamples = min(sapply(repSeqObj$samNames,function(x) length(table(repSeqObj$withinSampleClusters[[x]]$clusters))))
    
    # find good k across samples
    sils <- c()
    combinedCentroidsdist <- dist(combinedCentroids)
    
    for(j in minKacrossSamples:maxKacrossSamples){
      
      centroidKmcls <- kmeans(combinedCentroids,maxKacrossSamples,iter.max = 50,nstart = 50)
      
      sicentroid <- cluster::silhouette(centroidKmcls$cluster,combinedCentroidsdist)
      avgSilWidth <- summary(sicentroid)$avg.width
      
      nClusSilWidthAboveAverage <- sum(summary(sicentroid)$clus.avg.widths > avgSilWidth)/j # between 0-1
      #compressionGain <- 1-j/nrow(seqmersResampled) # between 0-1
      compressionGain <-j/nrow(combinedCentroids)
      avgSil <- (avgSilWidth + 1)/2 # add 1 to avg silhouette and divide by 2 to make it between 0 and 1
      
      #kscore <- (nClusSilWidthAboveAverage + compressionGain + avgSil)/3 # between 0 and 1
      
      kscore <- (nClusSilWidthAboveAverage + avgSil)/2
      sils <- c(sils,avgSilWidth)
      
    }  
    
    
    
    goodK = minKacrossSamples + which.max(sils) - 1
    centroidKmcls <- kmeans(combinedCentroids,goodK,iter.max = 100,nstart = 50)
    
    }else if(matchingMethod == "hc"){
      
      centroidKmcls <- list()
      dcalc1= dist.matrix(combinedCentroids, method=distMethod, convert=T, as.dist=TRUE)
      
      hc<-hclust(dcalc1,method="complete")
      cls <- dynamicTreeCut::cutreeDynamic(hc,distM=as.matrix(dcalc1),minClusterSize=1,verbose=0) 
      
      centroidKmcls$cluster <- cls
      names(centroidKmcls$cluster) <- rownames(combinedCentroids)
      
      goodK <- length(unique(centroidKmcls$cluster))
    
    }
    
    matchTable = matrix(data = NA, nrow = goodK, ncol = length(repSeqObj$samNames),dimnames = list(1:goodK,repSeqObj$samNames))
    
    for(i in 1:length(table(centroidKmcls$cluster))){
      tempCentroidClustMem <- centroidKmcls$cluster[centroidKmcls$cluster==i]
      availSamples <- sapply(names(tempCentroidClustMem),function(x) strsplit(x,"_")[[1]][1])
      matchingClustersInSamples <- sapply(names(tempCentroidClustMem),function(x) strsplit(x,"_")[[1]][2])
      # merge within sample clusters if they are matching the same other cluster in another sample
      # so merging happens if matching clusters across samples involve more than one cluster per sample at least in one case
      removedLables <- c()
      if(sum(table(availSamples) > 1) > 0){
        for(name in names(which(table(availSamples)>1))){
          clustersToBeMerged <- matchingClustersInSamples[availSamples == name]
          idxTobeRemoved <- which(availSamples == name)[-1]
          
          repSeqObj$withinSampleClusters[[name]]$clusters[repSeqObj$withinSampleClusters[[name]]$clusters %in% clustersToBeMerged] <- as.character(clustersToBeMerged[1])
          
          updatedCentroid <- getCenters(repSeqObj$withinSampleClusters[[name]]$seqmers[repSeqObj$withinSampleClusters[[name]]$clusters %in% clustersToBeMerged[1],],
                                        repSeqObj$withinSampleClusters[[name]]$clusters[repSeqObj$withinSampleClusters[[name]]$clusters %in% clustersToBeMerged[1]])
          
          repSeqObj$withinSampleClusters[[name]]$clusterCenters[as.character(clustersToBeMerged[1]),] <- updatedCentroid
          removedLables <- c(removedLables,idxTobeRemoved)
        }
      }
    
      
      if(length(removedLables) > 0){
        matchingClustersInSamples <- matchingClustersInSamples[- removedLables]
        availSamples <- availSamples[- removedLables]
      }
      #availSamples <- unique(availSamples)
      matchTable[i, match(availSamples, repSeqObj$samNames)] <- as.numeric(matchingClustersInSamples)
    }
    
  }
  
  return(matchTable)
  
}



getMatchingCluster <- function(combinedSams){
  dtocluster<- as.matrix(dist(combinedSams,diag = T,upper=T))[-1,1]
  imin <- names(which(dtocluster==min(dtocluster))) # imin holds cluster number in s2 that is closest to cluster i from samClusterCentroids.
  return(imin)
}





# differential abundance testing for matching clusters across sample groups/conditions ............................................................

# TO DO : following two functions are not for now used for DA testing. They are supposed to do the testing for each subrepetoire by
# performing permutation of cluster assignment labels in the withinSampleCluster assignments. This is not the same as the findDAClusters function which 
# instead implements the comparison based on abundance information gathered for matching clusters and then using standard statistical tests used for differential expression analysis. 

# This would probably be good when there are small number of samples, which is likely in Repseq studies (e.g upto 3 samples)

compareClusterAbundancesPaired <- function(sam1,sam2,clusMatchTable,s1cls,s2cls){
  
  clusterFcs=c()
  clsFCsObservedPermuted = NULL
  
  # comparison of clonal abundance between matching clusters in sam1 and sam2. 
  # This can be done two ways : c2/c1 where c2 and c1 are clone sizes relative to the total repertoire (sum(cloneSizes in cluster)/sum(totalclonSizes) ) 
  # or second way is c2/c1 where c2 and c1 are average clonesize in cluster relative to the total average in the matching clusters being compared..
  # (averageClonesize per cluster) / (averageCloneSize per totalrepertoire)
  
  # 
  
  for (ck in 1:nrow(clusMatchTable)){
    
    consensusK <- clusMatchTable[ck,]
    
    sam1selected <- sam1[s1cls$clusters==consensusK[1],]
    sam2selected <- sam2[s2cls$clusters==consensusK[2],]
    
    
    sam1CountPerClone <- sum(sam1$COUNT) / nrow(sam1)
    sam2CountPerClone <- sum(sam2$COUNT) / nrow(sam2)
    
    sam1SelectedCountPerClone <- sum(sam1selected$COUNT) / nrow(sam1selected)  
    sam2SelectedCountPerClone <- sum(sam2selected$COUNT) / nrow(sam2selected)
    
    sam1SelectedVsTotal <- sam1SelectedCountPerClone/sam1CountPerClone
    sam2SelectedVsTotal <- sam2SelectedCountPerClone/sam2CountPerClone
    
    sam1vssam2Selected <- round(sam2SelectedVsTotal/sam1SelectedVsTotal,digits = 2)
    
    clusterFcs <- c(clusterFcs,sam1vssam2Selected)
    
    
  }
  
  sam1new = sam1cls
  sam2new = sam2cls
  clsFCsObservedPermuted = rbind(clsFCsObservedPermuted,clusterFcs)
  
  for (s in 1:100){
    
    # shuffle cluster assignments and calculate fold change values
    
    sam1new$clusters <- sample(as.numeric(sam1cls$clusters),length(sam1cls$clusters),replace=F)
    sam2new$clusters <- sample(as.numeric(sam2cls$clusters),length(sam2cls$clusters),replace=F)
    
    pClusterFc <- getClusterFoldChanges(sam1,sam2,clusMatchTable,sam1new,sam2new)
    clsFCsObservedPermuted <- rbind(clsFCsObservedPermuted,pClusterFc)
    
  }
  
  # calculate pvalues from permutated values
  #print(head(clsFCs))
  
  ps <- c()
  for(i in 1:ncol(clsFCsObservedPermuted)){
    
    obs <- clsFCsObservedPermuted[1,i]
    otherValues <- clsFCsObservedPermuted[-1,i]
    pval <- sum(abs(otherValues) >= abs(obs)) / nrow(clsFCsObservedPermuted)
    
    #print(pval)
    ps <- c(ps,pval)
  }
  
  #print(clsFCsObservedPermuted)
  return(list(observedClusterFC=clusterFcs,Pval=ps))
  
}


compareClusterAbundancesUnPaired <- function(repSeqObj){
  
  clusterFcs=c()
  clsFCsObservedPermuted = NULL
  numberOfGroups = length(levels(factor(repSeqObj$group)))
  
  if("clusterMatchTable" %in% names(repSeqObj)){
    
    clusMatchTable = repSeqObj$clusterMatchTable
    
  }else
  {
    stop("Cluster match table not found !.") 
    
  }
  
  # comparison of clonal abundance between matching clusters in sam1 and sam2. 
  # This can be done two ways : c2/c1 where c2 and c1 are clone sizes relative to the total repertoire (sum(cloneSizes in cluster)/sum(totalclonSizes) ) 
  # or second way is c2/c1 where c2 and c1 are average clonesize in cluster relative to the total average in the matching clusters being compared..
  # (averageClonesize per cluster) / (averageCloneSize per totalrepertoire)
  
  
  for (ck in 1:nrow(clusMatchTable)){
    
    consensusK <- clusMatchTable[,]
    
    sam1selected <- sam1[s1cls$clusters==consensusK[1],]
    sam2selected <- sam2[s2cls$clusters==consensusK[2],]
    
    
    sam1CountPerClone <- sum(sam1$COUNT) / nrow(sam1)
    sam2CountPerClone <- sum(sam2$COUNT) / nrow(sam2)
    
    sam1SelectedCountPerClone <- sum(sam1selected$COUNT) / nrow(sam1selected)  
    sam2SelectedCountPerClone <- sum(sam2selected$COUNT) / nrow(sam2selected)
    
    sam1SelectedVsTotal <- sam1SelectedCountPerClone/sam1CountPerClone
    sam2SelectedVsTotal <- sam2SelectedCountPerClone/sam2CountPerClone
    
    sam1vssam2Selected <- round(sam2SelectedVsTotal/sam1SelectedVsTotal,digits = 2)
    
    clusterFcs <- c(clusterFcs,sam1vssam2Selected)
    
    
  }
  
  sam1new = sam1cls
  sam2new = sam2cls
  clsFCsObservedPermuted = rbind(clsFCsObservedPermuted,clusterFcs)
  
  for (s in 1:100){
    
    # shuffle cluster assignments and calculate fold change values
    
    sam1new$clusters <- sample(as.numeric(sam1cls$clusters),length(sam1cls$clusters),replace=F)
    sam2new$clusters <- sample(as.numeric(sam2cls$clusters),length(sam2cls$clusters),replace=F)
    
    pClusterFc <- getClusterFoldChanges(sam1,sam2,clusMatchTable,sam1new,sam2new)
    clsFCsObservedPermuted <- rbind(clsFCsObservedPermuted,pClusterFc)
    
  }
  
  # calculate pvalues from permutated values

  ps <- c()
  for(i in 1:ncol(clsFCsObservedPermuted)){
    
    obs <- clsFCsObservedPermuted[1,i]
    otherValues <- clsFCsObservedPermuted[-1,i]
    pval <- sum(abs(otherValues) >= abs(obs)) / nrow(clsFCsObservedPermuted)
    
    #print(pval)
    ps <- c(ps,pval)
  }
  
  #print(clsFCsObservedPermuted)
  return(list(observedClusterFC=clusterFcs,Pval=ps))
  
}




# work on the following to get differentially abundant clusters of clontoypes


getClusterAbundancesTable<- function(repSeqObj){
  
  
  
  if("clusterMatchTable" %in% names(repSeqObj)){
    
    clusMatchTable = repSeqObj$clusterMatchTable
    
  }else{
    
    stop("Cluster match table not found !.") 
    
  }
  
  
  
  cAbundanceTable=NULL 
  cNumberOfClones=NULL 
  cRelativeCloneSizeFrequencyTable=NULL #sum(cloneSizes in cluster) / sum(totalclonSizes) 
  cRelativeAverageCloneSizeTable=NULL #(averageClonesize in  cluster) / (averageCloneSize in total repertoire)
  
  
  
  for(subRep in 1:nrow(clusMatchTable)){
    
    subRepIndex <- clusMatchTable[subRep,]
    
    subRepertoireAbundance = c()
    subRepertoireNumberOfClonotypes = c()
    subRepertoireRelativeCloneSizeFreq = c()
    subRepertoireRelativeAverageCloneSize = c()
    
    
    for(s in names(subRepIndex)){
      sampleSelected = repSeqObj$scaledSampleData[[s]][repSeqObj$withinSampleClusters[[s]]$clusters == subRepIndex[s],]
      
      # record subrepertoire abundance
      subRepertoireAbundance = c(subRepertoireAbundance,sum(sampleSelected$COUNT))
      
      # record number of clonotypes in subrepertoire
      if(is.na(subRepIndex[s])){
        subRepertoireNumberOfClonotypes = c(subRepertoireNumberOfClonotypes,NA)
      }else{
        subRepertoireNumberOfClonotypes = c(subRepertoireNumberOfClonotypes,nrow(sampleSelected))
      }
      
      # record subrepertoire relative clone size frequency
      subRepertoireRelativeCloneSizeFreq = c(subRepertoireRelativeCloneSizeFreq,sum(sampleSelected$COUNT)/sum(repSeqObj$scaledSampleData[[s]]$COUNT))
      
      # record subrepertoire relative average clone size
      
      samSelectedCountPerClone = sum(sampleSelected$COUNT) / nrow(sampleSelected)  
      samCountPerClone = sum(repSeqObj$scaledSampleData[[s]]$COUNT) / nrow(repSeqObj$scaledSampleData[[s]])
      samSelectedVsTotal = samSelectedCountPerClone/samCountPerClone
      
      subRepertoireRelativeAverageCloneSize = c(subRepertoireRelativeAverageCloneSize,samSelectedVsTotal)
      
    }
    
    cAbundanceTable = rbind(cAbundanceTable,subRepertoireAbundance)
    cNumberOfClones = rbind(cNumberOfClones,subRepertoireNumberOfClonotypes)
    cRelativeCloneSizeFrequencyTable = rbind(cRelativeCloneSizeFrequencyTable,subRepertoireRelativeCloneSizeFreq)
    cRelativeAverageCloneSizeTable = rbind(cRelativeAverageCloneSizeTable,subRepertoireRelativeAverageCloneSize)
    
  }
  
  # add new tables to samples object
  # cAbundanceTable
  colnames(cAbundanceTable) = colnames(clusMatchTable)
  rownames(cAbundanceTable) = 1:nrow(clusMatchTable)
  
  repSeqObj <- addItemToObject(repSeqObj,cAbundanceTable,"cAbundanceTable")
  
  # cNumberOfClones
  colnames(cNumberOfClones) = colnames(clusMatchTable)
  rownames(cNumberOfClones) = 1:nrow(clusMatchTable)
  
  repSeqObj <- addItemToObject(repSeqObj,cNumberOfClones,"cNumberOfClones")
  
  # cRelativeCloneSizeFrequencyTable
  colnames(cRelativeCloneSizeFrequencyTable) = colnames(clusMatchTable)
  rownames(cRelativeCloneSizeFrequencyTable) = 1:nrow(clusMatchTable)
  
  repSeqObj <- addItemToObject(repSeqObj,cRelativeCloneSizeFrequencyTable,"cRelativeAbundanceTable")
  
  # cRelativeAverageCloneSizeTable
  colnames(cRelativeAverageCloneSizeTable) = colnames(clusMatchTable)
  rownames(cRelativeAverageCloneSizeTable) = 1:nrow(clusMatchTable)
  
  repSeqObj <- addItemToObject(repSeqObj,cRelativeAverageCloneSizeTable,"cRelativeAverageCloneSizeTable")
  
  
  
  return(repSeqObj)
  
}



#' Find differentially abundant subrepertoires
#' 
#'
findDAClusters <- function(repSeqObj,abundanceType=c("cAbundance","cRelAbundance","cRelCloneSize"),testType=c("t.test", "wilcox.test", "RankProd"), paired=F, ...){
  
  
  # testing implemented only for paired and unpaired two group comparison for now.
  
  
  grpLevels = levels(factor(repSeqObj$group))
  numberOfGroups = length(grpLevels)
  
  # get cluster abundance table using the cluster match table
  
  repSeqObj = getClusterAbundancesTable(repSeqObj)
  
  
  # Perform DA based on cluster abundance table (table type to use is selected by user)
  if(!abundanceType %in% c("cAbundance","cRelAbundance","cRelCloneSize"))
    stop("Abundance type for differential abundance testing not given. Please enter one of cAbundance, cRelAbundance or cRelCloneSize")
    
  if(abundanceType=="cAbundance"){selectedAbundanceTable = "cAbundanceTable"}
  if(abundanceType=="cRelAbundance"){selectedAbundanceTable = "cRelativeAbundanceTable"}
  if(abundanceType=="cRelCloneSize"){selectedAbundanceTable = "cRelativeAverageCloneSizeTable"}
  
  # subselect subrepertoires that will be used for testing based on their existence
  
  # first remove subrepertoires that exist only in one or 0 samples per group
  if(paired==T){
    tooFewSubReps <- apply(repSeqObj[[selectedAbundanceTable]],1,function(x) sum(!is.na(x[repSeqObj$group == grpLevels[2]] / x[repSeqObj$group == grpLevels[1]]) ) > 2)
  }else{
    tooFewSubReps <- apply(repSeqObj[[selectedAbundanceTable]],1,function(x) sum(!is.na(x[repSeqObj$group == grpLevels[1]])) > 2 & sum(!is.na(x[repSeqObj$group == grpLevels[2]])) > 2)
  }
  
  selectedSubrepertoiresTable <- repSeqObj[[selectedAbundanceTable]][tooFewSubReps,,drop=F]
  repSeqObj = addItemToObject(repSeqObj,selectedSubrepertoiresTable,"cSelectedSubRepertoireTable")
  selectedAbundanceTable = "cSelectedSubRepertoireTable"
  
  
  # Choose test type: for now t.test and Rank based Mann-whitney/wilcox, and RankProd implemented. 
 
  if(!testType %in% c("t.test", "wilcox.test", "RankProd")){
    
    cat("Warning: Proper test type not given. RankProd test will be used.")
    
    testType ="RankProd"
    
   }else{  
    
    if(testType[1]=="RankProd"){
      
      if(paired){
        # do one class RankProd analysis after getting the ratios of class 1 to 2 (since we have paired samples)
        d = repSeqObj[[selectedAbundanceTable]][,repSeqObj$group == grpLevels[1]]/repSeqObj[[selectedAbundanceTable]][,repSeqObj$group == grpLevels[2]]
        tocl=rep(1,ncol(d))
      }else{
        # do one class RankProd analysis after getting the ratios of class 1 to 2 (since we have paired samples)
        d = repSeqObj[[selectedAbundanceTable]]
        tocl=c(rep(0,sum(repSeqObj$group == grpLevels[1])),rep(1,sum(repSeqObj$group == grpLevels[2])))
      }
      
      # quantile normalization for rank prod only
      d_n=as.data.frame(normalize.quantiles(as.matrix(d)))
      colnames(d_n)<-colnames(d)
      rownames(d_n)<-rownames(d)
      d<-d_n
      
      sink("rpout");
      Rp.out <- RP(d,cl=tocl,logged=F,plot=F,rand=123,na.rm=T);
      sink(NULL);     
      
      topSubReps=topGene(Rp.out,gene.names=rownames(d),cutoff=1,method="pval",logged=F) 
      
        if(nrow(topSubReps$Table1) > 0){
          #print(head(topSubReps$Table1))
          # sometimes RP misses some row for some reason, for now we just skip those
          midx <- match(rownames(repSeqObj[[selectedAbundanceTable]]),rownames(topSubReps$Table1))
          midx <- midx[!is.na(midx)]
          repSeqObj[[selectedAbundanceTable]] <- repSeqObj[[selectedAbundanceTable]][which((rownames(repSeqObj[[selectedAbundanceTable]]) %in% rownames(topSubReps$Table1))),]
          
          pval = topSubReps$Table1[midx,colnames(topSubReps$Table1)=="P.value"]
        }else{#cat("No differentially abundant cluster (using the RankProd test)")
          }
      
      }else{
        
        if(testType[1]=="t.test"){selectedTest = t.test}
        if(testType[1]=="wilcox.test"){selectedTest = wilcox.test}
        pval=apply(repSeqObj[[selectedAbundanceTable]],1,function(x)selectedTest(x[repSeqObj$group == grpLevels[1]],x[repSeqObj$group == grpLevels[2]],paired=paired,na.rm=T,...)$p.value)
      
      }
    

    }

  
  #Rp.out <- RP(tempd,cl=c(0,0,0,0,1,1,1,1),logged=F,plot=F,rand=123,na.rm=T);
  #topSubReps=topGene(Rp.out,gene.names=rownames(tempd),cutoff=1,method="pval",logged=F) 
  
  #pval = topSubReps$Table1[match(rownames(repSeqObj[["cSelectedSubRepertoireTable"]]),rownames(topSubReps$Table1)),colnames(topSubReps$Table1)=="P.value"]
  
  
  
  # Multiple-testing correction using BH
  
  padj = round(p.adjust(pval, "BH"), 3)
  
  
  # Ratio of mean cAbundances 
  ratioAbundanceMeans = log2(rowMeans(repSeqObj[[selectedAbundanceTable]][,repSeqObj$group == grpLevels[2],drop=F],na.rm=T)/rowMeans(repSeqObj[[selectedAbundanceTable]][,repSeqObj$group == grpLevels[1],drop=F],na.rm=T))
  if(testType[1]=="RankProd")
    ratioAbundanceMeans = log2(1/topSubReps$Table1[rownames(repSeqObj[[selectedAbundanceTable]]),colnames(topSubReps$Table1)=="FC:(class1/class2)"])
  
  cDaResult = data.frame(logFC=ratioAbundanceMeans,pval,padj)
  
  repSeqObj = addItemToObject(repSeqObj,cDaResult,"cDaResult")
  
  #   if(paired==T & numberOfGroups==2){
  #     noSamples = length(repSeqObj$group)
  #     noPairs = noSamples/2
  #     observedFCs = NULL
  #     pvalues = NULL
  #     
  #     for(i in 1:noPairs){
  #       sam1=repSeqObj$scaledSampleData[[i]]
  #       sam1cls=repSeqObj$withinSampleClusters[[i]]
  #       
  #       sam2=repSeqObj$scaledSampleData[[i+noPairs]]
  #       sam2cls=repSeqObj$withinSampleClusters[[i+noPairs]]
  #       
  #       fcComparisonResult=compareClusterAbundances(sam1,sam2,repSeqObj$clusterMatchTable[,c(i,i+noPairs)],sam1cls,sam2cls)
  #       observedFCs =cbind(observedFCs,fcComparisonResult$observedClusterFC)
  #       pvalues=cbind(pvalues,fcComparisonResult$Pval)
  #        
  #     }
  #        
  #   }
  #   else{ #unpaired comparison
  #   fcComparisonResult = compareClusterAbundancesUnPaired(repSeqObj) 
  #   }
  #   
  
  
  return(repSeqObj)
  
  
  
}



# extract subrepertoires of interest ........................................................................

#' Extract DA subrepertoire and CDR3s contained in them
#' 
#'
extractDASubRepertoire <- function(repSeqObj,cutoff=0.1,method="pvalue"){
  
  # choose method
  if(method=="pvalue"){
    colToCheck="pval"
  }else{colToCheck="padj"}
  
  
  # select subrepertoires that showed significant difference in abundance
  
  # we are selecting now subrepertoires that show significant difference in abundance between groups, and enrichment in group2 (a positive fold change in log2)
  #significantSubRepertoires = rownames(repSeqObj$cDaResult)[which(repSeqObj$cDaResult[,colToCheck] < cutoff & repSeqObj$cDaResult$logFC > 0)]
  
  # both enrichment and de-enrichment
  significantSubRepertoires = rownames(repSeqObj$cDaResult)[which(repSeqObj$cDaResult[,colToCheck] < cutoff)]
  
  if(!length(significantSubRepertoires) >=1){
    #cat("No differentially abundant subrepertoires found with given cutoff!\n"); 
    return(repSeqObj)
    }
  
  # DA clonotypes
  DAClonotypes = NULL
  
  for(subR in significantSubRepertoires){
    
    sigSubRepIndex <- repSeqObj$clusterMatchTable[rownames(repSeqObj$clusterMatchTable) == subR,]
    sigSubRepIndex <- sigSubRepIndex[!is.na(sigSubRepIndex)]
    
    for(s in names(sigSubRepIndex)){
      sampleSelected = repSeqObj$scaledSampleData[[s]][repSeqObj$withinSampleClusters[[s]]$clusters == sigSubRepIndex[s],]
      datemp = cbind(sampleSelected$AMINOACID,rep(subR,length(sampleSelected$AMINOACID)))
      DAClonotypes = rbind(DAClonotypes,datemp)
    }  
  }
  
  
  # DA clonotypes in samples
  
  DAClonotypeAbundanceMatrix = NULL
  
  for(s in names(repSeqObj$scaledSampleData)){
    
    DAClonotypeAbundanceMatrix = cbind(DAClonotypeAbundanceMatrix,repSeqObj$scaledSampleData[[s]][match(DAClonotypes[,1],repSeqObj$scaledSampleData[[s]]$AMINOACID),c("COUNT")])
    
  }
  
  
  
  colnames(DAClonotypes) <- c("AMINOACID","Subrepertoire")
  
  rownames(DAClonotypeAbundanceMatrix) = DAClonotypes[,1]
  colnames(DAClonotypeAbundanceMatrix) = names(repSeqObj$scaledSampleData)
  DAClonotypeAbundanceMatrix[is.na(DAClonotypeAbundanceMatrix)] <- 0
  
  
  repSeqObj = addItemToObject(repSeqObj,significantSubRepertoires,"cDaClusters")
  repSeqObj = addItemToObject(repSeqObj,DAClonotypes,"cDaClonotypesWithCluster")
  repSeqObj = addItemToObject(repSeqObj,DAClonotypes[,1],"cDaClonotypes")
  repSeqObj = addItemToObject(repSeqObj,DAClonotypeAbundanceMatrix,"cDaClonotypeAbundanceMatrix")
  
  
return(repSeqObj) 
  
}


#' Extract any subrepertoire and CDR3s contained in them
#' 
#'
extractSubRepertoire <- function(repSeqObj,subReps=NULL){
  
  if(is.null(subReps))
    stop("No subrepertoires given.")
  
  if(!item.exists(repSeqObj,"clusterMatchTable")){
    stop("clusterMatchTable was not detected.")
  }
  
  
  # DA clonotypes
  subRepClonotypes = c()
  
  for(subR in subReps){
    if(!subR %in% rownames(repSeqObj$clusterMatchTable))
      next
    sigSubRepIndex <- repSeqObj$clusterMatchTable[rownames(repSeqObj$clusterMatchTable) == subR,]
    sigSubRepIndex <- sigSubRepIndex[!is.na(sigSubRepIndex)]
    
    for(s in names(sigSubRepIndex)){
      sampleSelected = repSeqObj$scaledSampleData[[s]][repSeqObj$withinSampleClusters[[s]]$clusters == sigSubRepIndex[s],]
      subRepClonotypes = c(subRepClonotypes,sampleSelected$AMINOACID)
    }  
  }
  
  
  # DA clonotypes in samples
  
  DAClonotypeAbundanceMatrix = NULL
  
  for(s in names(repSeqObj$scaledSampleData)){
    
    DAClonotypeAbundanceMatrix = cbind(DAClonotypeAbundanceMatrix,repSeqObj$scaledSampleData[[s]][match(subRepClonotypes,repSeqObj$scaledSampleData[[s]]$AMINOACID),c("COUNT")])
    
  }
  
  
  rownames(DAClonotypeAbundanceMatrix) = subRepClonotypes
  colnames(DAClonotypeAbundanceMatrix) = names(repSeqObj$scaledSampleData)
  DAClonotypeAbundanceMatrix[is.na(DAClonotypeAbundanceMatrix)] <- 0
  
  
  repSeqObj = addItemToObject(repSeqObj,subRepClonotypes,"cSubRepClonotypes")
  repSeqObj = addItemToObject(repSeqObj,DAClonotypeAbundanceMatrix,"subRepClonotypeAbundanceMatrix")
  
  
  return(repSeqObj) 
  
}



findDAClonotypesInDACluster<-function(clonotypes){
  clonotypes=as.vector(clonotypes)
    
  imp_clones_qiao <- read.table("Important clones qiao.txt", header=T, sep="\t",dec = ",")
  imp_clones_qiao=unique(as.vector(imp_clones_qiao[,1]))
  
  imp_clones_arnold1<- read.table("Important clones Arnold_CD4_gluten reactive.txt", header=T, sep="\t",dec = ",") 
  imp_clones_arnold1=unique(as.vector(imp_clones_arnold1[,1]))
  
  imp_clones_arnold2<- read.table("Important clones Arnold_CD8_isolated from challenged patients.txt", header=T, sep="\t",dec = ",") 
  imp_clones_arnold2=unique(as.vector(imp_clones_arnold2[,1]))
  
  imp_clones_arnold_total=c(as.vector(imp_clones_arnold1),as.vector(imp_clones_arnold2))
  
  imp_clones_petersen <- read.table("Important clones Petersen.txt", header=T, sep="\t",dec = ",")
  imp_clones_petersen=unique(as.vector(imp_clones_petersen[,8]))
  
  knownClones <- unique(c(imp_clones_qiao,imp_clones_arnold1,imp_clones_arnold2,imp_clones_petersen))
  
  found_in_knownclones=clonotypes[clonotypes %in% knownClones]
 
  return(found_in_knownclones)
  
}


# fishers exact test ranking:
compareAbundanceInPairedSamplesForRanking <- function(samObj,freqTable,pairs=NULL,paired=T){
  
  ### freqTable : an clone count table between sample groups, table has headers, headers 1 and 2 have the counts
  sam = round(freqTable)
  
  resultlist <- list()
  
  # comparison for every pair of samples

  calculateFisherStats <- function(k,kpair){
    pvals=c()
    ORs=c()
    CIs=c()
    inSample = c()
    NttoAA = c()
    
    for(i in 1:nrow(sam)){
      
      current <- c()
      othersSum <- c()
      
      nttoaa1 <- length(samObj$sampleData[[k]][samObj$sampleData[[k]]$AMINOACID %in% rownames(sam[i,,drop=F]),c("COUNT")])
      nttoaa2 <- length(samObj$sampleData[[kpair]][samObj$sampleData[[kpair]]$AMINOACID %in% rownames(sam[i,,drop=F]),c("COUNT")])
      
      # Nt to AA  of sample 2/sample 1
      if(nttoaa1==0)nttoaa1=1
      NttoAA <- c(NttoAA,nttoaa2/nttoaa1)
      
      # we add the counts at the AA level even if the differential abundance analysis was done at the NT level (clustering etc and cluster abundance was estimated at the Nt level)
      current[1] <- sum(samObj$sampleData[[k]][samObj$sampleData[[k]]$AMINOACID %in% rownames(sam[i,,drop=F]),c("COUNT")])
      current[2] <- sum(samObj$sampleData[[kpair]][samObj$sampleData[[kpair]]$AMINOACID %in% rownames(sam[i,,drop=F]),c("COUNT")])
      
      current[is.na(current)] <- 0
      
      if(current[1]==0 & current[2]==0){
        inSample <- c(inSample,F)
        othersSum[1] <- sum(samObj$sampleData[[k]][!samObj$sampleData[[k]]$AMINOACID %in% rownames(sam[i,,drop=F]),c("COUNT")])
        othersSum[2] <- sum(samObj$sampleData[[kpair]][!samObj$sampleData[[kpair]]$AMINOACID %in% rownames(sam[i,,drop=F]),c("COUNT")])
        
      }else if(current[1]==0 & current[2]!=0){
        inSample <- c(inSample,T)
        othersSum[1] <- sum(samObj$sampleData[[k]][,c("COUNT")])
        othersSum[2] <- sum(samObj$sampleData[[kpair]][!samObj$sampleData[[kpair]]$AMINOACID %in% rownames(sam[i,,drop=F]),c("COUNT")])
        
      }else if(current[1]!=0 & current[2]==0){
        inSample <- c(inSample,T)
        othersSum[1] <- sum(samObj$sampleData[[k]][!samObj$sampleData[[k]]$AMINOACID %in% rownames(sam[i,,drop=F]),c("COUNT")])
        othersSum[2] <- sum(samObj$sampleData[[kpair]][,c("COUNT")])
        
      }else{
        inSample <- c(inSample,T)
        othersSum[1] <- sum(samObj$sampleData[[k]][!samObj$sampleData[[k]]$AMINOACID %in% rownames(sam[i,,drop=F]),c("COUNT")])
        othersSum[2] <- sum(samObj$sampleData[[kpair]][!samObj$sampleData[[kpair]]$AMINOACID %in% rownames(sam[i,,drop=F]),c("COUNT")])
      }
      
      
      dataForComparison=matrix(c(as.numeric(current),as.numeric(othersSum)),nrow = 2,byrow=T)
      dataForComparison[is.na(dataForComparison)] <- 0
      
      dataForComparison = dataForComparison + 1
      dataForComparison <- round(dataForComparison)
      
      fresult=fisher.test(dataForComparison[,c(2,1)])
      pvals=c(pvals,fresult$p.value)
      ORs=c(ORs,as.numeric(fresult$estimate))
      CIs=c(CIs,paste(fresult$conf[1],fresult$conf[2],sep="-"))
    }
    
    res = cbind(inSample*1,NttoAA,as.numeric(pvals),ORs,CIs)
    
    res
  }
    
    
    
    #res=compareAbundanceFexactForRanking(sam[,c(k,kpair)])
  if(paired==T){
    rnds = length(pairs)/2
    for (k in 1:rnds){
      kpair = k + rnds 
      res <- calculateFisherStats(k,kpair)
      colnames(res) <- c("inSample","NTtoAA",paste(colnames(sam[,c(k,kpair)])[1],colnames(sam[,c(k,kpair)])[2],"pvalue",sep="_"),paste(colnames(sam[,c(k,kpair)])[1],colnames(sam[,c(k,kpair)])[2],"OR",sep="_"),paste(colnames(sam[,c(k,kpair)])[1],colnames(sam[,c(k,kpair)])[2],"CI",sep="_"))
      
      resultlist[[k]] <- res
      
    } #fishers exact test done for all clones, all samples
  }else{
    for(g1 in which(pairs==levels(factor(pairs))[1])){
      
      resTemp <- list()
      for(g2 in which(pairs==levels(factor(pairs))[2])){
        k=g1
        kpair=g2
        res <- calculateFisherStats(k,kpair)
        colnames(res) <- c("inSample","NTtoAA",paste(colnames(sam[,c(k,kpair)])[1],colnames(sam[,c(k,kpair)])[2],"pvalue",sep="_"),paste(colnames(sam[,c(k,kpair)])[1],colnames(sam[,c(k,kpair)])[2],"OR",sep="_"),paste(colnames(sam[,c(k,kpair)])[1],colnames(sam[,c(k,kpair)])[2],"CI",sep="_"))
        resTemp[[length(resTemp) + 1]] <- res
        
      }
      
       nSamplesWithTheCloneTemp <- c()
       NttoasInsamplesTemp <- c()
       pvalsTemp <- c()
       orsTemp <- c()
       
        for(i in 1:nrow(sam)){
        
        nSamplesWithTheCloneTemp <- c(nSamplesWithTheCloneTemp,round(mean(as.numeric(sapply(resTemp,function(x) x[i,1])))))
        NttoasInsamplesTemp  <- c(NttoasInsamplesTemp,mean(as.numeric(sapply(resTemp,function(x) x[i,2]))))
        pvalsTemp  <- c(pvalsTemp,mean(as.numeric(sapply(resTemp,function(x) x[i,3]))))
        orsTemp  <- c(orsTemp,mean(as.numeric(sapply(resTemp,function(x) x[i,4]))))
        }
      
        resultlist[[g1]] <- cbind(nSamplesWithTheCloneTemp,NttoasInsamplesTemp,pvalsTemp,orsTemp)
    }
    
  }
  
  
  
  # extract average pval and OR for clone across all samples
  avPvals <- c()
  avORs <- c()
  avNTtoAAs <- c()
  nOfSampesForClones <- c()
  
  for(i in 1:nrow(sam)){
    nSamplesWithTheClone <- as.numeric(sapply(resultlist,function(x) x[i,1]))
    NttoasInsamples <- as.numeric(sapply(resultlist,function(x) x[i,2]))
    pvals <- as.numeric(sapply(resultlist,function(x) x[i,3]))
    ors <- as.numeric(sapply(resultlist,function(x) x[i,4]))
    mpval <- mean(pvals[nSamplesWithTheClone==1])
    mors <- mean(ors[nSamplesWithTheClone==1])
    mNTtoAA <- mean(NttoasInsamples[nSamplesWithTheClone==1])
    
    avPvals<-c(avPvals,mpval)
    avORs<-c(avORs,mors)
    avNTtoAAs<-c(avNTtoAAs,mNTtoAA)
    
    # we give more weight to clones that appear in more number of samples after treatment than before
   
    #nOfSampesForClonesG2 <- c(nOfSampesForClones,sum(as.numeric(sapply(resultlist,function(x) x[i,1]))[pairs==1]))
    #nOfSampesForClonesG1 <- c(nOfSampesForClones,sum(as.numeric(sapply(resultlist,function(x) x[i,1]))[pairs==0]))
    
    nOfSampesForClones <- c(nOfSampesForClones,sum(sam[i,pairs==1] > 0) - sum(sam[i,pairs==0] > 0))
  }
  
  toRet <- data.frame(avNTtoAAs,avPvals,avORs,nOfSampesForClones)
  rownames(toRet) <- rownames(sam)
  
  return(toRet)
  
}



# Run DA analysis function :  .........................................................................

#' This function performs differential abundance analysis between groups of TCRB CDR3 samples (repseq data) to identify differentially abundant (DA) CDR3s. 
#' 
#' @description Function performs clustering based differential abundance analysis of CDR3 sequences in two sample groups with repeat resampling strategy. It first
#' performs within sample unsupervised clustering using subsequence frequency based unsupervised clustering, matches the clusters to their closest match across samples, and performs differential abundance testing at the level of matching
#' clusters to identify differentially abundant condition associated CDR3 sequences
#' @param repSeqObj is an object containing all repertoire sample data
#' @param clusterby boolean; subsequence type to consider, either NT (nucleotide) or AA (amino acid), default is NT
#' @param kmerWidth subsequence width to use, default is 4 for NT, and 3 for AA clusterby
#' @param paired boolean; whether to perform paired analysis for matched datasets,default is true. 
#' @param clusterDaPcutoff sub-repertoire level differential abundance testing cut off, default is 0.1. This works well for our test cases. 
#' @param positionWt boolean; whether to use positional weights for kmer frequencies, default is false
#' @param distMethod the distance method to be used for distance calculation between CDR3 feature vectors, use "euclidean" for nt 4-mer, and "cosine" for aa 3-mer feature vectors
#' @param useDynamicTreeCut boolean; default true, uses Dynamic Tree cut algorithm to cut clustering dendrograms. if false, findOptimalK will be used to find optimal k
#' @param matchingMethod matching method to match cluster centroids from all samples to identify subrepertoires; default is km (kmeans). If hc, hierarchical clustering will be used with 
#' dynamic tree cut to define clusters, if og an in house algorithm will be used that matches each cluster centroid in first sample to their closest centroids in all samples.
#' @param repeatResample boolean; perform repeat resampling, default is true. If false, all repertoire dataset will be used for analysis without downsampling.
#' @param nRepeats number of repeat resample runs to perform if repeatResample is true, default is 10
#' @param resampleSize the downsampling size in the repeat resample runs. default is 5000 
#' @param useProb boolean; if true, probabilistic sampling is performed for downsampling with most frequenty CDR3s being more likely to be resampled. If false, all CDR3s have equal chance of being resampled. Default is true.
#' @param returnAll boolean; if true, the function returns a list whose first element is a data frame with information on the candidate differentially abundant CDR3s and second element is the intermediate results for all repeat resamples. 
#' If false, the intermediate results fro all runs are not returned, only a data frame with candidate differentially abundant CDR3s and their filterig results is returned.
#' @param nRR the number of permutations to perform in the ranking step of candidate DA CDR3s to determine statistical significance.
#' @return a data frame with all candidate DA CDR3s if returnAll is false, a list with data frame of candidate DA CDR3s and all intermediate results if returnAll is true.
#' 
#' @examples results <- runDaAnalysis(repObj,clusterby="NT",kmerWidth=4,paired=T,clusterDaPcutoff=0.1,positionWt = F,distMethod="euclidean",matchingMethod="km",nRepeats=2,resampleSize=1000,useProb=T,returnAll=T,nRR=1000)
#' 
#' @export
# 
runDaAnalysis <- function(repSeqObj,clusterby="NT",kmerWidth=4,paired=T,clusterDaPcutoff=0.1,positionWt = F,distMethod=c("euclidean","cosine"),useDynamicTreeCut=T,matchingMethod="km",repeatResample=T,nRepeats=10,resampleSize=5000,useProb=T,returnAll=T,nRR=1000){
  
  clusterby <- toupper(clusterby)
  if(!clusterby %in% c("NT","AA"))
    stop("clusterby is not given. Given one of NT or AA.")
  
  if(!matchingMethod %in% c("hc","km","og"))
    stop("matchingMethod is not given. Given one of hc, km or og.")
  
  distMethod <- tolower(distMethod)
  if(!distMethod %in% c("euclidean","cosine"))
    stop("distMethod not given. Given one of euclidean or cosine. For better result use euclidean when clustering by NT, and cosine when clustering by AA.")
  
  registerDoParallel(cores=detectCores())  
  
  cDaClonotypesList <- list()
  
  if(repeatResample==T){
    nRepeats = nRepeats
  }else{
    nRepeats = 1
  }
  
  samObjResamples = foreach(i=1:nRepeats,.export=c("DNAStringSet","oligonucleotideFrequency","as","sparseMatrix","dist.matrix","cutreeDynamic","silhouette","summary")) %dopar% {
    
    samObj <- scaleSamples(repSeqObj,totalReads=1e5,resample=repeatResample,resampleSize=resampleSize,useProb=useProb)
    
    if(useDynamicTreeCut==T){
      samObjWithClusters = findOptimalClusters(samObj,k=NULL,clusterby,kmerWidth,positionWt,distMethod,useDynamicTreeCut=useDynamicTreeCut,matchingMethod=matchingMethod)
    }else{
      k = findOptimalK(sam2Obj,nSamEval=4,clusterby,kmerWidth=kmerWidth,posWt=positionWt,distMethod)
      samObjWithClusters = findOptimalClusters(samObj,k,clusterby,kmerWidth,positionWt,distMethod,useDynamicTreeCut=useDynamicTreeCut,matchingMethod=matchingMethod)
    }
    
    # Find DA clusters 
    samObj = findDAClusters(samObjWithClusters,abundanceType="cAbundance",testType="t.test",paired=paired) # first option
    #samObj = findDAClusters(samObjWithClusters,abundanceType="cAbundance",testType="RankProd",paired=T)
    #samObj = findDAClusters(samObjWithClusters,abundanceType="cRelCloneSize",testType="RankProd",paired=T) # second best option
    #samObj = findDAClusters(samObjWithClusters,abundanceType="cRelCloneSize",testType="t.test",paired=F)
    
    #Extract DA clones and association features (Vgenes etc)
    

    samObjWithDas = extractDASubRepertoire(samObj,cutoff=clusterDaPcutoff,method="pvalue")
    
    samObj = samObjWithDas
    
    samObj[4:length(samObj)]
    
  }
  
  # collect DA clones from all subsample runs
  
  for(i in 1:length(samObjResamples)){
    
    if(item.exists(samObjResamples[[i]],"cDaClusters")){
      cDaClonotypesList[[i]] <- samObjResamples[[i]]$cDaClonotypes
      
    }
    
  }
  
  
  allCandidateClones <- unlist(cDaClonotypesList)
  
  if(length(allCandidateClones) > 0 ){
    
    allnRepeats <- nRepeats
    nRepeatsWithHits <- length(cDaClonotypesList)
    
    cat("Candidate DA CDR3s detected from ",nRepeatsWithHits,"out of",allnRepeats,"repeat resamples.","\n Number of candidate CDR3 clonotypes detected as differentially abundant before filtering: ",length(allCandidateClones),"\n")
  }else{
    cat("No differentially abundant clonotypes detected.\n")
    
    if(returnAll==T){
      toRet <- list(cDaClonotypesList,samObjResamples)
      names(toRet) <- c("DaClonotypes","nRepeatResults")
      
    }else{
      toRet <- cDaClonotypesList
    }
    
    return(toRet)
  }
  
  
  ### filtering DA candidates :
  
  # DA clonotype ranking
  allCandidateClones <- unlist(cDaClonotypesList)
  CandidateClFreq <- sort(table(allCandidateClones),decreasing=T)
 
  commDaClones <- data.frame(hitCOUNT=as.numeric(CandidateClFreq))
  rownames(commDaClones) <- names(CandidateClFreq)
  
  #colnames(commDaClones) <- "hitCOUNT"
    
  
  DAClonotypeAbundanceMatrix2 = NULL
  
  for(s in names(repSeqObj$sampleData)){
    
    DAClonotypeAbundanceMatrix2 = cbind(DAClonotypeAbundanceMatrix2,repSeqObj$sampleData[[s]][match(rownames(commDaClones),repSeqObj$sampleData[[s]]$AMINOACID),c("COUNT")])
    
  }
  
  rownames(DAClonotypeAbundanceMatrix2) = rownames(commDaClones)
  colnames(DAClonotypeAbundanceMatrix2) = names(repSeqObj$sampleData)
  DAClonotypeAbundanceMatrix2[is.na(DAClonotypeAbundanceMatrix2)] <- 0
  
  # # select DA clones that show at least one incidence of increased abundance across the samples
  # selectedDAClonesIdx <- apply(DAClonotypeAbundanceMatrix2,1,function(x) sum((x[repSeqObj$group == unique(repSeqObj$group)[2]]/x[repSeqObj$group == unique(repSeqObj$group)[1]]) < 1,na.rm=T) == 0)
  # 
  # DAClonotypeAbundanceMatrix2_selected <- DAClonotypeAbundanceMatrix2[selectedDAClonesIdx,]
  # DAClonotypeAbundanceMatrix2_selected <- round(DAClonotypeAbundanceMatrix2_selected,digits=2)
  # 
  
  # decide enriched and de-enriched using the rank, with out selecting for increasing clones as commented above
  
  DAClonotypeAbundanceMatrix2_selected <- round(DAClonotypeAbundanceMatrix2)
  
  
  # ranking: highest count given top rank for repeat resamples
  selected_commDaClones <- commDaClones[match(rownames(DAClonotypeAbundanceMatrix2_selected),rownames(commDaClones)),,drop=F]
  resampleRank = as.numeric(factor(-selected_commDaClones[,1]))
  
  
  
  #... ranking the DA clonotypes: here use RF and fisher's exact based approach for each clone, write this !!
  # using randomForest
  clonecountM_class=as.data.frame(t(DAClonotypeAbundanceMatrix2_selected))
  classes=factor(repSeqObj$group)
  #clonecountM_class = data.frame(clonecountM_class,classes);
  
  nzaf <- ifelse(ncol(clonecountM_class) > 5000,5000,ncol(clonecountM_class))
  clones.rf <- randomForest(clonecountM_class,classes,ntree=nzaf,importance=TRUE)
  

  varimp=importance(clones.rf,type=1)
  rfRank <- as.numeric(factor(-varimp[,1]))
  
  
  fishExactRes <- compareAbundanceInPairedSamplesForRanking(repSeqObj,DAClonotypeAbundanceMatrix2_selected,repSeqObj$group,paired)
  
  fishExactRes$ntaaRank <- as.numeric(factor(-fishExactRes[,1]))
  fishExactRes$pvalRank <- as.numeric(factor(fishExactRes[,2]))
  fishExactRes$orRank <- as.numeric(factor(-fishExactRes[,3])) # taking log2 to get symmetric OR values
  fishExactRes$nSamRank <- as.numeric(factor(-fishExactRes[,4]))
  
  daClonotypesWithRank <- data.frame(DAClonotypeAbundanceMatrix2_selected,resampleRank,rfRank,ntaaRank=fishExactRes$ntaaRank,fpval=fishExactRes[,2],fPvalRank=fishExactRes$pvalRank,fOr=fishExactRes[,3],fOrRank=fishExactRes$orRank,nSamRank=fishExactRes$nSamRank)
  
  # combine rankings
  scaleRank <- function(r){
    return((r-min(r)) / (max(r)-min(r)))
  }
  
  r1 <- scaleRank(daClonotypesWithRank$rfRank)
  r2 <- scaleRank(daClonotypesWithRank$resampleRank)
  r3 <- scaleRank(daClonotypesWithRank$fPvalRank)
  r4 <- scaleRank(daClonotypesWithRank$fOrRank)
  r5 <- scaleRank(daClonotypesWithRank$ntaaRank)
  r6 <- scaleRank(daClonotypesWithRank$nSamRank)
  
  # rps is between 0 and 6, small rps means high enrichment, high rps means low enrichment
  rps <- r1 + r2 + r3 + r4 + r5 + r6
  rpRanks <- as.numeric(factor(rps))
  
  daClonotypesWithRank$rpRanks <- rpRanks
  
  #shuffle and calculate rpRanks
  permutedRps <- NULL
  tempdForRandomization <- daClonotypesWithRank[,c("resampleRank","rfRank","ntaaRank","fPvalRank","fOrRank","nSamRank")]
  
  # for(i in 1:nRR){
  #   tempdForRandomization <- tempdForRandomization[sample(nrow(tempdForRandomization)),]
  #   tempdForRandomization <- tempdForRandomization[,sample(ncol(tempdForRandomization))]
  #   
  #   r1 <- scaleRank(tempdForRandomization$rfRank)
  #   r2 <- scaleRank(tempdForRandomization$resampleRank)
  #   r3 <- scaleRank(tempdForRandomization$fPvalRank)
  #   r4 <- scaleRank(tempdForRandomization$fOrRank)
  #   r5 <- scaleRank(tempdForRandomization$ntaaRank)
  #   r6 <- scaleRank(tempdForRandomization$nSamRank)
  #   
  #   rpst <-  r1 + r2 + r3 + r4 + r5 + r6
  #   permutedRps <- cbind(permutedRps,rpst)
  #   
  # }
  #
  
  # permutedEnPval <- c()
  # for(i in 1:length(rps)){
  #   
  #   permutedEnPval <- c(permutedEnPval,sum(permutedRps[i,] <= rps[i]) / length(permutedRps[i,]) )
  #   
  # }
  # 
  
  # permutedDeEnPval <- c()
  # for(i in 1:length(rps)){
  #   
  #   permutedDeEnPval <- c(permutedDeEnPval,sum(permutedRps[i,] >= rps[i]) / length(permutedRps[i,]) )
  #   
  # }
  # 
  # 
  
  # Improved calculation of the permutation as well as p-values
  workWithPermMatrices <- function(x){
    tempPermTable <- as.data.frame(x)
    
    r1 <- scaleRank(tempPermTable$rfRank)
    r2 <- scaleRank(tempPermTable$resampleRank)
    r3 <- scaleRank(tempPermTable$fPvalRank)
    r4 <- scaleRank(tempPermTable$fOrRank)
    r5 <- scaleRank(tempPermTable$ntaaRank)
    r6 <- scaleRank(tempPermTable$nSamRank)
    
    rpst <-  r1 + r2 + r3 + r4 + r5 + r6
    return(rpst)
  }
  
  
  permutedRps <- lapply(1:nRR,function(x) workWithPermMatrices(tempdForRandomization[sample(nrow(tempdForRandomization)),sample(ncol(tempdForRandomization))]))
  
  
  permutedEnPval <- c()
  permutedDeEnPval <- c()
  
  for(i in 1:length(rps)){
    collectedPermutedRPST <- sapply(permutedRps,function(y) y[i])
    permutedEnPval <- c(permutedEnPval,(sum(collectedPermutedRPST <= rps[i]) + 1) / (length(collectedPermutedRPST) + 1) )
    permutedDeEnPval <- c(permutedDeEnPval,(sum(collectedPermutedRPST >= rps[i]) + 1) / (length(collectedPermutedRPST) + 1) )
    
  }
  
  
  
  
  permutedEnPval.adjusted=p.adjust(permutedEnPval, "BH") # pvals are adjusted
  permutedDeEnPval.adjusted=p.adjust(permutedDeEnPval, "BH") # pvals are adjusted
  
  #enrichment pvalues
  daClonotypesWithRank$permutedEnPval <- permutedEnPval
  daClonotypesWithRank$permutedEnPval.adjusted <- permutedEnPval.adjusted
  
  #deEnrichment pvalues
  daClonotypesWithRank$permutedDeEnPval <- permutedDeEnPval
  daClonotypesWithRank$permutedDeEnPval.adjusted <- permutedDeEnPval.adjusted
  
  DaClonotypesWithRank.pvalOrdered <- daClonotypesWithRank[order(daClonotypesWithRank$permutedEnPval,decreasing=F),]
  
  
  if(returnAll==T){
    toRet <- list(DaClonotypesWithRank.pvalOrdered,samObjResamples)
    names(toRet) <- c("DaClonotypes","nRepeatResults")
  }else{
    toRet <- DaClonotypesWithRank.pvalOrdered
  }
  
  return(toRet)
}




#' Extract top differentially abundant CDR3 sequences 
#' 
#' This function extracts condition associated CDR3s given the result of the runDaAnalysis function with a given p-value cutoff
#' 
#' @param candidateList the list of candidate DA CDR3s, the return value of runDaAnalysis function
#' @param enriched logical; true returns differentially enriched CDR3s, false returns differentially de-enriched CDR3s. Default is true
#' @param pValueCutoff the cutoff p-value, default is 0.05
#' @return function returns a data frame of the candidate CDR3s that have significant p-value, that is below the pValueCutoff
#' 
#' @examples TopDAClonotypes(results,enriched=T,pValueCutoff=0.05) # results is an object holding the result of running runDaAnalysis
#' 
#' @export
#' 
TopDAClonotypes <- function(candidateList,enriched=T,pValueCutoff=0.05){
 
  if(is.data.frame(candidateList)){
     if(nrow(candidateList) > 0){
      candidateCDR3s <- candidateList
     }else{
       stop("The list does not contain any candidate CDR3s.")
     }
  }else if(nrow(candidateList[[1]]) > 0){
    candidateCDR3s <- candidateList[[1]]
  }else{
    stop("The list does not contain any candidate CDR3s.")
  }
  
  if(enriched==T){
    daCDR3s <- candidateCDR3s[candidateCDR3s$permutedEnPval < pValueCutoff,]
  }else{
    daCDR3s <- candidateCDR3s[candidateCDR3s$permutedDeEnPval < pValueCutoff,]
    
  }
  
  daCDR3s 

}


#' V-gene usage analysis in differentially abundant CDR3s
#' 
#' This function compares observed V-gene frequencies to frequencies obtained from randomly sampled sets of CDR3s from pooled repertoires of all samples to evaluate bias in V-gene usage
#' 
#' @param repSamObj is a repseq data object containing all repertoire sample data 
#' @param DaClonotypes a vector of CDR3 amino acid sequences that are considered to be differentially abundant
#' @param n the number of repeat resamples to perform; each random resample samples CDR3s the same size as the DaClonotypes. Default is 20
#' @return function returns a data frame with the observed frequency, frequencies in random resamples, mean of the frequencies in the random resamples,
#' the p-value estimated from the null distribution of resampled frequencies, and percent increase in v-gene usage for every V-genes used in the DA DR3s
#' 
#' @export
#' 
computeVgeneEnrichement <- function(repSamObj,DaClonotypes,n=20){
  
  getVgenesInsamples <- function(aaClones){
    observedVgenes <- c()
    for(i in names(repSamObj$sampleData)){
      tempD <- repSamObj$sampleData[[i]][repSamObj$sampleData[[i]]$AMINOACID %in% aaClones,]
      observedVgenes <- c(observedVgenes,tempD$VGENENAME)
    }
    observedVgenes
  }
  
  observedVgenes <- getVgenesInsamples(DaClonotypes)
  observedVgenesFreq <- table(observedVgenes)/sum(table(observedVgenes))
  
  vgeneFreqObsAndSampled <- data.frame(observedVgeneFreq=as.numeric(observedVgenesFreq))
  rownames(vgeneFreqObsAndSampled) <- names(observedVgenesFreq)
  
  pooledRep <- getPooledSamples(repSamObj)
  
  for(i in 1:n){
    
    temp = pooledRep;
    randomDrawnClonesidx = sample(1:nrow(temp),length(DaClonotypes),replace=T,prob=pooledRep$FREQUENCYCOUNT)
    randomClones <- temp[randomDrawnClonesidx,]
    
    sampledVgenes <- getVgenesInsamples(randomClones$AMINOACID)
    sampledVgenesFreq <- table(sampledVgenes)/sum(table(sampledVgenes))
    
    sampledVgenesFreq <- sampledVgenesFreq[which(names(sampledVgenesFreq) %in% names(observedVgenesFreq))]
    
    notinsampled <- rep(0,sum(!names(observedVgenesFreq) %in% names(sampledVgenesFreq)))
    names(notinsampled) <- names(observedVgenesFreq)[!names(observedVgenesFreq) %in% names(sampledVgenesFreq)]
    
    sampledVgenesFreq <- c(sampledVgenesFreq,notinsampled)
    sampledVgenesFreq <- sampledVgenesFreq[sort(names(sampledVgenesFreq))]
    
    vgeneFreqObsAndSampled <- cbind(vgeneFreqObsAndSampled,as.numeric(sampledVgenesFreq))
    
  }
  
  
  meanVgeneFreqInRandomSamples <- apply(vgeneFreqObsAndSampled,1,function(x) mean(x[-1]))
  
  VfreqDiffPval <- apply(vgeneFreqObsAndSampled,1,function(x) sum(x[-1] >= x[1]) / length(x[-1]))
  perIncrease <- apply(vgeneFreqObsAndSampled,1,function(x) (x[1]-mean(x[-1]))/mean(x[-1]))
  
  vgeneFreqObsAndSampled$meanVgeneFreqInRandomSamples <- meanVgeneFreqInRandomSamples
  vgeneFreqObsAndSampled$VfreqDiffPval <- VfreqDiffPval
  vgeneFreqObsAndSampled$perIncrease <- perIncrease * 100
  
  rownames(vgeneFreqObsAndSampled) <- gsub("C","", rownames(vgeneFreqObsAndSampled))
  
  return(vgeneFreqObsAndSampled)
}





