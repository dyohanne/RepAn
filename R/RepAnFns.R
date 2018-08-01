# June 18, 2015
#
### ***************************************************************
### ******** Basic analysis for TCR repertoire data **************
### ***************************************************************
#require(ggplot2)

########################### External R packages ###########################

#library(seqRFLP) # write to fasta files
#library(Biostrings) #  for basic sequence manipulation
#library(permute)

loadPacks <- function(package.list = c("ggplot2","seqRFLP","stringr","permute","cluster","data.table","NbClust","doParallel","dendextend","kmer","ape","Matrix","wordspace","dynamicTreeCut","e1071","fclust","randomForest","preprocessCore","gplots")){
  # other probably needed packages that have been removed for now: "Rclusterpp"
  new.packages <-package.list[!(package.list %in% installed.packages()[,"Package"])]
  if(length(new.packages)) try(install.packages(new.packages))
  loadSuccess <- lapply(eval(package.list), require, character.only=TRUE,quietly=TRUE)
}


loadBioconductorPacks <- function(package.list = c("Biostrings","RankProd","preprocessCore","msa","ggseqlogo")){
  new.packages <-package.list[!(package.list %in% installed.packages()[,"Package"])]
  if(length(new.packages)){
    source("http://bioconductor.org/biocLite.R")
    biocLite(new.packages)
  }
  loadSuccess <- lapply(eval(package.list), require, character.only=TRUE,quietly=TRUE)
}



#loadPacks()
#loadBioconductorPacks()


########################### FUNCTIONS #########################


## read data files ####



readSample <- function(samplename,cloneStatus="p"){
  
  # read sample name and change colnames to uppercase
  
  sam <- read.table(samplename, header=T, sep="\t",dec = ".",stringsAsFactors=F)
  #sam <- as.data.frame(fread(samplename, header=T, sep="\t",stringsAsFactors=F))
  
  
  colnames(sam) <- toupper(colnames(sam))
  sam <- sam[,which(colnames(sam)=="NUCLEOTIDE"):length(colnames(sam))] # remove unwanted columns that come with some data e.g container etc

  # Extract cdr3 sequences and put on one column
  CDR3NT = substring(sam$NUCLEOTIDE,sam$VINDEX+1,sam$VINDEX + sam$CDR3LENGTH)
  sam <- data.frame(sam,CDR3NT,stringsAsFactors=F)
  
  # if there is a column named copy, change it to COUNT
  
  if(sum(colnames(sam) == "COPY") == 1){
    
    colnames(sam)[which(colnames(sam) == "COPY")] <- "UNNORMALIZEDCOUNT"
    colnames(sam)[which(colnames(sam) == "NORMALIZEDCOPY")] <- "COUNT"
    colnames(sam)[which(colnames(sam) == "NORMALIZEDFREQUENCY")] <- "FREQUENCYCOUNT"
    
    
  }
  
  # filter sample for only productive or unproductive clonotypes
  
  if(cloneStatus == "p"){
    sam = getProductive(sam)
  }else if(cloneStatus == "unp"){
    sam = getUnproductive(sam)
  }
  
  return(sam)
  
}

readSampleV2 <- function(samplename,cloneStatus="p"){
  
  # read sample name and change colnames to uppercase
  
  #sam <- read.table(samplename, header=T, sep="\t",dec = ".",stringsAsFactors=F)
  #sam <- as.data.frame(fread(samplename, header=T, sep="\t",stringsAsFactors=F,verbose=F,showProgress=F))
  
  sam <- read.table(samplename, header=T, sep="\t",dec = ".",stringsAsFactors=F)
  
  colnames(sam) <- toupper(colnames(sam))
  sam <- sam[,which(colnames(sam)=="NUCLEOTIDE"):length(colnames(sam))] # remove unwanted columns that come with some data e.g container etc
  
  # Extract cdr3 sequences and put on one column
  CDR3NT = substring(sam$NUCLEOTIDE,sam$VINDEX+1,sam$VINDEX + sam$CDR3LENGTH)
  sam <- data.frame(sam,CDR3NT,stringsAsFactors=F)
  
  # if there is a column named copy, change it to COUNT
    
  colnames(sam)[grep("TEMPLATES",colnames(sam))] <- "COUNT"
  colnames(sam)[grep("FREQUENCYCOUNT",colnames(sam))] <- "FREQUENCYCOUNT"
  
  #colnames(sam)[grep("AMINOACID",colnames(sam))] <- "AMINOACID"
  #colnames(sam)[grep("VGENENAME",colnames(sam))] <- "VGENENAME"
  #colnames(sam)[grep("SEQUENCESTATUS",colnames(sam))] <- "SEQUENCESTATUS"

  # filter sample for only productive or unproductive clonotypes
  
  if(cloneStatus == "p"){
    sam = getProductive(sam)
  }else if(cloneStatus == "unp"){
    sam = getUnproductive(sam)
  }
  
  return(sam)
  
}





renameCopyTOCount <- function(samplename){
  # read sample name and change colnames to uppercase
  for(rep in reps){
    sam <- get(rep)
    
    if(sum(colnames(sam) == "COPY") == 1){
      
      colnames(sam)[which(colnames(sam) == "COPY")] <- "UNNORMALIZEDCOUNT"
      colnames(sam)[which(colnames(sam) == "NORMALIZEDCOPY")] <- "COUNT"
      colnames(sam)[which(colnames(sam) == "NORMALIZEDFREQUENCY")] <- "FREQUENCYCOUNT"
      
      
    }
    
  }

}

extractCDR3 <- function(samplename){
  # read sample name and change colnames to uppercase
  for(rep in reps){
    p <- get(rep)
    CDR3NT = substring(p$NUCLEOTIDE,p$VINDEX,p$VINDEX+p$CDR3LENGTH)
    p <- data.frame(p,CDR3NT)
    
  }
  
}





#### Important Data related functions : eg. extract productive clones only etc. ####

getProductive <- function(rep) {
  p <- rep[rep$SEQUENCESTATUS=="Productive" | rep$SEQUENCESTATUS=="In",]
  return(p)
}

getUnproductive <- function(rep) {
  p <- rep[!(rep$SEQUENCESTATUS=="Productive" | rep$SEQUENCESTATUS=="In"),]
  return(p)
}


getOutOfFrame <- function(rep) {
  p <- rep[rep$SEQUENCESTATUS=="Out of frame" | rep$SEQUENCESTATUS=="Out",]
  return(p)
}

getHasStop <- function(rep) {
  p <- rep[rep$SEQUENCESTATUS=="Has stop" | rep$SEQUENCESTATUS=="Stop",]
  return(p)
}

getDefined <- function(rep) {
  p <- rep[rep$VGENENAME != "(undefined)" | rep$VGENENAME!= "unresolved",]
  return(p)
}


#### Getting all V and J genes in our data:  ####

getAllVgenes <- function(reps) 
{
  vgenes_GB=c()
  for(rep in reps){
    p <- get(rep)
    vgenes_GB=union(vgenes_GB,getUniqueVgenes(p))
  }
  
  return(sort(vgenes_GB))
}

getAllVfamilies <- function(reps) 
{
  vfam=c()
  for(rep in reps){
    p <- get(rep)
    vfam=union(vfam,getUniqueVfamilies(p))
  }
  
  return(sort(vfam))
}

getAllJgenes <- function(reps) 
{
  jgenes_GB=c()
  for(rep in reps){
    p <- get(rep)
    jgenes_GB=union(jgenes_GB,getUniqueJgenes(p))
  }
  
  return(sort(jgenes_GB))
}

getAllVJgenes <- function(reps) 
{
  jgenes=getAllJgenes(reps)
  vgenes=getAllVgenes(reps)
  ## All possible VJ combinations 
  vj <- expand.grid(vgenes,jgenes)
  colnames(vj) <- c("VGENENAME","JGENENAME")
  return(vj)
}


getUniqueVgenes <- function(rep) {
  p <- rep[rep$VGENENAME != "(undefined)" | rep$VGENENAME != "unresolved",]
  return(unique(p$VGENENAME)) # [-1] -1 inorder to drop the undefined Vgene type
}

getUniqueVfamilies <- function(rep) {
  p <- rep[rep$VFAMILYNAME != "(undefined)" | rep$VFAMILYNAME != "unresolved",]
  return(unique(p$VFAMILYNAME)) # [-1] -1 inorder to drop the undefined Vgene type
}

getUniqueVfamilies_PBMC <- function(rep) {
  p <- rep[rep$VFAMILYNAME != "(undefined)" | rep$VFAMILYNAME != "unresolved",]
  return(unique(p$VFAMILYNAME)) # [-1] -1 inorder to drop the undefined Vgene type
}

getUniqueJgenes <- function(rep) {
  p <- rep[rep$JGENENAME != "(undefined)" | rep$JGENENAME != "unresolved",]
  return(unique(p$JGENENAME))
}


#### V, J or VJ usage counts and proportions : 

# clone level
getVgeneCloneCountPMatrix<- function(reps,g){
  ### Returns a datamatrix whose first two colums show the V gene. The remaining columns show 
  ### The number of available clones using each v gene for each sample
  
  #reps = a vector of dataframe names. The data frames are TCR repertoire data. 
  
  m <- as.data.frame(g)
  colnames(m)  <- c("VGENENAME")
  for(rep in reps){
    print(rep)
    p <- get(rep)
    p <- getProductive(p)
    print (nrow(p))
    #p <- getDefined(p) 
    vfamAgg <- aggregate(AMINOACID~ VGENENAME ,data = p, FUN = length)
    vfamAgg[,2] <- vfamAgg[,2]/sum(vfamAgg[,2]) # get the proportion
    colnames(vfamAgg) <- c("VGENENAME",rep)  
    m <- merge(m,vfamAgg, by="VGENENAME",all.x=T) 
    
  }
  
 
  m[is.na(m)] <- 0
  
  rownames(m) <- m[,1]
  m <- m[,-1]
  
  return(m)
}

getVgeneCloneCountMatrix <- function(reps,g){
  ### Returns a datamatrix whose first two colums show the V gene. The remaining columns show 
  ### The number of available clones using each v gene for each sample
  
  #reps = a vector of dataframe names. The data frames are TCR repertoire data
  
  m <- as.data.frame(g)
  colnames(m)  <- c("VGENENAME")
  for(rep in reps){
    print(rep)
    p <- get(rep)
    p <- getProductive(p)
    print (nrow(p))
    #p <- getDefined(p) 
    vgeneAgg <- aggregate(AMINOACID~ VGENENAME ,data = p, FUN = length)
    colnames(vgeneAgg) <- c("VGENENAME",rep)  
    m <- merge(m,vgeneAgg, by="VGENENAME",all.x=T) 
    
  }
  m[is.na(m)] <- 0
  
  rownames(m) <- m[,1]
  m <- m[,-1]
  return(m)
}


getVgeneCloneCountPMatrixUnprod<- function(reps,g){
  ### Returns a datamatrix whose first two colums show the V gene. The remaining columns show 
  ### The number of available clones using each v gene for each sample
  
  #reps = a vector of dataframe names. The data frames are TCR repertoire data
  
  m <- as.data.frame(g)
  colnames(m)  <- c("VGENENAME")
  for(rep in reps){
    print(rep)
    p <- get(rep)
    p <- getUnproductive(p)
    print (nrow(p))
    #p <- getDefined(p) 
    vfamAgg <- aggregate(AMINOACID~ VGENENAME ,data = p, FUN = length)
    vfamAgg[,2] <- vfamAgg[,2]/sum(vfamAgg[,2]) # get the proportion
    colnames(vfamAgg) <- c("VGENENAME",paste(rep)) 
    m <- merge(m,vfamAgg, by="VGENENAME",all.x=T) 
    
  }
  m[is.na(m)] <- 0
  
  rownames(m) <- m[,1]
  m <- m[,-1]
  return(m)
}

getVgeneCloneCountMatrixUnprod <- function(reps,g){
  ### Returns a datamatrix whose first two colums show the V gene. The remaining columns show 
  ### The number of available clones using each v gene for each sample
  
  #reps = a vector of dataframe names. The data frames are TCR repertoire data
  
  m <- as.data.frame(g)
  colnames(m)  <- c("VGENENAME")
  for(rep in reps){
    print(rep)
    p <- get(rep)
    p <- getUnproductive(p)
    print (nrow(p))
    #p <- getDefined(p) 
    vgeneAgg <- aggregate(AMINOACID~ VGENENAME ,data = p, FUN = length)
    colnames(vgeneAgg) <- c("VGENENAME",rep)  
    m <- merge(m,vgeneAgg, by="VGENENAME",all.x=T) 
    
  }
  m[is.na(m)] <- 0
  rownames(m) <- m[,1]
  m <- m[,-1]
  return(m)
}


getVjCloneCountMatrix<- function(reps,vj){
  ### Returns a datamatrix whose first two colums show the V gene. The remaining columns show 
  ### The number of available clones using each v gene for each sample
  
  #reps = a vector of dataframe names. The data frames are TCR repertoire data
  
  m <- vj
  colnames(m)  <- c("VGENENAME","JGENENAME")
  for(rep in reps){
    print(rep)
    p <- get(rep)
    p <- getProductive(p)
    #p <- getDefined(p) 
    vjPairAgg <- aggregate(AMINOACID~ VGENENAME + JGENENAME   ,data = p, FUN = length)
    print(head(vjPairAgg))
    colnames(vjPairAgg) <- c("VGENENAME","JGENENAME",rep)  
    m <- merge(m,vjPairAgg, by=c("VGENENAME","JGENENAME"),all.x=T) 
    
  }
  m[is.na(m)] <- 0
  m[,1]<-paste(m[,1],m[,2])
  m<-m[,-2]
  
  rownames(m) <- m[,1]
  m <- m[,-1]
  return(m)
}

getVjCloneCountMatrixUnprod<- function(reps,vj){
  ### Returns a datamatrix whose first two colums show the V gene. The remaining columns show 
  ### The number of available clones using each v gene for each sample
  
  #reps = a vector of dataframe names. The data frames are TCR repertoire data
  
  m <- vj
  colnames(m)  <- c("VGENENAME","JGENENAME")
  
  for(rep in reps){
    print(rep)
    p <- get(rep)
    p <- getUnproductive(p)
    #p <- getDefined(p) 
    vjPairAgg <- aggregate(AMINOACID~ VGENENAME + JGENENAME   ,data = p, FUN = length)
    colnames(vjPairAgg) <- c("VGENENAME","JGENENAME",rep)  
    m <- merge(m,vjPairAgg, by=c("VGENENAME","JGENENAME"),all.x=T) 
    
  }
  m[is.na(m)] <- 0
  m[,1]<-paste(m[,1],m[,2])
  m<-m[,-2]
  
  rownames(m) <- m[,1]
  m <- m[,-1]
  return(m)
}

getVjCloneCountPMatrix<- function(reps,vj){
  ### Returns a datamatrix whose first two colums show the V gene. The remaining columns show 
  ### The number of available clones using each v gene for each sample
  
  #reps = a vector of dataframe names. The data frames are TCR repertoire data
  
  m <- vj
  colnames(m)  <- c("VGENENAME","JGENENAME")
  
  for(rep in reps){
    print(rep)
    p <- get(rep)
    p <- getProductive(p)
    #p <- getDefined(p) 
    vjPairAgg <- aggregate(AMINOACID~ VGENENAME + JGENENAME   ,data = p, FUN = length)
    vjPairAgg[,3] <- vjPairAgg[,3]/sum(vjPairAgg[,3]) # get the proportion
    colnames(vjPairAgg) <- c("VGENENAME","JGENENAME",rep)  
    m <- merge(m,vjPairAgg, by=c("VGENENAME","JGENENAME"),all.x=T) 
    
  }
  m[is.na(m)] <- 0
  m[,1]<-paste(m[,1],m[,2])
  m<-m[,-2]
  
  rownames(m) <- m[,1]
  m <- m[,-1]
  return(m)
}

getVjCloneCountPMatrixUnprod<- function(reps,vj){
  ### Returns a datamatrix whose first two colums show the V gene. The remaining columns show 
  ### The number of available clones using each v gene for each sample
  
  #reps = a vector of dataframe names. The data frames are TCR repertoire data
  
  m <- vj
  colnames(m)  <- c("VGENENAME","JGENENAME")
  
  for(rep in reps){
    print(rep)
    p <- get(rep)
    p <- getUnproductive(p)
    #p <- getDefined(p) 
    vjPairAgg <- aggregate(AMINOACID~ VGENENAME + JGENENAME   ,data = p, FUN = length)
    vjPairAgg[,3] <- vjPairAgg[,3]/sum(vjPairAgg[,3]) # get the proportion
    colnames(vjPairAgg) <- c("VGENENAME","JGENENAME",rep)  
    m <- merge(m,vjPairAgg, by=c("VGENENAME","JGENENAME"),all.x=T) 
    
  }
  m[is.na(m)] <- 0
  m[,1]<-paste(m[,1],m[,2])
  m<-m[,-2]
  
  rownames(m) <- m[,1]
  m <- m[,-1]
  return(m)
}


# read level

getVgeneReadCountPMatrix<- function(reps,g){
  ### Returns a datamatrix whose first two colums show the V gene. The remaining columns show 
  ### The number of available clones using each v gene for each sample
  
  #reps = a vector of dataframe names. The data frames are TCR repertoire data. 
  
  m <- as.data.frame(g)
  colnames(m)  <- c("VGENENAME")
  for(rep in reps){
    print(rep)
    p <- get(rep)
    p <- getProductive(p)
    print (nrow(p))
    #p <- getDefined(p) 
    vfamAgg <- aggregate(COUNT~ VGENENAME ,data = p, FUN = sum)
    vfamAgg[,2] <- vfamAgg[,2]/sum(vfamAgg[,2]) # get the proportion
    colnames(vfamAgg) <- c("VGENENAME",rep)  
    m <- merge(m,vfamAgg, by="VGENENAME",all.x=T) 
    
  }
  m[is.na(m)] <- 0
  
  rownames(m) <- m[,1]
  m <- m[,-1]
  return(m)
}

getVgeneReadCountMatrix <- function(reps,g){
  ### Returns a datamatrix whose first two colums show the V gene. The remaining columns show 
  ### The number of available clones using each v gene for each sample
  
  #reps = a vector of dataframe names. The data frames are TCR repertoire data
  
  m <- as.data.frame(g)
  colnames(m)  <- c("VGENENAME")
  for(rep in reps){
    print(rep)
    p <- get(rep)
    p <- getProductive(p)
    print (nrow(p))
    #p <- getDefined(p) 
    vgeneAgg <- aggregate(COUNT~ VGENENAME ,data = p, FUN = sum)
    colnames(vgeneAgg) <- c("VGENENAME",rep)  
    m <- merge(m,vgeneAgg, by="VGENENAME",all.x=T) 
    
  }
  m[is.na(m)] <- 0
  rownames(m) <- m[,1]
  m <- m[,-1]
  return(m)
}

getVgeneReadCountPMatrixUnprod<- function(reps,g){
  ### Returns a datamatrix whose first two colums show the V gene. The remaining columns show 
  ### The number of available clones using each v gene for each sample
  
  #reps = a vector of dataframe names. The data frames are TCR repertoire data
  
  m <- as.data.frame(g)
  colnames(m)  <- c("VGENENAME")
  for(rep in reps){
    print(rep)
    p <- get(rep)
    p <- getUnproductive(p)
    print (nrow(p))
    #p <- getDefined(p) 
    vfamAgg <- aggregate(COUNT~ VGENENAME ,data = p, FUN = sum)
    vfamAgg[,2] <- vfamAgg[,2]/sum(vfamAgg[,2]) # get the proportion
    colnames(vfamAgg) <- c("VGENENAME",paste(rep)) 
    m <- merge(m,vfamAgg, by="VGENENAME",all.x=T) 
    
  }
  m[is.na(m)] <- 0
  rownames(m) <- m[,1]
  m <- m[,-1]
  return(m)
}

getVgeneReadCountMatrixUnprod <- function(reps,g){
  ### Returns a datamatrix whose first two colums show the V gene. The remaining columns show 
  ### The number of available clones using each v gene for each sample
  
  #reps = a vector of dataframe names. The data frames are TCR repertoire data
  
  m <- as.data.frame(g)
  colnames(m)  <- c("VGENENAME")
  for(rep in reps){
    print(rep)
    p <- get(rep)
    p <- getUnproductive(p)
    print (nrow(p))
    #p <- getDefined(p) 
    vgeneAgg <- aggregate(COUNT~ VGENENAME ,data = p, FUN = sum)
    colnames(vgeneAgg) <- c("VGENENAME",rep)  
    m <- merge(m,vgeneAgg, by="VGENENAME",all.x=T) 
    
  }
  m[is.na(m)] <- 0
  
  rownames(m) <- m[,1]
  m <- m[,-1]
  return(m)
}

getVjReadCountMatrix<- function(reps,vj){
  ### Returns a datamatrix whose first two colums show the V gene. The remaining columns show 
  ### The number of available clones using each v gene for each sample
  
  #reps = a vector of dataframe names. The data frames are TCR repertoire data
  
  m <- vj
  colnames(m)  <- c("VGENENAME","JGENENAME")
  
  for(rep in reps){
    print(rep)
    p <- get(rep)
    p <- getProductive(p)
    #p <- getDefined(p) 
    vjPairAgg <- aggregate(COUNT~ VGENENAME + JGENENAME   ,data = p, FUN = sum)
    colnames(vjPairAgg) <- c("VGENENAME","JGENENAME",rep)  
    m <- merge(m,vjPairAgg, by=c("VGENENAME","JGENENAME"),all.x=T) 
    
  }
  m[is.na(m)] <- 0
  m[,1]<-paste(m[,1],m[,2])
  m<-m[,-2]
  
  rownames(m) <- m[,1]
  m <- m[,-1]
  return(m)
}

getVjReadCountMatrixUnprod<- function(reps,vj){
  ### Returns a datamatrix whose first two colums show the V gene. The remaining columns show 
  ### The number of available clones using each v gene for each sample
  
  #reps = a vector of dataframe names. The data frames are TCR repertoire data
  
  m <- vj
  colnames(m)  <- c("VGENENAME","JGENENAME")
  
  for(rep in reps){
    print(rep)
    p <- get(rep)
    p <- getUnproductive(p)
    #p <- getDefined(p) 
    vjPairAgg <- aggregate(COUNT~ VGENENAME + JGENENAME   ,data = p, FUN = sum)
    colnames(vjPairAgg) <- c("VGENENAME","JGENENAME",rep)  
    m <- merge(m,vjPairAgg, by=c("VGENENAME","JGENENAME"),all.x=T) 
    
  }
  m[is.na(m)] <- 0
  m[,1]<-paste(m[,1],m[,2])
  m<-m[,-2]
  
  rownames(m) <- m[,1]
  m <- m[,-1]
  return(m)
}

getVjReadCountPMatrix<- function(reps,vj){
  ### Returns a datamatrix whose first two colums show the V gene. The remaining columns show 
  ### The number of available clones using each v gene for each sample
  
  #reps = a vector of dataframe names. The data frames are TCR repertoire data
  
  m <- vj
  colnames(m)  <- c("VGENENAME","JGENENAME")
  
  for(rep in reps){
    print(rep)
    p <- get(rep)
    p <- getProductive(p)
    #p <- getDefined(p) 
    vjPairAgg <- aggregate(COUNT~ VGENENAME + JGENENAME   ,data = p, FUN = sum)
    vjPairAgg[,3] <- vjPairAgg[,3]/sum(vjPairAgg[,3]) # get the proportion
    colnames(vjPairAgg) <- c("VGENENAME","JGENENAME",rep)  
    m <- merge(m,vjPairAgg, by=c("VGENENAME","JGENENAME"),all.x=T) 
    
  }
  m[is.na(m)] <- 0
  m[,1]<-paste(m[,1],m[,2])
  m<-m[,-2]
  
  rownames(m) <- m[,1]
  m <- m[,-1]
  return(m)
}

getVjReadCountPMatrixUnprod<- function(reps,vj){
  ### Returns a datamatrix whose first two colums show the V gene. The remaining columns show 
  ### The number of available clones using each v gene for each sample
  
  #reps = a vector of dataframe names. The data frames are TCR repertoire data
  
  m <- vj
  colnames(m)  <- c("VGENENAME","JGENENAME")
  
  for(rep in reps){
    print(rep)
    p <- get(rep)
    p <- getUnproductive(p)
    #p <- getDefined(p) 
    vjPairAgg <- aggregate(COUNT~ VGENENAME + JGENENAME   ,data = p, FUN = sum)
    vjPairAgg[,3] <- vjPairAgg[,3]/sum(vjPairAgg[,3]) # get the proportion
    colnames(vjPairAgg) <- c("VGENENAME","JGENENAME",rep)  
    m <- merge(m,vjPairAgg, by=c("VGENENAME","JGENENAME"),all.x=T) 
    
  }
  m[is.na(m)] <- 0
  m[,1]<-paste(m[,1],m[,2])
  m<-m[,-2]
  
  rownames(m) <- m[,1]
  m <- m[,-1]
  return(m)
}



# General Counting function for V, J, VJ genes for each sequence status  ####

getAggregateCount <- function(reps,countfor="V",seqStatus="IN",lev="A",countsorpro="prop"){

  ### Returns a datamatrix whose first two colums show the features The remaining columns show 
  ### The number of available clones for each feature :
  
  #reps - a vector of dataframe names. The data frames are TCR repertoire data.
  # countfor - the features to be counted, default V for VgeneNames, vf for V family, J for J gene, VJ for vj genes
  # seqStatus - for which sequences is the count to be done. default In for productive, Out for out of frame, Stop for stop codon containing sequences
  # lev - at which level is the counting to be done . default A for Amino acid, N for nucleotide level
  # countsorpro - return counts raw counts or proportions
  
  # setting work parameters
  countfor = toupper(countfor)
  seqStatus= toupper(seqStatus)
  lev = toupper(lev)
  countsorpro= toupper(countsorpro)
  
  print(countfor)
  
  if(countfor=="V"){
    
    
    features= getAllVgenes(reps)
    
    m <- as.data.frame(features)
    colnames(m)  <- c("VGENENAME")
    
    
      for(rep in reps){
        print(rep)
        p <- get(rep)
        
        if(seqStatus=="IN"){
          p <- getProductive(p)
        }else if(seqStatus=="OUT"){
          p <- getOutOfFrame(p)
        }else if(seqStatus=="STOP"){
          p <- getHasStop(p)
        }
        
        print (nrow(p))
     
        if(lev=="N"){
        dAgg <- aggregate(NUCLEOTIDE~ VGENENAME ,data = p, FUN = length)
        }else {
        dAgg <- aggregate(AMINOACID~ VGENENAME ,data = p, FUN = length) 
        }
        
        if(countsorpro=="PROP"){
        dAgg[,2] <- dAgg[,2]/sum(dAgg[,2]) # get the proportion
        }
        colnames(dAgg) <- c("VGENENAME",rep)  
        m <- merge(m,dAgg, by="VGENENAME",all.x=T) 
        
      }
      m[is.na(m)] <- 0
      
      rownames(m) <- m[,1]
      m <- m[,-1]
      
     
     
  }
  else if(countfor=="J"){
    features= getAllJgenes(reps)
    
    m <- as.data.frame(features)
    colnames(m)  <- c("JGENENAME")
    
    for(rep in reps){
      print(rep)
      p <- get(rep)
      
      if(seqStatus=="IN"){
        p <- getProductive(p)
      }else if(seqStatus=="OUT"){
        p <- getOutOfFrame(p)
      }else if(seqStatus=="STOP"){
        p <- getHasStop(p)
      }
      
      print (nrow(p))
      
      if(lev=="N"){
        dAgg <- aggregate(NUCLEOTIDE~ JGENENAME ,data = p, FUN = length)
      }else {
        dAgg <- aggregate(AMINOACID~ JGENENAME ,data = p, FUN = length) 
      }
      
      if(countsorpro=="PROP"){
      dAgg[,2] <- dAgg[,2]/sum(dAgg[,2]) # get the proportion
      }
      colnames(dAgg) <- c("JGENENAME",rep)  
      m <- merge(m,dAgg, by="JGENENAME",all.x=T) 
      
    }
    m[is.na(m)] <- 0
    
    rownames(m) <- m[,1]
    m <- m[,-1]
    
    
    
  }
  
  else if(countfor=="VF"){
    features= getAllVfamilies(reps)
    
    m <- as.data.frame(features)
    colnames(m)  <- c("VFAMILYNAME")
    
    for(rep in reps){
      print(rep)
      p <- get(rep)
      
      if(seqStatus=="IN"){
        p <- getProductive(p)
      }else if(seqStatus=="OUT"){
        p <- getOutOfFrame(p)
      }else if(seqStatus=="STOP"){
        p <- getHasStop(p)
      }
      
      print (nrow(p))
      
      if(lev=="N"){
        dAgg <- aggregate(NUCLEOTIDE~ VFAMILYNAME ,data = p, FUN = length)
      }else {
        dAgg <- aggregate(AMINOACID~ VFAMILYNAME ,data = p, FUN = length) 
      }
      
      if(countsorpro=="PROP"){
      dAgg[,2] <- dAgg[,2]/sum(dAgg[,2]) # get the proportion
      }
      colnames(dAgg) <- c("VFAMILYNAME",rep)  
      m <- merge(m,dAgg, by="VFAMILYNAME",all.x=T) 
      
    }
    m[is.na(m)] <- 0
    
    rownames(m) <- m[,1]
    m <- m[,-1]
    
     
  }
  
  else if(countfor=="VJ"){
    features= getAllVJgenes(reps)
    
    m <- features
    colnames(m)  <- c("VGENENAME","JGENENAME")
    
    for(rep in reps){
      print(rep)
      p <- get(rep)
      
      if(seqStatus=="IN"){
        p <- getProductive(p)
      }else if(seqStatus=="OUT"){
        p <- getOutOfFrame(p)
      }else if(seqStatus=="STOP"){
        p <- getHasStop(p)
      }
      
      print (nrow(p))
      
      if(lev=="N"){
        dAgg <- aggregate(NUCLEOTIDE~ VGENENAME + JGENENAME ,data = p, FUN = length)
      }else {
        dAgg <- aggregate(AMINOACID~ VGENENAME + JGENENAME ,data = p, FUN = length) 
      }
      
      if(countsorpro=="PROP"){
      dAgg[,3] <- dAgg[,3]/sum(dAgg[,3]) # get the proportion
      }
      colnames(dAgg) <- c("VGENENAME","JGENENAME",rep)  
      m <- merge(m,dAgg, by=c("VGENENAME","JGENENAME"),all.x=T) 
      
    }
    m[is.na(m)] <- 0
    m[,1]<-paste(m[,1],m[,2])
    m<-m[,-2]
    
    rownames(m) <- m[,1]
    m <- m[,-1]
    
  }
  
  
  return(m)
 
}


# get overlapping clonotype counts for every pair of samples: using top 5% or 10% etc #### 

getAbundantClonotypesbyperc<- function(rep,topCutoff=10){
  
  # Return the top abundant clonotypes using the topcutoff , default is top 5%
  
  print(nrow(rep))
  y=rep$COUNT/sum(rep$COUNT)
  y=cumsum(y) #cumulitive read size proportion at each clone, ordered from biggest to smallest
  
  cuttingPoints=topCutoff/100
  y=cut(y,seq(0,1,by=cuttingPoints)) # Cut the cumulitive proportion values in y by 10% each
  
  
  # We can split the whole data frame using the 5% split factor
  
  rep_subgroups=split(rep,y)
  
  
  return(rep_subgroups[[1]])
  
  
}

getAbundantClonotypesbynumber<- function(rep,topCutoff=10){
  
  # Return the top abundant clonotypes using the topcutoff , default is top 10
  # it assumes rep is already sorted from biggest clone to smallest
  
  print(nrow(rep))
  
  return(rep[1:topCutoff,])
  
  
}


overLappingClonotypes<- function(reps,cutoff=20,cutoffby="perc"){
  ### Returns a datamatrix of top topnumber clonotopes from each sample plus: overlap count or Jaccard index, for prod and non productive sequences, N insertion size
  
  
  # setting work parameters
  cutoffby = toupper(cutoffby)
 
  
  m <- as.data.frame(expand.grid(reps,reps))
  colnames(m)  <- c("Sample1","Sample2")
  
  toadd= data.frame()
  
 
  
  for(r in 1:nrow(m)){
    
    overlapValuesProd=c()
    overlapValuesUnprod=c()
    overlapValuesProdJI=c()
    overlapValuesUnprodJI=c()
    
    overlapValuesProdJI_abundancebased=c()
    overlapValuesUnprodJI_abundancebased=c()
    
    NinsertionProd1=c()
    NinsertionUnprod1=c()
    
    NinsertionProd2=c()
    NinsertionUnprod2=c()
    
    
    print(r)
    sam1=m[r,1]
    sam2=m[r,2]
    
    #print(sam1)
    #print(sam2)
  
    sam1 <- get(as.character(sam1))
    sam2 <- get(as.character(sam2))
    
    
    sam1Prod <- getProductive(sam1)
    sam1Unprod <- getUnproductive(sam1)
    
    sam2Prod <- getProductive(sam2)
    sam2Unprod <- getUnproductive(sam2)
   
 
    if(cutoffby=="PERC"){
      topclones_sam1Prod = getAbundantClonotypesbyperc(sam1Prod,cutoff)
      topclones_sam1Unprod = getAbundantClonotypesbyperc(sam1Unprod,cutoff)
      
      topclones_sam2Prod = getAbundantClonotypesbyperc(sam2Prod,cutoff)
      topclones_sam2Unprod = getAbundantClonotypesbyperc(sam2Unprod,cutoff)
      
      # commons
      commonClonesProd=length(intersect(topclones_sam1Prod$NUCLEOTIDE,topclones_sam2Prod$NUCLEOTIDE))
      unionProd=length(union(topclones_sam1Prod$NUCLEOTIDE,topclones_sam2Prod$NUCLEOTIDE))
      
      # commonality checked at the Nt level for non productive sequences, no AA given in the original data
      commonClonesUnprod=length(intersect(topclones_sam1Unprod$NUCLEOTIDE,topclones_sam2Unprod$NUCLEOTIDE))
      unionUnprod=length(union(topclones_sam1Unprod$NUCLEOTIDE,topclones_sam2Unprod$NUCLEOTIDE))
      
      
      overlapValuesProd=c(overlapValuesProd,as.numeric(commonClonesProd))
      overlapValuesUnprod=c(overlapValuesUnprod,as.numeric(commonClonesUnprod))
      
      overlapValuesProdJI=c(overlapValuesProdJI,formatC(commonClonesProd/unionProd,digits= 8,format="f"))
      overlapValuesUnprodJI = c(overlapValuesUnprodJI,formatC(commonClonesUnprod/unionUnprod,digits= 8,format="f"))
      
      # Jaccard abundance-based similarity index (UV/(U+V-UV)). U is total relative abundance in of shared sample 1, U total relative of shared in sample 2
      
      commonClonesProd_U_data= topclones_sam1Prod[topclones_sam1Prod$NUCLEOTIDE %in% intersect(topclones_sam1Prod$NUCLEOTIDE,topclones_sam2Prod$NUCLEOTIDE),]
      commonClonesProd_U = sum(commonClonesProd_U_data$COUNT)/sum(topclones_sam1Prod$COUNT)
      
      commonClonesProd_V_data= topclones_sam2Prod[topclones_sam2Prod$NUCLEOTIDE %in% intersect(topclones_sam1Prod$NUCLEOTIDE,topclones_sam2Prod$NUCLEOTIDE),]
      commonClonesProd_V = sum(commonClonesProd_V_data$COUNT)/sum(topclones_sam2Prod$COUNT)
      
      uv=commonClonesProd_U*commonClonesProd_V
      uplusv=commonClonesProd_U+commonClonesProd_V
      
      AbundanceBasedJIProd = format(uv/(uplusv-uv),digits=8,format="f")
      
      overlapValuesProdJI_abundancebased=c(overlapValuesProdJI_abundancebased,as.numeric(AbundanceBasedJIProd))
      
        # for unproductive sequences
      
      commonClonesUnProd_U_data= topclones_sam1Unprod[topclones_sam1Unprod$NUCLEOTIDE %in% intersect(topclones_sam1Unprod$NUCLEOTIDE,topclones_sam2Unprod$NUCLEOTIDE),]
      commonClonesUnProd_U = sum(commonClonesUnProd_U_data$COUNT)/sum(topclones_sam1Unprod$COUNT)
      
      commonClonesUnProd_V_data= topclones_sam2Unprod[topclones_sam2Unprod$NUCLEOTIDE %in% intersect(topclones_sam1Unprod$NUCLEOTIDE,topclones_sam2Unprod$NUCLEOTIDE),]
      commonClonesUnProd_V = sum(commonClonesUnProd_V_data$COUNT)/sum(topclones_sam2Unprod$COUNT)
      
      uv_UnProd=commonClonesUnProd_U * commonClonesUnProd_V
      uplusv_UnProd=commonClonesUnProd_U + commonClonesUnProd_V
      
      AbundanceBasedJIUnProd = formatC(uv_UnProd/(uplusv_UnProd-uv_UnProd),digits=8,format="f")
      
      overlapValuesUnprodJI_abundancebased=c(overlapValuesUnprodJI_abundancebased,as.numeric(AbundanceBasedJIUnProd))
     
      
 
      # Ninsertions  in common sequences
      
      commonClonotypes1= topclones_sam1Prod[topclones_sam1Prod$NUCLEOTIDE %in% intersect(topclones_sam1Prod$NUCLEOTIDE,topclones_sam2Prod$NUCLEOTIDE),]
      commonClonotypes2= topclones_sam2Prod[topclones_sam2Prod$NUCLEOTIDE %in% intersect(topclones_sam1Prod$NUCLEOTIDE,topclones_sam2Prod$NUCLEOTIDE),]
      
      commonClonotypes1Unprod= topclones_sam1Unprod[topclones_sam1Unprod$NUCLEOTIDE %in% intersect(topclones_sam1Unprod$NUCLEOTIDE,topclones_sam2Unprod$NUCLEOTIDE),]
      commonClonotypes2Unprod= topclones_sam2Unprod[topclones_sam2Unprod$NUCLEOTIDE %in% intersect(topclones_sam1Unprod$NUCLEOTIDE,topclones_sam2Unprod$NUCLEOTIDE),]
      
      
      NinsertionsP1= mean(commonClonotypes1$N2INSERTION + commonClonotypes1$N1INSERTION)
      NinsertionsP2 = mean(commonClonotypes2$N2INSERTION + commonClonotypes2$N1INSERTION)
      
      NinsertionsUP1= mean(commonClonotypes1Unprod$N2INSERTION + commonClonotypes1Unprod$N1INSERTION)
      NinsertionsUP2 = mean(commonClonotypes2Unprod$N2INSERTION + commonClonotypes2Unprod$N1INSERTION)
      
      
      NinsertionProd1=c(NinsertionProd1,NinsertionsP1)
      NinsertionProd2= c(NinsertionProd2,NinsertionsP2)
      
      NinsertionUnprod1=c(NinsertionUnprod1,NinsertionsUP1)
      NinsertionUnprod2=c(NinsertionUnprod2,NinsertionsUP2)
      
      
    }else{
      topclones_sam1Prod = getAbundantClonotypesbynumber(sam1Prod,cutoff)
      topclones_sam1Unprod = getAbundantClonotypesbynumber(sam1Unprod,cutoff)
      
      topclones_sam2Prod = getAbundantClonotypesbynumber(sam2Prod,cutoff)
      topclones_sam2Unprod = getAbundantClonotypesbynumber(sam2Unprod,cutoff)
      
      
      # commons
      # commons
      commonClonesProd=length(intersect(topclones_sam1Prod$NUCLEOTIDE,topclones_sam2Prod$NUCLEOTIDE))
      unionProd=length(union(topclones_sam1Prod$NUCLEOTIDE,topclones_sam2Prod$NUCLEOTIDE))
      
      # commonality checked at the Nt level for non productive sequences, no AA given in the original data
      commonClonesUnprod=length(intersect(topclones_sam1Unprod$NUCLEOTIDE,topclones_sam2Unprod$NUCLEOTIDE))
      unionUnprod=length(union(topclones_sam1Unprod$NUCLEOTIDE,topclones_sam2Unprod$NUCLEOTIDE))
      
      
      overlapValuesProd=c(overlapValuesProd,as.numeric(commonClonesProd))
      overlapValuesUnprod=c(overlapValuesUnprod,as.numeric(commonClonesUnprod))
      
      overlapValuesProdJI=c(overlapValuesProdJI,formatC(commonClonesProd/unionProd,digits= 8,format="f"))
      overlapValuesUnprodJI = c(overlapValuesUnprodJI,formatC(commonClonesUnprod/unionUnprod,digits= 8,format="f"))
      
      
      # Jaccard abundance-based similarity index (UV/(U+V-UV)). U is total relative abundance in of shared sample 1, U total relative of shared in sample 2
      
      commonClonesProd_U_data= topclones_sam1Prod[topclones_sam1Prod$NUCLEOTIDE %in% intersect(topclones_sam1Prod$NUCLEOTIDE,topclones_sam2Prod$NUCLEOTIDE),]
      commonClonesProd_U = sum(commonClonesProd_U_data$COUNT)/sum(topclones_sam1Prod$COUNT)
      
      commonClonesProd_V_data= topclones_sam2Prod[topclones_sam2Prod$NUCLEOTIDE %in% intersect(topclones_sam1Prod$NUCLEOTIDE,topclones_sam2Prod$NUCLEOTIDE),]
      commonClonesProd_V = sum(commonClonesProd_V_data$COUNT)/sum(topclones_sam2Prod$COUNT)
      
      uv=commonClonesProd_U*commonClonesProd_V
      uplusv=commonClonesProd_U+commonClonesProd_V
      
      AbundanceBasedJIProd = format(uv/(uplusv-uv),digits=8,format="f")
      
      overlapValuesProdJI_abundancebased=c(overlapValuesProdJI_abundancebased,as.numeric(AbundanceBasedJIProd))
      
      # for unproductive sequences
      
      commonClonesUnProd_U_data= topclones_sam1Unprod[topclones_sam1Unprod$NUCLEOTIDE %in% intersect(topclones_sam1Unprod$NUCLEOTIDE,topclones_sam2Unprod$NUCLEOTIDE),]
      commonClonesUnProd_U = sum(commonClonesUnProd_U_data$COUNT)/sum(topclones_sam1Unprod$COUNT)
      
      commonClonesUnProd_V_data= topclones_sam2Unprod[topclones_sam2Unprod$NUCLEOTIDE %in% intersect(topclones_sam1Unprod$NUCLEOTIDE,topclones_sam2Unprod$NUCLEOTIDE),]
      commonClonesUnProd_V = sum(commonClonesUnProd_V_data$COUNT)/sum(topclones_sam2Unprod$COUNT)
      
      uv_UnProd=commonClonesUnProd_U * commonClonesUnProd_V
      uplusv_UnProd=commonClonesUnProd_U + commonClonesUnProd_V
      
      AbundanceBasedJIUnProd = formatC(uv_UnProd/(uplusv_UnProd-uv_UnProd),digits=8,format="f")
      
      overlapValuesUnprodJI_abundancebased=c(overlapValuesUnprodJI_abundancebased,as.numeric(AbundanceBasedJIUnProd))
      
      
      
      # Ninsertions  in common sequences
      
      commonClonotypes1= topclones_sam1Prod[topclones_sam1Prod$NUCLEOTIDE %in% intersect(topclones_sam1Prod$NUCLEOTIDE,topclones_sam2Prod$NUCLEOTIDE),]
      commonClonotypes2= topclones_sam2Prod[topclones_sam2Prod$NUCLEOTIDE %in% intersect(topclones_sam1Prod$NUCLEOTIDE,topclones_sam2Prod$NUCLEOTIDE),]
      
      commonClonotypes1Unprod= topclones_sam1Unprod[topclones_sam1Unprod$NUCLEOTIDE %in% intersect(topclones_sam1Unprod$NUCLEOTIDE,topclones_sam2Unprod$NUCLEOTIDE),]
      commonClonotypes2Unprod= topclones_sam2Unprod[topclones_sam2Unprod$NUCLEOTIDE %in% intersect(topclones_sam1Unprod$NUCLEOTIDE,topclones_sam2Unprod$NUCLEOTIDE),]
      
      
      NinsertionsP1= mean(commonClonotypes1$N2INSERTION + commonClonotypes1$N1INSERTION)
      NinsertionsP2 = mean(commonClonotypes2$N2INSERTION + commonClonotypes2$N1INSERTION)
      
      NinsertionsUP1= mean(commonClonotypes1Unprod$N2INSERTION + commonClonotypes1Unprod$N1INSERTION)
      NinsertionsUP2 = mean(commonClonotypes2Unprod$N2INSERTION + commonClonotypes2Unprod$N1INSERTION)
      
      
      NinsertionProd1=c(NinsertionProd1,NinsertionsP1)
      NinsertionProd2= c(NinsertionProd2,NinsertionsP2)
      
      NinsertionUnprod1=c(NinsertionUnprod1,NinsertionsUP1)
      NinsertionUnprod2=c(NinsertionUnprod2,NinsertionsUP2)
        
      }
    
    #totalNumberOfProdClonotypes= nrow(topclones_sam1Prod) + nrow(topclones_sam2Prod)
    #totalNumberOfUnprodClonotypes= nrow(topclones_sam1Unprod) + nrow(topclones_sam2Unprod)

    
    if(nrow(toadd) > 0){
    newrow=data.frame(overlapValuesProd,overlapValuesUnprod,overlapValuesProdJI,overlapValuesUnprodJI,overlapValuesProdJI_abundancebased,overlapValuesUnprodJI_abundancebased,NinsertionProd1,NinsertionProd2,NinsertionUnprod1,NinsertionUnprod2)
    colnames(newrow) <- c("OverlapCountProductive","OverlapCountUnProductive","OverlapProd_JI","OverlapUnProd_JI","AbundanceBased_JI","AbundanceBased_JI","NinsertionProdSample1","NinsertionProdSample2","NinsertionUnProdSample1","NinsertionUnProdSample2")
    
    toadd=rbind(toadd,newrow)
    }else{
      toadd=data.frame(overlapValuesProd,overlapValuesUnprod,overlapValuesProdJI,overlapValuesUnprodJI,overlapValuesProdJI_abundancebased,overlapValuesUnprodJI_abundancebased,NinsertionProd1,NinsertionProd2,NinsertionUnprod1,NinsertionUnprod2)
      colnames(toadd) <- c("OverlapCountProductive","OverlapCountUnProductive","OverlapProd_JI","OverlapUnProd_JI","AbundanceBased_JI","AbundanceBased_JI","NinsertionProdSample1","NinsertionProdSample2","NinsertionUnProdSample1","NinsertionUnProdSample2")  
    }
    
  }
   
  

    
  print(nrow(toadd))
  resultTable=data.frame(m,toadd)        
  
  write.table(resultTable,file="Pairwise Overlap Counts.txt",row.names = T,sep = "\t")
  return(resultTable)
}

overLappingClonotypes2<- function(reps,cutoff=20,cutoffby="perc"){
  ### Returns a datamatrix of top topnumber clonotopes from each sample plus: overlap count or Jaccard index, for prod and non productive sequences, N insertion size
  
  
  # setting work parameters
  cutoffby = toupper(cutoffby)
  
  
  m <- as.data.frame(expand.grid(reps,reps))
  colnames(m)  <- c("Sample1","Sample2")
  
  toadd= data.frame()
  
  
  
  for(r in 1:nrow(m)){
    
    overlapValuesProd=c()
    overlapValuesUnprod=c()
    overlapValuesProdJI=c()
    overlapValuesUnprodJI=c()
    
    overlapValuesProdJI_abundancebased=c()
    overlapValuesUnprodJI_abundancebased=c()
    
    NinsertionProd1=c()
    NinsertionUnprod1=c()
    
    NinsertionProd2=c()
    NinsertionUnprod2=c()
    
    
    print(r)
    sam1=m[r,1]
    sam2=m[r,2]
    
    #print(sam1)
    #print(sam2)
    
    sam1 <- get(as.character(sam1))
    sam2 <- get(as.character(sam2))
    
    
    sam1Prod <- getProductive(sam1)
    sam1Unprod <- getUnproductive(sam1)
    
    sam2Prod <- getProductive(sam2)
    sam2Unprod <- getUnproductive(sam2)
    
    
    if(cutoffby=="PERC"){
      topclones_sam1Prod = getAbundantClonotypesbyperc(sam1Prod,cutoff)
      topclones_sam1Unprod = getAbundantClonotypesbyperc(sam1Unprod,cutoff)
      
      topclones_sam2Prod = getAbundantClonotypesbyperc(sam2Prod,cutoff)
      topclones_sam2Unprod = getAbundantClonotypesbyperc(sam2Unprod,cutoff)
      
      # commons
      commonClonesProd=length(intersect(topclones_sam1Prod$NUCLEOTIDE,topclones_sam2Prod$NUCLEOTIDE)) 
      unionProd=length(union(topclones_sam1Prod$NUCLEOTIDE,topclones_sam2Prod$NUCLEOTIDE)) 
      
      # commonality checked at the Nt level for non productive sequences, no AA given in the original data
      commonClonesUnprod=length(intersect(topclones_sam1Unprod$NUCLEOTIDE,topclones_sam2Unprod$NUCLEOTIDE))
      unionUnprod=length(union(topclones_sam1Unprod$NUCLEOTIDE,topclones_sam2Unprod$NUCLEOTIDE))
      
      if(unionUnprod==0){
        unionUnprod = 1
      }
      
      if(unionProd==0){
        unionProd = 1
      }
      
      overlapValuesProd=c(overlapValuesProd,as.numeric(commonClonesProd))
      overlapValuesUnprod=c(overlapValuesUnprod,as.numeric(commonClonesUnprod))
      
      overlapValuesProdJI=c(overlapValuesProdJI,commonClonesProd/unionProd)
      overlapValuesUnprodJI = c(overlapValuesUnprodJI,commonClonesUnprod/unionUnprod)
      
      # Jaccard abundance-based similarity index (UV/(U+V-UV)). U is total relative abundance in of shared sample 1, U total relative of shared in sample 2
      
      commonClonesProd_U_data= topclones_sam1Prod[topclones_sam1Prod$NUCLEOTIDE %in% intersect(topclones_sam1Prod$NUCLEOTIDE,topclones_sam2Prod$NUCLEOTIDE),]
      commonClonesProd_U = sum(commonClonesProd_U_data$COUNT)/sum(topclones_sam1Prod$COUNT)
      
      commonClonesProd_V_data= topclones_sam2Prod[topclones_sam2Prod$NUCLEOTIDE %in% intersect(topclones_sam1Prod$NUCLEOTIDE,topclones_sam2Prod$NUCLEOTIDE),]
      commonClonesProd_V = sum(commonClonesProd_V_data$COUNT)/sum(topclones_sam2Prod$COUNT)
      
      uv=commonClonesProd_U*commonClonesProd_V
      uplusv=commonClonesProd_U+commonClonesProd_V
      
      AbundanceBasedJIProd = uv/(uplusv-uv)
      
      overlapValuesProdJI_abundancebased=c(overlapValuesProdJI_abundancebased,as.numeric(AbundanceBasedJIProd))
      
      # for unproductive sequences
      
      commonClonesUnProd_U_data= topclones_sam1Unprod[topclones_sam1Unprod$NUCLEOTIDE %in% intersect(topclones_sam1Unprod$NUCLEOTIDE,topclones_sam2Unprod$NUCLEOTIDE),]
      commonClonesUnProd_U = sum(commonClonesUnProd_U_data$COUNT)/sum(topclones_sam1Unprod$COUNT)
      
      commonClonesUnProd_V_data= topclones_sam2Unprod[topclones_sam2Unprod$NUCLEOTIDE %in% intersect(topclones_sam1Unprod$NUCLEOTIDE,topclones_sam2Unprod$NUCLEOTIDE),]
      commonClonesUnProd_V = sum(commonClonesUnProd_V_data$COUNT)/sum(topclones_sam2Unprod$COUNT)
      
      uv_UnProd=commonClonesUnProd_U * commonClonesUnProd_V
      uplusv_UnProd=commonClonesUnProd_U + commonClonesUnProd_V
      
      AbundanceBasedJIUnProd = uv_UnProd/(uplusv_UnProd-uv_UnProd)
      
      overlapValuesUnprodJI_abundancebased=c(overlapValuesUnprodJI_abundancebased,as.numeric(AbundanceBasedJIUnProd))
      
      
      
      # Ninsertions  in common sequences
      
      commonClonotypes1= topclones_sam1Prod[topclones_sam1Prod$NUCLEOTIDE %in% intersect(topclones_sam1Prod$NUCLEOTIDE,topclones_sam2Prod$NUCLEOTIDE),]
      commonClonotypes2= topclones_sam2Prod[topclones_sam2Prod$NUCLEOTIDE %in% intersect(topclones_sam1Prod$NUCLEOTIDE,topclones_sam2Prod$NUCLEOTIDE),]
      
      commonClonotypes1Unprod= topclones_sam1Unprod[topclones_sam1Unprod$NUCLEOTIDE %in% intersect(topclones_sam1Unprod$NUCLEOTIDE,topclones_sam2Unprod$NUCLEOTIDE),]
      commonClonotypes2Unprod= topclones_sam2Unprod[topclones_sam2Unprod$NUCLEOTIDE %in% intersect(topclones_sam1Unprod$NUCLEOTIDE,topclones_sam2Unprod$NUCLEOTIDE),]
      
      
      NinsertionsP1= mean(commonClonotypes1$N2INSERTION + commonClonotypes1$N1INSERTION)
      NinsertionsP2 = mean(commonClonotypes2$N2INSERTION + commonClonotypes2$N1INSERTION)
      
      NinsertionsUP1= mean(commonClonotypes1Unprod$N2INSERTION + commonClonotypes1Unprod$N1INSERTION)
      NinsertionsUP2 = mean(commonClonotypes2Unprod$N2INSERTION + commonClonotypes2Unprod$N1INSERTION)
      
      
      NinsertionProd1=c(NinsertionProd1,NinsertionsP1)
      NinsertionProd2= c(NinsertionProd2,NinsertionsP2)
      
      NinsertionUnprod1=c(NinsertionUnprod1,NinsertionsUP1)
      NinsertionUnprod2=c(NinsertionUnprod2,NinsertionsUP2)
      
      
    }else{
      topclones_sam1Prod = getAbundantClonotypesbynumber(sam1Prod,cutoff)
      topclones_sam1Unprod = getAbundantClonotypesbynumber(sam1Unprod,cutoff)
      
      topclones_sam2Prod = getAbundantClonotypesbynumber(sam2Prod,cutoff)
      topclones_sam2Unprod = getAbundantClonotypesbynumber(sam2Unprod,cutoff)
      
      
      # commons
      # commons
      commonClonesProd=length(intersect(topclones_sam1Prod$NUCLEOTIDE,topclones_sam2Prod$NUCLEOTIDE))
      unionProd=length(union(topclones_sam1Prod$NUCLEOTIDE,topclones_sam2Prod$NUCLEOTIDE))
      
      # commonality checked at the Nt level for non productive sequences, no AA given in the original data
      commonClonesUnprod=length(intersect(topclones_sam1Unprod$NUCLEOTIDE,topclones_sam2Unprod$NUCLEOTIDE))
      unionUnprod=length(union(topclones_sam1Unprod$NUCLEOTIDE,topclones_sam2Unprod$NUCLEOTIDE))
      
      if(unionUnprod==0){
        unionUnprod = 1
      }
      
      if(unionProd==0){
        unionProd = 1
      }
      
      overlapValuesProd=c(overlapValuesProd,as.numeric(commonClonesProd))
      overlapValuesUnprod=c(overlapValuesUnprod,as.numeric(commonClonesUnprod))
      
      overlapValuesProdJI=c(overlapValuesProdJI,commonClonesProd/unionProd)
      overlapValuesUnprodJI = c(overlapValuesUnprodJI,commonClonesUnprod/unionUnprod)
      
      
      # Jaccard abundance-based similarity index (UV/(U+V-UV)). U is total relative abundance in of shared sample 1, U total relative of shared in sample 2
      
      commonClonesProd_U_data= topclones_sam1Prod[topclones_sam1Prod$NUCLEOTIDE %in% intersect(topclones_sam1Prod$NUCLEOTIDE,topclones_sam2Prod$NUCLEOTIDE),]
      commonClonesProd_U = sum(commonClonesProd_U_data$COUNT)/sum(topclones_sam1Prod$COUNT)
      
      commonClonesProd_V_data= topclones_sam2Prod[topclones_sam2Prod$NUCLEOTIDE %in% intersect(topclones_sam1Prod$NUCLEOTIDE,topclones_sam2Prod$NUCLEOTIDE),]
      commonClonesProd_V = sum(commonClonesProd_V_data$COUNT)/sum(topclones_sam2Prod$COUNT)
      
      uv=commonClonesProd_U*commonClonesProd_V
      uplusv=commonClonesProd_U+commonClonesProd_V
      
      AbundanceBasedJIProd = uv/(uplusv-uv)
      
      overlapValuesProdJI_abundancebased=c(overlapValuesProdJI_abundancebased,as.numeric(AbundanceBasedJIProd))
      
      # for unproductive sequences
      
      commonClonesUnProd_U_data= topclones_sam1Unprod[topclones_sam1Unprod$NUCLEOTIDE %in% intersect(topclones_sam1Unprod$NUCLEOTIDE,topclones_sam2Unprod$NUCLEOTIDE),]
      commonClonesUnProd_U = sum(commonClonesUnProd_U_data$COUNT)/sum(topclones_sam1Unprod$COUNT)
      
      commonClonesUnProd_V_data= topclones_sam2Unprod[topclones_sam2Unprod$NUCLEOTIDE %in% intersect(topclones_sam1Unprod$NUCLEOTIDE,topclones_sam2Unprod$NUCLEOTIDE),]
      commonClonesUnProd_V = sum(commonClonesUnProd_V_data$COUNT)/sum(topclones_sam2Unprod$COUNT)
      
      uv_UnProd=commonClonesUnProd_U * commonClonesUnProd_V
      uplusv_UnProd=commonClonesUnProd_U + commonClonesUnProd_V
      
      AbundanceBasedJIUnProd = uv_UnProd/(uplusv_UnProd-uv_UnProd)
      
      overlapValuesUnprodJI_abundancebased=c(overlapValuesUnprodJI_abundancebased,as.numeric(AbundanceBasedJIUnProd))
      
      
      
      # Ninsertions  in common sequences
      
      commonClonotypes1= topclones_sam1Prod[topclones_sam1Prod$NUCLEOTIDE %in% intersect(topclones_sam1Prod$NUCLEOTIDE,topclones_sam2Prod$NUCLEOTIDE),]
      commonClonotypes2= topclones_sam2Prod[topclones_sam2Prod$NUCLEOTIDE %in% intersect(topclones_sam1Prod$NUCLEOTIDE,topclones_sam2Prod$NUCLEOTIDE),]
      
      commonClonotypes1Unprod= topclones_sam1Unprod[topclones_sam1Unprod$NUCLEOTIDE %in% intersect(topclones_sam1Unprod$NUCLEOTIDE,topclones_sam2Unprod$NUCLEOTIDE),]
      commonClonotypes2Unprod= topclones_sam2Unprod[topclones_sam2Unprod$NUCLEOTIDE %in% intersect(topclones_sam1Unprod$NUCLEOTIDE,topclones_sam2Unprod$NUCLEOTIDE),]
      
      
      NinsertionsP1= mean(commonClonotypes1$N2INSERTION + commonClonotypes1$N1INSERTION)
      NinsertionsP2 = mean(commonClonotypes2$N2INSERTION + commonClonotypes2$N1INSERTION)
      
      NinsertionsUP1= mean(commonClonotypes1Unprod$N2INSERTION + commonClonotypes1Unprod$N1INSERTION)
      NinsertionsUP2 = mean(commonClonotypes2Unprod$N2INSERTION + commonClonotypes2Unprod$N1INSERTION)
      
      
      NinsertionProd1=c(NinsertionProd1,NinsertionsP1)
      NinsertionProd2= c(NinsertionProd2,NinsertionsP2)
      
      NinsertionUnprod1=c(NinsertionUnprod1,NinsertionsUP1)
      NinsertionUnprod2=c(NinsertionUnprod2,NinsertionsUP2)
      
    }
    
    #totalNumberOfProdClonotypes= nrow(topclones_sam1Prod) + nrow(topclones_sam2Prod)
    #totalNumberOfUnprodClonotypes= nrow(topclones_sam1Unprod) + nrow(topclones_sam2Unprod)
    
    
    if(nrow(toadd) > 0){
      newrow=data.frame(overlapValuesProd,overlapValuesUnprod,overlapValuesProdJI,overlapValuesUnprodJI,overlapValuesProdJI_abundancebased,overlapValuesUnprodJI_abundancebased,NinsertionProd1,NinsertionProd2,NinsertionUnprod1,NinsertionUnprod2)
      colnames(newrow) <- c("OverlapCountProductive","OverlapCountUnProductive","OverlapProd_JI","OverlapUnProd_JI","AbundanceBased_JI","AbundanceBased_JI","NinsertionProdSample1","NinsertionProdSample2","NinsertionUnProdSample1","NinsertionUnProdSample2")
      
      toadd=rbind(toadd,newrow)
    }else{
      toadd=data.frame(overlapValuesProd,overlapValuesUnprod,overlapValuesProdJI,overlapValuesUnprodJI,overlapValuesProdJI_abundancebased,overlapValuesUnprodJI_abundancebased,NinsertionProd1,NinsertionProd2,NinsertionUnprod1,NinsertionUnprod2)
      colnames(toadd) <- c("OverlapCountProductive","OverlapCountUnProductive","OverlapProd_JI","OverlapUnProd_JI","AbundanceBased_JI","AbundanceBased_JI","NinsertionProdSample1","NinsertionProdSample2","NinsertionUnProdSample1","NinsertionUnProdSample2")  
    }
    
  }
  
  
  
  
  print(nrow(toadd))
  resultTable=data.frame(m,toadd)        
  
  # write pdf dendrogram and heatmaps
  
  # select the sample columns 1 and 2, and the JI column 5
  forProdselected = resultTable[,c(1,2,5)]
  forProdselected = forProdselected[order(forProdselected$Sample1),]
  
  l = length(unique(forProdselected$Sample1))
  mat = matrix(as.numeric(forProdselected[,3]),nrow=l,ncol=l,byrow=T)
  # replace NaN values with 0
  mat[is.nan(mat)] = 0
  colnames(mat) <- unique(forProdselected$Sample1)
  rownames(mat) <- unique(forProdselected$Sample1)
  
  
  mat = signif(mat,digits=2)
  hc = hclust(as.dist(log10(1-mat)))
  
  flname = paste("ProdJI_dendrogram.pdf")
  pdf(file=flname, useDingbats = FALSE)
  
  # very simple dendrogram
  plot(hc)
  
  dev.off()
  
  # for unprod
  # select the sample columns 1 and 2, and the JI column 5
  forUnprodselected = resultTable[,c(1,2,6)]
  forUnprodselected = forUnprodselected[order(forUnprodselected$Sample1),]
  
  l = length(unique(forUnprodselected$Sample1))
  mat = matrix(as.numeric(forUnprodselected[,3]),nrow=l,ncol=l,byrow=T)
  # replace NaN values with 0
  mat[is.nan(mat)] = 0
  colnames(mat) <- unique(forUnprodselected$Sample1)
  rownames(mat) <- unique(forUnprodselected$Sample1)
  
  
  mat = signif(mat,digits=2)
  hc = hclust(as.dist(log10(1-mat)))
  
  
  flname = paste("UnprodJI_dendrogram.pdf")
  pdf(file=flname, useDingbats = FALSE)
  
  # very simple dendrogram
  plot(hc)
  
  dev.off()
  
  write.table(resultTable,file="Pairwise Overlap Counts.txt",row.names = T,sep = "\t")
  return(resultTable)
}


# JI calculation for 'has stop' unproductive Nucleotides

overLappingClonotypesHasStop <- function(reps,cutoff=20,cutoffby="perc"){
  ### Returns a datamatrix of top topnumber clonotopes from each sample plus: overlap count or Jaccard index, for prod and non productive sequences, N insertion size
  
  ### Total read counts for the common clones of both samples are given for the productive/hasStop versions the JI analysis. This is implemented when the analysis is 
  ### done for by "perc". 
  
  
  # setting work parameters
  cutoffby = toupper(cutoffby)
  
  
  m <- as.data.frame(expand.grid(reps,reps))
  colnames(m)  <- c("Sample1","Sample2")
  
  toadd= data.frame()
  
  
  
  for(r in 1:nrow(m)){
    
    overlapValuesProd=c()
    overlapValuesUnprod=c()
    
    overlapValuesTotalSam1Prod=c()
    overlapValuesTotalSam2Prod=c()
    
    overlapValuesTotalSam1Unprod=c()
    overlapValuesTotalSam2Unprod=c()
    
    overlapValuesProdJI=c()
    overlapValuesUnprodJI=c()
    
    overlapValuesProdJI_abundancebased=c()
    overlapValuesUnprodJI_abundancebased=c()
    
    NinsertionProd1=c()
    NinsertionUnprod1=c()
    
    NinsertionProd2=c()
    NinsertionUnprod2=c()
    
    
    print(r)
    sam1=m[r,1]
    sam2=m[r,2]
    
    #print(sam1)
    #print(sam2)
    
    sam1 <- get(as.character(sam1))
    sam2 <- get(as.character(sam2))
    
    
    sam1Prod <- getProductive(sam1)
    sam1Unprod <- getHasStop(sam1)
    
    sam2Prod <- getProductive(sam2)
    sam2Unprod <- getHasStop(sam2)
    
    
    if(cutoffby=="PERC"){
      topclones_sam1Prod = getAbundantClonotypesbyperc(sam1Prod,cutoff)
      topclones_sam1Unprod = getAbundantClonotypesbyperc(sam1Unprod,cutoff)
      
      topclones_sam2Prod = getAbundantClonotypesbyperc(sam2Prod,cutoff)
      topclones_sam2Unprod = getAbundantClonotypesbyperc(sam2Unprod,cutoff)
      
      # commons
      commonClonesProd=length(intersect(topclones_sam1Prod$NUCLEOTIDE,topclones_sam2Prod$NUCLEOTIDE)) 
      unionProd=length(union(topclones_sam1Prod$NUCLEOTIDE,topclones_sam2Prod$NUCLEOTIDE)) 
      
      # commonality checked at the Nt level for non productive sequences, no AA given in the original data
      commonClonesUnprod=length(intersect(topclones_sam1Unprod$NUCLEOTIDE,topclones_sam2Unprod$NUCLEOTIDE))
      unionUnprod=length(union(topclones_sam1Unprod$NUCLEOTIDE,topclones_sam2Unprod$NUCLEOTIDE))
      
      
      
      
      if(unionUnprod==0){
        unionUnprod = 1
      }
      
      if(unionProd==0){
        unionProd = 1
      }
      
      overlapValuesProd=c(overlapValuesProd,as.numeric(commonClonesProd))
      overlapValuesUnprod=c(overlapValuesUnprod,as.numeric(commonClonesUnprod))
      
      
      
      overlapValuesProdJI=c(overlapValuesProdJI,commonClonesProd/unionProd)
      overlapValuesUnprodJI = c(overlapValuesUnprodJI,commonClonesUnprod/unionUnprod)
      
      # Jaccard abundance-based similarity index (UV/(U+V-UV)). U is total relative abundance in of shared sample 1, U total relative of shared in sample 2
      
      commonClonesProd_U_data= topclones_sam1Prod[topclones_sam1Prod$NUCLEOTIDE %in% intersect(topclones_sam1Prod$NUCLEOTIDE,topclones_sam2Prod$NUCLEOTIDE),]
      commonClonesProd_U = sum(commonClonesProd_U_data$COUNT)/sum(topclones_sam1Prod$COUNT)
      
      commonClonesProd_U_abundance = sum(commonClonesProd_U_data$COUNT)
      
      commonClonesProd_V_data= topclones_sam2Prod[topclones_sam2Prod$NUCLEOTIDE %in% intersect(topclones_sam1Prod$NUCLEOTIDE,topclones_sam2Prod$NUCLEOTIDE),]
      commonClonesProd_V = sum(commonClonesProd_V_data$COUNT)/sum(topclones_sam2Prod$COUNT)
      
      commonClonesProd_V_abundance = sum(commonClonesProd_V_data$COUNT)
      
      overlapValuesTotalSam1Prod=c(overlapValuesTotalSam1Prod,as.numeric(commonClonesProd_U_abundance))
      overlapValuesTotalSam2Prod=c(overlapValuesTotalSam2Prod,as.numeric(commonClonesProd_V_abundance))
      
      uv=commonClonesProd_U*commonClonesProd_V
      uplusv=commonClonesProd_U+commonClonesProd_V
      
      AbundanceBasedJIProd = uv/(uplusv-uv)
      
      overlapValuesProdJI_abundancebased=c(overlapValuesProdJI_abundancebased,as.numeric(AbundanceBasedJIProd))
      
      # for unproductive sequences
      
      commonClonesUnProd_U_data= topclones_sam1Unprod[topclones_sam1Unprod$NUCLEOTIDE %in% intersect(topclones_sam1Unprod$NUCLEOTIDE,topclones_sam2Unprod$NUCLEOTIDE),]
      commonClonesUnProd_U = sum(commonClonesUnProd_U_data$COUNT)/sum(topclones_sam1Unprod$COUNT)
      
      commonClonesUnProd_U_abundance = sum(commonClonesUnProd_U_data$COUNT)
      
      commonClonesUnProd_V_data= topclones_sam2Unprod[topclones_sam2Unprod$NUCLEOTIDE %in% intersect(topclones_sam1Unprod$NUCLEOTIDE,topclones_sam2Unprod$NUCLEOTIDE),]
      commonClonesUnProd_V = sum(commonClonesUnProd_V_data$COUNT)/sum(topclones_sam2Unprod$COUNT)
      
      commonClonesUnProd_V_abundance = sum(commonClonesUnProd_V_data$COUNT)
      
      
      overlapValuesTotalSam1Unprod=c(overlapValuesTotalSam1Unprod,as.numeric(commonClonesUnProd_U_abundance))
      overlapValuesTotalSam2Unprod=c(overlapValuesTotalSam2Unprod,as.numeric(commonClonesUnProd_V_abundance))
      
      
      uv_UnProd=commonClonesUnProd_U * commonClonesUnProd_V
      uplusv_UnProd=commonClonesUnProd_U + commonClonesUnProd_V
      
      AbundanceBasedJIUnProd = uv_UnProd/(uplusv_UnProd-uv_UnProd)
      
      overlapValuesUnprodJI_abundancebased=c(overlapValuesUnprodJI_abundancebased,as.numeric(AbundanceBasedJIUnProd))
      
      
      
      # Ninsertions  in common sequences
      
      commonClonotypes1= topclones_sam1Prod[topclones_sam1Prod$NUCLEOTIDE %in% intersect(topclones_sam1Prod$NUCLEOTIDE,topclones_sam2Prod$NUCLEOTIDE),]
      commonClonotypes2= topclones_sam2Prod[topclones_sam2Prod$NUCLEOTIDE %in% intersect(topclones_sam1Prod$NUCLEOTIDE,topclones_sam2Prod$NUCLEOTIDE),]
      
      commonClonotypes1Unprod= topclones_sam1Unprod[topclones_sam1Unprod$NUCLEOTIDE %in% intersect(topclones_sam1Unprod$NUCLEOTIDE,topclones_sam2Unprod$NUCLEOTIDE),]
      commonClonotypes2Unprod= topclones_sam2Unprod[topclones_sam2Unprod$NUCLEOTIDE %in% intersect(topclones_sam1Unprod$NUCLEOTIDE,topclones_sam2Unprod$NUCLEOTIDE),]
      
      
      NinsertionsP1= mean(commonClonotypes1$N2INSERTION + commonClonotypes1$N1INSERTION)
      NinsertionsP2 = mean(commonClonotypes2$N2INSERTION + commonClonotypes2$N1INSERTION)
      
      NinsertionsUP1= mean(commonClonotypes1Unprod$N2INSERTION + commonClonotypes1Unprod$N1INSERTION)
      NinsertionsUP2 = mean(commonClonotypes2Unprod$N2INSERTION + commonClonotypes2Unprod$N1INSERTION)
      
      
      NinsertionProd1=c(NinsertionProd1,NinsertionsP1)
      NinsertionProd2= c(NinsertionProd2,NinsertionsP2)
      
      NinsertionUnprod1=c(NinsertionUnprod1,NinsertionsUP1)
      NinsertionUnprod2=c(NinsertionUnprod2,NinsertionsUP2)
      
      
    }else{
      topclones_sam1Prod = getAbundantClonotypesbynumber(sam1Prod,cutoff)
      topclones_sam1Unprod = getAbundantClonotypesbynumber(sam1Unprod,cutoff)
      
      topclones_sam2Prod = getAbundantClonotypesbynumber(sam2Prod,cutoff)
      topclones_sam2Unprod = getAbundantClonotypesbynumber(sam2Unprod,cutoff)
      
      
      # commons
      # commons
      commonClonesProd=length(intersect(topclones_sam1Prod$NUCLEOTIDE,topclones_sam2Prod$NUCLEOTIDE))
      unionProd=length(union(topclones_sam1Prod$NUCLEOTIDE,topclones_sam2Prod$NUCLEOTIDE))
      
      # commonality checked at the Nt level for non productive sequences, no AA given in the original data
      commonClonesUnprod=length(intersect(topclones_sam1Unprod$NUCLEOTIDE,topclones_sam2Unprod$NUCLEOTIDE))
      unionUnprod=length(union(topclones_sam1Unprod$NUCLEOTIDE,topclones_sam2Unprod$NUCLEOTIDE))
      
      if(unionUnprod==0){
        unionUnprod = 1
      }
      
      if(unionProd==0){
        unionProd = 1
      }
      
      overlapValuesProd=c(overlapValuesProd,as.numeric(commonClonesProd))
      overlapValuesUnprod=c(overlapValuesUnprod,as.numeric(commonClonesUnprod))
      
      overlapValuesProdJI=c(overlapValuesProdJI,commonClonesProd/unionProd)
      overlapValuesUnprodJI = c(overlapValuesUnprodJI,commonClonesUnprod/unionUnprod)
      
      
      # Jaccard abundance-based similarity index (UV/(U+V-UV)). U is total relative abundance in of shared sample 1, U total relative of shared in sample 2
      
      commonClonesProd_U_data= topclones_sam1Prod[topclones_sam1Prod$NUCLEOTIDE %in% intersect(topclones_sam1Prod$NUCLEOTIDE,topclones_sam2Prod$NUCLEOTIDE),]
      commonClonesProd_U = sum(commonClonesProd_U_data$COUNT)/sum(topclones_sam1Prod$COUNT)
      
      commonClonesProd_V_data= topclones_sam2Prod[topclones_sam2Prod$NUCLEOTIDE %in% intersect(topclones_sam1Prod$NUCLEOTIDE,topclones_sam2Prod$NUCLEOTIDE),]
      commonClonesProd_V = sum(commonClonesProd_V_data$COUNT)/sum(topclones_sam2Prod$COUNT)
      
      uv=commonClonesProd_U*commonClonesProd_V
      uplusv=commonClonesProd_U+commonClonesProd_V
      
      AbundanceBasedJIProd = uv/(uplusv-uv)
      
      overlapValuesProdJI_abundancebased=c(overlapValuesProdJI_abundancebased,as.numeric(AbundanceBasedJIProd))
      
      # for unproductive sequences
      
      commonClonesUnProd_U_data= topclones_sam1Unprod[topclones_sam1Unprod$NUCLEOTIDE %in% intersect(topclones_sam1Unprod$NUCLEOTIDE,topclones_sam2Unprod$NUCLEOTIDE),]
      commonClonesUnProd_U = sum(commonClonesUnProd_U_data$COUNT)/sum(topclones_sam1Unprod$COUNT)
      
      commonClonesUnProd_V_data= topclones_sam2Unprod[topclones_sam2Unprod$NUCLEOTIDE %in% intersect(topclones_sam1Unprod$NUCLEOTIDE,topclones_sam2Unprod$NUCLEOTIDE),]
      commonClonesUnProd_V = sum(commonClonesUnProd_V_data$COUNT)/sum(topclones_sam2Unprod$COUNT)
      
      uv_UnProd=commonClonesUnProd_U * commonClonesUnProd_V
      uplusv_UnProd=commonClonesUnProd_U + commonClonesUnProd_V
      
      AbundanceBasedJIUnProd = uv_UnProd/(uplusv_UnProd-uv_UnProd)
      
      overlapValuesUnprodJI_abundancebased=c(overlapValuesUnprodJI_abundancebased,as.numeric(AbundanceBasedJIUnProd))
      
      
      
      # Ninsertions  in common sequences
      
      commonClonotypes1= topclones_sam1Prod[topclones_sam1Prod$NUCLEOTIDE %in% intersect(topclones_sam1Prod$NUCLEOTIDE,topclones_sam2Prod$NUCLEOTIDE),]
      commonClonotypes2= topclones_sam2Prod[topclones_sam2Prod$NUCLEOTIDE %in% intersect(topclones_sam1Prod$NUCLEOTIDE,topclones_sam2Prod$NUCLEOTIDE),]
      
      commonClonotypes1Unprod= topclones_sam1Unprod[topclones_sam1Unprod$NUCLEOTIDE %in% intersect(topclones_sam1Unprod$NUCLEOTIDE,topclones_sam2Unprod$NUCLEOTIDE),]
      commonClonotypes2Unprod= topclones_sam2Unprod[topclones_sam2Unprod$NUCLEOTIDE %in% intersect(topclones_sam1Unprod$NUCLEOTIDE,topclones_sam2Unprod$NUCLEOTIDE),]
      
      
      NinsertionsP1= mean(commonClonotypes1$N2INSERTION + commonClonotypes1$N1INSERTION)
      NinsertionsP2 = mean(commonClonotypes2$N2INSERTION + commonClonotypes2$N1INSERTION)
      
      NinsertionsUP1= mean(commonClonotypes1Unprod$N2INSERTION + commonClonotypes1Unprod$N1INSERTION)
      NinsertionsUP2 = mean(commonClonotypes2Unprod$N2INSERTION + commonClonotypes2Unprod$N1INSERTION)
      
      
      NinsertionProd1=c(NinsertionProd1,NinsertionsP1)
      NinsertionProd2= c(NinsertionProd2,NinsertionsP2)
      
      NinsertionUnprod1=c(NinsertionUnprod1,NinsertionsUP1)
      NinsertionUnprod2=c(NinsertionUnprod2,NinsertionsUP2)
      
    }
    
    #totalNumberOfProdClonotypes= nrow(topclones_sam1Prod) + nrow(topclones_sam2Prod)
    #totalNumberOfUnprodClonotypes= nrow(topclones_sam1Unprod) + nrow(topclones_sam2Unprod)
    
    
    if(nrow(toadd) > 0){
      newrow=data.frame(overlapValuesProd,overlapValuesTotalSam1Prod,overlapValuesTotalSam2Prod,overlapValuesUnprod,overlapValuesTotalSam1Unprod,overlapValuesTotalSam2Unprod,overlapValuesProdJI,overlapValuesUnprodJI,overlapValuesProdJI_abundancebased,overlapValuesUnprodJI_abundancebased,NinsertionProd1,NinsertionProd2,NinsertionUnprod1,NinsertionUnprod2)
      colnames(newrow) <- c("OverlapCountProductive","OverlapCountTotalProductiveSam1","OverlapCountTotalProductiveSam2","OverlapCountHasStop","OverlapCountTotalHasStopSam1","OverlapCountTotalHasStopSam2","OverlapProd_JI","OverlapHasStop_JI","AbundanceBasedProd_JI","AbundanceBasedHasStop_JI","NinsertionProdSample1","NinsertionProdSample2","NinsertionHasStopSample1","NinsertionHasStopSample2")
      
      toadd=rbind(toadd,newrow)
    }else{
      toadd=data.frame(overlapValuesProd,overlapValuesTotalSam1Prod,overlapValuesTotalSam2Prod,overlapValuesUnprod,overlapValuesTotalSam1Unprod,overlapValuesTotalSam2Unprod,overlapValuesProdJI,overlapValuesUnprodJI,overlapValuesProdJI_abundancebased,overlapValuesUnprodJI_abundancebased,NinsertionProd1,NinsertionProd2,NinsertionUnprod1,NinsertionUnprod2)
      colnames(toadd) <- c("OverlapCountProductive","OverlapCountTotalProductiveSam1","OverlapCountTotalProductiveSam2","OverlapCountHasStop","OverlapCountTotalHasStopSam1","OverlapCountTotalHasStopSam2","OverlapProd_JI","OverlapHasStop_JI","AbundanceBasedProd_JI","AbundanceBasedHasStop_JI","NinsertionProdSample1","NinsertionProdSample2","NinsertionHasStopSample1","NinsertionHasStopSample2")
    }
    
  }
  
  
  
  
  print(nrow(toadd))
  resultTable=data.frame(m,toadd)        
  
  
  write.table(resultTable,file="Pairwise Overlap Counts.txt",row.names = T,sep = "\t")
  return(resultTable)
  
  
}


# Amico acid based Jaccard index calculation

overLappingClonotypesAA<- function(reps,cutoff=20,cutoffby="perc"){
  ### Returns a datamatrix of top topnumber clonotopes (Amino acid level) from each sample plus: overlap count or Jaccard index, for prod and non productive sequences, N insertion size
  
  
  # setting work parameters
  cutoffby = toupper(cutoffby)
  
  
  m <- as.data.frame(expand.grid(reps,reps))
  colnames(m)  <- c("Sample1","Sample2")
  
  toadd= data.frame()
  
  
  
  for(r in 1:nrow(m)){
    
    overlapValuesProd=c()
    overlapValuesUnprod=c()
    overlapValuesProdJI=c()
    overlapValuesUnprodJI=c()
    
    overlapValuesProdJI_abundancebased=c()
    overlapValuesUnprodJI_abundancebased=c()
    
    NinsertionProd1=c()
    NinsertionUnprod1=c()
    
    NinsertionProd2=c()
    NinsertionUnprod2=c()
    
    
    print(r)
    sam1=m[r,1]
    sam2=m[r,2]
    
    #print(sam1)
    #print(sam2)
    
    sam1 <- get(as.character(sam1))
    sam2 <- get(as.character(sam2))
    
    
    sam1Prod <- getProductive(sam1)
    sam1Unprod <- getUnproductive(sam1)
    
    sam2Prod <- getProductive(sam2)
    sam2Unprod <- getUnproductive(sam2)
    
    
    if(cutoffby=="PERC"){
      topclones_sam1Prod = getAbundantClonotypesbyperc(sam1Prod,cutoff)
      topclones_sam1Unprod = getAbundantClonotypesbyperc(sam1Unprod,cutoff)
      
      topclones_sam2Prod = getAbundantClonotypesbyperc(sam2Prod,cutoff)
      topclones_sam2Unprod = getAbundantClonotypesbyperc(sam2Unprod,cutoff)
      
      # commons
      commonClonesProd=length(intersect(topclones_sam1Prod$AMINOACID,topclones_sam2Prod$AMINOACID)) 
      unionProd=length(union(topclones_sam1Prod$AMINOACID,topclones_sam2Prod$AMINOACID)) 
      
      # commonality checked at the Nt level for non productive sequences, no AA given in the original data
      commonClonesUnprod=length(intersect(topclones_sam1Unprod$AMINOACID,topclones_sam2Unprod$AMINOACID))
      unionUnprod=length(union(topclones_sam1Unprod$AMINOACID,topclones_sam2Unprod$AMINOACID))
      
      if(unionUnprod==0){
        unionUnprod = 1
      }
      
      if(unionProd==0){
        unionProd = 1
      }
      
      overlapValuesProd=c(overlapValuesProd,as.numeric(commonClonesProd))
      overlapValuesUnprod=c(overlapValuesUnprod,as.numeric(commonClonesUnprod))
      
      overlapValuesProdJI=c(overlapValuesProdJI,commonClonesProd/unionProd)
      overlapValuesUnprodJI = c(overlapValuesUnprodJI,commonClonesUnprod/unionUnprod)
      
      # Jaccard abundance-based similarity index (UV/(U+V-UV)). U is total relative abundance in of shared sample 1, U total relative of shared in sample 2
      
      commonClonesProd_U_data= topclones_sam1Prod[topclones_sam1Prod$AMINOACID %in% intersect(topclones_sam1Prod$AMINOACID,topclones_sam2Prod$AMINOACID),]
      commonClonesProd_U = sum(commonClonesProd_U_data$COUNT)/sum(topclones_sam1Prod$COUNT)
      
      commonClonesProd_V_data= topclones_sam2Prod[topclones_sam2Prod$AMINOACID %in% intersect(topclones_sam1Prod$AMINOACID,topclones_sam2Prod$AMINOACID),]
      commonClonesProd_V = sum(commonClonesProd_V_data$COUNT)/sum(topclones_sam2Prod$COUNT)
      
      uv=commonClonesProd_U*commonClonesProd_V
      uplusv=commonClonesProd_U+commonClonesProd_V
      
      AbundanceBasedJIProd = uv/(uplusv-uv)
      
      overlapValuesProdJI_abundancebased=c(overlapValuesProdJI_abundancebased,as.numeric(AbundanceBasedJIProd))
      
      # for unproductive sequences
      
      commonClonesUnProd_U_data= topclones_sam1Unprod[topclones_sam1Unprod$AMINOACID %in% intersect(topclones_sam1Unprod$AMINOACID,topclones_sam2Unprod$AMINOACID),]
      commonClonesUnProd_U = sum(commonClonesUnProd_U_data$COUNT)/sum(topclones_sam1Unprod$COUNT)
      
      commonClonesUnProd_V_data= topclones_sam2Unprod[topclones_sam2Unprod$AMINOACID %in% intersect(topclones_sam1Unprod$AMINOACID,topclones_sam2Unprod$AMINOACID),]
      commonClonesUnProd_V = sum(commonClonesUnProd_V_data$COUNT)/sum(topclones_sam2Unprod$COUNT)
      
      uv_UnProd=commonClonesUnProd_U * commonClonesUnProd_V
      uplusv_UnProd=commonClonesUnProd_U + commonClonesUnProd_V
      
      AbundanceBasedJIUnProd = uv_UnProd/(uplusv_UnProd-uv_UnProd)
      
      overlapValuesUnprodJI_abundancebased=c(overlapValuesUnprodJI_abundancebased,as.numeric(AbundanceBasedJIUnProd))
      
      
      
      # Ninsertions  in common sequences
      
      commonClonotypes1= topclones_sam1Prod[topclones_sam1Prod$AMINOACID %in% intersect(topclones_sam1Prod$AMINOACID,topclones_sam2Prod$AMINOACID),]
      commonClonotypes2= topclones_sam2Prod[topclones_sam2Prod$AMINOACID %in% intersect(topclones_sam1Prod$AMINOACID,topclones_sam2Prod$AMINOACID),]
      
      commonClonotypes1Unprod= topclones_sam1Unprod[topclones_sam1Unprod$AMINOACID %in% intersect(topclones_sam1Unprod$AMINOACID,topclones_sam2Unprod$AMINOACID),]
      commonClonotypes2Unprod= topclones_sam2Unprod[topclones_sam2Unprod$AMINOACID %in% intersect(topclones_sam1Unprod$AMINOACID,topclones_sam2Unprod$AMINOACID),]
      
      
      NinsertionsP1= mean(commonClonotypes1$N2INSERTION + commonClonotypes1$N1INSERTION)
      NinsertionsP2 = mean(commonClonotypes2$N2INSERTION + commonClonotypes2$N1INSERTION)
      
      NinsertionsUP1= mean(commonClonotypes1Unprod$N2INSERTION + commonClonotypes1Unprod$N1INSERTION)
      NinsertionsUP2 = mean(commonClonotypes2Unprod$N2INSERTION + commonClonotypes2Unprod$N1INSERTION)
      
      
      NinsertionProd1=c(NinsertionProd1,NinsertionsP1)
      NinsertionProd2= c(NinsertionProd2,NinsertionsP2)
      
      NinsertionUnprod1=c(NinsertionUnprod1,NinsertionsUP1)
      NinsertionUnprod2=c(NinsertionUnprod2,NinsertionsUP2)
      
      
    }else{
      topclones_sam1Prod = getAbundantClonotypesbynumber(sam1Prod,cutoff)
      topclones_sam1Unprod = getAbundantClonotypesbynumber(sam1Unprod,cutoff)
      
      topclones_sam2Prod = getAbundantClonotypesbynumber(sam2Prod,cutoff)
      topclones_sam2Unprod = getAbundantClonotypesbynumber(sam2Unprod,cutoff)
      
      
      # commons
      # commons
      commonClonesProd=length(intersect(topclones_sam1Prod$AMINOACID,topclones_sam2Prod$AMINOACID))
      unionProd=length(union(topclones_sam1Prod$AMINOACID,topclones_sam2Prod$AMINOACID))
      
      # commonality checked at the Nt level for non productive sequences, no AA given in the original data
      commonClonesUnprod=length(intersect(topclones_sam1Unprod$AMINOACID,topclones_sam2Unprod$AMINOACID))
      unionUnprod=length(union(topclones_sam1Unprod$AMINOACID,topclones_sam2Unprod$AMINOACID))
      
      if(unionUnprod==0){
        unionUnprod = 1
      }
      
      if(unionProd==0){
        unionProd = 1
      }
      
      overlapValuesProd=c(overlapValuesProd,as.numeric(commonClonesProd))
      overlapValuesUnprod=c(overlapValuesUnprod,as.numeric(commonClonesUnprod))
      
      overlapValuesProdJI=c(overlapValuesProdJI,commonClonesProd/unionProd)
      overlapValuesUnprodJI = c(overlapValuesUnprodJI,commonClonesUnprod/unionUnprod)
      
      
      # Jaccard abundance-based similarity index (UV/(U+V-UV)). U is total relative abundance in of shared sample 1, U total relative of shared in sample 2
      
      commonClonesProd_U_data= topclones_sam1Prod[topclones_sam1Prod$AMINOACID %in% intersect(topclones_sam1Prod$AMINOACID,topclones_sam2Prod$AMINOACID),]
      commonClonesProd_U = sum(commonClonesProd_U_data$COUNT)/sum(topclones_sam1Prod$COUNT)
      
      commonClonesProd_V_data= topclones_sam2Prod[topclones_sam2Prod$AMINOACID %in% intersect(topclones_sam1Prod$AMINOACID,topclones_sam2Prod$AMINOACID),]
      commonClonesProd_V = sum(commonClonesProd_V_data$COUNT)/sum(topclones_sam2Prod$COUNT)
      
      uv=commonClonesProd_U*commonClonesProd_V
      uplusv=commonClonesProd_U+commonClonesProd_V
      
      AbundanceBasedJIProd = uv/(uplusv-uv)
      
      overlapValuesProdJI_abundancebased=c(overlapValuesProdJI_abundancebased,as.numeric(AbundanceBasedJIProd))
      
      # for unproductive sequences
      
      commonClonesUnProd_U_data= topclones_sam1Unprod[topclones_sam1Unprod$AMINOACID %in% intersect(topclones_sam1Unprod$AMINOACID,topclones_sam2Unprod$AMINOACID),]
      commonClonesUnProd_U = sum(commonClonesUnProd_U_data$COUNT)/sum(topclones_sam1Unprod$COUNT)
      
      commonClonesUnProd_V_data= topclones_sam2Unprod[topclones_sam2Unprod$AMINOACID %in% intersect(topclones_sam1Unprod$AMINOACID,topclones_sam2Unprod$AMINOACID),]
      commonClonesUnProd_V = sum(commonClonesUnProd_V_data$COUNT)/sum(topclones_sam2Unprod$COUNT)
      
      uv_UnProd=commonClonesUnProd_U * commonClonesUnProd_V
      uplusv_UnProd=commonClonesUnProd_U + commonClonesUnProd_V
      
      AbundanceBasedJIUnProd = uv_UnProd/(uplusv_UnProd-uv_UnProd)
      
      overlapValuesUnprodJI_abundancebased=c(overlapValuesUnprodJI_abundancebased,as.numeric(AbundanceBasedJIUnProd))
      
      
      
      # Ninsertions  in common sequences
      
      commonClonotypes1= topclones_sam1Prod[topclones_sam1Prod$AMINOACID %in% intersect(topclones_sam1Prod$AMINOACID,topclones_sam2Prod$AMINOACID),]
      commonClonotypes2= topclones_sam2Prod[topclones_sam2Prod$AMINOACID %in% intersect(topclones_sam1Prod$AMINOACID,topclones_sam2Prod$AMINOACID),]
      
      commonClonotypes1Unprod= topclones_sam1Unprod[topclones_sam1Unprod$AMINOACID %in% intersect(topclones_sam1Unprod$AMINOACID,topclones_sam2Unprod$AMINOACID),]
      commonClonotypes2Unprod= topclones_sam2Unprod[topclones_sam2Unprod$AMINOACID %in% intersect(topclones_sam1Unprod$AMINOACID,topclones_sam2Unprod$AMINOACID),]
      
      
      NinsertionsP1= mean(commonClonotypes1$N2INSERTION + commonClonotypes1$N1INSERTION)
      NinsertionsP2 = mean(commonClonotypes2$N2INSERTION + commonClonotypes2$N1INSERTION)
      
      NinsertionsUP1= mean(commonClonotypes1Unprod$N2INSERTION + commonClonotypes1Unprod$N1INSERTION)
      NinsertionsUP2 = mean(commonClonotypes2Unprod$N2INSERTION + commonClonotypes2Unprod$N1INSERTION)
      
      
      NinsertionProd1=c(NinsertionProd1,NinsertionsP1)
      NinsertionProd2= c(NinsertionProd2,NinsertionsP2)
      
      NinsertionUnprod1=c(NinsertionUnprod1,NinsertionsUP1)
      NinsertionUnprod2=c(NinsertionUnprod2,NinsertionsUP2)
      
    }
    
    #totalNumberOfProdClonotypes= nrow(topclones_sam1Prod) + nrow(topclones_sam2Prod)
    #totalNumberOfUnprodClonotypes= nrow(topclones_sam1Unprod) + nrow(topclones_sam2Unprod)
    
    
    if(nrow(toadd) > 0){
      newrow=data.frame(overlapValuesProd,overlapValuesUnprod,overlapValuesProdJI,overlapValuesUnprodJI,overlapValuesProdJI_abundancebased,overlapValuesUnprodJI_abundancebased,NinsertionProd1,NinsertionProd2,NinsertionUnprod1,NinsertionUnprod2)
      colnames(newrow) <- c("OverlapCountProductive","OverlapCountUnProductive","OverlapProd_JI","OverlapUnProd_JI","AbundanceBased_JI","AbundanceBased_JI","NinsertionProdSample1","NinsertionProdSample2","NinsertionUnProdSample1","NinsertionUnProdSample2")
      
      toadd=rbind(toadd,newrow)
    }else{
      toadd=data.frame(overlapValuesProd,overlapValuesUnprod,overlapValuesProdJI,overlapValuesUnprodJI,overlapValuesProdJI_abundancebased,overlapValuesUnprodJI_abundancebased,NinsertionProd1,NinsertionProd2,NinsertionUnprod1,NinsertionUnprod2)
      colnames(toadd) <- c("OverlapCountProductive","OverlapCountUnProductive","OverlapProd_JI","OverlapUnProd_JI","AbundanceBased_JI","AbundanceBased_JI","NinsertionProdSample1","NinsertionProdSample2","NinsertionUnProdSample1","NinsertionUnProdSample2")  
    }
    
  }
  
  
  
  
  print(nrow(toadd))
  resultTable=data.frame(m,toadd)        
  
  # write pdf dendrogram and heatmaps
  
  # select the sample columns 1 and 2, and the JI column 5
  forProdselected = resultTable[,c(1,2,5)]
  forProdselected = forProdselected[order(forProdselected$Sample1),]
  
  l = length(unique(forProdselected$Sample1))
  mat = matrix(as.numeric(forProdselected[,3]),nrow=l,ncol=l,byrow=T)
  # replace NaN values with 0
  mat[is.nan(mat)] = 0
  colnames(mat) <- unique(forProdselected$Sample1)
  rownames(mat) <- unique(forProdselected$Sample1)
  
  
  mat = signif(mat,digits=2)
  hc = hclust(as.dist(log10(1-mat)))
  
  flname = paste("ProdJI_dendrogram.pdf")
  pdf(file=flname, useDingbats = FALSE)
  
  # very simple dendrogram
  plot(hc)
  
  dev.off()
  
  # for unprod
  # select the sample columns 1 and 2, and the JI column 5
  forUnprodselected = resultTable[,c(1,2,6)]
  forUnprodselected = forUnprodselected[order(forUnprodselected$Sample1),]
  
  l = length(unique(forUnprodselected$Sample1))
  mat = matrix(as.numeric(forUnprodselected[,3]),nrow=l,ncol=l,byrow=T)
  # replace NaN values with 0
  mat[is.nan(mat)] = 0
  colnames(mat) <- unique(forUnprodselected$Sample1)
  rownames(mat) <- unique(forUnprodselected$Sample1)
  
  
  mat = signif(mat,digits=2)
  hc = hclust(as.dist(log10(1-mat)))
  
  
  flname = paste("UnprodJI_dendrogram.pdf")
  pdf(file=flname, useDingbats = FALSE)
  
  # very simple dendrogram
  plot(hc)
  
  dev.off()
  
  write.table(resultTable,file="Pairwise AA Overlap Counts.txt",row.names = T,sep = "\t")
  return(resultTable)
}


# Overlap sequence calculation using various parameters

overLappingClonotypesByChoice <- function(reps,cutoff=100,cutoffby="perc",seqType="AA",unProdStatus="unprod"){
  ### Returns a datamatrix of top topnumber clonotopes from each sample plus: overlap count or Jaccard index, for prod and non productive sequences, N insertion size
  
  ### Total read counts for the common clones of both samples are given for the productive/ chosen unproductive types.
  ### cutoffby : Perc (percentage of total repertoire as indicated by cutoff) or noperc (select top clones as many as indicated by cutoff)
  ### seqType:  can be by Nt (nucleotide) or AA (AminoAcid)
  ### unProdStatus : which unproductive clones to consider Out or Stop or Unprod (union of Out and Stop)
  
  
  # setting work parameters
  #---------------------------------------------------------
  
  #reps <- samNames
  
  cutoffby = toupper(cutoffby)
  
  m <- as.data.frame(combinations(n = length(reps), r = 2, v = reps, repeats.allowed = TRUE))
  colnames(m)  <- c("Sample1","Sample2")
  
  
  
  toadd= data.frame()
  
  #change parameter values to uppercase to accept Nt, nt etc for example
  unProdStatus <- toupper(unProdStatus)
  seqType <- toupper(seqType)
  
  if(seqType=="AA"){
    seqType <- "AMINOACID"
  }else{
    seqType <- "NUCLEOTIDE"
  }
  
  #---------------------------------------------------------
  
  
  for(r in 1:nrow(m)){
    
    overlapValuesProd=c()
    overlapValuesUnprod=c()
    
    overlapValuesTotalSam1Prod=c()
    overlapValuesTotalSam2Prod=c()
    
    overlapValuesTotalSam1Unprod=c()
    overlapValuesTotalSam2Unprod=c()
    
    overlapValuesProdJI=c()
    overlapValuesUnprodJI=c()
    
    # the following two hold abundance based JI using : http://chao.stat.nthu.edu.tw/wordpress/paper/2006_Biometrics_62_P361.pdf
    
    overlapValuesProdJI_abundancebased=c()
    overlapValuesUnprodJI_abundancebased=c()
    
    # the following two hold abundance based JI calculation using (total size intersecting clones)/(total size of union clones)
    overlapValuesProdJI_abundancebased2=c()
    overlapValuesUnprodJI_abundancebased2=c()
    
    NinsertionProd1=c()
    NinsertionUnprod1=c()
    
    NinsertionProd2=c()
    NinsertionUnprod2=c()
    
    
    print(r)
    sam1=m[r,1]
    sam2=m[r,2]
    
    samplesCompared = paste(sam1,sam2,sep="_")
    #print(sam1)
    #print(sam2)
    
    sam1 <- get(as.character(sam1))
    sam2 <- get(as.character(sam2))
    
    #---------------------------------------------------------
    # productive and unproductive selection for each sample
    
    
    sam1Prod <- getProductive(sam1)
    
    if(unProdStatus=="STOP"){
      sam1Unprod <- getHasStop(sam1)
    }else if(unProdStatus=="OUT"){
      sam1Unprod <- getOutOfFrame(sam1)
    }else{
      sam1Unprod <- getUnproductive(sam1)
    }
    
    
    sam2Prod <- getProductive(sam2)
    
    if(unProdStatus=="STOP"){
      sam2Unprod <- getHasStop(sam2)
    }else if(unProdStatus=="OUT"){
      sam2Unprod <- getOutOfFrame(sam2)
    }else{
      sam2Unprod <- getUnproductive(sam2)
    }
    
    
    
    #---------------------------------------------------------
    # selection of sequence type index for each sample
    
    indSam1 <- which(colnames(sam1) == seqType)
    indSam2 <- which(colnames(sam2) == seqType)
    
    
    
    #---------------------------------------------------------
    # selecting top clonotypes for which to do the overlap counts, JI calculation etc. Either by percentage or actual number of top clonotypes
    
    if(cutoffby=="PERC"){
      topclones_sam1Prod = getAbundantClonotypesbyperc(sam1Prod,cutoff)
      topclones_sam1Unprod = getAbundantClonotypesbyperc(sam1Unprod,cutoff)
      
      topclones_sam2Prod = getAbundantClonotypesbyperc(sam2Prod,cutoff)
      topclones_sam2Unprod = getAbundantClonotypesbyperc(sam2Unprod,cutoff)
    }else{
      topclones_sam1Prod = getAbundantClonotypesbynumber(sam1Prod,cutoff)
      topclones_sam1Unprod = getAbundantClonotypesbynumber(sam1Unprod,cutoff)
      
      topclones_sam2Prod = getAbundantClonotypesbynumber(sam2Prod,cutoff)
      topclones_sam2Unprod = getAbundantClonotypesbynumber(sam2Unprod,cutoff)
    }
    
    #---------------------------------------------------------
    # Getting number of overlapping clonotypes, Jaccard index, average insertion size etc for the two samples being compared
    
    # commons for sequence type selected
    commonClonesProd=length(intersect(topclones_sam1Prod[,indSam1],topclones_sam2Prod[,indSam2])) 
    unionProd=length(union(topclones_sam1Prod[,indSam1],topclones_sam2Prod[,indSam2])) 
    
    # commonality checked at the Nt level for non productive sequences, no AA given in the original data
    commonClonesUnprod=length(intersect(topclones_sam1Unprod[,indSam1],topclones_sam2Unprod[,indSam2]))
    unionUnprod=length(union(topclones_sam1Unprod[,indSam1],topclones_sam2Unprod[,indSam2]))
    
    
    
    
    if(unionUnprod==0){
      unionUnprod = 1
    }
    
    if(unionProd==0){
      unionProd = 1
    }
    
    overlapValuesProd=c(overlapValuesProd,as.numeric(commonClonesProd))
    overlapValuesUnprod=c(overlapValuesUnprod,as.numeric(commonClonesUnprod))
    
    
    
    overlapValuesProdJI=c(overlapValuesProdJI,commonClonesProd/unionProd)
    overlapValuesUnprodJI = c(overlapValuesUnprodJI,commonClonesUnprod/unionUnprod)
    
    # Jaccard abundance-based similarity index (UV/(U+V-UV)). U is total relative abundance in of shared sample 1, U total relative of shared in sample 2
    # taken from : http://chao.stat.nthu.edu.tw/wordpress/paper/2006_Biometrics_62_P361.pdf
    
    
    commonClonesProd_U_data= topclones_sam1Prod[topclones_sam1Prod[,indSam1] %in% intersect(topclones_sam1Prod[,indSam1],topclones_sam2Prod[,indSam2]),]
    commonClonesProd_U = sum(commonClonesProd_U_data$COUNT)/sum(topclones_sam1Prod$COUNT)
    
    commonClonesProd_U_abundance = sum(commonClonesProd_U_data$COUNT)
    
    commonClonesProd_V_data= topclones_sam2Prod[topclones_sam2Prod[,indSam2] %in% intersect(topclones_sam1Prod[,indSam1],topclones_sam2Prod[,indSam2]),]
    commonClonesProd_V = sum(commonClonesProd_V_data$COUNT)/sum(topclones_sam2Prod$COUNT)
    
    commonClonesProd_V_abundance = sum(commonClonesProd_V_data$COUNT)
    
    overlapValuesTotalSam1Prod=c(overlapValuesTotalSam1Prod,as.numeric(commonClonesProd_U_abundance))
    overlapValuesTotalSam2Prod=c(overlapValuesTotalSam2Prod,as.numeric(commonClonesProd_V_abundance))
    
    uv=commonClonesProd_U*commonClonesProd_V
    uplusv=commonClonesProd_U+commonClonesProd_V
    
    AbundanceBasedJIProd = uv/(uplusv-uv)
    
    overlapValuesProdJI_abundancebased=c(overlapValuesProdJI_abundancebased,as.numeric(AbundanceBasedJIProd))
    
    # for abundance based JI type 2
    totalSam1Sam2Prod = sum(topclones_sam1Prod$COUNT) + sum(topclones_sam2Prod$COUNT)                                                                                        
    AbundanceBasedJIProd2 = (commonClonesProd_U_abundance + commonClonesProd_V_abundance) / totalSam1Sam2Prod
    overlapValuesProdJI_abundancebased2 = c(overlapValuesProdJI_abundancebased2,as.numeric(AbundanceBasedJIProd2))
    
    
    # for unproductive sequences
    
    commonClonesUnProd_U_data= topclones_sam1Unprod[topclones_sam1Unprod[,indSam1] %in% intersect(topclones_sam1Unprod[,indSam1],topclones_sam2Unprod[,indSam2]),]
    commonClonesUnProd_U = sum(commonClonesUnProd_U_data$COUNT)/sum(topclones_sam1Unprod$COUNT)
    
    commonClonesUnProd_U_abundance = sum(commonClonesUnProd_U_data$COUNT)
    
    commonClonesUnProd_V_data= topclones_sam2Unprod[topclones_sam2Unprod[,indSam2] %in% intersect(topclones_sam1Unprod[,indSam1],topclones_sam2Unprod[,indSam2]),]
    commonClonesUnProd_V = sum(commonClonesUnProd_V_data$COUNT)/sum(topclones_sam2Unprod$COUNT)
    
    commonClonesUnProd_V_abundance = sum(commonClonesUnProd_V_data$COUNT)
    
    
    overlapValuesTotalSam1Unprod=c(overlapValuesTotalSam1Unprod,as.numeric(commonClonesUnProd_U_abundance))
    overlapValuesTotalSam2Unprod=c(overlapValuesTotalSam2Unprod,as.numeric(commonClonesUnProd_V_abundance))
    
    
    uv_UnProd=commonClonesUnProd_U * commonClonesUnProd_V
    uplusv_UnProd=commonClonesUnProd_U + commonClonesUnProd_V
    
    AbundanceBasedJIUnProd = uv_UnProd/(uplusv_UnProd-uv_UnProd)
    
    overlapValuesUnprodJI_abundancebased=c(overlapValuesUnprodJI_abundancebased,as.numeric(AbundanceBasedJIUnProd))
    
    # for abundance based JI type 2
    totalSam1Sam2Unprod = sum(topclones_sam1Unprod$COUNT) + sum(topclones_sam2Unprod$COUNT)                                                                                          
    AbundanceBasedJIUnProd2 = (commonClonesUnProd_U_abundance + commonClonesUnProd_V_abundance) / totalSam1Sam2Unprod
    overlapValuesUnprodJI_abundancebased2 = c(overlapValuesUnprodJI_abundancebased2,as.numeric(AbundanceBasedJIUnProd2))
    
    
    # Ninsertions  in common sequences
    
    commonClonotypes1= topclones_sam1Prod[topclones_sam1Prod[,indSam1] %in% intersect(topclones_sam1Prod[,indSam1],topclones_sam2Prod[,indSam1]),]
    commonClonotypes2= topclones_sam2Prod[topclones_sam2Prod[,indSam2] %in% intersect(topclones_sam1Prod[,indSam1],topclones_sam2Prod[,indSam2]),]
    
    commonClonotypes1Unprod= topclones_sam1Unprod[topclones_sam1Unprod[,indSam1] %in% intersect(topclones_sam1Unprod[,indSam1],topclones_sam2Unprod[,indSam2]),]
    commonClonotypes2Unprod= topclones_sam2Unprod[topclones_sam2Unprod[,indSam2] %in% intersect(topclones_sam1Unprod[,indSam1],topclones_sam2Unprod[,indSam2]),]
    
    
    NinsertionsP1= mean(commonClonotypes1$N2INSERTION + commonClonotypes1$N1INSERTION)
    NinsertionsP2 = mean(commonClonotypes2$N2INSERTION + commonClonotypes2$N1INSERTION)
    
    NinsertionsUP1= mean(commonClonotypes1Unprod$N2INSERTION + commonClonotypes1Unprod$N1INSERTION)
    NinsertionsUP2 = mean(commonClonotypes2Unprod$N2INSERTION + commonClonotypes2Unprod$N1INSERTION)
    
    
    NinsertionProd1=c(NinsertionProd1,NinsertionsP1)
    NinsertionProd2= c(NinsertionProd2,NinsertionsP2)
    
    NinsertionUnprod1=c(NinsertionUnprod1,NinsertionsUP1)
    NinsertionUnprod2=c(NinsertionUnprod2,NinsertionsUP2)
    
    
    # write some common clonotypes of productive and unproductive repertoires
    # select columns to merge and write 
    
    selectedCols <- c("NUCLEOTIDE","AMINOACID","COUNT","VFAMILYNAME","VGENENAME","JGENENAME","CDR3LENGTH","SEQUENCESTATUS")
    selectedColsIdx <- which(colnames(commonClonotypes1) %in% selectedCols)
    
    # write for productive commons
    
    resultDir = paste("Pairwise_Overlap_results_","Productive_Unproductive(",unProdStatus,")",sep="")
    
    dir.create(file.path(resultDir),showWarnings = FALSE)
    
    
    commonsProdForWriting = head(merge(commonClonotypes1[,selectedColsIdx],commonClonotypes2[,selectedColsIdx],by="NUCLEOTIDE"),20)
    commonsUnProdForWriting = merge(commonClonotypes1Unprod[,selectedColsIdx],commonClonotypes2Unprod[,selectedColsIdx],by="NUCLEOTIDE")
    
    
    write.table(commonsProdForWriting,file=paste(resultDir,"/",samplesCompared,"_Common_productive_sequences_top20.txt",sep=""),row.names = F,sep = "\t")
    write.table(commonsUnProdForWriting,file=paste(resultDir,"/",samplesCompared,"_Common_Unproductive_sequences_all.txt",sep=""),row.names = F,sep = "\t")
    
    
    
    #totalNumberOfProdClonotypes= nrow(topclones_sam1Prod) + nrow(topclones_sam2Prod)
    #totalNumberOfUnprodClonotypes= nrow(topclones_sam1Unprod) + nrow(topclones_sam2Unprod)
    
    
    #unproductiveOverlapCountName = paste("OverlapCountUnProductive(",unProdStatus,")",sep="")
    
    if(nrow(toadd) > 0){
      newrow=data.frame(overlapValuesProd,overlapValuesTotalSam1Prod,overlapValuesTotalSam2Prod,overlapValuesUnprod,overlapValuesTotalSam1Unprod,overlapValuesTotalSam2Unprod,overlapValuesProdJI,overlapValuesUnprodJI,overlapValuesProdJI_abundancebased,overlapValuesUnprodJI_abundancebased,overlapValuesProdJI_abundancebased2,overlapValuesUnprodJI_abundancebased2,NinsertionProd1,NinsertionProd2,NinsertionUnprod1,NinsertionUnprod2)
      colnames(newrow) <- c("OverlapCountProductive","OverlapCountTotalProductiveSam1","OverlapCountTotalProductiveSam2","OverlapCountUnProductive","OverlapCountTotalUnProductiveSam1","OverlapCountTotalUnProductiveSam2","OverlapProd_JI","OverlapUnProductive_JI","AbundanceBasedProd_JI","AbundanceBasedUnProductive_JI","AbundanceBasedProd_JI2","AbundanceBasedUnProductive_JI2","NinsertionProdSample1","NinsertionProdSample2","NinsertionUnProductiveSample1","NinsertionUnProductiveSample2")
      
      toadd=rbind(toadd,newrow)
    }else{
      toadd=data.frame(overlapValuesProd,overlapValuesTotalSam1Prod,overlapValuesTotalSam2Prod,overlapValuesUnprod,overlapValuesTotalSam1Unprod,overlapValuesTotalSam2Unprod,overlapValuesProdJI,overlapValuesUnprodJI,overlapValuesProdJI_abundancebased,overlapValuesUnprodJI_abundancebased,overlapValuesProdJI_abundancebased2,overlapValuesUnprodJI_abundancebased2,NinsertionProd1,NinsertionProd2,NinsertionUnprod1,NinsertionUnprod2)
      colnames(toadd) <- c("OverlapCountProductive","OverlapCountTotalProductiveSam1","OverlapCountTotalProductiveSam2","OverlapCountUnProductive","OverlapCountTotalUnProductiveSam1","OverlapCountTotalUnProductiveSam2","OverlapProd_JI","OverlapUnProductive_JI","AbundanceBasedProd_JI","AbundanceBasedUnProductive_JI","AbundanceBasedProd_JI2","AbundanceBasedUnProductive_JI2","NinsertionProdSample1","NinsertionProdSample2","NinsertionUnProductiveSample1","NinsertionUnProductiveSample2")
    }
    
  }
  
  
  
  
  
  #head(toadd)
  #print(nrow(toadd))
  resultTable=data.frame(m,toadd)        
  
  outputFileName = paste(resultDir,"/","Pairwise_Overlap_Counts_",seqType,"_Productive_Unproductive(",unProdStatus,").txt",sep="")
  
  write.table(resultTable,file=outputFileName,row.names = T,sep = "\t")
  return(resultTable)
  
  
}




#### Making PCA for V or VJ usage ####

getPCAfordata <- function(dd,groups,scaling=T,oprefix="O"){
  # Assumes that the observations (samples) are on columns and variables on rows. We want to make the variables
  # on columns so it first transposes d.
  dd=t(dd)
  #head(d)
  
  pr<-prcomp(dd,scale=scaling)
  print(summary(pr)) # proportions of the total variance explained by each component
  
  #samples names
  sampleNames=rownames(dd)
  
  # get variance explained by each PC
  s=summary(pr)$importance[2,]
  
  # grouping of samples
  #groups=groups
  
  # Write and display PCA result
  
  flname = paste(oprefix,"PCA.pdf")
  pdf(file=flname, useDingbats = FALSE)
  
  tobedrawn=data.frame(pr$x,group=factor(groups),snames=rownames(dd))
  write.table(pr$x,file=paste(oprefix,"_PCArotatedData.txt"),row.names = T,col.names=T,sep = "\t")
  
  pcaP=ggplot(tobedrawn, aes(x=PC1, y=PC2,shape=group,colour=snames)) + xlab(paste("PC1 (",format(round(s[1]*100),nsmall=2),"%)")) + ylab(paste("PC2 (",format(round(s[2]*100),nsmall=2),"%)")) + geom_point(size=6) + scale_fill_brewer(palette="Set1") 
  pcaP = pcaP + theme_bw() + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank()) + theme(axis.line = element_line(colour = "black")) 
  
  
  
  
  print(pcaP)  
  #png(file=paste(oprefix,"PCA.png"),height=700, width=650)
  
  #Plot using ggplot, much more flexible
  #ggplot(as.data.frame(pr$x), aes(x=PC1, y=PC2)) + xlab(paste("PC1 (",format(round(s[1]*100),nsmall=2),"%)")) + ylab(paste("PC2 (",format(round(s[2]*100),nsmall=2),"%)")) + geom_point(size=6,aes(colour = sampleNames,shape=groups)) + scale_fill_brewer(palette="Set1") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
                                                                                                                                                                                                                                                                              
  
  dev.off()
  
  #ggsave(flname, plot = last_plot(), device = "pdf")
  
  #print(pcaP)
  
}

getPCAfordata1 <- function(dd,groups,oprefix="O"){
  # Assumes that the observations (samples) are on columns and variables on rows. We want to make the variables
  # on columns so it first transposes d.
  dd=t(dd)
  #head(d)
  
  pr<-prcomp(dd)
  print(summary(pr)) # proportions of the total variance explained by each component
  
  #samples names
  sampleNames=rownames(dd)
  
  # get variance explained by each PC
  s=summary(pr)$importance[2,]
  
  # grouping of samples
  #groups=groups
  
  # Write and display PCA result
  flname = paste(oprefix,"PCA.pdf")
  pdf(file=flname, useDingbats = FALSE)
  
  
  tobedrawn=data.frame(pr$x,group=factor(groups),snames=rownames(dd))
  write.table(pr$x,file=paste(oprefix,"_PCArotatedData.txt"),row.names = T,col.names=T,sep = "\t")
  
  pcaP=ggplot(tobedrawn, aes(x=PC1, y=PC2,shape=group,colour=snames)) + xlab(paste("PC1 (",format(round(s[1]*100),nsmall=2),"%)")) + ylab(paste("PC2 (",format(round(s[2]*100),nsmall=2),"%)")) + geom_point(size=6) + scale_fill_brewer(palette="Set1")
  pcaP = pcaP + theme_bw() + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank()) + theme(axis.line = element_line(colour = "black")) 
  
  
  #png(file=paste(oprefix,"PCA.png"),height=700, width=650)
  
  #Plot using ggplot, much more flexible
  #ggplot(as.data.frame(pr$x), aes(x=PC1, y=PC2)) + xlab(paste("PC1 (",format(round(s[1]*100),nsmall=2),"%)")) + ylab(paste("PC2 (",format(round(s[2]*100),nsmall=2),"%)")) + geom_point(size=6,aes(colour = sampleNames,shape=groups)) + scale_fill_brewer(palette="Set1") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
  print(pcaP)                                                                                                                                                                                                                                                                              
  
  dev.off()
  
  
  
  #print(pcaP)
  
}


getPCAfordata2 <- function(dd,groups,oprefix="O"){
  # Assumes that the observations (samples) are on columns and variables on rows. We want to make the variables
  # on columns so it first transposes d.
  dd=t(dd)
  #head(d)
  
  pr<-prcomp(dd)
  print(summary(pr)) # proportions of the total variance explained by each component
  
  #samples names
  sampleNames=rownames(dd)
  
  # get variance explained by each PC
  s=summary(pr)$importance[2,]
  
  # grouping of samples
  #groups=groups
  
  # Write and display PCA result
  flname = paste(oprefix,"PCA.pdf")
  pdf(file=flname, useDingbats = FALSE)
  
  tobedrawn=data.frame(pr$x,repType=factor(groups),snames=rownames(dd))
  write.table(pr$x,file=paste(oprefix,"_PCArotatedData.txt"),row.names = T,col.names=T,sep = "\t")
  
  #pcaP=ggplot(tobedrawn, aes(x=PC1, y=PC2,shape=repType,colour=repType,label=snames))
  # pcaP=ggplot(tobedrawn, aes(x=PC1, y=PC2,shape=repType,colour=repType,label=snames))
  pcaP=ggplot(tobedrawn, aes(x=PC1, y=PC2,label=snames)) + 
    xlab(paste("PC1 (",format(round(s[1]*100),nsmall=2),"%)")) + 
    ylab(paste("PC2 (",format(round(s[2]*100),nsmall=2),"%)")) + 
    scale_fill_brewer(palette="Set1")
  
  pcaP = pcaP + geom_text(size=3,check_overlap = TRUE,aes(colour = repType, fontface = "bold"),position=position_jitter())
  
  pcaP = pcaP + theme_bw() + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank()) + theme(axis.line = element_line(colour = "black")) 
  
  
  
  #png(file=paste(oprefix,"PCA.png"),height=700, width=650)
  
  #Plot using ggplot, much more flexible
  #ggplot(as.data.frame(pr$x), aes(x=PC1, y=PC2)) + xlab(paste("PC1 (",format(round(s[1]*100),nsmall=2),"%)")) + ylab(paste("PC2 (",format(round(s[2]*100),nsmall=2),"%)")) + geom_point(size=6,aes(colour = sampleNames,shape=groups)) + scale_fill_brewer(palette="Set1") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
  print(pcaP)                                                                                                                                                                                                                                                                              
  
  dev.off()
  
  
  
  #print(pcaP)
  
}

getPCAfordata3 <- function(dd,groups,oprefix="O"){
  # Assumes that the observations (samples) are on columns and variables on rows. We want to make the variables
  # on columns so it first transposes d.
  dd=t(dd)
  #head(d)
  
  pr<-prcomp(dd)
  print(summary(pr)) # proportions of the total variance explained by each component
  
  #samples names
  sampleNames=rownames(dd)
  
  # get variance explained by each PC
  s=summary(pr)$importance[2,]
  
  # grouping of samples
  #groups=groups
  
  # Write and display PCA result
  flname = paste(oprefix,"PCA.pdf")
  pdf(file=flname, useDingbats = FALSE)
  
  
  tobedrawn=data.frame(pr$x,repType=factor(groups),snames=rownames(dd))
  write.table(pr$x,file=paste(oprefix,"_PCArotatedData.txt"),row.names = T,col.names=T,sep = "\t")
  
  pcaP=ggplot(tobedrawn, aes(x=PC1, y=PC2,shape=repType,colour=repType)) + xlab(paste("PC1 (",format(round(s[1]*100),nsmall=2),"%)")) + ylab(paste("PC2 (",format(round(s[2]*100),nsmall=2),"%)")) + geom_point(size=6) + scale_fill_brewer(palette="Set1") 
  pcaP = pcaP + theme_bw() + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank()) + theme(axis.line = element_line(colour = "black")) 
  
  #png(file=paste(oprefix,"PCA.png"),height=700, width=650)
  
  #Plot using ggplot, much more flexible
  #ggplot(as.data.frame(pr$x), aes(x=PC1, y=PC2)) + xlab(paste("PC1 (",format(round(s[1]*100),nsmall=2),"%)")) + ylab(paste("PC2 (",format(round(s[2]*100),nsmall=2),"%)")) + geom_point(size=6,aes(colour = sampleNames,shape=groups)) + scale_fill_brewer(palette="Set1") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
  print(pcaP)                                                                                                                                                                                                                                                                              
  
  dev.off()
  
  
  
  #print(pcaP)
  
}


getPCAfordata4 <- function(dd,groups,scaling=T,oprefix="O"){
  # Assumes that the observations (samples) are on columns and variables on rows. We want to make the variables
  # on columns so it first transposes d.
  dd=t(dd)
  #head(d)
  
  pr<-prcomp(dd,scale=scaling)
  print(summary(pr)) # proportions of the total variance explained by each component
  
  #samples names
  sampleNames=rownames(dd)
  
  # get variance explained by each PC
  s=summary(pr)$importance[2,]
  
  # grouping of samples
  #groups=groups
  
  # Write and display PCA result
  
  flname = paste(oprefix,"PCA.pdf")
  pdf(file=flname, useDingbats = FALSE)
  
  tobedrawn=data.frame(pr$x,group=factor(groups),snames=rownames(dd))
  write.table(pr$x,file=paste(oprefix,"_PCArotatedData.txt"),row.names = T,col.names=T,sep = "\t")
  
  pcaP=ggplot(tobedrawn, aes(x=PC1, y=PC2,colour=group)) + xlab(paste("PC1 (",format(round(s[1]*100),nsmall=2),"%)")) + ylab(paste("PC2 (",format(round(s[2]*100),nsmall=2),"%)")) + geom_point(size=4) + scale_fill_brewer(palette="Set1") 
  pcaP = pcaP + theme_bw() + theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank()) + theme(axis.line = element_line(colour = "black")) 
  
  
  
  
  print(pcaP)  
  #png(file=paste(oprefix,"PCA.png"),height=700, width=650)
  
  #Plot using ggplot, much more flexible
  #ggplot(as.data.frame(pr$x), aes(x=PC1, y=PC2)) + xlab(paste("PC1 (",format(round(s[1]*100),nsmall=2),"%)")) + ylab(paste("PC2 (",format(round(s[2]*100),nsmall=2),"%)")) + geom_point(size=6,aes(colour = sampleNames,shape=groups)) + scale_fill_brewer(palette="Set1") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) 
  
  
  dev.off()
  
  #ggsave(flname, plot = last_plot(), device = "pdf")
  
  #print(pcaP)
  
}


#### Fisher's exact test for V,J or VJ usage counts ####

compareUsageProportions<-function(freqTable,pcutoff=0.05){
  ### freqTable : an Feature frequency table between sample groups, table has headers, headers 2 and 3 have the counts
  
  sam=freqTable
  #hdnames= colnames(sam)
  
  #hlaAlleles <- as.data.frame(sam[,1],stringsAsFactors=FALSE )
  #colnames(hlaAlleles)  <- hdnames[1]
  
  #signficanceFound=T # a flag to check if significant comparison is found in each round
  #topsignificantAllele=F # a flag to hold the allele showing the most significant difference between sample groups in each round, will be removed during next round comparison
  
  resultTable <- sam
  pvals=c()
  ORs=c()
  CIs=c()
  
  print(paste("There are",ncol(sam),"samples in the data. Comparison is done for every possible pair of samples..."))
  
  print(paste("Number of features to compare : ",nrow(sam)))
  
  toploop=ncol(sam)-1
  inloop=ncol(sam)
  
  for (j in 1:toploop){
    
    for (k in j+1:inloop){
      #print(inloop)
      #print(j)
      #print(k)
      if (k > inloop){break}
      res=compareUsageInTwoSamples(sam[,c(j,k)])
      colnames(res) <- c(paste(colnames(sam[,c(j,k)])[1],colnames(sam[,c(j,k)])[2],"padjusted",sep="_"),paste(colnames(sam[,c(j,k)])[1],colnames(sam[,c(j,k)])[2],"OR",sep="_"),paste(colnames(sam[,c(j,k)])[1],colnames(sam[,c(j,k)])[2],"CI",sep="_"))
      resultTable <- cbind(resultTable,res) 
    }
    
  }
  
  print("Result written to file....")
  write.table(resultTable,file="Muke_Usage_fishers_comparison_result.txt",row.names = T,sep = "\t")
  return(resultTable)
  
  
}


compareUsageInTwoSamples<-function(freqTable,pcutoff=0.05){
  sam=freqTable
  pvals=c()
  ORs=c()
  CIs=c()
  
  print(head(sam))
  
  for(i in 1:nrow(sam)){
    
    current<-sam[i,]
    othersSum<-colSums(sam[-i,])
    #print(othersSum)
    
    dataForComparison=matrix(c(as.numeric(current),as.numeric(othersSum)),nrow = 2,byrow=T)
    #colnames(dataForComparison)<-c(hdnames[2], hdnames[3])
    #rownames(dataForComparison)<-c(as.character(sam[i,i]), "OtherAlleles")
    
    fresult=fisher.test(dataForComparison)
    pvals=c(pvals,fresult$p.value)
    ORs=c(ORs,as.numeric(fresult$estimate))
    CIs=c(CIs,paste(fresult$conf[1],fresult$conf[2],sep="-"))
    
    #pvals=p.adjust(pvals, "BH") # pvals are adjusted
    
    #cnames=c(paste(colnames(sam)[selectedCols],"pval",sep="") ,paste(colnames(sam)[selectedCols],"OR",sep=""),paste(colnames(sam)[selectedCols],"CI",sep=""))
    
    
  }
  
  pvals=p.adjust(pvals, "BH") # pvals are adjusted
  res = cbind(pvals,ORs,CIs)
  
  return(res)  
}


compareAbundanceInTwoSamplesForShared<-function(freqTable,pcutoff=0.05){
  sam=freqTable
  pvals=c()
  ORs=c()
  CIs=c()
  
  print(head(sam))
  
  sharedClonesIndices = apply(sam,1,function(x) (sum(x>0) == 2))
  
  sharedCloneIndices = which(sharedClonesIndices)
  numOfShared = length(sharedCloneIndices)
  
    
  for(i in 1:numOfShared){
    
    current<-sam[sharedCloneIndices[i],]
    othersSum<-colSums(sam[-sharedCloneIndices[i],])
    #print(othersSum)
    
    dataForComparison=matrix(c(as.numeric(current),as.numeric(othersSum)),nrow = 2,byrow=T)
    #colnames(dataForComparison)<-c(hdnames[2], hdnames[3])
    #rownames(dataForComparison)<-c(as.character(sam[i,i]), "OtherAlleles")
    
    fresult=fisher.test(dataForComparison)
    pvals=c(pvals,fresult$p.value)
    ORs=c(ORs,as.numeric(fresult$estimate))
    CIs=c(CIs,paste(fresult$conf[1],fresult$conf[2],sep="-"))
    
    #pvals=p.adjust(pvals, "BH") # pvals are adjusted
    
    #cnames=c(paste(colnames(sam)[selectedCols],"pval",sep="") ,paste(colnames(sam)[selectedCols],"OR",sep=""),paste(colnames(sam)[selectedCols],"CI",sep=""))
    
    
  }
  
  pvals=p.adjust(pvals, "BH") # pvals are adjusted
  res = cbind(as.numeric(pvals),ORs,CIs)
  
  rownames(res) = rownames(sam[sharedCloneIndices,])
  
  return(res)  
}


compareUsageProportionsPaired<-function(freqTable,pcutoff=0.05,pairs){
  ### freqTable : an Feature frequency table between sample groups, table has headers, headers 2 and 3 have the counts
  
  sam=freqTable
  #hdnames= colnames(sam)
  
  #hlaAlleles <- as.data.frame(sam[,1],stringsAsFactors=FALSE )
  #colnames(hlaAlleles)  <- hdnames[1]
  
  #signficanceFound=T # a flag to check if significant comparison is found in each round
  #topsignificantAllele=F # a flag to hold the allele showing the most significant difference between sample groups in each round, will be removed during next round comparison
  
  resultTable <- sam
  pvals=c()
  ORs=c()
  CIs=c()
  
  
  
  print(paste("There are",ncol(sam),"samples in the data. Comparison is done for pair of samples..."))
  
  print(paste("Number of features to compare : ",nrow(sam)))
  
  rnds = length(pairs)/2
  for (k in 1:rnds){
      kpair = k + rnds
      res=compareUsageInTwoSamples(sam[,c(k,kpair)])
      colnames(res) <- c(paste(colnames(sam[,c(k,kpair)])[1],colnames(sam[,c(k,kpair)])[2],"padjusted",sep="_"),paste(colnames(sam[,c(k,kpair)])[1],colnames(sam[,c(k,kpair)])[2],"OR",sep="_"),paste(colnames(sam[,c(k,kpair)])[1],colnames(sam[,c(k,kpair)])[2],"CI",sep="_"))
      resultTable <- cbind(resultTable,res) 
    }
    
  
  
  print("Result written to file....")
  write.table(resultTable,file="Usage_fishers_comparison_result_forPairSamples.txt",row.names = T,sep = "\t")
  return(resultTable)
  
  
}


compareAbundanceInPairedSamples<-function(freqTable,pcutoff=0.01,pairs,writeResult=F,resultDir="matcheDA"){
  ### freqTable : an clone count table between sample groups, table has headers, headers 1 and 2 have the counts
  
  sam=freqTable
  
  resultlist <- list()
  pvals=c()
  ORs=c()
  CIs=c()
  
  
  
  print(paste("There are",ncol(sam),"samples in the data. Comparison is done for pair of samples..."))
  
  print(paste("Number of features to compare : ",nrow(sam)))
  
  
  # dir to write results to
    
  dir.create(file.path(resultDir),showWarnings = FALSE)
   
  
  # comparison for every pair of samples
  
  rnds = length(pairs)/2
  for (k in 1:rnds){
    kpair = k + rnds
    res=compareAbundanceInTwoSamplesForShared(sam[,c(k,kpair)])
    colnames(res) <- c(paste(colnames(sam[,c(k,kpair)])[1],colnames(sam[,c(k,kpair)])[2],"padjusted",sep="_"),paste(colnames(sam[,c(k,kpair)])[1],colnames(sam[,c(k,kpair)])[2],"OR",sep="_"),paste(colnames(sam[,c(k,kpair)])[1],colnames(sam[,c(k,kpair)])[2],"CI",sep="_"))
    
    print(str(res[,1]))
    
    
    siglist = as.numeric(res[,1]) < pcutoff
    
    resSig = res[siglist,]
    
    print(head(siglist))
    
    
    resultlist[[colnames(sam[,c(k,kpair)])[1]]] <- cbind(sam[rownames(resSig),c(k,kpair)],resSig)
    
    if(writeResult==T){
    print("Result written to file....")
    write.table(res,file=paste(resultDir,"/",colnames(sam[,c(k,kpair)])[1],"_",colnames(sam[,c(k,kpair)])[2],"_Abundance_fishers_comparison_result_forPairSamples.txt",sep=""),row.names = T,sep = "\t")
    }
   
    
      
  
  }
  
  return(resultlist)
  
  
}




#### composition based comparison of repertoires using 4-mer nt usage frequencies ####

getKmerNtFrequencies <- function(reps,kMer=4,asProb=F,vGenes=NULL,seqStatus="prod"){
  
  kmerFreqTable = NULL
  
  for(rep in reps){
    print(rep)
    p <- get(rep)
    
    if(seqStatus=="prod"){
      p <- getProductive(p)
    }else {
      p <- getUnproductive(p)
    }
    
  
    if(!is.null(vGenes)){    
      p <- p[p$VGENENAME %in% vGenes, ]
    }  
    
    print(nrow(p))
    #print(p[,c(8,22)])
    
    seqs <- unique(p$NUCLEOTIDE) # select unique AA clonotypes
    seqSet <- DNAStringSet(seqs)
    seq_mers <- oligonucleotideFrequency(seqSet,width=kMer,as.prob=asProb,simplify.as="collapsed")
    kmerFreqTable = rbind(kmerFreqTable,seq_mers)
    
  }
  
  
  rownames(kmerFreqTable) <- reps
  
  kmerFreqTable <- t(kmerFreqTable)
  
  print("Writing result to file....")
  
  write.table(kmerFreqTable,file="4merFrequencyTable.txt",row.names = T,col.names=T,sep = "\t")
  
  return(kmerFreqTable)
  
}

getKmerNtFrequenciesVJ <- function(reps,kMer=4,asProb=F,vjGenes=NULL,seqStatus="prod"){
  
  kmerFreqTable = NULL
  
  for(rep in reps){
    print(rep)
    p <- get(rep)
    
    if(seqStatus=="prod"){
      p <- getProductive(p)
    }else {
      p <- getUnproductive(p)
    }
    
    tempd = NULL
    if(!is.null(vjGenes)){
      
      res=strsplit(vjGenes," ")
      for(i in 1:length(res)){
        
        tempd <- rbind(tempd,p[(p$VGENENAME==res[[i]][1]) & (p$JGENENAME==res[[i]][2]),])
      }
    
      p <- tempd
      
    }  
   
  print(nrow(p))
  #print(p[,c(8,22)])
    
  seqs <- unique(p$NUCLEOTIDE) # select unique AA clonotypes
  seqSet <- DNAStringSet(seqs)
  seq_mers <- oligonucleotideFrequency(seqSet,width=kMer,as.prob=asProb,simplify.as="collapsed")
  kmerFreqTable = rbind(kmerFreqTable,seq_mers)
  
  }
  
 
  rownames(kmerFreqTable) <- reps
  
  kmerFreqTable <- t(kmerFreqTable)
  
  print("Writing result to file....")
  
  write.table(kmerFreqTable,file="4merFrequencyTable.txt",row.names = T,col.names=T,sep = "\t")
  
  return(kmerFreqTable)
  
}

getClonalKmerNtFrequencies <- function(reps,selectionType="random",nClones=100,kMer=4,asProb=F,vGenes=NULL,seqStatus="prod"){
  
  kmerFreqTable = NULL
  repnames = c()
  
  for(rep in reps){
    print(rep)
    p <- get(rep)
    
    if(seqStatus=="prod"){
      p <- getProductive(p)
    }else {
      p <- getUnproductive(p)
    }
    
    #select clones with given list of Vgenes
    
    if(!is.null(vGenes)){ 
      p <- p[p$VGENENAME %in% vGenes,]
    }  
    
    # then select 100 clonotypes
    
    if(selectionType=="random"){ 
    sampling.index = sample(1:nrow(p),nClones,replace=F)  
    }else{
      # just pick the top nClones
      sampling.index = 1:nClones  
    }
    
    p <- p[sampling.index,]
    
   
    
    print(nrow(p))
    
    seqs <- unique(p$NUCLEOTIDE) # select unique AA clonotypes
    seqSet <- DNAStringSet(seqs)
    seq_mers <- oligonucleotideFrequency(seqSet,width=kMer,as.prob=asProb)
    kmerFreqTable = rbind(kmerFreqTable,seq_mers)
    
    repnames = c(repnames,rep(rep,nClones))
  }
  
  rownames(kmerFreqTable) <- repnames  
  kmerFreqTable <- t(kmerFreqTable)
  
  print("Writing result to file....")
  
  write.table(kmerFreqTable,file="clonal4merFrequencyTable.txt",row.names = T,col.names=T,sep = "\t")
  
  return(list(freqTable = kmerFreqTable,sams = repnames))
  
}




