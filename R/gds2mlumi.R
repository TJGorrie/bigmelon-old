gds2mlumi <- function(gds, i, j){
 ###
  # Converts gds.class objects into methylumiset objects. 
  # Includes indexing features.
  # Args:
  #    x: gds object - made within bigmelon (ideally)
  #    i: row index, either numeric, logical or character
  #    j: column index, either numeric, logical or character
 ###
 ###
  # TODO: Error Handling
 ### 
  # Creating Assay Data n.b. dimensions must be equal. 
  history.submitted = as.character(Sys.time())
  x <- gds
  if("NBeads"%in%ls.gdsn(x)){
    aDat <- assayDataNew(betas        = x[i, j, node = "betas", name = T, drop = F],
                         pvals        = x[i, j, node = "pvals", name = T, drop = F],
                         NBeads       = x[i, j, node = "NBeads", name = T, drop = F],
                         methylated   = x[i, j, node = "methylated", name = T, drop = F],
                         unmethylated = x[i, j, node = "unmethylated", name = T, drop = F])
  } else {
    aDat <- assayDataNew(betas        = x[i, j, node = "betas", name = T, drop = F],
                         pvals        = x[i, j, node = "pvals", name = T, drop = F],
                         methylated   = x[i, j, node = "methylated", name = T, drop = F],
                         unmethylated = x[i, j, node = "unmethylated", name = T, drop = F])
  }

  # Creating empty MethylumiSet object.
  x.lumi = new("MethyLumiSet", assayData=aDat)

  # pData
  pdat <- pData(x)
  rownames(pdat) <- colnames(x)
  pData(x.lumi) <- pdat[j, , drop = F]
  
  # fData
  fdat <- fData(x)
  rownames(fdat) <- rownames(x)
  fData(x.lumi) <- fdat[i, , drop = F]

  # QC
  if(length(grep("QC", ls.gdsn(x), ignore.case = T))>1){ # Probably a better method
    qcm <- QCmethylated(x)
    qcu <- QCunmethylated(x)
    colnames(qcm) <- colnames(qcu) <- colnames(x)
    rownames(qcm) <- rownames(qcu) <- QCrownames(x)

    qc <- new("MethyLumiQC", assayData = assayDataNew(methylated   = qcm[ , j, drop = F],
                                                      unmethylated = qcu[ , j, drop = F])
          )
  
    x.lumi@QC <- qc
  }

  #  x.lumi@protocolData <- protocolData(NChannelSet)
  #  x.lumi@annotation <- annotation(NChannelSet)
  #  x.lumi@QC@annotation <- annotation(NChannelSet)
  # fvarLabels(x.lumi) <- possibleLabels[1:ncol(fdat)]
  # fvarMetadata(x.lumi)[,1] <- possibleMetadata[1:ncol(fdat)]
  # pval.detect(x.lumi) <- pval # default value

  history.finished <- as.character(Sys.time())
  history.command <- "Converted to methylumi with gds2mlumi (bigmelon)"
  x.lumi@history <- rbind(getHistory(x), 
                          data.frame(submitted = history.submitted, 
                                      finished = history.finished,    
                                       command = history.command))

  return(x.lumi)
}

gds2mset <- function(gds,i,j,anno=NULL){
# Gds to MethylSet method
 ###
  # Converts gds.class objects into methylumiset objects. 
  # Includes indexing features.
  # Args:
  #    x: gds object - made within bigmelon (ideally)
  #    i: row index, either numeric, logical or character
  #    j: column index, either numeric, logical or character
  #    anno: Annotation file to link to MethylSet
 ###
 ###
  # TODO:
  # > Error Handling
 ### 
  x <- gds
  if(!is.null(anno)){ 
    if(!anno%in%c("27k", "450k", "epic", "unknown")){
      stop("anno needs to be either: \'27k\', \'450k\', or \'epic\' or \'unknown\'")
    }
  }
  M <- x[i = i, j = j,   "methylated", name = T, drop = F]
  U <- x[i = i, j = j, "unmethylated", name = T, drop = F]
  pd <- pData(x)[j, , drop = F]
  rownames(pd) <- colnames(x)[j]
  pd <- annotatedDataFrameFrom(object = as.matrix(pd), byrow = T)

  if(!is.null(anno)){
    if(anno=="27k"){ anno <- c("IlluminaHumanMethylation27k", "ilmn12.hg19")
    } else if(anno=="450k"){ anno <- c("IlluminaHumanMethylation450k", "ilmn12.hg19")
    } else if(anno=="epic"){ anno <- c("IlluminaHumanMethylationEPIC", "ilmn12.hg19")
    } else if(anno=="unknown"){ anno <- c("Unknown", "Unknown")
    }
  }

  # Guess Array Type - will not get correct array if performed on subset.
  if(is.null(anno)){
    nr <- nrow(fData(x))
    # Will guess array type based on number of rows, will mess up on subsets!
    if(nr > 50000 & nr < 500000){ anno <- c("IlluminaHumanMethylation450k", "ilmn12.hg19")
    } else if(nr >= 500000){ anno <- c("IlluminaHumanMethylationEPIC", "ilmn12.hg19")
    } else if(nr <=50000){ anno <- c("IlluminaHumanMethylation27k", "ilmn12.hg19")
    }
  }
  names(anno) <- c("array", "annotation")
  out <- MethylSet(Meth = M, Unmeth = U, phenoData = pd, annotation = anno)

  # Finicky, fix later important.
  out@preprocessMethod <- c(rg.norm="Converted from gdsfmt to MethylSet (bigmelon)",
                            minfi = as.character(packageVersion("minfi")),
                            manifest = NA #packageVersion(getManifest(anno))
                            )
  out

}










