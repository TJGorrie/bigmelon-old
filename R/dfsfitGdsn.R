###
 # Performs the functions from dfsfit.R and dfs2.R within 'wateRmelon'
 # used in many flavours of quantile normalization methods.
 #
###
 # Used by:
 # 'dasen.gds' | location: dasenGdsn.R 
 #
###
 # Arguments :
 # gds       : Gds object
 # targetnode: Node within gdsobject where intensities exist (gdsn.class)
 # newnode   : Name of new node to store output within object (character)
 # roco      : "Location" where barcodes are stored in gds object. String Manipulated if supplied
 #             Usually found in 'pData/position', Default = NULL
 # onetwo    : "Location" where Probe Design is stored OR vector containing Probe Design types.
 # 
###
 # NOTE: This step is considerably slower than the wateRmelon method.
 #
###

dfsfit.gdsn <- function(gds,
                        targetnode,
                        newnode,
                        roco,
                        onetwo
                        ){ # {{{


  # Converting supplied 'nodes' into environment
  if(length(onetwo) == 1)           onetwo <- read.gdsn(index.gdsn(gds, onetwo))
  if(class(onetwo) == 'gdsn.class') onetwo <- read.gdsn(onetwo)
  if(length(roco) == 1)             roco <- read.gdsn(index.gdsn(gds, roco))
  if(class(roco) == 'gdsn.class')   roco <- read.gdsn(roco)
  
  datnod <- targetnode

  dim <- objdesp.gdsn(datnod)$dim

  # Replicating: dfs2.R 'apply(mn, 2, dfs2, onetwo)'

  # Use apply.gdsn?
  mdf <- apply.gdsn(datnod,
             margin = 2,
             FUN = function(val, onetwo){
    one <- density(val[onetwo == 'I'], na.rm = TRUE, n = 2^15, from = 0, to = 5000)
    two <- density(val[onetwo == 'II'], na.rm = TRUE, n = 2^15, from = 0, to = 5000)
    one$x[which.max(one$y)] - two$x[which.max(two$y)] 
             },
             as.is = "double",
             onetwo = onetwo)

#  mdf <- list()
#  for(x in 1:500){                            
#    val <- readex.gdsn(datnod, sel = list(NULL, x)) 
#    one <- density(val[onetwo=='I'], na.rm=TRUE, n = 2^15, from = 0, to = 5000)
#    two <- density(val[onetwo=='II'],na.rm=TRUE, n = 2^15, from = 0, to = 5000)
#    mdf[x] <- one$x[which.max(one$y)] - two$x[which.max(two$y)]
#  }
#  mdf <- unlist(mdf)

  # Replicating: dfsfit.R 
  if (! is.null(roco) ) {
     scol <- as.numeric(substr(roco,6,6))
     srow <- as.numeric(substr(roco,3,3))
     fit  <- try(  lm(mdf ~ srow + scol ), silent=TRUE) 
      if (! inherits (fit, "try-error") ) { mdf   <- fit$fitted.values}
      else { message ('Sentrix position model failed, skipping') }
  }

  # Creating newnode:
  n.t <- add.gdsn(gds, newnode, storage = "float64", 
                  valdim = c(dim[1],0), val = NULL, replace=TRUE)
  
  for(x in 1:dim[2]){
    val <- readex.gdsn(datnod, sel = list(NULL, x)) 
    # otcor <- matrix(rep(mdf, sum(onetwo=='I')), byrow=T, nrow=sum(onetwo=='I'))
    # Not needed since working with single column at a time!
#    val[onetwo=='I'] <- val[onetwo=='I'] - rep(mdf[x], sum(onetwo=='I'))
    val[onetwo=='I'] <- val[onetwo=='I'] - rep(mdf[x], sum(onetwo=='I'))
    # Commit to New Node.
    append.gdsn(n.t, val)
  }

} # }}}

#TGS OK
