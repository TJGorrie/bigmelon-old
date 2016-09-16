###
 # pwod.gdsn
 # gdsn.class "method" for wateRmelon::pwod
 # Functions identically to wateRmelon::pwod although takes
 # a gdsn.class object as input. Writes to a temporary file locally
 # akin to how all of the bigmelon "memory-efficient" methods work. 
 #
### 
 # Args:
 # node: gdsn.class object typically containing betas.
 # mul : Number of interquartile ranges to determine outlying probes.
 #       Default = 4 to ensure only very obvious outliers are removed.
 #
###
 # Output:
 #  Overwrites specified node with identical matrix that has NAs where
 #  outlying probes were found. 
 #
###
 # Note: 
 # gds object must be opened with `readonly=F` for this to function correctly.
 #
###

pwod.gdsn <- function(node, mul = 4){ # {{{
  # create temporary file in working directory
  f <- createfn.gds("temp.gds", allow.duplicate = TRUE) 
  dim <- objdesp.gdsn(node)$dim

  # Create new node in "temp.gds"
  n.t <- add.gdsn(node = f, name = "pwod", storage = "float64", 
                  valdim = c(dim[1], 0), val = NULL, replace = TRUE)

  # pwod.R outputs to created node.
  apply.gdsn(node,
             margin = 1,
             FUN = function(x,y){
                   # pwod.R compute row by row
                   quan <- fivenum(x)
                   iqr <- quan[4] - quan[2]
                   bounds <- c(quan[4] + (iqr * y),
                   quan[2] - (iqr * y)) 
                   # Upper Bound is [1], Lower Bound is [2]
                   d <- x > bounds[1] | x < bounds[2]
                   x[d] <- NA 
                   x
                   },
             y = mul,
             as.is = "gdsnode",
             target.node = n.t)

  # Create new node that replaces original node.
#  h.t <- add.gdsn(getfolder.gdsn(node),
#                  name = objdesp.gdsn(node)$name,
#                  valdim = c(dim[1], 0),
#                  val = NULL,
#                  storage = "float64", # Link this to trait possibly
#                  replace = TRUE)
  # Append new array col by col.
#  for(i in 1:dim[2]){
#  val <- readex.gdsn(index.gdsn(f, "pwod"), sel = list(NULL,i))
#  append.gdsn(node = h.t, val = val)
#  }
  # Testing
  copyto.gdsn(node = getfolder.gdsn(node), source = index.gdsn(f, "pwod"), name = paste0("pwod_", objdesp.gdsn(node)$name))

  # Close + Delete temp file
  closefn.gds(f)
  unlink("temp.gds")

} # }}}
