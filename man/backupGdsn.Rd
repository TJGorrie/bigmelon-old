%backup.gdsn <- function(gds=NULL, node){
%  # Quick function to quickly copy specific nodes into a new-folder 
%  # called "backup". Really for personal use - not recommended for standard workflow.
%  if(is.null(gds)) gds <- getfolder.gdsn(node)
%  if(!("backup"%in%ls.gdsn(gds))) addfolder.gdsn(gds, "backup")
%  copyto.gdsn(index.gdsn(gds, "backup"), index.gdsn(gds, node))
%}

\name{backup.gdsn}
\alias{backup.gdsn}

\title{
Copy gds node to a backup folder within gds object
}

\description{
Quick function that will copy designated gdsn.class node within a gds object to a 'backup'
folder. If 'backup' folder does not exist, this is created. 
}

\usage{
backup.gdsn(gds = NULL, node)
}

\arguments{
\item{gds}{
If NULL, function will call \code{\link[gdsfmt]{getfolder.gdsn}} to find the root node. 
Otherwise, user can specify a separate gds.class object to copy the specified node to.
}

\item{node}{
gdsn.class object (a gds node)
}
}

\author{
Tyler Gorrie-Stone <tgorri@essex.ac.uk>
}

\seealso{\code{\link[gdsfmt]{copyto.gdsn}}}

\examples{
data(melon)
e <- es2gds(melon, "melon.gds")
nod <- index.gdsn(e, "betas")
backup.gdsn(gds = NULL, node = nod)

closefn.gds(e)
unlink("melon.gds")
}

\keyword{backup}
