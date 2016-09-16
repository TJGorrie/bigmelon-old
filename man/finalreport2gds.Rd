\name{finalreport2gds}
\alias{finalreport2gds}
\title{
Read finalreport files and convert to genomic data structure files
}
\description{
Function to easily load Illumina methylation data into a genomic data structure (GDS) file.
}
\usage{
finalreport2gds(finalreport, gds, ...)
}
\arguments{
  \item{finalreport}{
A filename of the text file exported from GenomeStudio
}
  \item{gds}{
The filename for the gds file to be created
}
  \item{\dots}{
Additional arguments passed to \code{\link[methylumi]{methylumiR}}
}
}
\details{
Creates a .gds file.
}
\value{
A gds.class object
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Tyler Gorrie-Stone
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
# finalreport2gds("finalreport.txt", gds="finalreport.gds")
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
