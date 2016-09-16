\name{gds2mlumi}
\alias{gds2mlumi}
\alias{gds2mset}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert Genomic Data Structure file to Methylumiset or Methylset object.
}
\description{
Convert a Genomic Data Structure object back into a methylumi object, with subsetting features.
}
\usage{
gds2mlumi(gds, i, j)
gds2mset(gds, i, j, anno)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gds}{
a gds object
}
  \item{i}{
Index of rows
}
  \item{j}{
Index of Columns
}
  \item{anno}{
If NULL, function will attempt to guess the annotation to be used. Otherwise can be
specified with either "27k", "450k", "epic" or "unknown".
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A methylumi object
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
data(melon)
e <- es2gds(melon, "melon.gds")
gds2mlumi(e)
closefn.gds(e)
unlink("melon.gds")
  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
