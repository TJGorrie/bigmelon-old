\name{pwod.gdsn}
\alias{pwod.gdsn}
\alias{pwod, gds.class}
\alias{pwod, gdsn.class}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Probe-Wise Outlier Detection
}
\description{
'P'robe-'W'ise 'O'utlier 'D'etection via interquartile ranges
}
\usage{
pwod.gdsn(node, mul = 4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{node}{
gdsn.class object containing array to be filtered
}
  \item{mul}{
Number of interquartile ranges used to determine outlying
          probes. Default is 4 to ensure only very obvious outliers are
          removed.
}
}
\details{
Detects outlying probes across arrays in methylumi and minfi
     objects. Outliers are probable low MAF/SNP heterozygotes. 

}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
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

\seealso{
\code{\link[wateRmelon]{pwod}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
data(melon)
e <- es2gds(melon, "melon.gds")
pwod(e)
# Or can use pwod(betas(e))
closefn.gds(e)
unlink("melon.gds")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
