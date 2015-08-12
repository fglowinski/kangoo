\name{analyzeInclIndex}
\alias{analyzeInclIndex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
analyzeInclIndex
}
\description{
This function takes the output of the ABI SDS progamm and calculates the alternative isoform fraction (AIF), which describes the proportion of the alternative isoform compared to the contitutive expression as defined by the conservative primer pair.

Combines the loadExp and calculateInclIndex function into one call.}
\usage{
analyzeInclIndex(plate.setup = "splicing", plate.format = "96-well")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{plate.setup}{
Specify the sample setup on the plate. Default is "splicing".
}
  \item{plate.format}{
Specify the plate format of the RT-PCR cycler. Default is "96-well".
}
}
\details{
This function outputs 3 plots:
1. plot of the raw CT values for quality control
2. plot of the AIF without p-value
3. plot of the AIF including the p-value for each gene

Also it outputs a table describing the results and showing the dAIF. This output is also appended to the general output file in the e-book.
In addition a small txt file is written giving the session info.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
