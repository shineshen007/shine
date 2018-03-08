\name{StaAnalysis}
\alias{StaAnalysis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{StaAnalysis
%%  ~~function to do ... ~~
}
\description{a function can generate PCA,PLSDA,heatmap,s plot of foldchange and volcano plot, also can calculate vip value
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
StaAnalysis(data = NULL, sample.info = NULL,group =c("case","control"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{a dataframe include name,mz,rt and isotope columns,the rest of all are sample and QC columns.
%%     ~~Describe \code{data} here~~
}
  \item{sample.info}{a dataframe include sample.name,injection.order,class,batch and group columns.
%%     ~~Describe \code{sample.info} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
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
\author{Shen Shine
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
##---- Be sure the format of data and sample.info is correct!! ----
#load the demo data
data(data, package = "Shine")
data(sample.information, package = "Shine")

##create a folder for demo
dir.create("demo")
setwd("demo")

# export the demo data as csv
write.csv(data, "data.csv", row.names = FALSE)
write.csv(sample.info, "sample.info.csv", row.names = FALSE)

# Analysis process
StaAnalysis(data = "data.csv",sample.info = "sample.info.csv",group =c("S","M"))

