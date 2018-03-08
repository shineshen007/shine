\name{PathwayBarplot}
\alias{PathwayBarplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{PathwayBarplot
%%  ~~function to do ... ~~
}
\description{a function can generate Pathway Barplot
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
PathwayBarplot(data = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{a dataframe include pathway and pvalue columns.
%%     ~~Describe \code{data} here~~
}
\item{}{
%%     ~~Describe \code{data} here~~
}
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
##---- Be sure the format of data is correct!! ----
#load the demo data
data(pathway_results, package = "Shine")

##create a folder for demo
dir.create("demo")
setwd("demo")

# export the demo data as csv
write.csv(pathway_results, "data.csv", row.names = FALSE)


# Analysis process
PathwayBarplot(data=data)

