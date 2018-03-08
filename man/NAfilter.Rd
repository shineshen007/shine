\name{NAfilter}
\alias{NAfilter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{NAfilter
%%  ~~function to do ... ~~
}
\description{a function filter NA in compound.name column
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
NAfilter(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{a dataframe must include compound.name column.
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
##---- Be sure the format of data is correct, compound.name column must be included in your data!! ----
#load the demo data
data(unique.data, package = "Shine")

##create a folder for demo
dir.create("demo")
setwd("demo")

# export the demo data as csv
write.csv(unique.data, "data.csv", row.names = FALSE)

# Analysis process
NAfilter(data)


