\name{get.file}
\alias{get.file}
\title{INTERNAL: Get data files}
\description{Internal function to get required design candidate, data, design and any other 
  files from either \code{R} sources or as tab-delimited files stored in the working directory.
}

\usage{
get.file(dat.in)
}

\arguments{
  \item{dat.in}{The input data to be read referencing either an \code{R} data frame or a
  string referencing a tab-delimited file stored in the working directory.}
}

\details{
  This internal function obtains data from different sources and formats them as data frames
  for use in other functions.
  
  If \code{dat.in} is already a data frame in \code{R}, the data frame is returned. 
  No other types of \code{R} objects are permitted. 
  
  If data are stored in a tab-delimited file in the working directory, \code{dat.in} 
  references the filename as a character string enclosed in quotes. 
  The tab-delimited file must be readable using 
  \code{read.table(}\emph{filename}\code{, header=TRUE, sep="\t")}.
} 

\value{
  The input data (\code{dat.in}) coded as a data frame.
}

\seealso{\code{dcm.design.cand}}

\examples{
## INTERNAL USE ONLY
}