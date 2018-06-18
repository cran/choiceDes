\name{write.tab}
\alias{write.tab}
\title{INTERNAL: Write a data frame as tab-delimited file}
\description{Internal function that acts as an alias to \code{write.table}, appending extra
  arguments.
}

\usage{
write.tab(x, f)
}

\arguments{
  \item{x}{A data frame object in \code{R}.}
  \item{f}{A character string, usually ending in "*.txt", indicating the name of the file to
    be generated.}
}

\details{
  Writes a data frame to the file indicated by \code{f}, using \code{write.table}
  and appending the following arguments: \code{row.names=FALSE}, \code{col.names=TRUE},
  \code{quote=FALSE}, and \code{sep="\t"}.
  
  Called internally by \code{dcm.design}, \code{dcm.design.cand}, and
  \code{tradeoff.des}.
} 

\value{
  Does not return any value.
}

\seealso{\code{dcm.design}, \code{dcm.design.cand}, \code{tradeoff.des}}

\examples{
## INTERNAL USE ONLY
}