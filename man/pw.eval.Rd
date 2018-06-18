\name{pw.eval}
\alias{pw.eval}
\title{INTERNAL: Evaluate two-way frequencies}
\description{Internal function to evaluate whether the two-way (pairwise) frequencies
  of items in a matrix are balanced.
}

\usage{
pw.eval(items, shown, drows, des)
}

\arguments{
  \item{items}{The total number of items shown in a tradeoff exercise.}
  \item{shown}{The number of items shown in each tradeoff task.}
  \item{drows}{The number of rows in the tradeoff design matrix.}
  \item{des}{A matrix consisting of the levels-coded tradeoff design.}
}

\details{
  Evaluates the two-way (pairwise) frequencies of items in a tradeoff design matrix
  and returns those frequencies as well as the off-diagonal mean and standard deviation of the
  frequencies.
  
  This function is called internally by \code{tradeoff.des}.
} 

\value{
  \item{tbl}{A matrix of the two-way frequencies (pairs of items) in the tradeoff design.}
  \item{od.mean}{The mean of the off-diagonal elements in \code{tbl}.}
  \item{od.stdv}{The standard deviation of the off-diagonal elements in \code{tbl}.}
}

\seealso{\code{tradeoff.des}}

\examples{
## INTERNAL USE ONLY
}