\name{optBlockC}
\alias{optBlockC}
\title{Optimal design blocking}
\description{INTERNAL: Simplified wrapper for blocking of experimental designs using \code{optBlock} in the 
  \emph{AlgDesign} package.
}

\usage{
optBlockC(withinData, blocksizes, nRepeats=5)
}

\arguments{
  \item{withinData}{A matrix or data frame describing the variables. Data types allowed include: factors, 
  effects-coded designs, or dummy-coded designs. Other data types will lead to errors or unbalanced designs. 
  If the columns are not named, they will be assumed to have the names X1, X2, etc. The number of rows in 
  \code{withinData} must be at least as large as the sum of the number of terms plus the number of blocks.}
  \item{blocksizes}{A vector giving the block sizes for each block. The length of \code{blocksizes} specifies 
  the number of blocks.}
  \item{nRepeats}{The number of times the entire process is repeated.}
}

\details{
  Simplified wrapper for \code{optBlock} that optimizes blocks on a pre-existing design or a set of factors using 
  the \code{D} criterion. Does not permit whole plot factors to interact with within plot factors. See 
  \code{optBlock} for additional details.
} 

\value{
  \item{D}{\emph{det(M)^(1/k)}, where \emph{det(M)} is the determinant of the normalized dispersion matrix 
  \emph{M}, or \emph{m=X'X/N}, where each row of \emph{X} has the appropriate block mean subtracted.}
  \item{diagonality}{\emph{|M|/P^(1/k)}, where \emph{P} is the product of the diagonal elements of 
  \emph{M}.}
  \item{Blocks}{A list of blocks, labeled B1, B2, etc.}
  \item{design}{A data frame. The design with blocks in stacked order.}
  \item{rows}{Numeric row numers of the design rows corresponding to \code{withinData} rows.}
}

\references{
  Wheeler, R.E. (2004). optBlock. \emph{AlgDesign.} The R project for statistical computing.
    (http://www.r-project.org).
}


\examples{
##INTERNAL USE ONLY
}