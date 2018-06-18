\name{optFederovC}
\alias{optFederovC}
\title{Optimal design}
\description{INTERNAL: Simplified wrapper for calculating exact algorithmic designs using Federov's exchange
 algorithm. Based on \code{optFederov} in the \emph{AlgDesign} package.
}

\usage{
optFederovC(modelData, nTrials, nRepeats=5)
}

\arguments{
  \item{modelData}{The candidate list. A matrix or data frame describing the variables. If a matrix is input 
    and the columns are not named, they will be assigned names X1, X2, etc. Permitted data types include 
	factors or levels- or effects-coded designs.}
  \item{nTrials}{The number of trials in the final design.}
  \item{nRepeats}{The number of times the whole process is repeated.}
}

\details{
  Generates exact algorithmic designs using Federov's exchange algorithm, and optimizing the \code{D} criterion. 
  See \code{optFederov} for algorithmic details. A vignette is also available by typing \code{vignette("AlgDesign").}
  
  Input data, i.e., \code{modelData}, must be of a form that \code{model.matrix(~., modelData} results in an 
  effects-coded design or candidate set.
} 

\value{
  \item{D}{The \emph{k}th root of the generalized variance: \emph{det(M)^(1/k)}, where \emph{det(M)} is the 
    determinant of the normalized dispersion matrix, or \emph{m=Z'Z/n}, where \emph{Z=X[rows,]}.}
  \item{A}{The average coefficient variance: \emph{trace(M')/k}, where \emph{M'} is the inverse of \emph{M}.}
  \item{Ge}{The minimax normalized variance over \emph{X}, expressed as an efficiency with respect to the 
    optimal approximate theory design. It is defined as \emph{k/max(d)}, where \emph{max(d)} is the 
	maximum normalized variance over \emph{X}, or the maximum of \emph{x'(M')x}, over all rows \emph{x'} 
	of \emph{X}.}
  \item{Dea}{A lower bound on \code{D} efficiency for approximate theory designs. It is equal to 
    \emph{exp(1-1/Ge)}.}
  \item{design}{The design.}
  \item{rows}{A numerical vector of the design row numbers from \code{modelData}.}
}

\references{
  Wheeler, R.E. (2004). optFederov. \emph{AlgDesign.} The R project for statistical computing.
    (http://www.r-project.org).
}


\examples{
##INTERNAL USE ONLY
}