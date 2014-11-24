\name{dcm.design.effcy}
\alias{dcm.design.effcy}
\title{INTERNAL: Calculate design efficiencies}
\description{Internal function to calculate mathematical efficiencies of designs.
}

\usage{
dcm.design.effcy(des)
}

\arguments{
  \item{des}{An effects-coded design to be evaluated.}
}

\details{
  This function calculates overall \code{D}-efficiency, the variance-covariance matrix,
  and standard deviations for each parameter from an effects coded design.
  
  Called internally by \code{dcm.design} and \code{dcm.design.cand}.
} 

\value{
  \item{D}{Overall \code{D}-efficiency: \emph{det(M)}^\emph{(}-1\emph{/k)}, where \emph{det(M)} 
    is the determinant of the dispersion matrix \emph{X'X}, and \emph{k} is the number of 
	parameters.
  }
  \item{V}{Variance-covariance matrix derived from \emph{M}.}
  \item{s}{Parameter standard deviations: \emph{sqrt(diag(V))}.}
}

\seealso{\code{dcm.design}, \code{dcm.design.cand}}

\examples{
des <- dcm.design(c(3,3,4,3), 10, 8, 3)$effects$design
eff <- dcm.design.effcy(des)
}