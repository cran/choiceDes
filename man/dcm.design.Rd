\name{dcm.design}
\alias{dcm.design}
\title{Optimal fractional factorial design}
\description{Generate an optimal fractional factorial design given
  vectors of factor lengths.
}

\usage{
dcm.design(cand, nb, sets, alts, fname=NULL, Rd=20, print=TRUE)
}

\arguments{
  \item{cand}{A vector of factor lengths, or a list containing vectors of factor lengths if a
    combinatorial design is needed.}
  \item{nb}{The number of blocks or versions in the final design.}
  \item{sets}{The number of choice sets in \emph{each version} of the final design.}
  \item{alts}{The number of alternatives in each choice set.}
  \item{fname}{A character string, usually ending in ".txt", indiciating the name of the file containing 
    the levels-coded design.}
  \item{Rd}{The number of repeats used by the initial design and blocking algorithms. See arg 
    \code{nRepeats} in \code{optFederov} and \code{optBlock} for additional details}
  \item{print}{Boolean indicating whether there is output to the console during execution.}
}

\details{
  This function generates balanced and blocked choice sets from one or more specified full factorial candidate 
  set(s) using a modified Federov (1972) algorithm. See \code{optFederov} in \emph{AlgDesign} (Wheeler, 2004) 
  for a more complete description of the algorithm. Starting points are chosen randomly (as opposed to by 
  nullification) and may be seeded using \code{set.seed}. The \code{D} criterion is used for optimization.
  
  See \code{optBlock} for a description of the blocking method used.
  
  If \code{fname} is specified in the call a tab-delimited plain-text file is generated in the working 
  directory containing the levels-coded design.
  
  Large problems will complete faster by setting \code{Rd} to a smaller value. However, this may come at the
  expense of a more efficient design.
} 

\value{
  \item{levels}{A data frame consisting of the levels-coded design with blocks stacked in order. Variables for
    card, version and task are appended.}
  \item{effects}{A list of the effects-coded, blocked design and diagnostics. See \code{optBlock} for
    additional details.}
  \item{d.eff}{A list containing \code{D} efficiency, the variance-covariance matrix, and parameter
    standard deviations from the effects-coded design. See \code{dcm.design.effcy} for additional details.}
}

\references{
  Federov, V.V. (1972). \emph{Theory of optimal experiments.} Academic Press, New York.
  
  Wheeler, R.E. (2004). \emph{AlgDesign.} The R project for statistical computing.
    (http://www.r-project.org).
}

\seealso{\code{optFederov}, \code{optBlock}}

\examples{
## Example 1:
## design from a single candidate set
levs1 <- c(3,3,5,4)
des <- dcm.design(levs1, 10, 6, 2)

## Example 2:
## combinatorial design from more than one candidate set
levs2 <- list(c(3,3), c(5,4))
des <- dcm.design(levs2, 10, 6, 2)
}