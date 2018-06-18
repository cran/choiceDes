\name{dcm.design.cand}
\alias{dcm.design.cand}
\title{Optimal fractional factorial design}
\description{Generate an optimal restricted fractional factorial design given
  a pre-generated candidate set.
}

\usage{
dcm.design.cand(cand, nb, sets, alts, fname=NULL, Rd=20, print=TRUE)
}

\arguments{
  \item{cand}{A data frame of columns representing factors in the design OR a tab-delimited file readable 
    using \code{read.table(}\emph{filename}\code{)}.\cr
	If \code{cand} is not a data frame, an external file is assumed.}
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
  Generates balanced and blocked choice sets from a \emph{pre-generated} candidate set. 
  Typical use will involve (1) generating a full factorial candidate set (see \code{gen.factorial}), 
  (2) manipulating levels as desired (e.g., adding restrictions) and,
  (3) using the manipulated set as input into the function.
  
  Design optimization and blocking use the same algorithms as those in \code{dcm.design}.

  If \code{fname} is not \code{NULL} a tab-delimited plain-text file is generated in the working 
  directory containing the levels-coded design.  
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

\seealso{\code{dcm.design}, \code{optFederov}, \code{optBlock}}

\examples{
## generate full factorial candidate set
cand <- gen.factorial(c(3,3,4), factors="all")

## restrict the candidate set so that level 3 in the first factor 
## cannot occur with level 1 in the second factor
remove.rows <- which(cand[,1] == 3 & cand[,2] == 1)
cand.restr <- cand[-remove.rows,]

## generate the design from the restricted candidate set
## and check that no design rows violate the restriction
des <- dcm.design.cand(cand.restr, 10, 6, 2)
which(des$levels[,4] == 3 & des$levels[,5] == 1)
}