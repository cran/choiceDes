\name{tradeoff.des}
\alias{tradeoff.des}
\title{MaxDiff and other tradeoff designs}
\description{Generate a design to be used for MaxDiff and related tradeoff exercises.
}

\usage{
tradeoff.des(items, shown, vers, tasks, fname=NULL, Rd=20, Rc=NULL, print=TRUE)
}

\arguments{
  \item{items}{The total number of items in the tradeoff exercise.}
  \item{shown}{The number of items shown in each tradeoff task.}
  \item{vers}{The number of blocks or versions in the final design.}
  \item{tasks}{The number of tradeoff tasks in \emph{each version} of the final design.}
  \item{fname}{A character string, usually ending in ".txt", indicating the name of the file containing 
    the tradeoff design.}
  \item{Rd}{The number of iterations in the design and blocking processes.}
  \item{Rc}{The number of iterations in the item by column position optimization routine.}
  \item{print}{Boolean indicating whether there is output to the console during execution.}
}

\details{
  This function replicates the functionality of Sawtooth Software MaxDiff Designer for designing MaxDiff and
  related tradeoff tasks.
  
  A modified Federov (1972) algorithm is applied to a factor equal in length to the number of items to optimize 
  the BIB design at \code{vers x tasks} rows and \code{shown} columns.
  
  The optimized design is evaluated for one-way frequencies (equal representation of each item across all 
  versions and column positions). Designs are also optimized for two-way or pairwise balance across all tasks. 
  Column position balance is the more computationally intensive process. The number of iterations required for
  this step is determined by the \code{Rc} argument which defaults to the larger of 1,000 or 10 \code{x} the 
  number of rows in the design. Large design problems may require a larger number of iterations to achieve 
  optimal column position balance.
  
  Once an optimal design has been found, it is blocked into versions using \code{optBlock} to ensure equal
  representation of items \emph{within} each version. See Wheeler (2004) for a more complete description of the
  modified Federov and blocking algorithms used in optimizing these designs.
} 

\value{
  \item{design}{A matrix consisting of the optimized design and additional variables for card, version and task.}
  \item{balance}{Tables of one-way item frequencies, two-way (pairwise) item frequencies, and item frequencies by 
	column position. Means and standard deviations are calculated from all elements of the one-way and column
	position tables, and from the off-diagonal elements of the two-way (pairwise) table.}
  \item{Rc.crit}{The criterion that is minimized to achieve column position balance is output as a vector
	(\code{crit.vec}) along with the number of iterations executed since the last change in this criterion 
	(\code{crit.stable}). If \code{crit.stable} is large or is a large proportion of the total number of iterations 
	(\code{Rc}), the solution is stable in terms of column position balance. If \code{crit.stable} is small, the 
	solution is likely unstable and column position balance could be improved by increasing \code{Rc}. See also
	\code{cp.scree} which produces a line plot of \code{crit.vec} as a function of \code{Rc}. }
  \item{time.elapsed}{The time required for the function to execute.}
}

\references{
  Federov, V.V. (1972). \emph{Theory of optimal experiments.} Academic Press, New York.
  
  Wheeler, R.E. (2004). \emph{AlgDesign.} The R project for statistical computing.
    (http://www.r-project.org).
}


\examples{
## Example 1:
## typical MaxDiff design with 12 items
des <- tradeoff.des(12, 4, 10, 9)

## Example 2:
## typical paired comparisons design with 14 items
des <- tradeoff.des(14, 2, 6, 14)
}