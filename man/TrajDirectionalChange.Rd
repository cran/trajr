% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/straightness.R
\name{TrajDirectionalChange}
\alias{TrajDirectionalChange}
\title{Directional change (DC)}
\usage{
TrajDirectionalChange(trj, nFrames = 1)
}
\arguments{
\item{trj}{Track to calculate DC for.}

\item{nFrames}{Frame delta to process: if 1, every frame is processed, if 2,
every 2nd frame is processed, and so on. Default is 1.}
}
\value{
The directional change (DC) in degrees between every pair of
  consecutive segments in the trajectory, i.e. if \code{nFrames} is 1, the
  returned vector will have length \code{nrow(trj) - 2}.
}
\description{
Calculates the time variation of directional change (DC) of a trajectory
\emph{sensu} Kitamura & Imafuku (2015). Directional change is defined as the
angular change (in degrees) between two steps in the trajectory, divided by
the time difference between the two steps.
}
\details{
This function returns the DC for each pair of consecutive steps. Kitamura &
Imafuku (2015) used the mean and the standard deviation of DC for portions of
trajectories as index values of nonlinearity and irregularity respectively.
}
\examples{
set.seed(42)
trj <- TrajGenerate()
SD = mean(TrajDirectionalChange(trj))
SDDC = sd(TrajDirectionalChange(trj))

}
\references{
Kitamura, T., & Imafuku, M. (2015). Behavioural mimicry in flight
  path of Batesian intraspecific polymorphic butterfly Papilio polytes.
  Proceedings of the Royal Society B: Biological Sciences, 282(1809).
  doi:10.1098/rspb.2015.0483
}
