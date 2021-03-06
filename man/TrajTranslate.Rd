% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trajectory.R
\name{TrajTranslate}
\alias{TrajTranslate}
\title{Translate a trajectory}
\usage{
TrajTranslate(trj, dx, dy, dt = 0)
}
\arguments{
\item{trj}{The Trajectory to be translated.}

\item{dx}{Delta x.}

\item{dy}{Delta y.}

\item{dt}{Delta time.}
}
\value{
A new trajectory which is a translated version of the input trajectory.
}
\description{
Shifts an entire trajectory by the specified delta x and y.
}
\examples{
# Shift a trajectory so that its origin is (10, 15).
# Note that trajectories created by TrajGenerate always start at (0, 0)
set.seed(42)
trj <- TrajGenerate()
trj <- TrajTranslate(trj, 10, 15)

# Translate a trajectory so its origin (0, 0) and it starts at time 0
trj <- TrajTranslate(trj, -trj$x[1], -trj$y[1], -trj$time[1])

}
