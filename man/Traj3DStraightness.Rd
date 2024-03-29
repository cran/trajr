% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/3D.R
\name{Traj3DStraightness}
\alias{Traj3DStraightness}
\title{Straightness of a 3D Trajectory}
\usage{
Traj3DStraightness(trj3d)
}
\arguments{
\item{trj3d}{3-dimensional trajectory to calculate straightness of.}
}
\value{
The straightness index of \code{trj}, which is a value between 0
  (infinitely tortuous) to 1 (a straight line).
}
\description{
Calculates the straightness index of a 3-dimensional trajectory, \eqn{D / L},
where \code{D} is the beeline distance between the first and last points in
the trajectory,and \code{L} is the path length travelled (Batschelet, 1981).
Benhamou (2004) considers the straightness index to be a reliable measure of
the efficiency of a directed walk, but inapplicable to random trajectories.
The straightness index of a random walk tends towards zero as the number of
steps increases, hence should only be used to compare the tortuosity of
random walks consisting of a similar number of steps.
}
\details{
The straightness index is also known as the net-to-gross displacement ratio.
According to Batschelet (1981), this value (termed \emph{d}) is an
approximation of \emph{r}, which is the length of the mean vector of turning
angles of a constant step-length trajectory (see
\code{\link{TrajMeanVectorOfTurningAngles}} and
\code{\link{TrajRediscretize}} for creating a constant step-length
trajectory).
}
\references{
Batschelet, E. (1981). Circular statistics in biology. ACADEMIC PRESS, 111
FIFTH AVE., NEW YORK, NY 10003, 1981, 388.

Benhamou, S. (2004). How to reliably estimate the tortuosity of an animal's
path. Journal of Theoretical Biology, 229(2), 209-220.
doi:10.1016/j.jtbi.2004.03.016
}
\seealso{
\code{\link{Traj3DFromCoords}}, \code{\link{Traj3DDistance}} for
  trajectory distance (or displacement), and \code{\link{Traj3DLength}} for
  trajectory path length, \code{\link{Traj3DStraightness}} for the
  straightness of a 2D trajectory.
}
