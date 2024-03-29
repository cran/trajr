# Trajectory basic query functions

# Substitute for the R base function `Arg(z)` that treats the argument of zero
# as undefined. The base R function returns 0.
TrajArg <- function(z) {
  ifelse(Mod(z) == 0,
         NA,
         Arg(z))
}

# ---- Trajectory query ----

#' Trajectory frames-per-second
#'
#' Returns the frames-per-second recorded for this trajectory.
#'
#' @param trj Trajectory to query
#'
#' @export
TrajGetFPS <- function(trj) { attr(trj, .TRAJ_FPS) }

#' Trajectory number of coordinates
#'
#' Returns the number of coordinates recorded for this trajectory, i.e. 1 more
#' than the number of steps.
#'
#' @param trj Trajectory to query
#'
#' @export
TrajGetNCoords <- function(trj) { attr(trj, .TRAJ_NFRAMES) }

#' Trajectory spatial units
#'
#' Returns the spatial units specified for a scaled trajectory.
#'
#' @param trj Trajectory to query
#'
#' @seealso \code{\link{TrajScale}}, \code{\link{TrajGetTimeUnits}}.
#'
#' @export
TrajGetUnits <- function(trj) { attr(trj, .TRAJ_UNITS) }

#' Trajectory temporal units
#'
#' Returns the temporal units specified for a scaled trajectory.
#'
#' @param trj Trajectory to query
#'
#' @seealso \code{\link{TrajFromCoords}}, \code{\link{TrajGetUnits}}.
#'
#' @export
TrajGetTimeUnits <- function(trj) {
  attr(trj, .TRAJ_TIME_UNITS)
}

# ---- Trajectory analysis ----

#' Trajectory step lengths
#'
#' Returns the lengths of each step in a trajectory.
#'
#' @param trj Trajectory to query.
#'
#' @seealso \code{\link{TrajLength}}
#'
#' @export
TrajStepLengths <- function(trj) {
  Mod(utils::tail(trj$displacement, -1)) # First displacement is not a step (also usually 0)
}

#' Trajectory distance
#'
#' Calculates the distance between the start and end of a trajectory (or a
#' portion of a trajectory). Also called the diffusion distance, net distance,
#' displacement, or bee-line from start to finish.
#'
#' @param trj Trajectory whose distance is to be calculated.
#' @param startIndex Index of the starting point.
#' @param endIndex Index of the ending point.
#'
#' @return Numeric distance from the start to the end of the trajectory.
#'
#' @export
TrajDistance <- function(trj, startIndex = 1, endIndex = nrow(trj)) {
  Mod(diff(trj$polar[c(startIndex, endIndex)]))
}

#' Trajectory length
#'
#' Calculates the cumulative length of a trajectory (or a portion of a
#' trajectory), which is the total distance travelled along the trajectory.
#'
#' @param trj Trajectory whose length is to be calculated.
#' @param startIndex Index of the starting point.
#' @param endIndex Index of the ending point.
#'
#' @return Numeric length of the trajectory.
#'
#' @seealso \code{\link{TrajStepLengths}}
#'
#' @export
TrajLength <- function(trj, startIndex = 1, endIndex = nrow(trj)) {
  sum(Mod(diff(trj$polar[startIndex:endIndex])))
}

#' Trajectory duration
#'
#' Calculates the temporal duration of a trajectory (or a portion of a
#' trajectory).
#'
#' @param trj Trajectory whose duration is to be calculated.
#' @param startIndex Index of the starting point.
#' @param endIndex Index of the ending point.
#'
#' @return Numeric duration of the trajectory, in time units.
#'
#' @seealso \code{\link{TrajGetTimeUnits}}
#'
#' @export
TrajDuration <- function(trj, startIndex = 1, endIndex = nrow(trj)) {
  diff(trj$displacementTime[c(startIndex, endIndex)])
}

#' Trajectory mean velocity
#'
#' Calculates the mean or net velocity of a trajectory (or a portion of a
#' trajectory). This is the velocity from the start point to the end point,
#' ignoring the path that was taken.
#'
#' @param trj Trajectory whose duration is to be calculated.
#' @param startIndex Index of the starting point.
#' @param endIndex Index of the ending point.
#'
#' @return Numeric duration of the trajectory, in time units.
#'
#' @seealso \code{\link{TrajGetTimeUnits}}
#'
#' @export
TrajMeanVelocity <- function(trj, startIndex = 1, endIndex = nrow(trj)) {
  d <- (trj[endIndex, c("x", "y")] - trj[startIndex, c("x", "y")]) / TrajDuration(trj, startIndex, endIndex)
  complex(real = d[1], imaginary = d[2])
}

#' Turning angles of a Trajectory
#'
#' Calculates the step angles (in radians) of each segment, either relative to
#' the previous segment or relative to the specified compass direction.
#'
#' Note that since turning angles are circular quantities, i.e. 360° == 0°, it
#' is incorrect to treat them as linear quantities. In particular, do not
#' calculate arithmetic means or standard deviations of turning angles. See
#' Batschelet, (1981) for a detailed explanation and techniques for dealing with
#' circular quantities.
#'
#' The turning angle before and after every zero-length segment will be
#' \code{NA}, since the angle of a zero-length segment is undefined. This
#' behaviour began in \code{trajr} version 1.5.0 (or development version
#' 1.4.0.9000). Prior to this fix, the angle of a zero-length segment was
#' assumed to be 0, which led to incorrect turning angles being returned. One
#' approach to dealing with zero-length segments is to simply remove them from
#' the trajectory. See \code{\link{TrajFromTrjPoints}} for a means to achieve this.
#'
#' @param trj the trajectory whose angles are to be calculated.
#' @param lag Angles between every lag'th segment are calculated. Only applies
#'   to non-directed walks, i.e. \code{compass.direction} is \code{NULL}.
#' @param compass.direction If not \code{NULL}, step angles are calculated
#'   relative to this angle (in radians), otherwise they are calculated relative
#'   to the previous step angle.
#'
#' @return Step angles in radians, normalised so that \code{-pi < angle <= pi}.
#'   If \code{compass.direction} is \code{NULL} (the default), the returned
#'   vector will have length \code{nrow(trj) - 2}, i.e. one angle for every pair
#'   of adjacent segments. If \code{compass.direction} is not \code{NULL}, the
#'   returned vector will have length \code{nrow(trj) - 1}, i.e. one angle for
#'   every segment.
#'
#' @seealso \code{\link{TrajStepLengths}},
#'   \code{\link{TrajMeanVectorOfTurningAngles}},
#'   \code{\link{TrajFromTrjPoints}}
#'
#' @references
#'
#' Batschelet, E. (1981). Circular statistics in biology. ACADEMIC PRESS, 111
#' FIFTH AVE., NEW YORK, NY 10003, 1981, 388.
#'
#' @export
TrajAngles <- function(trj, lag = 1, compass.direction = NULL) {
  if (is.null(compass.direction)) {
    angles <- diff(TrajArg(trj$displacement[2:nrow(trj)]), lag)
  } else {
    angles <- TrajArg(trj$displacement[2:nrow(trj)]) - compass.direction
  }

  # Normalise so that -pi < angle <= pi
  ii <- which(angles <= -pi)
  angles[ii] <- angles[ii] + 2 * pi
  ii <- which(angles > pi)
  angles[ii] <- angles[ii] - 2 * pi

  angles
}

#' Trajectory expected square displacement
#'
#' Calculates the expected square displacement for a trajectory assuming it is a
#' correlated random walk, using the formula in Kareiva & Shigesada, (1983).
#'
#' Note that Cheung, Zhang, Stricker, and Srinivasan (2007) define an
#' alternative formulation for expected maximum displacement, Emax (see
#' \code{\link{TrajEmax}}).
#'
#' @param trj A Trajectory.
#' @param n Number of steps to calculate.
#' @param eqn1 If \code{TRUE}, calculate using equation 1, otherwise using
#'   equation 2. Equation 2 applies when the mean of turning angles is 0,
#'   i.e.turns are unbiased.
#' @param compass.direction If not \code{NULL}, step angles are calculated
#'   relative to this angle (in radians), otherwise they are calculated relative
#'   to the previous step angle.
#'
#' @examples
#' set.seed(1)
#' # A random walk
#' trj <- TrajGenerate(200)
#' smoothed <- TrajSmoothSG(trj)
#'
#' # Calculate actual squared displacement at all points along the trajectory
#' sd2 <- sapply(2:nrow(smoothed), function(n) TrajDistance(smoothed, 1, n) ^ 2)
#' # Calculate expected squared displacement
#' ed2_1 <- sapply(2:nrow(smoothed), function(n) TrajExpectedSquareDisplacement(smoothed, n, TRUE))
#' ed2_2 <- sapply(2:nrow(smoothed), function(n) TrajExpectedSquareDisplacement(smoothed, n, FALSE))
#'
#' # Plot expected against actual. According to Kareiva & Shigesada, (1983), if actual
#' # (approximately) matches expected, the trajectory is probably a correlated random walk
#' par(mar = c(5, 5, 0.1, 0.1) + .1)
#' plot(2:nrow(smoothed), sd2, type = 'l', pch = 16, cex = .2, lwd = 2,
#'      xlab = 'Number of consecutive moves',
#'      ylab = expression('Squared displacement, ' * R[n]^2))
#' lines(2:nrow(smoothed), ed2_1, col = "grey", lwd = 2)
#' lines(2:nrow(smoothed), ed2_2, col = "pink", lwd = 2)
#'
#' legend("bottomright",
#'        c(expression("Actual displacement"^2),
#'          expression("Expected displacement"^2 * " (eqn 1)"),
#'          expression("Expected displacement"^2 * " (eqn 2)")),
#'        col = c('black', 'grey', 'pink'), lwd = 2,
#'        inset = c(0.01, 0.02))
#'
#' @seealso \code{\link{TrajEmax}}
#'
#' @references
#'
#' Cheung, A., Zhang, S., Stricker, C., & Srinivasan, M. V. (2007). Animal
#' navigation: the difficulty of moving in a straight line. Biological
#' Cybernetics, 97(1), 47-61. doi:10.1007/s00422-007-0158-0
#'
#' Kareiva, P. M., & Shigesada, N. (1983). Analyzing insect movement as a
#' correlated random walk. Oecologia, 56(2), 234-238. doi:10.1007/bf00379695
#'
#' @export
TrajExpectedSquareDisplacement <- function(trj, n = nrow(trj), eqn1 = TRUE, compass.direction = NULL) {
  sl <- TrajStepLengths(trj)
  ta <- TrajAngles(trj, compass.direction = compass.direction)
  l <- mean(sl)
  l2 <- mean(sl ^ 2)
  c <- mean(cos(ta), na.rm = TRUE)
  s <- mean(sin(ta), na.rm = TRUE)
  s2 <- s^2

  if (eqn1) {
    # Eqn 1
    alpha <- atan2(s, c)
    gamma <- ((1 - c)^2 - s2) * cos((n + 1) * alpha) - 2 * s * (1 - c) * sin((n + 1) * alpha)
    esd <- n * l2 + 2 * l^2 * ((c - c^2 - s2) * n  - c) / ((1 - c)^2 + s2) +
      2 * l^2 * ((2 * s2 + (c + s2) ^ ((n + 1) / 2)) / ((1 - c)^2 + s2)^2) * gamma
    # There seems to be a bug in the expression - for very straight trajectories,
    # value is negative although it seems to have a reasonable absolute value
    abs(esd)
  } else {
    # Eqn 2
    n * l2 + 2 * l^2 * c / (1 - c) * (n - (1 - c^n) / (1 - c))
  }
}
