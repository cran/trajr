## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(trajr)

## ---- echo = FALSE------------------------------------------------------------
# NOTE that I have to do some confusing chunk duplication so that I can use functions before printing their definitions

## ----defFns, echo = FALSE-----------------------------------------------------
# Some colours
DARKBLUE <- "#102050"
MIDRED <- "#f82010"
MIDBLUE <- "#2040b8"
LIGHTBLUE <- "#60a0ff"

# Return the coordinates of a point in a trajectory as vector. By default, it's the end point.
trjPt <- function(trj, rowIdx = nrow(trj)) {
  as.numeric(trj[rowIdx, c(1, 2)])
}

# Rotate pt around origin by angle
rotatePt <- function(origin, pt, angle) {
  rm <- matrix(c(cos(angle), sin(angle), -sin(angle), cos(angle)), ncol = 2)
  npt <- as.numeric(t(pt) - origin)
  rm %*% npt + origin
}

# Generate a dataframe of points that lie along an arc
arcPts <- function(radius = 1, startAngle = 0, endAngle = 2 * pi, cx = 0, cy = 0, numPts = 10) {
  angles <- seq(startAngle, endAngle, length.out = numPts)
  data.frame(x = cx + radius * cos(angles), y = cy + radius * sin(angles))
}


## ----rotate, echo = FALSE-----------------------------------------------------
# Rotates trj around origin so that its starting point lies inside the boundary.
# Uses a brute force algorithm to find the minimal rotation: simply tests lots
# of potential rotations.
#
# @param origin The trajectory is rotated around this point.
# @param trj The trajectory to be rotated.
# @param boundary The region to stay inside.
# @param inc Angular increment (in radians) of points to test. The first point
#   tested is inc, then -inc, then 2 * inc, 2 * -inc and so on.
RotateToDeflect <- function(origin, trj, boundary, inc = pi / 90) {
  pt2 <- trjPt(trj, 1) # Starting point of trj

  # Find a rotation such that pt2 is inside the boundary
  angle <- findRotation(origin, pt2, inc)

  # Now rotate the whole trajectory around the origin point
  TrajRotate(trj, angle, origin = origin, relative = FALSE)
}

# This is the algorithm to find the minimal rotation angle. Simply generates
# lots of angles, then tests them until a suitable angle is found
findRotation <- function(pt1, pt2, inc) {
  for (alpha in seq(inc, pi, by = inc)) {
    # Rotate pt2 around pt1 by +- alpha
    npt2 <- rotatePt(pt1, pt2, alpha)
    # point.in.polygon returns 1 if the point is inside
    if (sp::point.in.polygon(npt2[1], npt2[2], boundary$x, boundary$y) == 1)
      return(alpha)
    npt2 <- rotatePt(pt1, pt2, -alpha)
    if (sp::point.in.polygon(npt2[1], npt2[2], boundary$x, boundary$y) == 1)
      return(-alpha)
  }
  stop("Cannot find suitable rotation")
}


## ----constrain, echo = FALSE--------------------------------------------------
# Constrains a trajectory to the inside of a boundary, using a simple model of
# behaviour which is: don't walk through walls.
ConstrainTrajectory <- function(trj, boundary) {
  # Start adjusting the trajectory so it doesn't cross any walls.
  # Find the first crossing, and split into 2 parts
  l <- TrajSplitAtFirstCrossing(trj, boundary)
  
  # Now l is a list containing 2 trajectories (which we earlier referred to as t1 & t2).
  # Build up a list of rotated parts as we go
  parts <- list(l[[1]])
  # When l has length 1, the remainder of the trajectory lies inside the boundary
  while (length(l) == 2) {
    
    # Rotate the section t2 about the last point in the previous section
    t2 <- RotateToDeflect(trjPt(l[[1]]), l[[2]], boundary)
    
    # Work out where the trajectory now crosses the boundary
    l <- TrajSplitAtFirstCrossing(t2, boundary)
    
    # Save the rotated section that's inside
    parts[[length(parts) + 1]] <- l[[1]]
  }
  
  # Put all the parts back together into a single trajectory
  TrajMerge(parts)
}

## ----echo=FALSE, fig.cap="_Figure 1. Steps to constrain a trajectory following our rules_", fig.height=7----
par(mfrow = c(3, 2), mar = c(5, 4, 1, 2) + .1)
boundary <- data.frame(x = c(-10, 10, 10, -10), y = c(-10, -10, 10, 10))
set.seed(1)
trj <- TrajGenerate(n = 6, stepLength = 4, angularErrorSd = .4)
xlim <- range(c(trj$x, boundary$x))
ylim <- range(c(trj$y, boundary$y))

# Step 1
plot(trj, xlim = xlim, ylim = ylim, col = DARKBLUE)
polygon(boundary, border = "brown", lwd = 2)
adj <- -.2
line <- -1
mtext("A", line = line, adj = adj)

# Step 2
l <- TrajSplitAtFirstCrossing(trj, boundary)
plot(trj, xlim = xlim, ylim = ylim, col = DARKBLUE)
polygon(boundary, border = "brown", lwd = 2)
pI <- trjPt(l[[1]])
pO <- trjPt(l[[2]], 1)
points(pI[1], pI[2], pch = 4)
text(pI[1], pI[2], "i", pos = 1, font = 4)
points(pO[1], pO[2], pch = 4)
text(pO[1], pO[2], "o", pos = 1, font = 4)
mtext("B", line = line, adj = adj)

# Step 3
plot(l[[1]], xlim = xlim, ylim = ylim, col = MIDRED, lwd = 2)
lines(l[[2]], col = MIDBLUE)
polygon(boundary, border = "brown", lwd = 2, lwd = 2)
text(4, -.8, "t1", pos = 1, font = 4)
text(18, -2.5, "t2", pos = 1, font = 4)
mtext("C", line = line, adj = adj)

# Step 4
angle <- findRotation(pI, pO, pi / 10)

# Step 5
plot(l[[1]], xlim = xlim, ylim = ylim, col = MIDRED, lwd = 2)
polygon(boundary, border = "brown", lwd = 2, lwd = 2)
x <- sapply(seq(0, angle, length.out = 5), function(alpha) lines(TrajRotate(l[[2]], alpha, origin = pI, relative = FALSE), col = LIGHTBLUE, start.pt.cex = .6))
mtext("D", line = line, adj = adj)

# Step 6
plot(l[[1]], xlim = range(boundary$x), ylim = range(boundary$y), col = MIDRED, lwd = 2)
polygon(boundary, border = "brown", lwd = 2, lwd = 2)
parts <- list(l[[1]])
while (length(l) == 2) {
  t2 <- RotateToDeflect(trjPt(l[[1]]), l[[2]], boundary)
  l <- TrajSplitAtFirstCrossing(t2, boundary)
  parts[[length(parts) + 1]] <- l[[1]]
  lines(l[[1]], start.pt.cex = 1.2)
}
mtext("E", line = line, adj = adj)

# Step 7
plot(NULL, xlim = range(boundary$x), ylim = range(boundary$y), asp = 1)
polygon(boundary, border = "brown", lwd = 2, lwd = 2)
lines(TrajMerge(parts), col = DARKBLUE, lwd = 2)
mtext("F", line = line, adj = adj)


## ----defFns, eval = FALSE-----------------------------------------------------
#  # Some colours
#  DARKBLUE <- "#102050"
#  MIDRED <- "#f82010"
#  MIDBLUE <- "#2040b8"
#  LIGHTBLUE <- "#60a0ff"
#  
#  # Return the coordinates of a point in a trajectory as vector. By default, it's the end point.
#  trjPt <- function(trj, rowIdx = nrow(trj)) {
#    as.numeric(trj[rowIdx, c(1, 2)])
#  }
#  
#  # Rotate pt around origin by angle
#  rotatePt <- function(origin, pt, angle) {
#    rm <- matrix(c(cos(angle), sin(angle), -sin(angle), cos(angle)), ncol = 2)
#    npt <- as.numeric(t(pt) - origin)
#    rm %*% npt + origin
#  }
#  
#  # Generate a dataframe of points that lie along an arc
#  arcPts <- function(radius = 1, startAngle = 0, endAngle = 2 * pi, cx = 0, cy = 0, numPts = 10) {
#    angles <- seq(startAngle, endAngle, length.out = numPts)
#    data.frame(x = cx + radius * cos(angles), y = cy + radius * sin(angles))
#  }
#  

## ----rotate, eval = FALSE-----------------------------------------------------
#  # Rotates trj around origin so that its starting point lies inside the boundary.
#  # Uses a brute force algorithm to find the minimal rotation: simply tests lots
#  # of potential rotations.
#  #
#  # @param origin The trajectory is rotated around this point.
#  # @param trj The trajectory to be rotated.
#  # @param boundary The region to stay inside.
#  # @param inc Angular increment (in radians) of points to test. The first point
#  #   tested is inc, then -inc, then 2 * inc, 2 * -inc and so on.
#  RotateToDeflect <- function(origin, trj, boundary, inc = pi / 90) {
#    pt2 <- trjPt(trj, 1) # Starting point of trj
#  
#    # Find a rotation such that pt2 is inside the boundary
#    angle <- findRotation(origin, pt2, inc)
#  
#    # Now rotate the whole trajectory around the origin point
#    TrajRotate(trj, angle, origin = origin, relative = FALSE)
#  }
#  
#  # This is the algorithm to find the minimal rotation angle. Simply generates
#  # lots of angles, then tests them until a suitable angle is found
#  findRotation <- function(pt1, pt2, inc) {
#    for (alpha in seq(inc, pi, by = inc)) {
#      # Rotate pt2 around pt1 by +- alpha
#      npt2 <- rotatePt(pt1, pt2, alpha)
#      # point.in.polygon returns 1 if the point is inside
#      if (sp::point.in.polygon(npt2[1], npt2[2], boundary$x, boundary$y) == 1)
#        return(alpha)
#      npt2 <- rotatePt(pt1, pt2, -alpha)
#      if (sp::point.in.polygon(npt2[1], npt2[2], boundary$x, boundary$y) == 1)
#        return(-alpha)
#    }
#    stop("Cannot find suitable rotation")
#  }
#  

## ----constrain, eval = FALSE--------------------------------------------------
#  # Constrains a trajectory to the inside of a boundary, using a simple model of
#  # behaviour which is: don't walk through walls.
#  ConstrainTrajectory <- function(trj, boundary) {
#    # Start adjusting the trajectory so it doesn't cross any walls.
#    # Find the first crossing, and split into 2 parts
#    l <- TrajSplitAtFirstCrossing(trj, boundary)
#  
#    # Now l is a list containing 2 trajectories (which we earlier referred to as t1 & t2).
#    # Build up a list of rotated parts as we go
#    parts <- list(l[[1]])
#    # When l has length 1, the remainder of the trajectory lies inside the boundary
#    while (length(l) == 2) {
#  
#      # Rotate the section t2 about the last point in the previous section
#      t2 <- RotateToDeflect(trjPt(l[[1]]), l[[2]], boundary)
#  
#      # Work out where the trajectory now crosses the boundary
#      l <- TrajSplitAtFirstCrossing(t2, boundary)
#  
#      # Save the rotated section that's inside
#      parts[[length(parts) + 1]] <- l[[1]]
#    }
#  
#    # Put all the parts back together into a single trajectory
#    TrajMerge(parts)
#  }

## ---- fig.cap="_Figure 2. Trajectory constrained to a circular arena_"--------
# Circular arena
boundary <- arcPts(100, 0, 2 * pi, numPts = 60)

# Create a random trajectory
set.seed(1)
trj <- TrajGenerate(n = 5000, angularErrorSd = .14, spatialUnits = "mm")

# Constrain it to the arena
constrained <- ConstrainTrajectory(trj, boundary)

plot(constrained, col = DARKBLUE)
polygon(boundary, border = "brown", lwd = 2)


## ---- fig.cap="_Figure 3. Trajectory constrained to an hourglass arena_"------
# Build an hourglass-shaped arena similar to Creed & Miller, (1990)
hourglassArena <- function() {
  # Define lower-left quadrant shape
  c1 <- arcPts(30, pi, 1.5 * pi, -70, -30)
  c2 <- arcPts(30, 1.5 * pi, 1.9 * pi, -49, -30)
  c3 <- arcPts(20, .9 * pi, pi / 2, 0, -40)
  xs <- c(c1$x, c2$x, c3$x)
  ys <- c(c1$y, c2$y, c3$y)
  # Exploit the 2 axes of symmetry
  data.frame(x = c(xs, -rev(xs), -xs, rev(xs)),
             y = c(ys, rev(ys), -(ys), -rev(ys)))
}

boundary <- hourglassArena()

# Create a random trajectory
set.seed(2)
trj <- TrajGenerate(n = 10000, stepLength = 2, angularErrorSd = .1, spatialUnits = "mm")

# Constrain it to the arena
constrained <- ConstrainTrajectory(trj, boundary)

plot(constrained, col = DARKBLUE)
polygon(boundary, border = "brown", lwd = 2)


## ---- fig.cap="_Figure 4. Heatmap of trajectory constrained to an hourglass arena_"----
d <- MASS::kde2d(constrained$x, constrained$y, n = 300)
par(mar = c(3, 2, 1, 1) + .1)
image(d, asp = 1)
polygon(boundary, border = "blue", lwd = 2)


