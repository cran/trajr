library(trajr)

context("trajectory creation")

trjFromAnglesAndLengths <- function(angles, lengths) {
  coords <- c(0, cumsum(complex(modulus = lengths, argument = angles)))
  TrajFromCoords(data.frame(Re(coords), Im(coords)))
}

test_that("Trajectory creation", {
  csvFile <- "../testdata/096xypts.csv"
  expect_true(file.exists(csvFile))
  coords <- utils::read.csv(csvFile, stringsAsFactors = FALSE)
  expect_false(is.null(coords))
  trj <- TrajFromCoords(coords, fps = 1000)

  expect_false(is.null(trj))
  expect_equal(2030, nrow(trj))
  xRange <- c(997.31, 1541.549436)
  expect_equal(range(trj$x), xRange)
  yRange <- c(669.883810, 956.924828)
  expect_equal(range(trj$y), yRange)

  # Scaling
  scale <- 1 / 2500
  scaled <- TrajScale(trj, scale, "m")
  #plot(scaled)
  expect_false(is.null(scaled))
  expect_equal(nrow(trj), nrow(scaled))
  expect_equal(range(scaled$x), xRange * scale)
  expect_equal(range(scaled$y), yRange * scale)

  # Smoothing
  smoothed <- TrajSmoothSG(scaled, 3, 101)
  #plot(smoothed)
  expect_true(TrajLength(smoothed) < TrajLength(scaled))
  expect_true(abs(TrajDistance(smoothed) - TrajDistance(scaled)) < TrajDistance(scaled) / 10)

  # Derivatives
  derivs <- TrajDerivatives(smoothed)
  #plot(derivs$speed, type = 'l', col = 'red')
  #plot(derivs$acceleration, type = 'l')

  # Rediscretization
  rd <- TrajRediscretize(smoothed, .001)
  #plot(rd)

  expect_true(TrajStraightness(smoothed) < 1)
  expect_true(TrajStraightness(smoothed) > 0)

  corr <- TrajDirectionAutocorrelations(rd)
  # plot(corr, type='l')
  mn <- TrajDAFindFirstMinimum(corr, 10)
  # points(mn["deltaS"], mn["C"], pch = 16, col = "red", lwd = 2)
  # points(mn["deltaS"], mn["C"], col = "black", lwd = 2)
  mx <- TrajDAFindFirstMaximum(corr, 5)
  # points(mx["deltaS"], mx["C"], pch = 16, col = "green", lwd = 2)
  # points(mx["deltaS"], mx["C"], col = "black", lwd = 2)

})

test_that("Speed intervals", {

  # 1 Interval with no start and 1 stop
  set.seed(1)
  trj <- TrajGenerate(200, random = TRUE)
  slowerThan = NULL
  fasterThan = 120
  smoothed <- TrajSmoothSG(trj, 3, 101)
  intervals <- TrajSpeedIntervals(smoothed, slowerThan = slowerThan, fasterThan = fasterThan)
  #plot(intervals)
  expect_true(nrow(intervals) == 1)

  # 1 Interval with 1 start and no stop
  set.seed(2)
  trj <- TrajGenerate(200, random = TRUE)
  slowerThan = NULL
  fasterThan = 120
  smoothed <- TrajSmoothSG(trj, 3, 101)
  intervals <- TrajSpeedIntervals(smoothed, slowerThan = slowerThan, fasterThan = fasterThan)
  #plot(intervals)
  expect_true(nrow(intervals) == 1)

  # 0 intervals
  set.seed(3)
  trj <- TrajGenerate(200, random = TRUE)
  slowerThan = NULL
  fasterThan = 200
  smoothed <- TrajSmoothSG(trj, 3, 101)
  intervals <- TrajSpeedIntervals(smoothed, slowerThan = slowerThan, fasterThan = fasterThan)
  #plot(intervals)
  expect_true(nrow(intervals) == 0)

  # 3 intervals
  set.seed(4)
  trj <- TrajGenerate(200, random = TRUE)
  slowerThan = 150
  fasterThan = 90
  smoothed <- TrajSmoothSG(trj, 3, 101)
  intervals <- TrajSpeedIntervals(smoothed, slowerThan = slowerThan, fasterThan = fasterThan)
  #plot(intervals)
  expect_true(nrow(intervals) == 3)

  # 3 intervals
  set.seed(4)
  trj <- TrajGenerate(200, random = TRUE)
  slowerThan = 50
  fasterThan = NULL
  smoothed <- TrajSmoothSG(trj, 3, 101)
  intervals <- TrajSpeedIntervals(smoothed, slowerThan = slowerThan, fasterThan = fasterThan)
  #plot(intervals)
  expect_true(nrow(intervals) == 3)

  # 2 intervals
  set.seed(4)
  trj <- TrajGenerate(20, random = TRUE)
  slowerThan = 92
  fasterThan = NULL
  intervals <- TrajSpeedIntervals(trj, slowerThan = slowerThan, fasterThan = fasterThan)
  #plot(intervals)
  expect_true(nrow(intervals) == 2)

  # Interval wholly contained within a segment
  set.seed(4)
  trj <- TrajGenerate(10, random = TRUE)
  slowerThan = 110
  fasterThan = 107
  intervals <- TrajSpeedIntervals(trj, slowerThan = slowerThan, fasterThan = fasterThan)
  #plot(intervals)
  expect_true(nrow(intervals) == 0)

  set.seed(1)
  trj <- TrajGenerate(10, random = TRUE)
  slowerThan = NULL
  fasterThan = 110
  intervals <- TrajSpeedIntervals(trj, slowerThan = slowerThan, fasterThan = fasterThan)
  #plot(intervals)
  expect_true(nrow(intervals) == 2)

  slowerThan = 107
  fasterThan = NULL
  intervals <- TrajSpeedIntervals(trj, slowerThan = slowerThan, fasterThan = fasterThan)
  #plot(intervals)
  expect_true(nrow(intervals) == 3)
})

test_that("Emax", {
  set.seed(1)
  trj1 <- TrajGenerate(1000, angularErrorSd = .5)
  trj2 <- TrajGenerate(1000, angularErrorSd = .2)
  # trj2 (black) should be straighter than trj1 (red), hence Emax(trj1) < Emax(trj2)
  #plot(trj2, asp = NULL, xlim = range(c(trj1$x, trj2$x)), ylim = range(c(trj1$y, trj2$y)))
  #plot(trj1, col = "red", add = TRUE)

  expect_true(TrajEmax(trj1) < TrajEmax(trj2))
})

test_that("Sinuosity", {
  set.seed(1)
  trj1 <- TrajGenerate(1000, angularErrorSd = .5)
  trj2 <- TrajGenerate(1000, angularErrorSd = .2)
  # trj2 (black) should be straighter than trj1 (red), hence Sinuosity(trj1) > Sinuosity(trj2)
  #plot(trj2, asp = NULL, xlim = range(c(trj1$x, trj2$x)), ylim = range(c(trj1$y, trj2$y)))
  #plot(trj1, col = "red", add = TRUE)

  expect_true(TrajSinuosity(trj1) > TrajSinuosity(trj2))
})

test_that("Directional change", {

  # Test that directional change as implemented gives the same results as the equation in the book
  L <- c(1, 1, 1, 1)
  A <- c(pi / 4, 0, -pi / 8, pi / 6)
  trj <- trjFromAnglesAndLengths(A, L)
  #plot(trj, turning.angles = "random")

  .bookCalc <- function(trj) {
    # Lengths between consecutive points
    lengths1 <- Mod(diff(trj$polar))
    # Lengths between points 2 apart
    lengths2 <- Mod(trj$polar[3:nrow(trj)] - trj$polar[1:(nrow(trj) - 2)])
    # Times between points 2 apart
    times2 <- trj$displacementTime[3:nrow(trj)] - trj$displacementTime[1:(nrow(trj) - 2)]

    sapply(1:(nrow(trj) - 2), function(i) {
      a <- lengths1[i]
      b <- lengths1[i+1]
      c <- lengths2[i]
      t <- times2[i]
      (180 - (180 / pi * acos((a ^ 2 + b ^ 2 - c ^ 2) / (2 * a * b)))) / t
    })
  }
  expect_equal(TrajDirectionalChange(trj), .bookCalc(trj))

  set.seed(1)
  trj <- TrajGenerate()
  expect_equal(TrajDirectionalChange(trj), .bookCalc(trj))

  #microbenchmark(TrajDirectionalChange(trj), .bookCalc(trj), times = 1000)

})

test_that("Reverse", {
  set.seed(1)
  trj <- TrajGenerate()
  rv <- TrajReverse(trj)
  expect_equal(nrow(rv), nrow(trj))
  expect_equal(rv$polar[1], trj$polar[nrow(trj)])
  expect_equal(rv$polar[nrow(rv)], trj$polar[1])
  expect_equal(TrajLength(rv), TrajLength(trj))
  expect_equal(TrajEmax(rv), TrajEmax(trj))
})

test_that("Step lengths", {
  set.seed(1)
  nSteps <- 100
  nTrajs <- 4
  stepLength <- 1
  trjs <- lapply(1:nTrajs, TrajGenerate, n = nSteps, stepLength = stepLength)
  sl <- TrajsStepLengths(trjs)
  expect_equal(length(sl), nSteps * nTrajs)
  # Expect mean and median to be roughly equal to the specified step length
  expect_equal(mean(sl), stepLength, tolerance = 2e-2)
  expect_equal(median(sl), stepLength, tolerance = 2e-2)
})

test_that("Generate", {

  unifDist <- function(n) runif(n, -1, 1)

  set.seed(1)
  sd <- 0.5
  trj <- TrajGenerate(angularErrorSd = sd, linearErrorDist = unifDist)
  # Should NOT be able to reject the NULL hypothesis that turning angle errors are normally distributed
  expect_true(shapiro.test(TrajAngles(trj))$p.value > 0.05)
  expect_equal(sd(TrajAngles(trj)), sd, tolerance = 5e-2)
  # Should be able to reject the NULL hypothesis that linear errors are normally distributed
  expect_true(shapiro.test(TrajStepLengths(trj))$p.value <= 0.05)
  trj <- TrajGenerate(angularErrorDist = unifDist)
  # Should be able to reject the NULL hypothesis that turning angles are normally distributed
  expect_true(shapiro.test(TrajAngles(trj))$p.value <= 0.05)
  # Should NOT be able to reject the NULL hypothesis that linear errors are normally distributed
  expect_true(shapiro.test(TrajStepLengths(trj))$p.value > 0.05)
})

test_that("Smoothing", {
  set.seed(1)
  sd <- 0.5
  trj <- TrajGenerate(angularErrorSd = sd)
  smoothed <- TrajSmoothSG(trj, 3, 41)
  expect_true(TrajEmax(trj) < TrajEmax(smoothed))
  smoothed2 <- TrajSmoothSG(trj, 3, 101)
  expect_true(TrajEmax(smoothed) < TrajEmax(smoothed2))
})

test_that("Convenience", {
  # Reads a set of points from a file. The points come from multiple tracks
  # due to noise in the video conversion process.
  # The longest track is the one we are interested in
  #
  # Value - data frame with values x & y, and an attribute "numFrames" which records the number of frames in the source video
  .MreadPoints <- function(file, ...) {
    points <- read.csv(file, comment.char = '#')

    # Save the number of frames in the file in case the track doesn't extend until the end
    maxFrame <- max(points$Frame)

    # Save number of frames
    attr(points, 'numFrames') <- maxFrame

    points
  }

  tracks <- rbind(
    data.frame(file = "3527.csv", species = "Zodariid2 sp1", category = "spider"),
    data.frame(file = "3530.csv", species = "Daerlac nigricans", category = "mimic bug"),
    data.frame(file = "3534.csv", species = "Daerlac nigricans", category = "mimic bug"),
    data.frame(file = "3537.csv", species = "Myrmarachne erythrocephala", category = "mimic spider"),
    data.frame(file = NA, species = "", category = ""),
    data.frame(file = "3542.csv", species = "Polyrhachis sp1", category = "ant"),
    data.frame(file = "3543.csv", species = "Polyrhachis sp1", category = "ant"),
    data.frame(file = "3548.csv", species = "Crematogaster sp1", category = "ant"),
    data.frame(file = NA, species = "", category = ""),
    stringsAsFactors = FALSE
  )
  csvStruct <- list(x = "x", y = "y", time = "Time")

  # Expect to fail with a message when there are blank
  expect_error(TrajsBuild(tracks$file, scale = .220 / 720, spatialUnits = "m", timeUnits = "s", csvStruct = csvStruct, rootDir = "..", csvReadFn = .MreadPoints),
               "Trajectory input file name is blank or NULL.*")
  tracks <- na.omit(tracks)
  trjs <- TrajsBuild(tracks$file, scale = .220 / 720, spatialUnits = "m", timeUnits = "s", csvStruct = csvStruct, rootDir = "..", csvReadFn = .MreadPoints)

  expect_equal(length(trjs), nrow(tracks))
  expect_equal(TrajGetUnits(trjs[[2]]), "m")
  expect_equal(TrajGetTimeUnits(trjs[[2]]), "s")

  # Define a function which calculates some statistics
  # of interest for a single trajectory
  characteriseTrajectory <- function(trj) {
    # Measures of speed
    derivs <- TrajDerivatives(trj)
    mean_speed <- mean(derivs$speed)
    sd_speed <- sd(derivs$speed)

    # Measures of straightness
    sinuosity <- TrajSinuosity(trj)
    resampled <- TrajRediscretize(trj, .001)
    Emax <- TrajEmax(resampled)

    # Periodicity
    corr <- TrajDirectionAutocorrelations(resampled, round(nrow(resampled) / 4))
    first_min <- TrajDAFindFirstMinimum(corr)

    # Return a list with all of the statistics for this trajectory
    list(mean_speed = mean_speed,
         sd_speed = sd_speed,
         sinuosity = sinuosity,
         Emax = Emax,
         first_min_deltaS = first_min[1],
         first_min_C = first_min[2])
  }

  stats <- TrajsMergeStats(trjs, characteriseTrajectory)

})

test_that("Convenience-multi", {

  # Test building multiple trajectories from each "file"
  readFn <- function(filename, ...) {
    # Return 2 very different trajectories
    t1 <- TrajGenerate(50)
    t2 <- TrajGenerate(20, random = FALSE)
    list(t1 = t1[, c('x', 'y', 'time')], t2 = t2[, c('x', 'y', 'time')])
  }

  trjs <- TrajsBuild(c("one", "two"), csvReadFn = readFn, smoothN = 11)

  expect_equal(length(trjs), 4)
})

test_that("Sinuosity", {
  set.seed(1)
  for(aa in seq(0, 2, by = .1)) {
    trj <- TrajGenerate(angularErrorSd = aa)
    # Don't expect equal, just close
    expect_equal(TrajSinuosity(trj), TrajSinuosity2(trj), tolerance = 0.2)
  }
})