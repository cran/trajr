## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(trajr)

x = seq(0, 6, .5)
# For simplicity, these functions assume frequency and amplitude of oscillations in speed, i.e. it's not parameterised
disp <- function(times) {
  ft <- floor(times)
  f <- times - ft
  ft * .5 + ifelse(f > 0.5, f - 0.5, 0)
}
fadedCol <- "#d3d3c0"
addSpeed <- function() lines(seq(0, 6, .5), c(rep(c(0, 1), 6), 0), type = 's', col = fadedCol)
tFromT <- function(times) TrajFromCoords(data.frame(disp(times), rep(0, length(times)), times), timeCol = 3)
plotSpeedAndDisp <- function(times, cex = 1) {
  t <- tFromT(times)
  f <- 1 / mean(diff(times))
  
  plot(x, disp(x), type = "l", col = fadedCol, xlim = range(c(0, t$time)), xlab = "Time", ylab = "Displacement", main = sprintf("Sampled Displacement (f = %g Hz)", f))
  lines(t$time, cumsum(t$displacement))
  points(t$time, cumsum(t$displacement), pch = 16, cex = cex)
  
  plot(TrajSpeedIntervals(t, slowerThan = 0.1), xlim = range(t$time), ylim = c(0, 1), main = sprintf("Derived speed (f = %g Hz)", f))
  addSpeed()
  td <- TrajDerivatives(t)
  points(td$speedTimes + t$time[1], td$speed, pch = 16, cex = cex)
}


## ----echo=FALSE, fig.height=4, fig.width=4.5----------------------------------
plot(x, disp(x), type = "l", xlab = "Time", ylab = "Displacement", main = "Actual Displacement")
plot(seq(0, 6, .5), c(rep(c(0, 1), 6), 0), type = 's', xlab = "Time", ylab = "Speed", main = "Actual speed")

## ----echo=FALSE, fig.height=4, fig.width=4.5----------------------------------
plotSpeedAndDisp(seq(0, 6, 1 / 2))

## ----echo=FALSE, fig.height=4, fig.width=4.5----------------------------------
times <- seq(0, 6, 1 / 2) + .25
plotSpeedAndDisp(times)

## ----echo=FALSE, fig.height=4, fig.width=4.5----------------------------------
times <- seq(0, 6, 1 / 6) + .25
plotSpeedAndDisp(times)

## ----echo=FALSE, fig.height=4, fig.width=4.5----------------------------------
times <- seq(0, 6, 1 / 12) + .25
plotSpeedAndDisp(times, cex = .5)

