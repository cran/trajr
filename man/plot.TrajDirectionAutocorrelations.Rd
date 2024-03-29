% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dirn_autocorr.R
\name{plot.TrajDirectionAutocorrelations}
\alias{plot.TrajDirectionAutocorrelations}
\title{Plot method for direction autocorrelation}
\usage{
\method{plot}{TrajDirectionAutocorrelations}(
  x,
  firstMinWindowSize = 10,
  type = "l",
  ylab = expression("C(" * Delta * s * ")"),
  xlab = expression(Delta * s),
  ...
)
}
\arguments{
\item{x}{Trajectory to be plotted.}

\item{firstMinWindowSize}{If not NULL, specifies a window size used to
calculate the first local minimum, which is then plotted as a point.}

\item{type, xlab, ylab}{Defaults for plotting.}

\item{...}{Additional arguments passed to \code{\link{plot}}.}
}
\description{
The \code{plot} method for \code{\link{TrajDirectionAutocorrelations}} objects.
Plots the direction autocorrelation function as returned by a call to
\code{TrajDirectionAutocorrelations}, with a optional dot at the first
local minimum.
}
