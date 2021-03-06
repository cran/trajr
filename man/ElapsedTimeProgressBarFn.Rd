% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/progressbar.R
\name{ElapsedTimeProgressBarFn}
\alias{ElapsedTimeProgressBarFn}
\title{A general purpose progress bar that reports elapsed time rather than number of items}
\usage{
ElapsedTimeProgressBarFn(numItems, reportFn)
}
\arguments{
\item{numItems}{Number of items to be processed}

\item{reportFn}{A function used to report changing progress}
}
\value{
A function which should be called for each item as it is processed.
}
\description{
A general purpose progress bar that reports elapsed time rather than number of items
}
