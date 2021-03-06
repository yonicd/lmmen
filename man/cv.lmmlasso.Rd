% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cv.lmmlasso.R
\name{cv.lmmlasso}
\alias{cv.lmmlasso}
\title{Cross Validation for lmmlasso package}
\usage{
cv.lmmlasso(dat, lambda = seq(0, 500, 5), ...)
}
\arguments{
\item{dat}{matrix, containing y,X,Z and subject variables}

\item{lambda}{numeric, path of positive regularization parameter, Default: seq(0, 500, 5)}

\item{...}{parameters to pass to lmmlasso}
}
\value{
lmmlasso fit object
}
\description{
Cross Validation for lmmlasso package as shown in example xxx
}
\examples{
 \dontrun{cv.lmmlasso(initialize_example(seed = 1))}
}
\seealso{
\code{\link[lmmlasso]{lmmlasso}}
}
