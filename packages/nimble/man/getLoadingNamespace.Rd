% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in Rtmph3UnEb/file63ef76b18d3d
\name{getLoadingNamespace}
\alias{getLoadingNamespace}
\title{return the namespace in which a nimbleFunction is being loaded}
\usage{
getLoadingNamespace()
}
\description{
\code{nimbleFunction} constructs and evals a reference class definition.  When \code{nimbleFunction} is used in package source code, this can lead to problems finding things due to namespace issues.  Using \code{where = getLoadingNamespace()} in a \code{nimbleFunction} in package source code should solve this problem.
}
\details{
\code{nimbleFunction}s defined in the NIMBLE source code use \code{where = getLoadingNamespace()}.  Please let the NIMBLE developers know if you encounter problems with this.
}

