% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in Rtmph3UnEb/file63ef3330b755
\name{buildEnsembleKF}
\alias{buildEnsembleKF}
\title{Create an Ensemble Kalman filter algorithm to sample from latent states.}
\usage{
buildEnsembleKF(model, nodes, control = list())
}
\arguments{
\item{model}{A NIMBLE model object, typically representing a state space model or a hidden Markov model}

\item{nodes}{A character vector specifying the latent model nodes which the Ensemble Kalman filter will estimate.}

\item{control}{A list specifying different control options for the particle filter.  Options are described in the details section below.}
}
\description{
Create an Ensemble Kalman filter algorithm for a given NIMBLE state space model.
}
\details{
The \code{control()} list option is described in detail below:
\describe{
 \item{saveAll}{Indicates whether to save state samples for all time points (TRUE), or only for the most recent time point (FALSE)}
\item{timeIndex}{An integer used to manually specify which dimension of the latent state variable indexes time.
 Only needs to be set if the number of time points is less than or equal to the size of the latent state at each time point.}
}

Runs an Ensemble Kalman filter to estimate a latent state given observations at each time point.  The ensemble Kalman filter
is a Monte Carlo approximation to a Kalman filter that can be used when the model's transition euqations do not follow a normal distribution.
Latent states (x[t]) and observations (y[t]) can be scalars or vectors at each time point,
and sizes of observations can vary from time point to time point.
In the BUGS model, the observations (y[t]) must be equal to some (possibly nonlinear) deterministic function
of the latent state (x[t]) plus an additive error term.  Currently only normal and multivariate normal
error terms are supported.
The transition from x[t] to x[t+1] does not have to be normal or linear.  Output from the posterior distribution of the latent
states is stored in \code{mvSamples}.
}
\examples{
\dontrun{
model <- nimbleModel(code = ...)
my_ENKFF <- buildEnsembleKF(model, 'x')
Cmodel <- compileNimble(model)
Cmy_ENKF <- compileNimble(my_ENKF, project = model)
Cmy_ENKF$run(m = 100000)
ENKF_X <- as.matrix(Cmy_ENKF$mvSamples, 'x')
hist(ENKF_X)
}
}
\author{
Nicholas Michaud
}
\references{
Houtekamer, P.L., and H.L. Mitchell. (1998). Data assimilation using an ensemble Kalman filter technique. \emph{Monthly Weather Review}, 126(3), 796-811.
}

