% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in Rtmph3UnEb/file63ef24c37a5a
\name{MCMCsuite}
\alias{MCMCsuite}
\title{Executes multiple MCMC algorithms and organizes results.}
\usage{
MCMCsuite(code, constants = list(), data = list(), inits = list(),
  monitors = character(), niter = 10000, burnin = 2000, thin = 1,
  summaryStats = c("mean", "median", "sd", "CI95_low", "CI95_upp"),
  calculateEfficiency = FALSE, MCMCs = "nimble", MCMCdefs = list(),
  winbugs_directory = "C:/WinBUGS14", winbugs_program = "WinBUGS",
  openbugs_directory = "C:/OpenBUGS323", openbugs_program = "OpenBUGS",
  stan_model = "", stan_inits = NULL, stan_data = NULL,
  stanNameMaps = list(), makePlot = TRUE, savePlot = TRUE,
  plotName = "MCMCsuite", setSeed = TRUE,
  check = getNimbleOption("checkModel"), debug = FALSE)
}
\arguments{
\item{code}{The quoted code expression representing the model, such as the return value from a call to \code{nimbleCode}).
No default value, this is a required argument.}

\item{constants}{A named list giving values of constants for the model.
This is the same as the \code{constants} argument which would be passed to \code{nimbleModel}.
Default value is list().}

\item{data}{A named list giving the data values for the model.
This is the same as the \code{data} argument which would be passed to \code{nimbleModel} or \code{model$setData}.
Default value is \code{list()}.}

\item{inits}{A named list giving the initial values for the model.
This is the same as the \code{inits} argument which would be passed to \code{nimbleModel} or \code{model$setInits}.
Default value is \code{list()}.}

\item{monitors}{A character vector giving the node names or variable names to monitor.
The samples corresponding to these nodes will be stored in the output samples, will have summary statistics calculated, and density and trace plots generated.
Default value is all top-level stochastic nodes of the model.}

\item{niter}{Number of MCMC iterations to run.
This applies to all MCMC algorithms in the suite.
Default value is 10,000.}

\item{burnin}{Number of initial, post-thinning, MCMC iterations to discard.
Default value is 2,000.}

\item{thin}{Thinning interval for the MCMC samples.
This applies to all MCMC algorithms in the suite.  The thinning occurs prior to the burnin samples being discarded.
Default value is 1.}

\item{summaryStats}{A character vector, specifying the summary statistics to calculate on the MCMC samples.
Each element may be the character name of an exisiting R function (possibly user-defined) which acts on a numeric vector and returns a scalar (e.g., \code{mean} or \code{sd},
or a character string which when parsed and evaluted will define such a function (e.g., \code{function(x) mean(sqrt(x))}).
Default value is \code{c('mean', 'median', 'sd', 'CI95_low', 'CI95_upp')}, where the final two elements are functions which calculate the limits of a 95 percent Bayesian credible interval.}

\item{calculateEfficiency}{A logical, specifying whether to calculate the efficiency for each MCMC algorithm.  Efficiency is defined as the effective sample size (ESS) of each model parameter divided by the algorithm runtime (in seconds).  Default is FALSE.}

\item{MCMCs}{A character vector specifying the MCMC algorithms to run.
\code{'winbugs'} specifies WinBUGS;
\code{'openbugs'} specifies OpenBUGS;
\code{'jags'} specifies JAGS;
\code{'stan'} specifies Stan; in this case, must also provide the \code{'stan_model'} argument;
\code{'nimble'} specifies NIMBLE's default MCMC algorithm;
\code{'nimble_noConj'} specifies NIMBLE's default MCMC algorithm without the use of any conjugate Gibbs sampling;
\code{'nimble_RW'} specifies NIMBLE MCMC algorithm using only random walk Metropolis-Hastings (\code{'RW'}) samplers;
\code{'nimble_slice'} specifies NIMBLE MCMC algorithm using only slice (\code{'slice'}) samplers;
\code{'autoBlock'} specifies NIMBLE MCMC algorithm with block sampling of dynamically determined parameter groups attempting to maximize sampling efficiency;
Anything else will be interpreted as NIMBLE MCMC algorithms, and must have associated entries in the MCMCdefs argument.
Default value is \code{'nimble'}, which specifies NIMBLE's default MCMC algorithm.}

\item{MCMCdefs}{A named list of MCMC definitions.  The names of list elements should corespond to any custom MCMC algorithms specified in the \code{MCMCs} argument.
The list elements should be quoted expressions, enclosed in {} braces.  When executed, the internal code must return an MCMC configuration object,
specifying the corresponding MCMC algorithm; in particular, setting the appropriate samplers.  The code may assume existance of the R model object \code{Rmodel},
and must *return* the MCMC configuration object.  Therefore, the final line of such a code block would frequently be a standalone \code{MCMCconf}, to return this object.}

\item{winbugs_directory}{A character string giving the directory of the executable WinBUGS program for the WinBUGS MCMC.
This argument will be passed directly to the bugs(...) call, from the R2WinBUGS library.
Default value is \code{'C:/WinBUGS14'}.}

\item{winbugs_program}{A character string giving the name of the WinBUGS program, for the WinBUGS MCMC.
This argument will be passed directly to the bugs(...) call, from the R2WinBUGS library.
Default value is \code{'WinBUGS'}.}

\item{openbugs_directory}{A character string giving the directory of the executable OpenBUGS program for the OpenBUGS MCMC.
This argument will be passed directly to the bugs(...) call, from the R2WinBUGS library.
Default value is \code{'C:/OpenBUGS323'}.}

\item{openbugs_program}{A character string giving the name of the OpenBUGS program, for the OpenBUGS MCMC.
This argument will be passed directly to the bugs(...) call, from the R2WinBUGS library.
Default value is \code{'OpenBUGS'}.}

\item{stan_model}{A character string specifying the location and name of the model file (\code{'modelName.stan'}) for use with the Stan MCMC program.
This argument must include the \code{'.stan'} extension, and must be provided whenever the \code{MCMCs} argument includes \code{'stan'}.}

\item{stan_inits}{A character string specifying the location and name of the inits file (\code{'modelName.init.R'}) for use with the Stan MCMC program.
This argument must include the \code{'.init.R'} extension, and must be provided whenever the \code{MCMCs} argument includes \code{'stan'}.
If omitted, it will attempt to locate an inits file in the same directory as the Stan model file.}

\item{stan_data}{A character string specifying the location and name of the data file (in the form \code{'modelName.data.R'}) for use with the Stan MCMC program.
This argument must include the \code{'.data.R'} extension, and must be provided whenever the \code{MCMCs} argument includes \code{'stan'}.
If omitted, it will attempt to locate a data file in the same directory as the Stan model file.}

\item{stanNameMaps}{A list specifying name mappings between Stan and WinBUGS/OpenBUGS.
The syntax for list elements is list(BUGS_PARAM_NAME = list(StanSourceName = 'STAN_PARAM_NAME', transform = function(x) TRANSFORMATION_FUNCTION(x))).
The transformation is optional.}

\item{makePlot}{Logical argument, specifying whether to generate the trace plots and posterior density plots, for each monitored node.
Default value is \code{TRUE}.}

\item{savePlot}{Logical argument, specifying whether to save the trace plots and density plots.
Plots will be saved into the current working directory.
Only used when \code{makePlot == TRUE}.
Default value is \code{TRUE}.}

\item{plotName}{Character string, giving the file name for saving the trace plots and density plots.
Only used when \code{makePlot == TRUE} and \code{savePlot == TRUE}.
Default value is \code{'MCMCsuite'}.}

\item{setSeed}{Logical argument, specifying whether to set.seed(0) prior to MCMC sampling.
Default value is \code{TRUE}.}

\item{check}{Logical argument, specifying whether to check the model object for missing or invalid values.  Default is given by the NIMBLE option 'checkModel', see help on \code{nimbleOptions} for details.}

\item{debug}{Logical argument, specifying whether to enter a \code{browser()} at the onset of executing each MCMC algrithm.
For use in debugging individual MCMC algorithms, if necessary.
Default value is FALSE.}

\item{...}{For internal use only}
}
\value{
Returns a named list containing elements:
samples: A 3-dimensional array containing samples from each MCMC algorithm.
summary: A 3-dimensional array containing summary statistics for each variable and algorithm.
timing: A numeric vector containing timing information.
efficiency: Minimum and mean sampling efficiencies for each algorithm (only provided if option calculateEfficiency = TRUE).
See the NIMBLE User Manual for more information about the organization of the return object.
}
\description{
Creates, runs, and organizes output from a suite of MCMC algorithms, all applied to the same model, data, and initial values.
This can include WinBUGS, OpenBUGS, JAGS and Stan MCMCs, as well as NIMBLE MCMC algorithms.
Trace plots and density plots for the MCMC samples may also be generated and saved.
}
\details{
Creates and runs an MCMC Suite.
By default, this will execute the specified MCMCs, record all samples, generate summary statistics, and create and save trace plots and posterior density plots.
This default behavior can ben altered via a variety of arguments.
Following execution of the MCMC algorithms, returns a named list containing \code{samples}, \code{summary}, and \code{timing} elements.
See the NIMBLE User Manual for more information about the organization of the return object.
}
\examples{
\dontrun{
code <- nimbleCode({
    mu ~ dnorm(0, 1)
    x ~ dnorm(mu, 1)
})
output <- MCMCsuite(code,
                    data = list(x=3),
                    inits = list(mu=0),
                    niter = 10000,
                    monitors = 'mu',
                    MCMCs = c('nimble', 'nimble_RW'),
                    summaryStats = c('mean', 'sd', 'max', 'function(x) max(abs(x))'),
                    makePlot = FALSE)
}
}
\author{
Daniel Turek
}

