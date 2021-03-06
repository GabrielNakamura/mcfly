% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcfly_function.R
\name{mcfly}
\alias{mcfly}
\title{Title Estimating the influence of stabilizing selection on species distribution}
\usage{
mcfly(
  comm,
  phylo,
  envir,
  xy.coords,
  occurrence = TRUE,
  entropy.order = 1,
  niche.breadth = 10,
  m = 0.5,
  n.timestep = 50,
  OU.alpha = c("uniform", "half-life"),
  W.r.prior = FALSE,
  summary.stat = "correlation",
  tol = 0.2,
  sample.size.posterior = 240,
  max.sample.size.prior = 2400,
  HPD = 0.9,
  return.comm = FALSE,
  return.w.priors = FALSE,
  return.alpha.priors = TRUE,
  parallel = NULL,
  scenario.ID = "mcfly",
  output.dir.path = "delorean"
)
}
\arguments{
\item{comm}{Matrix containing occurrences or abundances of species in sites. Species in columns and sites in rows.}

\item{phylo}{Newick object containing the phylogenetic relationship among species.}

\item{envir}{A one column matrix containing environmental variable for each community}

\item{xy.coords}{A two column matrix containing the coordinates of each community}

\item{occurrence}{Logical argument (TRUE or FALSE) indicating if community matrix must be transformed to presence/absence}

\item{entropy.order}{Numeric value indicating the scale of Rényi diversity, as accepted by \code{\link{renyi}}. Default is 1}

\item{niche.breadth}{Numeric value indicating the width of niche of species in the metacommunity, as accepted by \code{\link{metasim}}. Default is 10}

\item{m}{Numeric value indicating the immigration rate at each site, reported as Hubbel´s m. This is the same parameter accepted by \code{\link{metasim}}.}

\item{n.timestep}{Numeric value indicating the number of timesteps used in the simulation of metacommunities,
this is the same argument used in \code{\link{metasim}}. Default is 50, it is not recommended the use of lower values.}

\item{OU.alpha}{Character indicating the type of prior that will be used in ABC model. The options were "uniform" for a uniform sample of
alpha values and "half-life" for a prior of alpha values represented as being half-life values, calculated as being log().}

\item{W.r.prior}{Logical (TRUE or FALSE) indicating if the the W.r parameter would be a single value (FALSE) with value of 0, indicating a panmictic metacommunity
or follow a prior distribution (TRUE) of values calculated as being the slopes of dispersal kernel indicating the contribution of species from neighboring patches
to the local immigrant pool.}

\item{summary.stat}{Character indicating the type of summary statistic that will be used in ABC model. Default is "correlation", that is calculated
as the correlation between the diversity values calculated according to the Rényi scale defined in entropy.order argument. Another option is "dimensionality"
but it is not implemented yet.}

\item{tol}{Numeric value that defines the tolerance value (calculated as 1 - correlation) used in ABC model to assemble the posterior distribution. Default is 0.2.}

\item{sample.size.posterior}{Numeric value that defines the minimum size of the posterior distribution. Default is 240.}

\item{max.sample.size.prior}{Numeric value that defines the maximum size of the posterior distribution. Default is 2400.}

\item{HPD}{Numeric value indicating the probability mass for the Highest Density Interval for the posterior
probability distribution obtained in ACB model. This is the same value used in \code{\link{hdi}}. Default is 0.9.}

\item{return.comm}{Logical (TRUE/FALSE), indicating if the simulated metacommunities must be returned in the output. Default is FALSE.}

\item{return.w.priors}{Logical (TRUE/FALSE), indicating if the prior distribution of W.r values used in ABC model must be returned in the output.
Default is FALSE}

\item{return.alpha.priors}{Logical (TRUE/FALSE), indicating if the the prior distribution of alpha values must be returned in the output. Default is FALSE.}

\item{parallel}{Numerical value indicating the numbers of cores that must be used in the parallel computation. Default is NULL, indicating that the
calculations of ABC model will not be parallelized.}

\item{scenario.ID}{Character indicating the name of the simulation scenario. The same as used in \code{\link{metasim}}. Default is "mcfly".}

\item{output.dir.path}{Character indicating the name of directory to save simulations results and metadata used in \code{\link{metasim}}. Default is "delorean".}
}
\value{

}
\description{
mcfly function to estimate the influence of stabilizing niche selection on species diversity across environmental gradients
}
\details{
This function estimate the influence of stabilizing niche selection on species diversity across environmental gradients by applying to
a occurrence matrix of species containing presence/absence or abundance an Approximate Bayesian Computation (ABC) framework. We used in ABC a individual
based-model from MCSim package
}
