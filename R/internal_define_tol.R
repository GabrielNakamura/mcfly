#' @title Define the tolerance value to be used in ABC model
#'
#' @description Auxiliary function to define tolerance value to be used in the mcfly function.
#'
#' @details This function provides a procedure to estimate the values of tolerance that
#'     can be used in [mcfly()] function, since the tolerance value will vary accordingly to
#'     different alpha values.
#'
#' @param envir Vector with values of environmental variable describing each community of metacommunity.
#' @param xy.coords Two column matrix with coordinates of communities in metacommunity.
#' @param parallel Numeric value indicating the number of cores to be used in parallel computation.
#' @param m Numeric value indicating the immigration rate at each site, reported as Hubbel´s m. This is the same parameter accepted by \code{\link{metasim}}.
#' @param W.r.prior Logical (TRUE or FALSE) indicating if the the W.r parameter would be a single value (FALSE) with value of 0, indicating a panmictic metacommunity
#'     or follow a prior distribution (TRUE) of values calculated as being the slopes of dispersal kernel indicating the contribution of species from neighboring patches
#'     to the local immigrant pool.
#' @param OU.alpha Character indicating the type of prior that will be used in ABC model. The options were "uniform" for a uniform sample of
#'     alpha values and "half-life" for a prior of alpha values represented as being half-life values, calculated as being log().
#' @param n.timestep Numeric value indicating the number of timesteps to simulate the metacommunity.
#'     same argument used in \code{\link{metasim}} function.
#' @param scenario.ID Character indicating the name of the simulation scenario. The same as used in \code{\link{metasim}}. Default is "mcfly".
#' @param output.dir.path Character indicating the name of directory to save simulations results and metadata used in \code{\link{metasim}}. Default is "delorean".
#' @param comm  community matrix
#' @param phylo phylo object
#' @param occurrence logical
#' @param entropy.order numeric
#' @param niche.breadth numeric
#' @param summary.stat character
#' @param max.sample.size.prior numeric
#' @param sample.size.posterior numeric
#' @param probs numeric
#'
#' @return numeric indicating tolerance value
#'
define_tolerance<-function(comm,phylo,envir,xy.coords,
                           occurrence,entropy.order,
                           niche.breadth=10,m=0.5,
                           n.timestep=50,
                           summary.stat="correlation",
                           W.r.prior,
                           OU.alpha=c("uniform","half-life"),
                           max.sample.size.prior = 100*parallel,
                           sample.size.posterior = 20*parallel,
                           parallel = NULL,scenario.ID="doc",
                           output.dir.path="einstein",probs=0.99){
  tol<-as.numeric(1-quantile(mcfly(comm=comm,phylo=phylo,envir=envir,xy.coords=xy.coords,
                                   occurrence=occurrence,entropy.order=entropy.order,
                                   niche.breadth=niche.breadth,m=m,
                                   n.timestep=n.timestep,summary.stat=summary.stat,
                                   return.comm=FALSE,return.w.priors=FALSE,
                                   return.alpha.priors = FALSE,tol=1,
                                   scenario.ID=scenario.ID,W.r.prior= W.r.prior,OU.alpha=OU.alpha,
                                   max.sample.size.prior = max.sample.size.prior,
                                   sample.size.posterior = sample.size.posterior,parallel = parallel,
                                   output.dir.path=output.dir.path)$Summary.Statistics,probs = probs))
  return(tol)
}
