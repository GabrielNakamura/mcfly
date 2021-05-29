#' Define tolerance values for ABC model in mcfly
#'
#' @param comm Matrix containing occurrences or abundances of species in sites. Species in columns and sites in rows.
#' @param phylo Newick object containing the phylogenetic relationship among species.
#' @param envir A one column matrix containing environmental variable for each community
#' @param xy.coords A two column matrix containing the coordinates of each community
#' @param occurrence Logical argument (TRUE or FALSE) indicating if community matrix must be transformed to presence/absence
#' @param entropy.order Numeric value indicating the scale of Rényi diversity, as accepted by \code{\link{renyi}}. Default is 1
#' @param niche.breadth Numeric value indicating the width of niche of species in the metacommunity, as accepted by \code{\link{metasim}}. Default is 10
#' @param m Numeric value indicating the immigration rate at each site, reported as Hubbel´s m. This is the same parameter accepted by \code{\link{metasim}}.
#' @param n.timestep Numeric value indicating the number of timesteps used in the simulation of metacommunities,
#'     this is the same argument used in \code{\link{metasim}}. Default is 50, it is not recommended the use of lower values.
#' @param summary.stat Character indicating the type of summary statistic that will be used in ABC model. Default is "correlation", that is calculated
#'     as the correlation between the diversity values calculated according to the Rényi scale defined in entropy.order argument.
#' @param OU.alpha Character indicating the type of prior that will be used in ABC model. The options were "uniform" for a uniform sample of
#'     alpha values and "half-life" for a prior of alpha values represented as being half-life values, calculated as being log().
#' @param max.sample.size.priorNumeric value that defines the maximum size of the posterior distribution. Default is 2400.
#' @param sample.size.posterior Numeric value that defines the minimum size of the posterior distribution. Default is 240.
#' @param parallel Numerical value indicating the numbers of cores that must be used in the parallel computation. Default is NULL, indicating that the
#'     calculations of ABC model will not be parallelized.
#' @param scenario.ID Character indicating the name of the simulation scenario. The same as used in \code{\link{metasim}}. Default is "doc".
#' @param output.dir.path Character indicating the name of directory to save simulations results and metadata used in \code{\link{metasim}}. Default is "einstein".
#' @param probs Numeric indicating the quantiles to be used in tolerance distribution values.
#'     Default is c(0.8,0.9,0.95,0.99)
#'
#' @return
#' @export
#'
#' @examples
define_tolerance<-function(comm,phylo,envir,xy.coords,
                      occurrence,entropy.order,
                      niche.breadth=10,m=0.5,
                      n.timestep=50,
                      summary.stat="correlation",
                      OU.alpha=c("uniform","half-life"),
                      max.sample.size.prior = 100*parallel,
                      sample.size.posterior = 20*parallel,
                      parallel = NULL,scenario.ID="doc",
                      output.dir.path="einstein",probs=c(0.8,0.9,0.95,0.99)){
 test<-mcfly(comm=comm,phylo=phylo,envir=envir,xy.coords=xy.coords,
                       occurrence=occurrence,entropy.order=entropy.order,
                       niche.breadth=niche.breadth,m=m,
                       n.timestep=n.timestep,summary.stat=summary.stat,
                       return.comm=FALSE,return.w.priors=FALSE,
                       return.alpha.priors = FALSE,tol=1,
                       scenario.ID=scenario.ID,W.r.prior= FALSE,OU.alpha=OU.alpha,
                       max.sample.size.prior = max.sample.size.prior,
                       sample.size.posterior = sample.size.posterior,parallel = parallel,
                       output.dir.path=output.dir.path)
 posterior<-test$Summary.Statistics
 tolerance.options<-1-quantile(test$Summary.Statistics,probs = probs)
 res<-list(Posterior=posterior,Tolerance=tolerance.options)
 return(res)
}
