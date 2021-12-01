#' Title Estimating the influence of stabilizing selection on species distribution
#'
#' @description mcfly function to estimate the influence of stabilizing niche selection on species diversity across environmental gradients
#'
#' @details This function estimate the influence of stabilizing niche selection on species diversity across environmental gradients by applying to
#' a occurrence matrix of species containing presence/absence or abundance an Approximate Bayesian Computation (ABC) framework. We used in ABC a individual
#' based-model from MCSim package
#'
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
#' @param OU.alpha Character indicating the type of prior that will be used in ABC model. The options were "uniform" for a uniform sample of
#'     alpha values and "half-life" for a prior of alpha values represented as being half-life values, calculated as being log().
#' @param W.r.prior Logical (TRUE or FALSE) indicating if the the W.r parameter would be a single value (FALSE) with value of 0, indicating a panmictic metacommunity
#'     or follow a prior distribution (TRUE) of values calculated as being the slopes of dispersal kernel indicating the contribution of species from neighboring patches
#'     to the local immigrant pool.
#' @param summary.stat Character indicating the type of summary statistic that will be used in ABC model. Default is "distance", that is calculated
#'     as the complement of the correlation between the diversity values calculated according to the Rényi scale defined in entropy.order argument. Another option is "dimensionality"
#'     but it is not implemented yet.
#' @param tol Numeric value that defines the tolerance value (calculated as 1 - correlation) used in ABC model to assemble the posterior distribution. Default is 0.2.
#' @param sample.size.posterior Numeric value that defines the minimum size of the posterior distribution. Default is 240.
#' @param max.sample.size.prior Numeric value that defines the maximum size of the posterior distribution. Default is 2400.
#' @param HPD Numeric value indicating the probability mass for the Highest Density Interval for the posterior
#'     probability distribution obtained in ACB model. This is the same value used in \code{\link{hdi}}. Default is 0.9.
#' @param return.comm Logical (TRUE/FALSE), indicating if the simulated metacommunities must be returned in the output. Default is FALSE.
#' @param return.w.priors  Logical (TRUE/FALSE), indicating if the prior distribution of W.r values used in ABC model must be returned in the output.
#'     Default is FALSE
#' @param return.alpha.priors Logical (TRUE/FALSE), indicating if the the prior distribution of alpha values must be returned in the output. Default is FALSE.
#' @param parallel Numerical value indicating the numbers of cores that must be used in the parallel computation. Default is NULL, indicating that the
#'     calculations of ABC model will not be parallelized.
#' @param scenario.ID Character indicating the name of the simulation scenario. The same as used in \code{\link{metasim}}. Default is "mcfly".
#' @param output.dir.path Character indicating the name of directory to save simulations results and metadata used in \code{\link{metasim}}. Default is "delorean".
#'
#' @import stats
#'
#' @return
#' @export
#'
#
mcfly <- function(comm, phylo, envir, xy.coords,
                     occurrence = TRUE, entropy.order = 1,
                     niche.breadth = 10,
                     m = 0.5,
                     n.timestep = 50,
                     OU.alpha=c("uniform","half-life"),
                     W.r.prior = FALSE,
                     summary.stat = "distance",
                     tol = 0.2,
                     sample.size.posterior = 240, max.sample.size.prior = 2400,
                     HPD = 0.9,
                     return.comm = FALSE,
                     return.w.priors = FALSE,
                     return.alpha.priors = TRUE,
                     parallel = NULL,
                     scenario.ID="mcfly",
                     output.dir.path = "delorean"){
  date.mat<-matrix(NA,2,1,dimnames=list(c("Started on","Finished on")," "))
  date.mat[1,] <- date()
  if(!sample.size.posterior%%1==0){
    stop("\n sample.size.posterior must be an integer")
  }
  if(!max.sample.size.prior%%1==0){
    stop("\n max.sample.size.prior must be an integer")
  }
  if (sample.size.posterior > max.sample.size.prior) {
    stop("\n max.sample.size.prior must be equal or higher than
         sample.size.posterior")
  }
  if (tol < 0 | tol>1) {
    stop("\n tol must vary between 0 and 1")
  }
  if (HPD < 0 | HPD>1) {
    stop("\n HPD must vary between 0 and 1")
  }
  SUMMARY.stat <- c("distance", "dimensionality")
  summary.stat <- pmatch(summary.stat, SUMMARY.stat)
  if (length(summary.stat) != 1) {
    stop("\n Only one argument is accepted in summary.stat")
  }
  if (is.na(summary.stat)) {
    stop("\n Invalid summary.stat")
  }

  if(any(!all(colnames(comm) %in% phylo$tip.label))){
    stop("there is at least one species in the metacommunity not present in the phylogeny")
  }

  if(!ape::is.ultrametric(phylo)){
    stop("phylo is not an ultrametric tree")
  }
  # forge alpha prior distribution
  DRoot.mat<-matrix(NA,1,3,dimnames=list("Value",c("Tree_age",
                                      "Minimum_alpha","Maximum_alpha")))
  DRoot.mat[,1]<-max(phytools::nodeHeights(phylo))
  DRoot.mat[,2]<-log(2)/(DRoot.mat[,1])
  DRoot.mat[,3]<-log(2)/(0.03333333*DRoot.mat[,1])
  if(OU.alpha=="uniform"){
    prior.alpha <- runif(10*max.sample.size.prior, min = DRoot.mat[,2],
                         max=DRoot.mat[,3])
    alpha.mode<-"uniform"
  }
  if(OU.alpha=="half-life"){
    prior.alpha<-log(2)/runif(10*max.sample.size.prior,
                          min=0.03333333*DRoot.mat[,1],max=DRoot.mat[,1])
    alpha.mode<-"half-life"
  }

  if(return.alpha.priors){
    prior.alpha.res <- prior.alpha
  } else{
    prior.alpha.res<- NULL
  }

  if(occurrence){
    comm <- ifelse(comm > 0, yes = 1, no = 0)
  }
  names.comm <- rownames(comm)

  if(length(entropy.order)>1){
    stop("entropy.order must have only one value")
  }
  # statistics of observed communities and tree --------------------------
  div <- vegan::renyi(comm, scales = entropy.order)
  if(inherits(envir, "matrix") | inherits(envir, "data.frame")){
    if(ncol(envir)>1){
      stop("envir must have only one column")
    }
    names.envir <- rownames(envir)
    envir <- envir[, 1, drop = TRUE]
    names(envir) <- names.envir
  }
  envir<-scales::rescale(envir,c(0,100))
  root.value <- mean(envir)
  niche.sigma<-sqrt(sd(envir))
  theta.val <- as.numeric(envir[which(div >= quantile(x=div, probs=0.95))])

  # Defines dispersal limitation
  if(W.r.prior){
    dist.xy <- scales::rescale(dist(xy.coords,diag=T,upper=T),c(0,1))
    r <- max(dist.xy*as.dist(ape::mst(dist.xy),diag=T,upper=T))
    dlk <- runif(n=max.sample.size.prior,min=0,max=1)
    prior.w <- round(-log(dlk)/(r^2),3)
  } else {
    prior.w <- 0
  }
  if(return.w.priors){
    prior.w.res <- prior.w
  } else{
     prior.w.res<- NULL
  }
  # define carrying capacity of sites
  if(occurrence){
    comm <- comm*100
  }
  JL <-rowSums(comm)
  JM <- sum(JL)
  # define species names
  spp.names<-colnames(comm)
  # define seeder for recruitment from regional species pool
  spp.freq<-ifelse(colSums(comm)>0,colSums(comm),50)

  # initiating simulated metacommunity with MCSim -----------------------------
  utils::capture.output({my.landscape <- MCSim::make.landscape(site.coords = xy.coords,
                                        Ef = envir,
                                        m = m,
                                        JM = JM,
                                        JL = JL)}) # landscape for metacommunity
                                                      #simulation
  print("I can guess you guys aren't ready for that yet...")
  # makeCluster
  newClusters <- FALSE
  if (is.numeric(parallel)) {
    if(!parallel%%1==0){
      stop("\n parallel must an integer")
    }
    n.cluster <- parallel
    parallel <- parallel::makeCluster(parallel, type = "PSOCK")
    newClusters <- TRUE
  }
  if (!inherits(parallel, "cluster")) {
    RES <- f.internal(k = 0,
                      sample.size.posterior = sample.size.posterior,
                      max.sample.size.prior = max.sample.size.prior,
                      prior.alpha = prior.alpha,
                      prior.w = prior.w,
                      theta.val = theta.val,
                      phylo = phylo,
                      niche.sigma = niche.sigma,
                      root.value = root.value,
                      my.landscape = my.landscape,
                      spp.freq = spp.freq,
                      JM = JM,
                      n.timestep = n.timestep,
                      niche.breadth = niche.breadth,
                      occurrence = occurrence,
                      spp.names = spp.names,
                      names.comm = names.comm,
                      entropy.order = entropy.order,
                      summary.stat = summary.stat,
                      div = div,
                      tol = tol,
                      return.comm = return.comm,
                      scenario.ID=scenario.ID,
                      output.dir.path = output.dir.path)
    # Total sample size
    # If NA in the last position total.sample.size is equal to
    #sample.size.posterior
    RES_prior <- RES$prior # priors
    RES <- RES$posterior # posteriors
    total.sample.size <- sapply(RES, function(x) ifelse(is.null(x), NA,
                                      x$sample.size))[sample.size.posterior]
    total.sample.size <- ifelse(is.na(total.sample.size),
                        yes = max.sample.size.prior, no = total.sample.size)
  }
  else {
    # Recalculate the sample.size.posterior and max.sample.size.prior to
    #redistribute between clusters
    sample.size.posterior <- ceiling(sample.size.posterior/n.cluster)
    max.sample.size.prior <- ceiling(max.sample.size.prior/n.cluster)
    RES <- parallel::parLapply(parallel, seq_len(n.cluster), fun = f.internal,
                               sample.size.posterior = sample.size.posterior,
                               max.sample.size.prior = max.sample.size.prior,
                               prior.alpha = prior.alpha,
                               prior.w = prior.w,
                               theta.val = theta.val,
                               phylo = phylo,
                               niche.sigma = niche.sigma,
                               root.value = root.value,
                               my.landscape = my.landscape,
                               JM = JM,
                               n.timestep = n.timestep,
                               spp.freq = spp.freq,
                               niche.breadth = niche.breadth,
                               occurrence = occurrence,
                               spp.names = spp.names,
                               names.comm = names.comm,
                               entropy.order = entropy.order,
                               summary.stat = summary.stat,
                               div = div,
                               tol = tol,
                               return.comm = return.comm,
                               scenario.ID=scenario.ID,
                               output.dir.path = output.dir.path)
    RES_prior <- do.call(rbind, lapply(RES, function(x) x$prior))
    RES <- unlist(lapply(RES, function(x) x$posterior),
                                recursive = FALSE)

    # Total sample size
    # If NA in each last position total.sample.size is equal to
    #sample.size.posterior
    last.set <- seq.int(from = sample.size.posterior, to =
                  sample.size.posterior*n.cluster, by = sample.size.posterior)
    total.sample.size <- sapply(RES, function(x) ifelse(is.null(x), NA,
                                                  x$sample.size))[last.set]
    total.sample.size <- sum(ifelse(is.na(total.sample.size), yes =
                                max.sample.size.prior, no = total.sample.size))
    # Calculate the effective sample.size.posterior and max.sample.size.prior
    sample.size.posterior <- sample.size.posterior*n.cluster
    max.sample.size.prior <- max.sample.size.prior*n.cluster
  }
  # stopCluster
  if (newClusters) {
    parallel::stopCluster(parallel)
  }
  # Organize results
  if(return.comm){
    COMM.sim <- lapply(RES, function(x) if(is.null(x)) NULL else x$comm.sim)
    COMM.sim <- COMM.sim[unlist(lapply(COMM.sim, function(x) !is.null(x)))]
  } else{
    COMM.sim <- NULL
  }
  w.simul.ent <- sapply(RES, function(x) ifelse(is.null(x), NA, x$w.simul.ent))
  alpha.simul.ent <- sapply(RES, function(x) ifelse(is.null(x), NA,
                                                    x$alpha.simul.ent))
  theta.simul.ent <- sapply(RES, function(x) ifelse(is.null(x), NA,
                                                    x$theta.simul.ent))
  cor.posterior.ent <- sapply(RES, function(x) ifelse(is.null(x), NA,
                                                      x$cor.posterior.ent))
  k.niche.simul <- sapply(RES, function(x) ifelse(is.null(x), NA,
                                                        x$k.niche.simul))
  # Remove NA
  theta.simul.ent<-theta.simul.ent[!is.na(theta.simul.ent)]
  cor.posterior.ent<-cor.posterior.ent[!is.na(cor.posterior.ent)]
  k.niche.simul<-k.niche.simul[!is.na(k.niche.simul)]
  posterior.dist.alpha <- alpha.simul.ent[!is.na(alpha.simul.ent)]
  if(W.r.prior){
    posterior.dist.w <- w.simul.ent[!is.na(w.simul.ent)]
  } else {
    posterior.dist.w<-0
  }
  n.tol <- sum(!is.na(alpha.simul.ent))
  if(n.tol>0){
    HPD.alpha <- HDInterval::hdi(object = posterior.dist.alpha, credMass =
                                                      HPD, allowSplit=TRUE)
      if(W.r.prior){
        HPD.w <- HDInterval::hdi(object = posterior.dist.w, credMass = HPD,
                                                            allowSplit=TRUE)
      } else {HPD.w<-NA
        }
  } else {
    HPD.alpha <- NA
    HPD.w <- NA
    }
  spp.mat<-matrix(NA,2,1,dimnames=list(c("Spp.phylogeny","Spp.metacommunity"),"N"))
  spp.mat[1,]<-length(phylo$tip.label)
  spp.mat[2,]<-ncol(comm)
  size.mat<-matrix(NA,4,1,dimnames=list(c("Maximum_prior","Total_prior",
                      "Nominal_posterior","Final_posterior"),"Sample_size"))
  size.mat[1,] <- max.sample.size.prior
  size.mat[2,] <- total.sample.size
  size.mat[3,] <- sample.size.posterior
  size.mat[4,] <- n.tol
  date.mat[2,] <- date()
  res.list <- list(Time.spent=date.mat,
                   COMM.sim = COMM.sim,
                   Species.Pools = spp.mat,
                   Sample_Attributes = size.mat,
                   Alpha_Limits=DRoot.mat,
                   Alpha.prior.mode=alpha.mode,
                   Alpha_Prior_Distribution = RES_prior[, 1][!is.na(RES_prior[, 1])],
                   W_Prior_Distribution = RES_prior[ , 2][!is.na(RES_prior[, 2])],
                   Theta = theta.simul.ent,
                   K_niche = k.niche.simul,
                   Summary.Statistics = cor.posterior.ent,
                   Alpha_Posterior_Distribution = posterior.dist.alpha,
                   HPD_Alpha = HPD.alpha,
                   W_Posterior_Distribution = posterior.dist.w,
                   HPD_w = HPD.w
                   )
  print("...but your kids are gonna love it!!!")
  return(res.list)
}

