#' Title Simulation framework to estimate the influence of niche adaptation rate and dispersal limitation on species diversity distribution across environmental gradients.
#'
#' @description We introduce mcfly, an R package to estimate the influence of niche adaptation rate and dispersal limitation on species diversity distribution
#' of a given phylogenetic lineage across environmental gradients. For this, mcfly adopts the individual-based metacommunity simulation framework algorithm
#' implemented in the package MCSim (See MCSim package for more information on individual based model) coupled to Approximate Bayesian Computation (ABC) to estimate the posterior distribution of (i) the adaptation rate parameter
#' of Ornstein-Uhlenbeck (OU) evolutionary model and (ii) the slope the dispersal kernel, i.e, the probability density function of the dispersal success from
#' a source to a sink assemblage)that maximizes the association between species diversity across a set of sites and a given environmental gradient.
#'
#' @details This function estimate the influence of niche selection on species diversity across environmental gradients by applying to
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
#' @param Hill.numbers Logical argument indicating if Hill numbers must be used instead of Renyi entropy. Default is FALSE, indicating that Renyi entropy is the default
#' @param niche.breadth Numeric value indicating the width of niche of species in the metacommunity, as accepted by \code{\link{metasim}}. Default is 10
#' @param m Numeric value indicating the immigration rate at each site, reported as Hubbel´s m. This is the same parameter accepted by \code{\link{metasim}}.
#' @param n.timestep Numeric value indicating the number of timesteps used in the simulation of metacommunities,
#'     this is the same argument used in \code{\link{metasim}}. Default is 50, it is not recommended the use of lower values.
#' @param OU.alpha Character indicating the type of prior that will be used in ABC model. The options were "uniform" for a uniform sample of
#'     alpha values and "half-life" for a prior of alpha values represented as being half-life values, calculated as being log().
#' @param W.r.prior Logical (TRUE or FALSE) indicating if the the W.r parameter would be a single value (FALSE) with value of 0, indicating a panmictic metacommunity
#'     or follow a prior distribution (TRUE) of values calculated as being the slopes of dispersal kernel indicating the contribution of species from neighboring patches
#'     to the local immigrant pool.
#' @param tol Numeric value that defines the tolerance value (calculated as 1 - correlation) used in ABC model to assemble the posterior distribution. Default is 0.2.
#' @param sample.size.posterior Numeric value that defines the minimum size of the posterior distribution. Default is 240.
#' @param max.sample.size.prior Numeric value that defines the maximum size of the posterior distribution. Default is 2400.
#' @param HPD Numeric value indicating the probability mass for the Highest Density Interval for the posterior
#'     probability distribution obtained in ACB model. This is the same value used in \code{\link{hdi}}. Default is 0.9.
#' @param return.comm Logical (TRUE/FALSE), indicating if the simulated metacommunities must be returned in the output. Default is FALSE.
#' @param plot.res Logical, indicates if a plot with posterior and prior distribution, as well as the HPD must be showed. Default is FALSE
#' @param parallel Numerical value indicating the numbers of cores that must be used in the parallel computation. Default is NULL, indicating that the
#'     calculations of ABC model will not be parallelized.
#' @param scenario.ID Character indicating the name of the simulation scenario. The same as used in \code{\link{metasim}}. Default is "mcfly".
#' @param output.dir.path Character indicating the name of directory to save simulations results and metadata used in \code{\link{metasim}}. Default is "delorean".
#'
#' @import stats
#' @import magrittr
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_histogram
#' @importFrom ggplot2 geom_density
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_area
#' @import patchwork
#' @importFrom dplyr mutate
#' @importFrom dplyr case_when
#' @importFrom tibble tibble
#'
#' @return A list containing the following objects:
#'     \item{Time.spent}{The amout of time spent to run the analysis.}
#'     \item{Simulated.Metacomm}{Simulated community matrix.}
#'     \item{Data.Attributes}{A matrix containing the number of species in the phylogeny,
#'     the number of species in the metacommunity, number of sites in the metacommunity,
#'     maximmun distance between two species in the phylogenetic tree, tree depth}
#'      \item{Sample.Attributes}{A matrix containing information of posterior distribution as the maximum size of the posterior,
#'          total sample sample size of prior, the type of posterior (if alpha or half-life) and total size of posterior}
#'      \item{Alpha.Limits}{A matrix containing the minimum and maximum alfa parameters in the posterior distribution}
#'      \item{Alpha.prior.mode}{A character indicating the type of prior used in ABC (alpha or half-life)}
#'      \item{Alpha.prior.mode}{A character indicating the type of prior used in ABC (alpha or half-life)}
#'      \item{Alpha_Prior_Distribution}{Numeric vector with alpha values used as the prior distribution}
#'      \item{W.Prior.Distribution}{Numeric vector with w values used as the prior distribution}
#'      \item{Theta}{Numeric values indicating theta parameters in the simulation}
#'      \item{K.niche}{Numeric vector with Blomberg's K statistic of phylogenetic signal calculated for simulated traits}
#'      \item{Distance.measure}{Numeric vector with correlations between observed and simulated diversity metrics}
#'      \item{Alpha.Posterior.Distribution}{Numeric vector containing alpha values estimated from ABC model (posterior distribution)}
#'      \item{HPD.Alpha}{Numeric vector indicating the range of alpha posterior distribution corresponding to the High Posterior Density in a probability density distribution}
#'      \item{W.Posterior.Distribution}{Numeric vector containing the w values estimated from ABC model (posterior distribution)}
#'      \item{HPD.w}{Numeric vector indicating the range of w posterior distribution corresponding to the High Posterior Density in a probability density distribution}
#'      \item{Coordinates}{A matrix with x and y coordinates of sites in a metacommunity}
#'      \item{Observed.Entropy}{A vector containing the observed values of diversity for each site in the metacommunity}
#'      \item{Mean.Simulated.Entropy}{A vector containing mean values of diversity for each site computed from simulated metacommunities}
#'
#' @export
#'
#

mcfly <- function(comm, phylo, envir, xy.coords,
                  occurrence = TRUE, entropy.order = 1,
                  Hill.numbers = FALSE,
                  niche.breadth = 10,
                  m = 0.5,
                  n.timestep = 50,
                  OU.alpha=c("uniform","half-life"),
                  W.r.prior = FALSE,
                  tol = 0.2,
                  sample.size.posterior = 240, max.sample.size.prior = 2400,
                  HPD = 0.9,
                  return.comm = FALSE,
                  plot.res = FALSE,
                  parallel = NULL,
                  scenario.ID="mcfly",
                  output.dir.path = "delorean"){
  date.mat <- matrix(NA,2,1,dimnames=list(c("Started on","Finished on")," "))
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
    prior.alpha <-log(2)/runif(10*max.sample.size.prior,
                               min=0.03333333*DRoot.mat[,1],max=DRoot.mat[,1])
    alpha.mode<-"half-life"
  }

  if(occurrence){
    comm <- ifelse(comm > 0, yes = 1, no = 0)
  }
  names.comm <- rownames(comm)

  if(length(entropy.order)>1){
    stop("entropy.order must have only one value")
  }
  # statistics of observed communities and tree --------------------------

  div <- vegan::renyi(comm, scales = entropy.order, hill = Hill.numbers)
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

  prior.w.res <- prior.w

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
                      Hill.numbers = Hill.numbers,
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
    Mean.Simulated.Entropy <- RES$mean.sim.entropy
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
                               Hill.numbers = Hill.numbers,
                               div = div,
                               tol = tol,
                               return.comm = return.comm,
                               scenario.ID=scenario.ID,
                               output.dir.path = output.dir.path)
    RES_prior <- do.call(rbind, lapply(RES, function(x) x$prior))
    RES_posterior <- unlist(lapply(RES, function(x) x$posterior),
                            recursive = FALSE)

    # Total sample size
    # If NA in each last position total.sample.size is equal to
    #sample.size.posterior
    last.set <- seq.int(from = sample.size.posterior, to =
                          sample.size.posterior*n.cluster, by = sample.size.posterior)
    total.sample.size <- sapply(RES_posterior, function(x) ifelse(is.null(x), NA,
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
  # Organize RES_posterior results
  if(return.comm){
    Simulated.Metacomm <- lapply(RES_posterior, function(x) if(is.null(x)) NULL else x$comm.sim)
    Simulated.Metacomm <- Simulated.Metacomm[unlist(lapply(Simulated.Metacomm, function(x) !is.null(x)))]
  } else{
    Simulated.Metacomm <- NULL
  }
  w.simul.ent <- sapply(RES_posterior, function(x) ifelse(is.null(x), NA, x$w.simul.ent))
  alpha.simul.ent <- sapply(RES_posterior, function(x) ifelse(is.null(x), NA,
                                                              x$alpha.simul.ent))
  theta.simul.ent <- sapply(RES_posterior, function(x) ifelse(is.null(x), NA,
                                                              x$theta.simul.ent))
  cor.posterior.ent <- sapply(RES_posterior, function(x) ifelse(is.null(x), NA,
                                                                x$cor.posterior.ent))
  k.niche.simul <- sapply(RES_posterior, function(x) ifelse(is.null(x), NA,
                                                            x$k.niche.simul))

  # Remove NA
  theta.simul.ent <- theta.simul.ent[!is.na(theta.simul.ent)]
  cor.posterior.ent <- cor.posterior.ent[!is.na(cor.posterior.ent)]
  k.niche.simul <- k.niche.simul[!is.na(k.niche.simul)]
  posterior.dist.alpha <- alpha.simul.ent[!is.na(alpha.simul.ent)]
  if(W.r.prior){
    posterior.dist.w <- w.simul.ent[!is.na(w.simul.ent)]
  } else {
    posterior.dist.w<-0
  }
  n.tol <- sum(!is.na(alpha.simul.ent))
  if(n.tol > 0){
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
  if(W.r.prior == TRUE){
    spp.mat <- matrix(NA, 5, 1,
                      dimnames = list(c("Spp.phylogeny", "Spp.metacommunity", "n.sites", "max.dist.mst", "age.max.phylo"),
                                      "N")
    )
    spp.mat[1, ] <- length(phylo$tip.label)
    spp.mat[2, ] <- ncol(comm)
    spp.mat[3, ] <- nrow(comm)
    spp.mat[4, ] <- r # minimun spanning tree
    spp.mat[5, ] <- max(phytools::nodeHeights(phylo))
  } else{
    spp.mat <- matrix(NA, 4, 1,
                      dimnames = list(c("Spp.phylogeny", "Spp.metacommunity", "n.sites", "age.max.phylo"),
                                      "N")
    )
    spp.mat[1, ] <- length(phylo$tip.label)
    spp.mat[2, ] <- ncol(comm)
    spp.mat[3, ] <- nrow(comm)
    spp.mat[4, ] <- max(phytools::nodeHeights(phylo))
  }
  size.mat<-matrix(NA,4,1,dimnames=list(c("Maximum_prior","Total_prior",
                                          "Nominal_posterior","Final_posterior"),"Sample_size"))
  size.mat[1,] <- max.sample.size.prior
  size.mat[2,] <- total.sample.size
  size.mat[3,] <- sample.size.posterior
  size.mat[4,] <- n.tol
  date.mat[2,] <- date()
  res.list <- list(Time.spent = date.mat,
                   Simulated.Metacomm = Simulated.Metacomm,
                   Data.Attributes = spp.mat,
                   Sample.Attributes = size.mat,
                   Alpha.Limits = DRoot.mat,
                   Alpha.prior.mode = alpha.mode,
                   Alpha_Prior_Distribution = RES_prior[, 1][!is.na(RES_prior[, 1])],
                   W.Prior.Distribution = RES_prior[ , 2][!is.na(RES_prior[, 2])],
                   Theta = theta.simul.ent,
                   K.niche = k.niche.simul,
                   Distance.measure = (1 - cor.posterior.ent),
                   Alpha.Posterior.Distribution = posterior.dist.alpha,
                   HPD.Alpha = HPD.alpha,
                   W.Posterior.Distribution = posterior.dist.w,
                   HPD.w = HPD.w,
                   Coordinates = xy.coords,
                   Observed.Entropy = div,
                   Mean.Simulated.Entropy = do.call(rbind, lapply(RES, function(x) x$mean.sim.entropy))
  )
  print("...but your kids are gonna love it!!!")
  if(plot.res == TRUE){
    if(W.r.prior == TRUE){
      df_res_alpha <- data.frame(alpha_values = c(res.list$Alpha_Prior_Distribution,
                                                  res.list$Alpha.Posterior.Distribution),
                                 distribution = c(rep("prior", length(res.list$Alpha_Prior_Distribution)),
                                                  rep("posterior", length(res.list$Alpha.Posterior.Distribution))
                                 )
      )
      df_res_w <- data.frame(w_values = c(res.list$W.Prior.Distribution,
                                          res.list$W.Posterior.Distribution),
                             distribution = c(rep("prior", length(res.list$W.Prior.Distribution)),
                                              rep("posterior", length(res.list$W.Posterior.Distribution))
                             )
      )
      density_alpha<- density(x=res.list$Alpha.Posterior.Distribution,
                              from=res.list$Alpha.Limits[2],
                              to=res.list$Alpha.Limits[3])
      hdi_09_alpha <- HDInterval::hdi(object = density_alpha, credMass =
                                        0.9, allowSplit = TRUE)
      ord_alpha <- order(apply(hdi_09_alpha, MARGIN = 1, FUN = diff), decreasing = T)


      data_density_alpha <- tibble(x = density_alpha$x, y = density_alpha$y) %>%
        mutate(variable = case_when(
          (x >= hdi_09_alpha[ord_alpha[1], 1] & x <= hdi_09_alpha[ord_alpha[1], 2]) ~ "On",
          (x >= hdi_09_alpha[ord_alpha[2], 1]  & x <= hdi_09_alpha[ord_alpha[2], 2]) ~ "Off",
          TRUE ~ NA_character_))

      alpha_plot <- ggplot(data = df_res_alpha, aes(x = alpha_values, fill =  distribution)) +
        geom_histogram(aes(y=..density..), alpha=0.5,
                       position="identity") +
        geom_density(data = df_res_alpha, aes(x = alpha_values, color = distribution), alpha=.2, show.legend = F) +
        geom_area(data = filter(data_density_alpha, variable == 'On'), aes(x = x, y = y), fill = 'red', alpha = .6) +
        geom_area(data = filter(data_density_alpha, variable == 'Off'), aes(x = x, y = y), fill = 'red', alpha = .6) +
        labs(y = "Density", x = "Alpha values" , fill = "Distribution") +
        theme(legend.position = "right", panel.background = element_rect(fill = "transparent"),
              plot.margin = unit(c(0.4, 0.4, 0.4, 0.4), "mm"),
              legend.title = element_text(family = "Times", color = "black", face = "bold", size = 12),
              legend.text = element_text(family = "Times", color = "black", size = 12),
              axis.text = element_text(family = "Times", color = "black", size = 12),
              axis.line = element_line(colour = "black"),
              plot.subtitle = element_text(family = "Arial",
                                           color = "black",
                                           size = 9,
                                           hjust = 0.5,
                                           margin = margin(b = 6)
              ))
      density_w <- density(x=res.list$W.Posterior.Distribution)
      hdi_09_w <- HDInterval::hdi(object = density_w, credMass =
                                    0.9, allowSplit = TRUE)
      ord_w <- order(apply(hdi_09_w, MARGIN = 1, FUN = diff), decreasing = T)

      data_density_w <- tibble(x = density_w$x, y = density_w$y) %>%
        mutate(variable = case_when(
          (x >= hdi_09_w[ord_w[1], 1] & x <= hdi_09_w[ord_w[1],2]) ~ "On",
          (x >= hdi_09_w[ord_w[2], 1]  & x <= hdi_09_w[ord_w[2],2]) ~ "Off",
          TRUE ~ NA_character_))


      W_plot <- ggplot(data = df_res_w, aes(x = w_values, fill =  distribution)) +
        geom_histogram(aes(y=..density..), alpha=0.5,
                       position="identity") +
        geom_density(data = df_res_w, aes(x = w_values, color = distribution), alpha=.2, show.legend = F) +
        geom_area(data = filter(data_density_w, variable == 'On'), aes(x = x, y = y), fill = 'red', alpha = .6) +
        geom_area(data = filter(data_density_w, variable == 'Off'), aes(x = x, y = y), fill = 'red', alpha = .6) +
        labs(y = "Density", x = "W values" , fill = "Distribution") +
        theme(legend.position = "right", panel.background = element_rect(fill = "transparent"),
              plot.margin = unit(c(0.4, 0.4, 0.4, 0.4), "mm"),
              legend.title = element_text(family = "Times", color = "black", face = "bold", size = 12),
              legend.text = element_text(family = "Times", color = "black", size = 12),
              axis.text = element_text(family = "Times", color = "black", size = 12),
              axis.line = element_line(colour = "black"),
              plot.subtitle = element_text(family = "Arial",
                                           color = "black",
                                           size = 9,
                                           hjust = 0.5,
                                           margin = margin(b = 6)
              ))
      alpha_plot + W_plot
    } else{
      alpha_plot
    }
  }
  return(res.list)
}
