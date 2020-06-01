#' @title Mcfly 
#' 
#' @description Mcfly function to estimate the influence of stabilizing niche selection on species diversity across environmental gradients 
#'
#' @importFrom MCSim make.landscape metasim
#' @importFrom vegan decostand renyi
#' @importFrom ape rTraitCont
#' @importFrom geiger fitContinuous
#' @importFrom reshape2 acast
#' @importFrom picante Kcalc
#' @importFrom parallel makeCluster clusterEvalQ parLapply stopCluster
#' @param comm Matrix containing occurrences or abundances of species in sites. Species in columns and sites in rows.
#' @param subset Logical. If TRUE, only the subset of species present in comm will be considered in community simulations. Note that niche position will be simulated for all species in the phylogeny.
#' @param occurrence Logical. If TRUE, comm will be transformed to presence/absence matrix.
#' @param env Numeric vector containing values of the environmental variable for the sites.
#' @param site.coords Geographic coordinates (longitude and latitude) of sites.
#' @param tree Newick object containing a phylogenetic tree.
#' @param OU.alpha Numeric vector containing the values of alpha parameter of Ornstein-Uhlenbeck model.
#' @param sigma Numeric value indicating the sigma parameter to be used in \code{\link[ape:rTraitCont]{rTraitCont}} function.
#' @param theta Numeric value indicating the theta parameter to be used in \code{\link{rTraitCont}}.
#' @param root.value Numeric value indicating the root.value parameter to be used in \code{\link{rTraitCont}}.
#' @param runs Numeric value indicating the number of simulations to be performed.
#' @param ncores Numeric value indicating the number of cores to perform parallel computation. Default ncores = NULL.
#' @param area.m2 A scalar or a vector. Area of each site to be used in metacommunity simulation. This only matters if immigration is defined as no. individuals / m2. Default value = 1.
#' @param m Immigration rate at each site reported as Hubbell’s m. Same as used in \code{\link{make.landscape}} of MCSim package. Default value = 0.5.
#' @param JM The number of individuals (abundance data) or species (presence/absence data) in the metacommunity. Same argument used in \code{\link{make.landscape}}.
#' @param JM.limit The maximum number of individuals in the metacommunity. Same argument used in \code{\link{make.landscape}}.
#' @param JL The number of individuals (abundance data) or species (presence/absence data) in each individual community. Same argument used in \code{\link{make.landscape}}.
#' @param nu The probability that a novel species will appear during a recruitment event. Same argument used in \code{\link{metasim}}. Default value = 0.
#' @param speciation.limit The speciation rate of simulated metacommunities. Same argument used in \code{\link{metasim}}. Default value = 0.
#' @param n.timestep Number of time steps (generations) in each simulation. Default value = 50.
#' @param W.r The slope of the dispersal kernel, indicating the intensity of dispersal of individuals among sites. Default = 0, indicating a panmitic metacommunity.
#' @param scenario.ID Provides the name of scenario that will be simulated.
#' @param sim.ID Provides a name for the simulation, which is saved along with parameter values in the output.dir.path.
#' @param output.dir.path Name of the folder that will be created in working directory of R session to store parameter values of simulations.
#' @param OU.alpha.v Numeric vector containing values of OU alpha parameter.
#' @param mat.env Matrix containing values for environmental gradient.
#' @param landscape Landscape where the simulations will go on. Object generated with \code{\link{make.landscape}} from MCSim package.
#' @param spp.freq Numeric vector of regional abundance, used in \code{\link{metasim}} from MCSim package.
#' @param reps Tricky argument to work in parallel process.
#' 
#' @return List, with the following components: 
#'     \item{Entropy}{Entropy for empirical communities.}
#'     \item{Predicted.entropy.1}{Mean predicted entropy of order 1 for simulated metacommunities.}
#'     \item{Predicted.entropy.2}{Mean predicted entropy of order 2 for simulated metacommunities.}
#'     \item{Predicted.entropy.12}{Mean predicted entropy of order 12 for simulated metacommunities.}
#'     \item{K.niche.position}{Mean values for K statistic calculated with bootstrap procedure for simulated niche position.}
#'     \item{Alpha.niche.position}{Mean values for alfa parameter of OU model calculated using bootstrap procedure for simulated niche position.}
#'     \item{Z0.niche.position}{Mean values of simulated niche position at the root of the phylogeny calculated using bootstrap procedure.}
#'     \item{AIC}{Numerical matrix containing AIC statistic for each OU alpha value and entropy order.}
#'     \item{W}{Numerical matrix containing Akaike weights (wi) from model selection, for each OU alpha value and entropy order.}
#'     \item{R2}{Numerical matrix containing R2 statistic from linear model relating the observed entropies to mean predicted entropy values, for each OU alpha value and entropy order.}
#'     
#' @export
Mcfly <- function(comm, subset, occurrence = FALSE, env, site.coords, tree,
                  OU.alpha, sigma,theta, root.value, runs,
                  ncores = NULL, area.m2 = 1, m = 0.5, JM = sum(comm), JM.limit = JM, JL = rowSums(comm),
                  nu = 0, speciation.limit = 0, n.timestep = 50, W.r = 0, scenario.ID = "species.sorting",
                  sim.ID = "data", output.dir.path= "OUTPUT_DATA"){
  
  mat.env <- as.matrix(env)
  if(dim(mat.env)[2] > 1){
    stop("\n Only one environmental gradient must be supplied in env \n")
  }
  if(!inherits(comm, what = "matrix")){
    stop("\n comm must be a matrix with presence/absence or species abundance in communities \n")
  }  
  if(dim(site.coords)[2] > 2){
    stop("\n site.coords must be a matrix with x y sites coordinates \n")
  }
  if(n.timestep < 50){
    warning("\n The number of timesteps used in the community simulation is lower than 50 \n")
  }
  rownames(mat.env) <- rownames(comm)
  obs.ent <- vegan::renyi(comm, scales = c(1, 2, 12))
  obs.ent.1 <- obs.ent$`1`
  obs.ent.2 <- obs.ent$`2`
  obs.ent.12 <- obs.ent$`12`
  landscape <- MCSim::make.landscape(site.coords = site.coords, Ef = env,
                                     area.m2 = area.m2, m = m, JM = JM, JL = JL)
  spp.freq <- colSums(comm)
  L <- matrix(NA, nrow(comm), ncol(comm), dimnames = list(rownames(comm), colnames(comm)))
  if (is.numeric(ncores)) {
    #######innitiating parallel procces######
    SET <- rep(1, runs)
    CL <- parallel::makeCluster(ncores, type = "PSOCK", setup_timeout = 0.5)
    parallel::clusterEvalQ(CL, library(MCSim))
    parallel::clusterEvalQ(CL, library(ape))
    parallel::clusterEvalQ(CL, library(geiger))
    res_run <- parallel::parLapply(CL, SET, mcfly_RUN, OU.alpha.v = OU.alpha, comm = comm, subset = subset, occurrence = occurrence, tree = tree,
                                   mat.env = mat.env,landscape = landscape, JM = JM, spp.freq = spp.freq, sigma = sigma,
                                   theta = theta, W.r= W.r, root.value = root.value, n.timestep = n.timestep, 
                                   scenario.ID = scenario.ID, sim.ID= sim.ID, output.dir.path = output.dir.path)
    parallel::stopCluster(CL)
  } else{
    res_run <- vector("list", runs)
    for(i in 1:runs){
      res_run[[i]] <- mcfly_RUN(OU.alpha.v = OU.alpha, comm = comm, subset = subset, occurrence = occurrence, tree = tree,
                                mat.env = mat.env,landscape = landscape, JM = JM, spp.freq = spp.freq, sigma = sigma,
                                theta = theta, W.r= W.r, root.value = root.value, n.timestep = n.timestep, 
                                scenario.ID = scenario.ID, sim.ID= sim.ID, output.dir.path = output.dir.path)
    }
  }
  ####mean values for each comm for each alpha value####
  sigsq.P <- matrix(unlist(lapply(res_run, function(x) x[[1]][,1])), 
                    nrow = length(OU.alpha), ncol = runs, byrow = FALSE, 
                    dimnames = list(OU.alpha, paste("run", 1:runs)))
  K.P <- matrix(unlist(lapply(res_run, function(x) x[[1]][,2])), 
                nrow = length(OU.alpha), ncol = runs, byrow = FALSE, 
                dimnames = list(OU.alpha, paste("run", 1:runs)))
  alpha.P <- matrix(unlist(lapply(res_run, function(x) x[[1]][,3])), 
                    nrow = length(OU.alpha), ncol = runs, byrow = FALSE, 
                    dimnames = list(OU.alpha, paste("run", 1:runs)))
  z0.P <- matrix(unlist(lapply(res_run, function(x) x[[1]][,4])), 
                 nrow = length(OU.alpha), ncol = runs, byrow = FALSE, 
                 dimnames = list(OU.alpha, paste("run", 1:runs)))
  mean.ent.1 <- matrix(NA, nrow = nrow(comm), ncol = length(OU.alpha), 
                       dimnames = list(rownames(comm), c(paste("alpha", OU.alpha, sep = "="))))
  mean.ent.2 <- matrix(NA, nrow = nrow(comm), ncol = length(OU.alpha),
                       dimnames = list(rownames(comm), c(paste("alpha", OU.alpha, sep = "="))))
  mean.ent.12 <- matrix(NA, nrow = nrow(comm), ncol = length(OU.alpha),
                        dimnames = list(rownames(comm), c(paste("alpha", OU.alpha, sep = "="))))
  for(i in 1:nrow(comm)){
    mean.ent.1[i,] <- apply(matrix(unlist(lapply(res_run, function(x) x[[2]][,i])), 
                                   nrow = runs, ncol = length(OU.alpha), byrow = TRUE), 2, mean)
    mean.ent.2[i,] <- apply(matrix(unlist(lapply(res_run, function(x) x[[3]][,i])), 
                                   nrow = runs, ncol = length(OU.alpha), byrow = TRUE), 2, mean)
    mean.ent.12[i,] <- apply(matrix(unlist(lapply(res_run, function(x) x[[4]][,i])),
                                    nrow = runs, ncol = length(OU.alpha), byrow = TRUE), 2, mean)
  }
  #####AIC selection#####
  AIC.mod <- matrix(NA, nrow = length(OU.alpha), ncol = 3,
                    dimnames = list(c(paste("alpha", OU.alpha, sep = "=")), paste(c("Ent.1", "Ent.2", "Ent.12"))))
  R2 <- matrix(NA, nrow = length(OU.alpha), ncol = 3,
               dimnames = list(c(paste("alpha", OU.alpha, sep = "=")), paste(c("Ent.1", "Ent.2", "Ent.12"))))
  for (k in 1:length(OU.alpha)){
    mod.temp1 <- stats::lm(obs.ent.1~mean.ent.1[,k])
    mod.temp2 <- stats::lm(obs.ent.2~mean.ent.2[,k])
    mod.temp3 <- stats::lm(obs.ent.12~mean.ent.12[,k])
    AIC.mod[k,1] <- stats::AIC(mod.temp1)
    AIC.mod[k,2] <- stats::AIC(mod.temp2)
    AIC.mod[k,3] <- stats::AIC(mod.temp3)
    R2[k,1] <- summary(mod.temp1)$r.squared
    R2[k,2] <- summary(mod.temp2)$r.squared
    R2[k,3] <- summary(mod.temp3)$r.squared
  }
  W <- matrix(NA, nrow = length(OU.alpha), ncol = 3, 
              dimnames = list(c(paste("alpha", OU.alpha, sep = "=")), paste(c("Ent.1", "Ent.2", "Ent.12"))))
  for (l in 1:ncol(AIC.mod)){
    W[,l] <- MuMIn::Weights(AIC.mod[,l])
  }
  K.mean <- matrix(NA, nrow = length(OU.alpha), ncol = runs,
                   dimnames = list(c(paste("alpha", OU.alpha, sep = "=")), 1:runs))
  alpha.mean <- matrix(NA, nrow = length(OU.alpha), ncol = runs,
                       dimnames = list(c(paste("alpha", OU.alpha, sep = "=")), 1:runs))
  sigsq.mean <- matrix(NA, nrow = length(OU.alpha), ncol = runs,
                       dimnames = list(c(paste("alpha", OU.alpha, sep = "=")), 1:runs))
  z0.mean <- matrix(NA, nrow = length(OU.alpha), ncol = runs,
                    dimnames = list(c(paste("alpha", OU.alpha, sep = "=")), 1:runs))
  for(j in 1:runs){
    K.mean[,j] <- apply(K.P[, sample(1:ncol(K.P), size = runs, replace = TRUE), drop = FALSE], 1, mean)
    alpha.mean[,j] <- apply(alpha.P[, sample(1:ncol(alpha.P), size = runs, replace = TRUE), drop = FALSE], 1, mean)
    sigsq.mean[,j] <- apply(sigsq.P[, sample(1:ncol(sigsq.P), size = runs, replace = TRUE), drop = FALSE], 1, mean)
    z0.mean[,j] <- apply(z0.P[, sample(1:ncol(z0.P), size = runs, replace = TRUE), drop = FALSE], 1, mean)
  }
  K.mean.boot<-matrix(NA, nrow=3, ncol = length(OU.alpha),
                      dimnames = list(c("lower", "mean", "upper"), c(paste("alpha", OU.alpha, sep = "="))))
  K.mean.boot[1,] <- apply(K.mean, 1, function(x) stats::quantile(x, probs = 0.025))
  K.mean.boot[3,] <- apply(K.mean, 1, function(x) stats::quantile(x, probs = 0.975))
  K.mean.boot[2,] <- apply(K.mean, 1, function(x) stats::quantile(x, probs = 0.5))
  alpha.mean.boot <- matrix(NA, nrow = 3, ncol = length(OU.alpha),
                            dimnames = list(c("lower", "mean", "upper"), c(paste("alpha", OU.alpha, sep = "="))))
  alpha.mean.boot[1,] <- apply(alpha.mean, 1,function(x) stats::quantile(x, probs = 0.025))
  alpha.mean.boot[3,] <- apply(alpha.mean, 1,function(x) stats::quantile(x, probs = 0.975))
  alpha.mean.boot[2,] <- apply(alpha.mean, 1,function(x) stats::quantile(x, probs = 0.5))
  sigsq.mean.boot <- matrix(NA, nrow = 3, ncol = length(OU.alpha),
                            dimnames = list(c("lower", "mean", "upper"), c(paste("alpha", OU.alpha, sep = "="))))
  sigsq.mean.boot[1,] <- apply(sigsq.mean, 1, function(x) stats::quantile(x, probs = 0.025))
  sigsq.mean.boot[3,] <- apply(sigsq.mean, 1, function(x) stats::quantile(x, probs = 0.975))
  sigsq.mean.boot[2,] <- apply(sigsq.mean, 1, function(x) stats::quantile(x, probs = 0.5))
  z0.mean.boot <- matrix(NA, nrow = 3, ncol = length(OU.alpha),
                         dimnames = list(c("lower", "mean", "upper"), c(paste("alpha", OU.alpha, sep = "="))))
  z0.mean.boot[1,] <- apply(z0.mean, 1, function(x) stats::quantile(x, probs = 0.025))
  z0.mean.boot[3,] <- apply(z0.mean, 1, function(x) stats::quantile(x, probs = 0.975))
  z0.mean.boot[2,] <- apply(z0.mean, 1, function(x) stats::quantile(x, probs = 0.5))
  Res <- list(Entropy = obs.ent, Predicted.entropy.1 = mean.ent.1, Predicted.entropy.2 = mean.ent.2, 
              Predicted.entropy.12 = mean.ent.12, K.niche.position = K.mean.boot, 
              Alpha.niche.position = alpha.mean.boot, Sigsq.niche.position = sigsq.mean.boot, 
              Z0.niche.position = z0.mean.boot, AIC = AIC.mod, W = W, R2 = R2)
  return(Res)
}
