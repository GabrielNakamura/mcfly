#' @title Define the tolerance value to be used in ABC model
#'
#' @description Auxiliary function to define tolerance value to be used in the mcfly function.
#'
#' @details This function provides a procedure to estimate the values of tolerance that 
#'     can be used in [mcfly()] function, since the tolerance value will vary accordingly to 
#'     different alpha values.
#'    
#' @param comm.obs.meta Matrix containing the observed metacommunity.
#' @param phylo.meta Newick object with phylogenetic relationship among species of metacommunity.
#' @param envir Vector with values of environmental variable describing each community of metacommunity.
#' @param xy.coords Two column matrix with coordinates of communities in metacommunity.
#' @param Nmeta Numeric indicating the number of metacommunities to be simulated in ABC model.
#' @param parallel Numeric value indicating the number of cores to be used in parallel computation.
#' @param m Numeric value indicating the immigration rate at each site, reported as Hubbel´s m. This is the same parameter accepted by \code{\link{metasim}}.
#' @param W.r.prior Logical (TRUE or FALSE) indicating if the the W.r parameter would be a single value (FALSE) with value of 0, indicating a panmictic metacommunity
#'     or follow a prior distribution (TRUE) of values calculated as being the slopes of dispersal kernel indicating the contribution of species from neighboring patches 
#'     to the local immigrant pool.
#' @param OU.alpha Character indicating the type of prior that will be used in ABC model. The options were "uniform" for a uniform sample of 
#'     alpha values and "half-life" for a prior of alpha values represented as being half-life values, calculated as being log().
#' @param tol Numeric value indicating the initial tolerance to be used in mcfly function. Default is 1
#' @param n.timestep Numeric value indicating the number of timesteps to simulate the metacommunity.
#'     same argument used in \code{\link{metasim}} function.
#' @param max.sample.prior Numeric value indicating the maximmum size of sampled prior distribution. Default is 100*parallel.
#' @param max.size.posterior Numeric value indicating the maximum size of posterior distribution. Default is 20*parallel.
#' @param scenario.ID Character indicating the name of the simulation scenario. The same as used in \code{\link{metasim}}. Default is "mcfly".
#' @param output.dir.path Character indicating the name of directory to save simulations results and metadata used in \code{\link{metasim}}. Default is "delorean".
#' 
#' @return A matrix containing the total sample size used to find the posterior distribution, the length of posterior distribution,
#'     the minimum alpha value that can be sample given the phylogenetic tree used,
#'     the maximum alpha value that can be sample given the phylogenetic tree used,
#'     the highest 1% of tolerance values found in ABC model and the alpha values sampled from 
#'     prior distribution.
#'     
#' @examples 
#' \dontrun{
#'     
#'      # read datasets 
#'      comm <- data("Furnariidae") # community matrix
#'      phylo <- data("phylo_Furnariidae") # phylogenetic hypothesis
#'      envir_furnariidae <- data("envir") 
#'      coords <- envir_furnariidae[, c(1, 2)] # site coordinates
#'      envir <- envir_furnariidae[, -c(1, 2)] # environmental variables
#'      # run define_tol function
#'      res_tol <- define_tol(comm.obs.meta = comm,
#'            phylo.meta = phylo, 
#'            envir = envir, 
#'            xy.coords = coords, 
#'            parallel = NULL
#'            )
#'      tol_val <- res_tol[, 5]
#'            
#'     }
#' @seealso \code{\link{mcfly}}
#' 
#' @export
#'
define_tol <- function(comm.obs.meta, 
                       phylo.meta, 
                       envir, 
                       xy.coords, 
                       Nmeta,
                       parallel,
                       m = 0.5, 
                       W.r.prior = FALSE,
                       OU.alpha = "uniform",
                       tol = 1,
                       n.timestep = 50,
                       max.sample.size.prior = 100 * parallel,
                       max.sample.size.posterior = 20 * parallel, 
                       scenario.ID= "mcfly",
                       output.dir.path = "meta_output"){
  sum.res<-matrix(NA, Nmeta, 19, dimnames=list(1:Nmeta,c("Total_Prior_Size",
                                                       "Final_Posterior_Size",
                                                       "Min_Alpha_Prior",
                                                       "Max_Alpha_Prior",
                                                       "Tol",
                                                       "Nspp_Metacommunity",
                                                       "Alpha_Metacommunity",
                                                       "K_L-Posterior",
                                                       "K_H-Posterior",
                                                       "Mean_Posterior_Alpha",
                                                       "SD_Posterior_Alpha",
                                                       "Median_Posterior_Alpha",
                                                       "Alpha_L-Posterior",
                                                       "Alpha_H-Posterior",
                                                       "Mean_Posterior_W",
                                                       "SD_Posterior_W",
                                                       "Median_Posterior_W",
                                                       "W_L-Posterior",
                                                       "W_H-Posterior")))
  for(i in 1:Nmeta){
    output.dir.path.1 <- paste(output.dir.path, i, sep = ".")
    sim.ID <- paste(scenario.ID,i,sep=".")
    DRoot <- suppressWarnings(as.numeric(adephylo::distRoot(phylo.meta,
                                                            method = "patristic")[1])
    )
    min.alpha <- log(2)/(2*DRoot)
    max.alpha <- log(2)/(0.03333333*DRoot)
    prior.alpha <- runif(10*max.sample.size.prior,
                         min=min(alpha.interval),max=max(alpha.interval))
    
    
    comm.names.meta <- rownames(as.matrix(xy.coords), FALSE, prefix = "comm")
    names(envir) <- comm.names.meta
    alpha.meta <- sample(prior.alpha, 1)
    sigma.meta <- sqrt(sd(envir))
    theta.meta <- mean(envir)
    root.value.meta <- mean(envir)
    JL.meta <- rep(1000, nrow(comm.obs.meta))
    JM.meta <- sum(JL.meta)
    my.landscape.meta <- MCSim::make.landscape(site.coords = xy.coords,
                                               Ef = envir,
                                               m= m, JM = JM.meta, JL = JL.meta)
    thresh <- TRUE
    while(thresh){
      niche.pos.meta <- ape::rTraitCont(phy= phylo.meta, model= "OU", 
                                        alpha= alpha.meta,
                                        sigma= sigma.meta, theta= theta.meta, 
                                        root.value= root.value.meta)
      fitOU.comm <- suppressWarnings(geiger::fitContinuous(phylo.meta,
                                                           niche.pos.meta,
                                                           model = "OU", bounds = list(alpha = c(min.alpha, max.alpha)),
                                                           ncores = parallel)
      )
      thresh <- fitOU.comm$opt$alpha < (alpha.meta - alpha.meta * 0.2) |
        fitOU.comm$opt$alpha > (alpha.meta + alpha.meta*0.2)
    }
    
    if(subset == TRUE){
      niche.pos.meta <- sample(niche.pos.meta, size= ncol(comm.obs.meta))
      niche.pos.meta[order(as.numeric(gsub("s" , "", names(niche.pos.meta)
      )
      )
      )
      ]
    } else {niche.pos.meta <- niche.pos.meta}
    
    spp.names.meta <- names(niche.pos.meta)
    spp.freq.meta <- rmultinom(1, size= JM.meta, 
                               prob = scales::rescale(rnorm(length(niche.pos.meta),
                                                            0,10), to = c(0.4, 1)
                               )
    )
    spp.gamma<- spp.freq.meta/JM.meta
    
    ####initiate training####
    output.dir.path.2 <- paste("einstein", output.dir.path.1, sep=".")
    training<-mcfly(comm = comm.obs.meta,
                    phylo = envir,
                    envir = env.meta,
                    xy.coords = xy.coords,
                    m= 0.5,
                    W.r.prior = W.r.prior,
                    OU.alpha = OU.alpha,
                    tol = 1,
                    max.sample.size.prior = max.sample.size.prior,
                    sample.size.posterior = sample.size.posterior,
                    parallel = parallel,
                    output.dir.path = output.dir.path.2)
    tol.i <- 1 - quantile(training$Summary.Statistics, probs = 0.99)
    output.dir.path.3 <- paste("delorean", output.dir.path.1, sep = ".")
    test <- mcfly(comm = comm.obs.meta,
                phylo = phylo.meta,
                envir = envir,
                xy.coords = xy.coords,
                m = 0.5,
                W.r.prior = W.r.prior,
                tol = tol.i,
                OU.alpha = OU.alpha,
                max.sample.size.prior = max.sample.size.prior,
                sample.size.posterior = sample.size.posterior,
                parallel = parallel,
                output.dir.path = output.dir.path.3)
    sum.res[i,1] <- test$Sample_Attributes[2,]
    sum.res[i,2] <- test$Sample_Attributes[4,]
    sum.res[i,3] <- min.alpha
    sum.res[i,4] <- max.alpha
    sum.res[i,5] <- tol.i
    sum.res[i,7] <- alpha.meta
  }
  
  return(sum.res)
}
