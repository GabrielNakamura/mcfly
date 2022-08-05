#' Internal auxiliar function to parallel process in mcfly
#'
#' @param k numeric, auxiliar parameter to paralell proccess
#' @param sample.size.posterior numeric indicating the size of posterior distribution
#' @param max.sample.size.prior numeric indicating the sie of prior distribution
#' @param prior.alpha character indication the type of prior alpha same as mcfly
#' @param prior.w logical indicating if the w must be estimated
#' @param theta.val numeric theta value
#' @param phylo phylogeny
#' @param niche.sigma numeric indicating the sigma for nich
#' @param root.value numeric indicating the mean value of traits
#' @param my.landscape landscape
#' @param JM total number of individuals in simulation process
#' @param n.timestep number of timesteps in simulation process
#' @param spp.freq frequency of species
#' @param niche.breadth numeric indicating the niche breath of species
#' @param occurrence logical
#' @param spp.names names of species
#' @param names.comm names for communities
#' @param entropy.order diversity value accordingly to entropy order
#' @param div diversity value
#' @param tol tolerance value
#' @param return.comm logical
#' @param scenario.ID character
#' @param output.dir.path character
#'
#' @return list with results to be used in mcfly function
#' @import stats
#'
#' @examples
f.internal <- function(k,
                       sample.size.posterior,
                       max.sample.size.prior,
                       prior.alpha,
                       prior.w,
                       theta.val,
                       phylo,
                       niche.sigma,
                       root.value,
                       my.landscape,
                       JM,
                       n.timestep,
                       spp.freq,
                       niche.breadth,
                       occurrence,
                       spp.names,
                       names.comm,
                       entropy.order,
                       div,
                       tol,
                       return.comm,
                       scenario.ID,
                       output.dir.path){
  # List of results
  RES <- vector("list", sample.size.posterior)
  scenario.ID=paste(scenario.ID,k,sep=".")
  output.dir.path<-paste(output.dir.path,k,sep = ".")
  # number of values in posterior distribution
  cont.size.ent <- 0
  total.sample.size <- 0
  RES.prior <- matrix(NA, nrow =  max.sample.size.prior, ncol = 2)
  mat.ent <- matrix(NA, nrow = length(names.comm), ncol = length(max.sample.size.prior),
                dimnames = list(names.comm,
                                paste("sim", 1:max.sample.size.prior, sep = ""))
                )
  for(i in 1:max.sample.size.prior){
    sim.ID<-paste(scenario.ID,i,sep=".")
    total.sample.size <- total.sample.size+1
    # sampling alpha value
    alpha.sim <- sample(x = prior.alpha, size = 1)
    # sampling W.r value
    W.r.sim <- sample(x = prior.w, size = 1)
    RES.prior[i, 1] <- alpha.sim
    RES.prior[i, 2] <- W.r.sim
    # sampling theta
    theta.sim <- sample(x = theta.val, size = 1)
    niche <- ape::rTraitCont(phy = phylo, model= "OU", alpha = alpha.sim,
                             sigma = niche.sigma, theta = theta.sim, root.value = root.value)
    # trait parameters and niche position----------
    k.niche <- as.numeric(picante::Kcalc(niche, phylo))
    niche.pos<-niche[spp.names]
    sim <- MCSim::metasim(landscape = my.landscape, nu = 0,
                          speciation.limit = 0, JM.limit = JM,
                          n.timestep = n.timestep,
                          W.r = W.r.sim, save.sim = FALSE, trait.Ef = niche.pos,
                          trait.Ef.sd = niche.breadth, gamma.abund = spp.freq,
                          taxa.list = spp.names,scenario.ID=scenario.ID,sim.ID=sim.ID,output.dir.path=output.dir.path)
    comm.out <- sim$J.long
    comm.frame <- comm.out[which(comm.out[,1]==n.timestep),]
    comm.sim <- tapply(comm.frame$count,list(comm.frame$site, comm.frame$spp),
                       sum)
    rownames(comm.sim) <- names.comm
    if(occurrence){
      comm.sim <- ifelse(comm.sim > 0, yes = 1, no = 0)
    }
    # entropy values----
    ent <- vegan::renyi(comm.sim, scales = entropy.order)
    # summary statistic correlation -----
      # correlation between observed and simulated entropy
      cor.obs.simul.ent <- suppressWarnings(cor(ent, div))
      # tolerance value of ABC
      tol.sim.obs.ent <- 1 - abs(cor.obs.simul.ent)
      if(is.na(tol.sim.obs.ent)){
        tol.sim.obs.ent <- 1
      }

    # posterior distribution-----
    if(tol.sim.obs.ent <= tol){
      cont.size.ent <- cont.size.ent + 1
      if(return.comm){
        RES[[cont.size.ent]]$comm.sim <- comm.sim
      }
      RES[[cont.size.ent]]$w.simul.ent <- W.r.sim
      RES[[cont.size.ent]]$alpha.simul.ent <- alpha.sim
      RES[[cont.size.ent]]$theta.simul.ent <- theta.sim
      RES[[cont.size.ent]]$cor.posterior.ent <- cor.obs.simul.ent
      RES[[cont.size.ent]]$k.niche.simul <- k.niche
      RES[[cont.size.ent]]$sample.size <- total.sample.size
      mat.ent[, cont.size.ent] <- vegan::renyi(comm.sim, scales = entropy.order)[, 1]
    }
    if(cont.size.ent == sample.size.posterior){
      break
    }
  }
  list_res <- vector(mode = "list", length = 3)
  list_res[[1]] <- RES
  list_res[[2]] <- RES.prior
  mean.sim.entropy <- apply(mat.ent, MARGIN = 1, function(x) mean(x, na.rm = TRUE))
  list_res[[3]] <- mean.sim.entropy
  names(list_res) <- c("posterior", "prior", "mean.sim.entropy")
  return(list_res)
}
