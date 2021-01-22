#'  
#' 
#' 
#' 
#' @importFrom stats cophenetic hclust cor
#' @param k only for parallel computation
#' @param sample.size.posterior size of posterior distribution
#' @param max.sample.size.prior maximum sample size for prior distribution
#' @param prior.alpha type of prior alpha, the same as in mcfly
#' @param prior.w type of prior w
#' @param theta.val theta
#' @param phylo phylogenetic tree in newick format
#' @param niche.sigma sigma
#' @param root.value root value
#' @param my.landscape simulated landscape
#' @param JM JM parameter, the same accepted by \code{\link{metasim}}
#' @param n.timestep number of timesteps for metacommunity simulation
#' @param spp.freq species frequency
#' @param niche.breadth niche breadth
#' @param occurrence if occurence or abundance
#' @param spp.names species names
#' @param names.comm communities names
#' @param entropy.order entropy order for Rényi series
#' @param summary.stat type of summary statistic
#' @param div diversity 
#' @param tol tolerance value
#' @param return.comm logical
#' @param scenario.ID name of simulation scenario
#' @param output.dir.path name of output directory
#'
#' @rdname mcfly
#' @include mcfly.R
#' @encoding UTF-8
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
                       summary.stat,
                       div,
                       tol,
                       return.comm,
                       return.w.priors,
                       return.alpha.priors,
                       scenario.ID,
                       output.dir.path){
  # List of results
  RES <- vector("list", sample.size.posterior)
  scenario.ID <- paste(scenario.ID, k, sep=".")
  output.dir.path <- paste(output.dir.path, k, sep = ".")
  # number of values in posterior distribution
  cont.size.ent <- 0  
  total.sample.size <- 0
  for(i in 1:max.sample.size.prior){
    sim.ID <- paste(scenario.ID, i, sep= ".")
    total.sample.size <- total.sample.size + 1
    # sampling alpha value
    alpha.sim <- sample(x = prior.alpha, size = 1) 
    # sampling W.r value
    W.r.sim <- sample(x = prior.w, size = 1) 
    # sampling theta
    theta.sim <- sample(x = theta.val, size = 1) 
    niche <- ape::rTraitCont(phy = phylo, model= "OU", alpha = alpha.sim,
                             sigma = niche.sigma, theta = theta.sim, root.value = root.value)
    # trait parameters and niche position----------
    k.niche <- as.numeric(picante::Kcalc(niche, phylo))
    niche.pos<- niche[spp.names]
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
    ent<- vegan::renyi(comm.sim, scales = entropy.order)
    # summary statistic correlation -----
    if(summary.stat == 1){
      # correlation between observed and simulated entropy
      cor.obs.simul.ent <- suppressWarnings(cor(ent, div))
      # tolerance value of ABC
      tol.sim.obs.ent <- 1 - abs(cor.obs.simul.ent)
      if(is.na(tol.sim.obs.ent)){
        tol.sim.obs.ent <- 1
      }
    }
    
    # summary statistic dimensionality (IV) -------
    if(summary.stat == 2){
      # calculating diversity metrics 
      # phylo metrics
      PDfaith_sim <-picante::pd(comm, phylo)$PD #phylo diversity
      mntd_sim <- picante::mntd(samp = comm, dis = cophenetic(phylo))
      PSV_sim <- picante::psv(samp = comm, tree = phylo, compute.var = TRUE, scale.vcv = TRUE)$PSVs
      DBPhylo_sim <- FD::dbFD(x = cophenetic(phylo), a = comm[,match(rownames(cophenetic(phylo)), colnames(comm))], 
                              calc.FRic = F, w.abun = FALSE, calc.FDiv = TRUE, calc.CWM = FALSE, calc.FGR = FALSE) 
      Peve_sim <- DBPhylo_sim$FEve
      Peve_sim[which(is.na(Peve_sim))] <- 0 
      # functional metrics
      distrait_sim <- vegan::vegdist(niche.pos, method = "euclidean")
      funct_dendro_sim <- ape::as.phylo(hclust(distrait_sim, method = "average"))
      FDfaith_sim <- picante::pd(comm, funct_dendro_sim)$PD #func diversity
      FEve_sim[which(is.na(FEve_sim))] <- 0
      FDiv_sim <- DBFunc_sim$FDiv #Functional divergence
      FDiv_sim[which(is.na(FDiv_sim))]<- 0
      # taxonomic metric
      rich <- rowSums(comm)
      # matrix M
      matrix_Msim <- as.matrix(data.frame(PDfaith_sim, mntd_sim, PSV_sim, Peve_sim, FDfaith_sim, FEve_sim, FDiv_sim, richness= rich_sim))
      IVs_sim_res <- ImportanceVal(matrix.M = matrix_Msim, scale = TRUE, method = "max", stopRule = TRUE)
      IVs_sim<- IVs_sim_res$IV.obs_stopRule
      cor_IV <- cor(IVs_sim, IVs_obs)
      tol_sim_obs_ent <- 1 - abs(cor_IV) # tolerance value of ABC
      if(is.na(tol_sim_obs_ent)){
        tol_sim_obs_ent <- 1
      }
    }
    
    # summary statistic dimensionality ----
    
    
    # posterior distribution-----
    if(tol.sim.obs.ent <= tol){
      cont.size.ent <- cont.size.ent + 1
      if(return.comm){
        RES[[cont.size.ent]]$comm.sim <- comm.sim
      }
      if(return.alpha.priors == T){
        RES[[cont.size.ent]]$alpha.prior <- alpha.sim
      }
      if(return.w.priors == TRUE){
        RES[[cont.size.ent]]$w.r.prior <- W.r.sim
      }
      RES[[cont.size.ent]]$w.simul.ent <- W.r.sim
      RES[[cont.size.ent]]$alpha.simul.ent <- alpha.sim
      RES[[cont.size.ent]]$theta.simul.ent <- theta.sim
      RES[[cont.size.ent]]$cor.posterior.ent <- cor.obs.simul.ent
      RES[[cont.size.ent]]$k.niche.simul <- k.niche
      RES[[cont.size.ent]]$sample.size <- total.sample.size
    } 
    if(cont.size.ent == sample.size.posterior){
      break
    }
  }
  return(RES)
}
