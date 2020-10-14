
comm = comm.obs
subset = subset 
occurrence = occurrence
env = env
site.coords = xy.coords 
tree = tree 
OU.alpha = alpha.values 
sigma = sqrt(sd(env)) 
theta = theta.sim 
root.value = mean(env) 
runs = 10
ncores = 2
n.timestep = n.timestep
W.r = W.r


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
      i= 1
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
