#' @rdname Mcfly
#' @include Mcfly.R
#' @encoding UTF-8

OU.alpha.v <- alpha.values
comm<- comm.obs

mcfly_RUN <- function(OU.alpha.v, comm, subset, occurrence, tree, mat.env, landscape, JM, n.timestep, spp.freq, W.r,
                      sigma, theta, root.value, scenario.ID, sim.ID, output.dir.path, reps = 1, size= 100){
  f.internal.run <- function(x, comm, subset, occurrence, tree, mat.env, landscape, JM, n.timestep, spp.freq, W.r,
                             sigma, theta, root.value, scenario.ID, sim.ID, output.dir.path){
    output.dir.path <- paste(output.dir.path, x, sep = ".")
    thresh <- TRUE
    
    while(thresh){
      
      P <- ape::rTraitCont(phy = tree, model = "OU", sigma = sigma, alpha = x,
                           theta = theta, root.value = root.value)
      fitOU <- geiger::fitContinuous(tree, P, model = "OU", bounds = list(alpha=c(min=exp(-500),max=exp(1))))
      k.niche <- picante::Kcalc(P, tree, checkdata = FALSE)
      alpha.niche <- fitOU$opt$alpha
      if(x==0){
        thresh <- alpha.niche>0.0001|k.niche>1.2|k.niche<0.8
      } else {
        thresh <- alpha.niche<(x - x*0.05)|alpha.niche>(x + x*0.05)
      }
    }
    if(subset){
      subset.spp <- colnames(comm)
      niche <- P[subset.spp]
    } else{
      niche <- P
    }
    sigsq.P <- fitOU$opt$sigsq
    K.P <- as.numeric(picante::Kcalc(P, tree, checkdata = FALSE))
    alpha.P <- alpha.niche
    z0.P <- fitOU$opt$z0
    t.P <- t(niche)
    NB <- matrix(NA, nrow(mat.env), ncol(t.P))
    for(l in 1:ncol(t.P)){
      for(k in 1:nrow(mat.env)){
        NB[k,l] <- sqrt((mat.env[k,]-(t.P[,l]))^2)
      }
    }
    sd.NB <- sqrt(max(NB)-apply(NB,2, stats::sd))
    meta <- MCSim::metasim(landscape = landscape, scenario.ID = scenario.ID,
                           sim.ID = sim.ID, alpha.fisher = 1, nu = 0, speciation.limit = 0,
                           JM.limit = JM, n.timestep = n.timestep, W.r = W.r, save.sim = FALSE,
                           output.dir.path = output.dir.path, trait.dispersal = NULL,
                           trait.dispersal.median = 1, trait.dispersal.range = 0,
                           trait.Ef = niche, trait.Ef.sd = sd.NB, gamma.abund = spp.freq,
                           J.t0 = NULL, taxa.list = colnames(comm))
    L.frame <- meta$J.long
    L.frame <- L.frame[which(L.frame[,1]==n.timestep),]
    L <- as.matrix(reshape2::acast(L.frame, site~spp, value.var = "count"))
    rownames(L) <- rownames(comm)
    col.order <- colnames(comm)
    L <- L[, col.order]
    if(occurrence){
      L <- vegan::decostand(L, method = "pa")
    }
    ent <- vegan::renyi(L, scales = c(1, 2, 12))
    ent.1 <- ent$`1`
    ent.2 <- ent$`2`
    ent.12 <- ent$`12`
    return(list(cbind(sigsq.P, K.P, alpha.P, z0.P), ent.1, ent.2, ent.12))
  }
  f.internal.org <- function(x){
    x <- lapply(lapply(x, function(x) x[1]), function(y) y[[1]])
    d <- cumsum(c(0, sapply(x, nrow)))
    res <- matrix(NA, sum(sapply(x, nrow)), sapply(x, ncol)[1])
    colnames(res) <- colnames(x[[1]])
    rownames(res) <- names(x)
    for(i in 1:length(x)){
      res[(d[i]+1):d[i+1],] <- x[[i]]
    }
    return(res)
  }
  
  samp_OUalpha.values <- sample(OU.alpha.v, size = size, replace = T) # sample from alpha values ranging from 0 to 1
  
  value <- lapply(samp_OUalpha.values, f.internal.run, 
                  comm = comm, subset = subset, occurrence = occurrence, tree = tree,
                  mat.env = mat.env, landscape = landscape, JM = JM, n.timestep = n.timestep, 
                  spp.freq = spp.freq, W.r= W.r, sigma = sigma, theta = theta.sim, root.value = mean(env),  
                  scenario.ID = scenario.ID, sim.ID= sim.ID, output.dir.path = output.dir.path)
  res.par <- f.internal.org(value)
  #ent.1.mat <- matrix(unlist(lapply(lapply(value, function(x) x[2]), function(y) y[[1]])), 
  #                    nrow = length(OU.alpha.v), ncol = nrow(comm),
  #                    dimnames = list(c(paste("alpha", OU.alpha.v, sep = "=")), rownames(comm)), byrow = TRUE)
  #ent.2.mat <- matrix(unlist(lapply(lapply(value, function(x) x[3]), function(y) y[[1]])), 
  #                    nrow = length(OU.alpha.v), ncol = nrow(comm),
  #                    dimnames = list(c(paste("alpha", OU.alpha.v, sep = "=")), rownames(comm)), byrow = TRUE)
  #ent.12.mat <- matrix(unlist(lapply(lapply(value, function(x) x[4]), function(y) y[[1]])), 
  #                     nrow = length(OU.alpha.v), ncol = nrow(comm),
  #                     dimnames = list(c(paste("alpha", OU.alpha.v, sep = "=")), rownames(comm)), byrow = TRUE)
  #return(list(res.par, ent.1.mat, ent.2.mat, ent.12.mat))
  list(res.par, value)
}
