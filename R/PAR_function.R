#' parallel_PAR
#' 
#' Auxiliar function to allow parallel computations
#'
#' @param OU.alpha.v Numeric vector containing values of alpha parameter to be tested
#' @param comm Matrix with species occurence or abundances, sampling units in rows and species in columns
#' @param subset Logical, if TRUE only a subset of species will be used in caculations
#' @param occurrence Logical, if TRUE community matrix will be converted to species occurence
#' @param tree Phylogenetic tree in newick format
#' @param mat.env Matrix containing values for environmental gradient
#' @param landscape Landscape where the simulations will go on. Object generated with \code{\link{make.landscape}} from MCSim package
#' @param JM Numeric value indicating the total number of individuals to include in \code{\link{make.landscape}} from MCSim package
#' @param n.timestep Numeric value indicating the number of generations in the simulation process performed with \code{\link{metasim}} from MCSim package
#' @param spp.freq Numeric vector of regional abundance, used in \code{\link{metasim}} from MCSim package
#' @param W.r Numeric indicating the dispersal kernel slope in \code{\link{metasim}} from MCSim package
#' @param sigma Numeric vector indicating the standard-dviation of niche position. Argument used in \code{\link{rTraitCont}}
#' @param theta Numeric value used in  \code{\link{rTraitCont}} to indicate the optimun for each branch
#' @param root.value Numeric value indicating the trait value at the root to be used in \code{\link{rTraitCont}}
#' @param scenario.ID Character setting the name for the simulation scenario. Used in \code{\link{metasim}}
#' @param sim.ID Character string indicating the name for a particular simulation. Used in \code{\link{rTraitCont}}
#' @param output.dir.path Character string indicating the name of directory to save simulation results. Used in \code{\link{rTraitCont}}
#' @param reps tricky argument to work in parallel process
#'
#' @export
#'
#' @return A list with parameters from simulated communities
#'
#' @examples
parallel_PAR<- function(OU.alpha.v,comm,subset,occurrence,tree,mat.env,landscape,JM,n.timestep,spp.freq,W.r,
                        sigma,theta,root.value,scenario.ID,sim.ID,output.dir.path, reps= 5){
  value<- lapply(OU.alpha.v, function(x){
    output.dir.path<-paste(output.dir.path, x, sep = ".")
    thresh<-TRUE
    while(thresh){
      P<-ape::rTraitCont(phy= tree, model= "OU", sigma= sigma,alpha= x,
                         theta=theta,root.value=root.value)
      fitOU<-geiger::fitContinuous(tree,P,model="OU",bounds=list(alpha=c(min=exp(-500),max=exp(1))),ncores=1)
      k.niche<-picante::Kcalc(P,tree,checkdata=F)
      alpha.niche<-fitOU$opt$alpha
      thresh<-
        if(x==0){
          alpha.niche>0.0001|k.niche>1.2|k.niche<0.8
        } else {alpha.niche<(x - x*0.05)|alpha.niche>(x + x*0.05)
        }
    }
    if(subset==TRUE){
      subset.spp<-colnames(comm)
      niche<-P[subset.spp]
    } else{niche<-P}
    
    sigsq.P<-fitOU$opt$sigsq
    K.P<-as.numeric(picante::Kcalc(P,tree,checkdata=F))
    alpha.P<-alpha.niche
    z0.P<-fitOU$opt$z0
    t.P<-as.matrix(t(niche))
    NB<-matrix(NA, nrow(mat.env), ncol(t.P))
    for(l in 1:ncol(t.P)){
      for(k in 1:nrow(mat.env)){
        NB[k,l]=sqrt((mat.env[k,]-(t.P[,l]))^2)
      }
    }
    sd.NB<-sqrt(max(NB)-apply(NB,2, stats::sd))
    meta<- MCSim::metasim(landscape=landscape,scenario.ID=scenario.ID,
                          sim.ID=sim.ID,alpha.fisher=1,nu=0,speciation.limit=0,
                          JM.limit=JM,n.timestep=n.timestep,W.r=W.r,save.sim=FALSE,
                          output.dir.path=output.dir.path,trait.dispersal=NULL,
                          trait.dispersal.median=1,trait.dispersal.range=0,
                          trait.Ef=niche,trait.Ef.sd=sd.NB,gamma.abund=spp.freq,
                          J.t0=NULL,taxa.list=colnames(comm))
    L.frame<-meta$J.long
    L.frame<-L.frame[which(L.frame[,1]==n.timestep),]
    L<-as.matrix(reshape2::acast(L.frame,site~spp,value.var="count"))
    rownames(L)<-rownames(comm)
    col.order<-colnames(comm)
    L<-L[,col.order]
    if(occurrence==TRUE){
      vegan::decostand(L,"pa")
    } else {L}
    ent<-vegan::renyi(L,scales=c(1,2,12))
    ent.1<-ent$`1`
    ent.2<-ent$`2`
    ent.12<-ent$`12`
    return(list(cbind(sigsq.P,K.P,alpha.P,z0.P),ent.1,ent.2,ent.12))
  }
  )
  res.par <- organize.list(value)
  ent.1.mat<-matrix(unlist(lapply(lapply(value, function(x) x[2]), function(y) y[[1]])), nrow= length(OU.alpha.v), ncol= nrow(comm),dimnames=
                      list(c(paste("alpha",OU.alpha.v,sep="=")),rownames(comm)), byrow= T)
  ent.2.mat<-matrix(unlist(lapply(lapply(value, function(x) x[3]), function(y) y[[1]])), nrow= length(OU.alpha.v), ncol= nrow(comm),dimnames=
                      list(c(paste("alpha",OU.alpha.v,sep="=")),rownames(comm)), byrow= T)
  ent.12.mat<-matrix(unlist(lapply(lapply(value, function(x) x[4]), function(y) y[[1]])), nrow= length(OU.alpha.v), ncol= nrow(comm),dimnames=
                       list(c(paste("alpha",OU.alpha.v,sep="=")),rownames(comm)), byrow= T)
  return(list(res.par, ent.1.mat, ent.2.mat, ent.12.mat))
}
