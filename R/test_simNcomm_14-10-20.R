Nmeta <- 100
Nspp <- 50
Nspp.phy <- 50
Nspp.comm <- 50
Ncomm <- 50
subset <- FALSE
occurrence <- FALSE
alpha.comm <- 1
noise <- 0 
n.timestep <- 50 
W.r <- 0
theta <- "mean"
sim.ID <- "OU1"
scenario.ID <- "OU1"
sigma<- 0.1
n.timestep<- 10
output.dir.path <- here::here()

sim.Ncomm<-function(scenario.ID,sim.ID,output.dir.path,Nmeta,Nspp.phy,Nspp.comm,Ncomm,subset=TRUE,occurrence=FALSE,
                    alpha.comm,alpha.values,theta="mean",noise=0,W.threshold=0.9,n.timestep,W.r=0,
                    runs,ncores){
  
  alpha.niche.pos<- matrix(NA, nrow= Nmeta, ncol= 1, dimnames= list(1:Nmeta, "Alpha"))
  k.niche.pos<-matrix(NA,nrow=Nmeta,ncol=1,dimnames=list(1:Nmeta,"K.value"))
  k.niche.pos.noise<-matrix(NA,nrow=Nmeta,ncol=1,dimnames=list(1:Nmeta,"K.value"))
  
  W.list<-list()
  R2.list<-list()
  for(m in 1:Nmeta){
    m <- 1
    tree<-geiger::sim.bdtree(b=0.1,d=0,stop="taxa",n=Nspp.phy,t=30,extinct=FALSE)
    env<-runif(Ncomm,0,100)
    env<-scales::rescale(env,c(1,100))
    mat.env<-as.matrix(env)
    JL<-sqrt(round(runif(length(env),9000,10000),0)^2)
    JM<-sum(JL)
    xy.coords<-data.frame(x=runif(length(env),1,100),y=runif(length(env),1,100))
    comm.names<-rownames(as.matrix(xy.coords),FALSE,prefix="comm")
    rownames(mat.env)<-comm.names
    my.landscape<-MCSim::fn.make.landscape(site.coords=xy.coords,Ef=env,
                                           area.m2=1,m=0.5,JM=JM,JL=JL)
    
    #### theta ####
    if(theta=="max"){
      theta<-max(env)
    }
    if(theta=="mean"){
      theta<-mean(env)
    }
    
    threshold<-TRUE
    while(threshold){
      P<-ape::rTraitCont(phy=tree,model="OU",alpha=alpha.comm,
                         sigma=sqrt(sd(env)),theta=theta,root.value=mean(env))
      fitOU.comm<-suppressWarnings(geiger::fitContinuous(tree,P,
                                                         model="OU",bounds=list(alpha=c(0,1)),ncores=4))
      alpha.meta<-fitOU.comm$opt$alpha
      k.niche<-picante::Kcalc(P,tree)
      threshold<-
        if(alpha.comm==0){
          alpha.meta>0.0001|k.niche>1.2|k.niche<0.8
        } else {alpha.meta<(alpha.comm-alpha.comm*0.05)|alpha.meta>(alpha.comm+alpha.comm*0.05)
        }
    }
    
    sigsq<-fitOU.comm$opt$sigsq
    niche.noise<-rnorm(Nspp.phy,0,sd=noise)
    niche.pos.noise<-P+niche.noise
    if(subset==TRUE){
      niche.pos<-sample(niche.pos.noise,size=Nspp.comm)
      niche.pos[order(as.numeric(gsub("s","",names(niche.pos))))]
    } else {niche.pos<-niche.pos.noise}
    
    spp.names<-names(niche.pos)
    tniche.pos<-as.matrix(t(niche.pos))
    alpha.niche.pos[m,]<-alpha.meta
    k.niche.pos[m,]<-picante::Kcalc(P,tree)
    k.niche.pos.noise[m,]<-picante::Kcalc(niche.pos.noise,tree)
    niche.desv<-matrix(1,length(comm.names),length(niche.pos),dimnames=
                         list(comm.names,names(niche.pos)))
    for(j in 1:ncol(tniche.pos)){
      for(i in 1:nrow(mat.env)){
        niche.desv[i,j]=sqrt((mat.env[i,]-tniche.pos[,j])^2)
      }
    }
    niche.bre<-sqrt(max(niche.desv)-(apply(niche.desv,2,sd)))
    spp.freq<-rmultinom(1,size=JM,prob=scales::rescale(rnorm(length(niche.pos),
                                                             0,10),to=c(0.3,1)))
    spp.gamma<-spp.freq/JM
    sim<- MCSim::fn.metaSIM(landscape= my.landscape, scenario.ID= scenario.ID,
                            sim.ID= sim.ID, alpha.fisher= 1, nu= 0, speciation.limit= 0, JM.limit= JM,
                            n.timestep = n.timestep, W.r= W.r, save.sim= FALSE, 
                            output.dir.path= output.dir.path, trait.dispersal= NULL, 
                            trait.dispersal.median= 1, trait.dispersal.range= 0, trait.Ef= niche.pos, 
                            trait.Ef.sd= niche.bre, gamma.abund= spp.gamma, J.t0= NULL, taxa.list= spp.names)
    comm.out<-sim$J.long
    comm.frame<-comm.out[which(comm.out[,1]==n.timestep),]
    comm.obs<-as.matrix(reshape2::acast(comm.frame,site~spp,value.var="count"))
    rownames(comm.obs)<-comm.names
    if(occurrence==TRUE){
      vegan::decostand(comm.obs,"pa")
    } else {comm.obs}
    div<- vegan::renyi(comm.obs, scales= 1)
    ED<- cbind(mat.env,div)
    #theta.sim<-ED[which.max(div),1]
    theta.sim <- sample(ED[which(ED[,2]>=quantile(ED[,2],prob=0.9)),1],1)
    
    #### inicio simulacao dos dados - modelo ABC #####
    
    
    #simcomm<- mcfly::Mcfly(comm = comm.obs, subset = subset, occurrence = occurrence, env = env,
    #                       site.coords = xy.coords, tree = tree, OU.alpha = alpha.values, 
    #                       sigma = sqrt(sd(env)), theta = theta.sim, root.value = mean(env), 
    #                       runs = runs, ncores = ncores, n.timestep = n.timestep, W.r = W.r)
    
    alpha.values <- runif(n = 10000, min = 0, max = 1) # alpha values to be sampled
    
    res_model <- mcfly_RUN(OU.alpha.v = alpha.values, comm = comm, subset = subset, occurrence = occurrence, tree = tree,
              mat.env = mat.env,landscape = my.landscape, JM = JM, spp.freq = spp.freq, sigma = sigma,
              theta = theta, W.r= W.r, root.value = mean(env), n.timestep = n.timestep, 
              scenario.ID = scenario.ID, sim.ID= sim.ID, output.dir.path = output.dir.path) # run simulation with sampled values from OU.alpha
    
    res_model_entropy<- res_model[[-1]] # only entropy values
    res_model_param<- res_model[[1]] # parameters values
    
    ## matrix with entropy values from simulation model
    #sqt_ent<- matrix(
    #  unlist(
    #    lapply(res_model_entropy, 
    #           function(x){
    #             lapply(x[-1], 
    #                    function(y){
    #                      cor(y, ED[,2])
    #               #sum ((y - mean(y)) ^ 2)
    #               }
    #               )
    #             }
    #           )
    #    ),
    #  nrow= runs, ncol= 3, dimnames= list(paste("run", 1:runs), c("ent1", "ent2", "ent5")
    #                                      ),
    #  byrow= T)
    
    cor_model <- unlist(lapply(lapply(res_model_entropy, 
           function(x){
             lapply(x[-1], 
                    function(y){
                      cor(y, ED[,2])
                      #sum ((y - mean(y)) ^ 2)
                    }
             )
           }
    ), function(z) z[[1]]))
    
    # parameters and correlation
    param_cor_data <- data.frame(res_model_param, cor_model)
    
    # extracting only alpha values to ensemble posterior distribution 
    posterior_alpha <- param_cor_data[which(param_cor_data$cor_model >= 0.6), "alpha.P"]
   
    #W.list[[m]]<-simcomm$W
    #R2.list[[m]]<-simcomm$R2
    #print(m)
  }
  # end Nmeta
  
  R2.ent.1.full<-matrix(NA,nrow=Nmeta,ncol=length(alpha.values),dimnames=list(1:Nmeta,paste("alpha",alpha.values,sep="=")))
  R2.ent.2.full<-matrix(NA,nrow=Nmeta,ncol=length(alpha.values),dimnames=list(1:Nmeta,paste("alpha",alpha.values,sep="=")))
  R2.ent.12.full<-matrix(NA,nrow=Nmeta,ncol=length(alpha.values),dimnames=list(1:Nmeta,paste("alpha",alpha.values,sep="=")))
  
  W.ent.1.full<-matrix(NA,nrow=Nmeta,ncol=length(alpha.values),dimnames=list(1:Nmeta,paste("alpha",alpha.values,sep="=")))
  W.ent.2.full<-matrix(NA,nrow=Nmeta,ncol=length(alpha.values),dimnames=list(1:Nmeta,paste("alpha",alpha.values,sep="=")))
  W.ent.12.full<-matrix(NA,nrow=Nmeta,ncol=length(alpha.values),dimnames=list(1:Nmeta,paste("alpha",alpha.values,sep="=")))
  
  for(m in 1:Nmeta){
    R2.ent.1.full[m,]<-R2.list[[m]][,1]
    R2.ent.2.full[m,]<-R2.list[[m]][,2]
    R2.ent.12.full[m,]<-R2.list[[m]][,3]
    
    
    W.ent.1.full[m,]<-W.list[[m]][,1]
    W.ent.2.full[m,]<-W.list[[m]][,2]
    W.ent.12.full[m,]<-W.list[[m]][,3]
  }
  
  clean.W.ent.1.full<-ifelse(W.ent.1.full>=W.threshold,1,0)
  clean.W.ent.2.full<-ifelse(W.ent.2.full>=W.threshold,1,0)
  clean.W.ent.12.full<-ifelse(W.ent.12.full>=W.threshold,1,0)
  
  
  power<-matrix(NA, nrow=3, ncol=length(alpha.values), dimnames=list(c("Ent.1","Ent.2","Ent.12"), 
                                                                     c(paste("alpha", alpha.values, sep="=")
                                                                     )
  )
  )
  power[1,]<-colSums(clean.W.ent.1.full)
  power[2,]<-colSums(clean.W.ent.2.full)
  power[3,]<-colSums(clean.W.ent.12.full)
  
  prop.power<-matrix(NA, nrow= 3, ncol=length(alpha.values), 
                     dimnames= list(c("Ent.1", "Ent.2", "Ent.12"), c(paste("alpha", alpha.values, sep="=")
                     )
                     )
  )
  prop.power[1,]<-power[1,]/rowSums(power)[1]
  prop.power[2,]<-power[2,]/rowSums(power)[2]
  prop.power[3,]<-power[3,]/rowSums(power)[3]
  
  R2.ent.1.mean<-matrix(NA,nrow=length(alpha.values),ncol=Nmeta,dimnames=list(c(paste("alpha",alpha.values,sep="=")),1:Nmeta))
  R2.ent.2.mean<-matrix(NA,nrow=length(alpha.values),ncol=Nmeta,dimnames=list(c(paste("alpha",alpha.values,sep="=")),1:Nmeta))
  R2.ent.12.mean<-matrix(NA,nrow=length(alpha.values),ncol=Nmeta,dimnames=list(c(paste("alpha",alpha.values,sep="=")),1:Nmeta))
  
  W.ent.1.mean<-matrix(NA,nrow=length(alpha.values),ncol=Nmeta,dimnames=list(c(paste("alpha",alpha.values,sep="=")),1:Nmeta))
  W.ent.2.mean<-matrix(NA,nrow=length(alpha.values),ncol=Nmeta,dimnames=list(c(paste("alpha",alpha.values,sep="=")),1:Nmeta))
  W.ent.12.mean<-matrix(NA,nrow=length(alpha.values),ncol=Nmeta,dimnames=list(c(paste("alpha",alpha.values,sep="=")),1:Nmeta))
  
  
  for(m in 1:Nmeta){
    R2.ent.1.mean[,m]<-apply(R2.ent.1.full[sample(1:nrow(R2.ent.1.full),size=Nmeta,replace=TRUE), drop = FALSE,],2,mean)
    R2.ent.2.mean[,m]<-apply(R2.ent.2.full[sample(1:nrow(R2.ent.2.full),size=Nmeta,replace=TRUE), drop = FALSE,],2,mean)
    R2.ent.12.mean[,m]<-apply(R2.ent.12.full[sample(1:nrow(R2.ent.12.full),size=Nmeta,replace=TRUE), drop = FALSE,],2,mean)
    
    W.ent.1.mean[,m]<-apply(W.ent.1.full[sample(1:nrow(W.ent.1.full),size=Nmeta,replace=TRUE), drop = FALSE,],2,mean)
    W.ent.2.mean[,m]<-apply(W.ent.2.full[sample(1:nrow(W.ent.2.full),size=Nmeta,replace=TRUE), drop = FALSE,],2,mean)
    W.ent.12.mean[,m]<-apply(W.ent.12.full[sample(1:nrow(W.ent.12.full),size=Nmeta,replace=TRUE), drop = FALSE,],2,mean)
  }
  
  mean.R2.ent.1.boot<-matrix(NA,nrow=3,ncol=length(alpha.values),dimnames=
                               list(c("lower","mean","upper"),c(paste("alpha",alpha.values,sep="="))))
  mean.R2.ent.2.boot<-matrix(NA,nrow=3,ncol=length(alpha.values),dimnames=
                               list(c("lower","mean","upper"),c(paste("alpha",alpha.values,sep="="))))
  mean.R2.ent.12.boot<-matrix(NA,nrow=3,ncol=length(alpha.values),dimnames=
                                list(c("lower","mean","upper"),c(paste("alpha",alpha.values,sep="="))))
  
  
  mean.R2.ent.1.boot[1,]<-apply(R2.ent.1.mean,1,function(x)quantile(x,probs = 0.025))
  mean.R2.ent.1.boot[2,]<-apply(R2.ent.1.mean,1,function(x)quantile(x,probs = 0.5))
  mean.R2.ent.1.boot[3,]<-apply(R2.ent.1.mean,1,function(x)quantile(x,probs = 0.975))
  
  mean.R2.ent.2.boot[1,]<-apply(R2.ent.2.mean,1,function(x)quantile(x,probs = 0.025))
  mean.R2.ent.2.boot[2,]<-apply(R2.ent.2.mean,1,function(x)quantile(x,probs = 0.5))
  mean.R2.ent.2.boot[3,]<-apply(R2.ent.2.mean,1,function(x)quantile(x,probs = 0.975))
  
  mean.R2.ent.12.boot[1,]<-apply(R2.ent.12.mean,1,function(x)quantile(x,probs = 0.025))
  mean.R2.ent.12.boot[2,]<-apply(R2.ent.12.mean,1,function(x)quantile(x,probs = 0.5))
  mean.R2.ent.12.boot[3,]<-apply(R2.ent.12.mean,1,function(x)quantile(x,probs = 0.975))
  
  
  mean.W.ent.1.boot<-matrix(NA,nrow=3,ncol=length(alpha.values),dimnames=
                              list(c("lower","mean","upper"),c(paste("alpha",alpha.values,sep="="))))
  mean.W.ent.2.boot<-matrix(NA,nrow=3,ncol=length(alpha.values),dimnames=
                              list(c("lower","mean","upper"),c(paste("alpha",alpha.values,sep="="))))
  mean.W.ent.12.boot<-matrix(NA,nrow=3,ncol=length(alpha.values),dimnames=
                               list(c("lower","mean","upper"),c(paste("alpha",alpha.values,sep="="))))
  
  
  mean.W.ent.1.boot[1,]<-apply(W.ent.1.mean,1,function(x)quantile(x,probs = 0.025))
  mean.W.ent.1.boot[2,]<-apply(W.ent.1.mean,1,function(x)quantile(x,probs = 0.5))
  mean.W.ent.1.boot[3,]<-apply(W.ent.1.mean,1,function(x)quantile(x,probs = 0.975))
  
  mean.W.ent.2.boot[1,]<-apply(W.ent.2.mean,1,function(x)quantile(x,probs = 0.025))
  mean.W.ent.2.boot[2,]<-apply(W.ent.2.mean,1,function(x)quantile(x,probs = 0.5))
  mean.W.ent.2.boot[3,]<-apply(W.ent.2.mean,1,function(x)quantile(x,probs = 0.975))
  
  mean.W.ent.12.boot[1,]<-apply(W.ent.12.mean,1,function(x)quantile(x,probs = 0.025))
  mean.W.ent.12.boot[2,]<-apply(W.ent.12.mean,1,function(x)quantile(x,probs = 0.5))
  mean.W.ent.12.boot[3,]<-apply(W.ent.12.mean,1,function(x)quantile(x,probs = 0.975))
  
  
  Res<-list(K.niche.position=simcomm$K.niche.position,
            Alpha.niche.position=simcomm$Alpha.niche.position,
            Sigsq.niche.position=simcomm$Sigsq.niche.position,
            Z0.niche.position=simcomm$Z0.niche.position,Alpha.metacomm=alpha.niche.pos,
            K.Metacomm=k.niche.pos,
            K.Metacomm.noise=k.niche.pos.noise,R2.ent.1.CI95=mean.R2.ent.1.boot,
            R2.ent.2.CI95=mean.R2.ent.2.boot,R2.ent.12.CI95=mean.R2.ent.12.boot,
            W.ent.1.CI95=mean.W.ent.1.boot,W.ent.2.CI95=mean.W.ent.2.boot,W.ent.12.CI95=mean.W.ent.12.boot,Power=power,Prop.Power=prop.power)
  return(Res)
}
