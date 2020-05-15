parallel::makeCluster()

n<- c(1, 2, 3)
CL<- cl <- parallel::makeCluster(2, setup_timeout = 0.5)


devtools::install_github("GabrielNakamura/mcfly", force = T, build_vignettes = T)
1
Yes
library(mcfly)
vignette("McFly_vignette")

data("Furnariidae")
data("phylo_Furnariidae")
data("envir")

sub.tree<- ape::subtrees(phylo_Furnariidae) #splitting internal nodes in the phylogeny
n.nodes<-length(sub.tree)
nspp.nodes<-matrix(NA,n.nodes,1)
for(i in 1:n.nodes){
  nspp.nodes[i,]<-sub.tree[[i]]$Ntip #number of species in each internal node
}

dim.comm.subset<-matrix(NA,nspp.nodes,2)
for(j in 1:n.nodes){
  sub.phy<-sub.tree[[j]]
  comm.sub.phy<-Furnariidae[,sub.phy$tip.label]
  zero.row<-which(rowSums(comm.sub.phy)==0)
  comm.subset<-comm.sub.phy[-zero.row,]
  dim.comm.subset[j,]<-dim(comm.subset) # number of sites showing occurrences of species belonging to each internal node
}


plot(sub.tree[[49]]) #Synallaxinae subfamily
sub.phy<- sub.tree[[49]]
comm.sub.phy<- Furnariidae[, sub.tree[[49]]$tip.label] #selecting only sites that present occurrences of species of Synallaxinae
zero.row<- which(rowSums(comm.sub.phy) == 0)
comm.subset<- comm.sub.phy[-zero.row, ] #community matrix with occurrences of Synalaxynae species
coords.subset<- as.matrix(envir[rownames(comm.subset), c(1,2)]) #extracting coordinates

env.subset<-envir[rownames(comm.subset),]
envir<- scales::rescale(as.matrix(env.subset$Lat),c(1,100)) #scaling the environmental variable to vary between 1 and 100
sigma<- sd(envir) #standard deviation of environmental variable
root.value<- mean(envir) #mean value
div<- vegan::renyi(comm.subset,scales=1) #diversity measure
mod<- lm(div~envir) #linear model
pred<- predict.lm(mod) #predicted relationship among diversity and the environmental variable
ED<- cbind(envir,pred)
theta<- ED[which.max(pred),1] #optimum niche value

OU.alpha<- c(0,0.05,0.25,1)


ncores<- parallel::detectCores() - 1

test.synallaxinae<- mcfly::Mcfly(comm= comm.subset, subset= FALSE, occurrence= TRUE, env= envir, site.coords= coords.subset, tree= sub.phy, OU.alpha= OU.alpha, sigma= sigma, theta= theta, root.value= root.value, runs= 100, ncores= 3, area.m2= 1, m= 0.5, JM= sum(comm.subset), JM.limit= JM,JL= rowSums(comm.subset),nu= 0, speciation.limit= 0, n.timestep= 50,W.r= 0, scenario.ID= "species.sorting", sim.ID= "data", output.dir.path= "OUTPUT_DATA")


