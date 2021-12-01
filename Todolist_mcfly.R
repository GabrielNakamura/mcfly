#To do list - mcfly

# 1. Corrigir Summary.stat -> Distance. Corrigir tb o help.
# 2. Inserir data.frame com coordenadas espaciais, obs.entropy e mean.sim.entropy no output
# 3. substituir Species.Pools por Data.Attributes,  contendo spp.phylogeny, spp.metacommunity, n.sites, max.dist.mst., tempo máximo da árvore (raiz tip)
# 4. calcular density functions e extrair valores
# 5. incluir módulo gráfico, incluindo density functions para os posteriores e a curva dlk x r. e incluir gráfico meia vida. calcular o hpd da meia vida com o log e depois converter


# Defines dispersal limitation in km
  dist.km<-as.dist(geodist::geodist(x=xy.coords,measure = "geodesic"),diag=T,upper=T)/1000
  dist.xy <- scales::rescale(dist.km,c(0,1))
  r <- max(dist.xy*as.dist(ape::mst(dist.xy),diag=T,upper=T))
  dlk <- runif(n=max.sample.size.prior,min=0,max=1)
  prior.w <- round(-log(dlk)/(r^2),3)
