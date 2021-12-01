r<-round(seq(0,1, length.out = 100),3)
w<-c(0,1,5,15,100)


dlk<-matrix(NA,length(r),length(w),dimnames=list(r,w))
for (i in 1:length(w)){
  for (p in 1:length(r)){
    dlk[p,i]<-round(exp(1)^-(w[i]*r[p]^2),3)
  }
}


2.71^-(0*0.5^2)
plot(r,dlk[,1])
plot(r,dlk[,2])
plot(r,dlk[,3])
plot(r,dlk[,4])
plot(r,dlk[,5])
