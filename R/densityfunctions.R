#compute density function
dens.a.prior<-density(x=test.phyllostomidae$Alpha_Prior_Distribution,from=test.phyllostomidae$Alpha_Limits[2],to=test.phyllostomidae$Alpha_Limits[3])
alpha.prior<-dens.a.prior$x
dens.alpha.prior<-dens.a.prior$y
plot(alpha.prior,dens.alpha.prior)

dens.a.post<-density(x=test.phyllostomidae$Alpha_Posterior_Distribution,from=test.phyllostomidae$Alpha_Limits[2],to=test.phyllostomidae$Alpha_Limits[3])
alpha.post<-dens.a.post$x
dens.alpha.post<-dens.a.post$y
plot(alpha.post,dens.alpha.post)
hdi.09.alpha.post<-HDInterval::hdi(dens.a.post,allowSplit=T,0.9)#uncorrected
plot(dens.a.post)

median.alpha.interval.1.uncorr<-quantile(test.phyllostomidae$Alpha_Posterior_Distribution[
  which(test.phyllostomidae$Alpha_Posterior_Distribution>=hdi.09.alpha.post[1,1]&
          test.phyllostomidae$Alpha_Posterior_Distribution<=hdi.09.alpha.post[1,2])],
  probs=0.5)
median.alpha.interval.2.uncorr<-quantile(test.phyllostomidae$Alpha_Posterior_Distribution[
  which(test.phyllostomidae$Alpha_Posterior_Distribution>=hdi.09.alpha.post[2,1]&
          test.phyllostomidae$Alpha_Posterior_Distribution<=hdi.09.alpha.post[2,2])],
  probs=0.5)


dens.hl<-density(x=log(log(2)/test.phyllostomidae$Alpha_Posterior_Distribution),
                 to=log(log(2)/test.phyllostomidae$Alpha_Limits[2]),
                 from=log(log(2)/test.phyllostomidae$Alpha_Limits[3]))
plot(dens.hl)
hdi.09.hl<-HDInterval::hdi(dens.hl,allowSplit=TRUE,0.9) 
b_b<-log(2)/exp(hdi.09.hl[1,1])
b_a<-log(2)/exp(hdi.09.hl[1,2])
a_b<-log(2)/exp(hdi.09.hl[2,1])
a_a<-log(2)/exp(hdi.09.hl[2,2])
corr_hdi.09.alpha<-matrix(c(a_a,b_a,a_b,b_b),2,2,
                          dimnames=list(c(1,2),c("begin","end")))#corrected