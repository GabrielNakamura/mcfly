#' Internal function
#'
#' Auxiliar function to organize the results from parallel computation
#'
#' @param x List object returned from \code{\link{parallel_PAR}}
#'
#' @return Matrix with results obtained from simulation process
#' 
#'
#' @examples
organize.list<-function(x){
  x<-lapply(lapply(x, function(x) x[1]), function(y) y[[1]])
  d<-cumsum(c(0,sapply(x,nrow)))
  res<-matrix(NA,sum(sapply(x,nrow)),sapply(x,ncol)[1])
  colnames(res)<-colnames(x[[1]])
  rownames(res)<-names(x)
  for(i in 1:length(x)){
    res[(d[i]+1):d[i+1],]<-x[[i]]
  }
  return(res)
}
