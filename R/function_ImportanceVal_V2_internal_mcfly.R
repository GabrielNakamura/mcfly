#' @title Importance Values
#'
#'
#' @importFrom stats prcomp setNames
#' @param matrix.M Matrix object containing diversity metric values for each community.
#' @param scale Logical argument. If TRUE (default) the diversity metrics in matrix M are standardized. 
#' @param method Character. The method of standardization to be used in matrix M.
#' @param stopRule Logical. If TRUE a criterion for select eigenvalues will be used.
#' 
ImportanceVal<- function(matrix.M, scale= TRUE, method= "max", stopRule= TRUE){
  if(is.matrix(matrix.M) == FALSE){
    matrix.M<- as.matrix(matrix.M)
    if(ncol(matrix.M)<3){
      stop("\n matrix M must be at least 3 components of diversity\n")
    }
    if(nrow(matrix.M)<3){
      stop("\n Matrix M must be at least 3 communities\n")
    }
  } 
  matrix.M.stand<-vegan::decostand(x = matrix.M, method = method, MARGIN = 2)[1:nrow(matrix.M),]
  if(scale == TRUE){
    metric.sqrt.corr<- (prcomp(x = matrix.M.stand, scale.= FALSE)$rotation ^ 2)
    prop.var<- summary(prcomp(x = matrix.M.stand, scale.= FALSE))$importance[2,]
    IVs.resul<- matrix(nrow= ncol(metric.sqrt.corr), ncol= ncol(matrix.M), dimnames= list(paste("PC",1:ncol(metric.sqrt.corr)), colnames(matrix.M)))
    for(i in 1:nrow(metric.sqrt.corr)){
      IVs.resul[,i]<- metric.sqrt.corr[i,] * as.matrix(prop.var)
    }
    IVs.proportion<-matrix(nrow= ncol(metric.sqrt.corr), ncol= ncol(matrix.M), dimnames= list(paste("PC",1:ncol(metric.sqrt.corr)), colnames(matrix.M)))
    for(i in 1:nrow(IVs.resul)){
      IVs.proportion[i,]<-IVs.resul[i,]/prop.var[i]
    }
    if(stopRule==TRUE){
      sig.eigen<-which(prcomp(matrix.M.stand, scale. = FALSE)$sdev^2>mean(prcomp(matrix.M.stand, scale. = FALSE)$sdev^2))
      IVs.resul.sig<-IVs.resul[sig.eigen,]
      IV.resul.sig<- setNames(list(IVs.resul.sig, prop.var, metric.sqrt.corr), c("IV.obs_stopRule", "Var.by.axis","Metrics_correlation"))
      return(IV.resul.sig)
    } else{
      IV.resul<- setNames(list(IVs.resul, prop.var, metric.sqrt.corr), c("IV.obs", "Var.by.axis","Metrics_correlation"))
      return(IV.resul)
    }
  }
  if(scale == FALSE){
    metric.sqrt.corr<- (prcomp(x = matrix.M, scale.= FALSE)$rotation ^ 2)
    prop.var<- summary(prcomp(x = matrix.M, scale.= FALSE))$importance[2,]
    names.matrix.IV<- list("IV.resu", colnames(matrix.M))
    IVs.resul<- matrix(nrow= ncol(metric.sqrt.corr), ncol= ncol(matrix.M), dimnames= names.matrix.IV)
    for(i in 1:nrow(metric.sqrt.corr)){
      IVs.resul[,i]<- metric.sqrt.corr[i,] * as.matrix(prop.var)
    }
    if(stopRule==TRUE){
      sig.eigen<-which(prcomp(matrix.M, scale. = FALSE)$sdev^2>mean(prcomp(matrix.M, scale. = FALSE)$sdev^2))
      IVs.resul.sig<-IVs.resul[sig.eigen,]
      IV.resul.sig<- setNames(list(IVs.resul.sig, prop.var, metric.sqrt.corr), c("IV.obs_stopRule", "Var.by.axis","Metrics_correlation"))
      return(IV.resul.sig)
    } else{
      IV.resul<- setNames(list(IVs.resul, prop.var, metric.sqrt.corr), c("IV.obs_stopRule", "Var.by.axis","Metrics_correlation"))
      return(IV.resul)
    }
  }
}