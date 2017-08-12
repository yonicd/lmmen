#' @title Cross Validation for lmmlasso package
#' @description Cross Validation for lmmlasso package as shown in example xxx
#' @param data matrix, containing y,X,Z and subject variables
#' @param lambda numeric, path of positive regularization parameter, Default: seq(0, 500, 5)
#' @param ... parameters to pass to lmmlasso
#' @return numeric, vector that includes values of estimated fixed, random and bic
#' @examples 
#' \donttest{
#'  cv.lmmlasso(initialize_example(seed=1))
#'  }
#' @rdname cv.lmmlasso
#' @importFrom lmmlasso lmmlasso
#' @importFrom utils capture.output
#' @export 
cv.lmmlasso<-function(data, lambda=seq(0, 500, 5), ...){
  max.iter=length(lambda)
  data=as.matrix(data)
  y=matrix(data[,grepl('^y',colnames(data))],ncol=1)
  X=cbind(rep(1,nrow(data)),data[,grepl('^X',colnames(data))])
  Z=cbind(rep(1,nrow(data)),data[,grepl('^Z',colnames(data))])
  grp=factor(row.names(data))
  if(!'pdMat'%in%names(match.call()[3])) pdMat="pdSym"
  BIC_vec<-BIC_DIFF<-BIC_DIFF_I<-Inf
  i=0

  for(i in 1:max.iter){

    utils::capture.output({
      suppressWarnings({object <- lmmlasso::lmmlasso(x=X,y=y,z=Z,grp=grp,lambda=lambda[i],pdMat = pdMat,...)})
    })
    out<-object
    BIC_vec<-c(BIC_vec,object$bic)
    if(i>1){
      BIC_DIFF=c(BIC_DIFF,BIC_vec[i]-BIC_vec[i-1])
      if(abs(BIC_DIFF[i])<1e-4) break 
    }
  }

  out
}
