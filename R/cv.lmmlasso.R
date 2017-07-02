#' @title Cross Validation for lmmlasso package
#' @description Cross Validation for lmmlasso package as shown in example xxx
#' @param data matrix, containing y,X,Z and subject variables
#' @param lambda numeric, path of positive regularization parameter, Default: seq(0, 500, 5)
#' @param ... parameters to pass to lmmlasso
#' @return numeric, vector that includes values of estimated fixed, random and bic
#' @examples 
#'  cv.lmmlasso(initialize_example(seed=1))
#' @rdname cv.lmmlasso
#' @importFrom lmmlasso lmmlasso
#' @export 

cv.lmmlasso<-function(data,lambda=seq(0,500,5),...){
  max.iter=length(lambda)
  data=as.matrix(data)
  y=matrix(data[,grepl('^y',colnames(data))],ncol=1)
  X=cbind(rep(1,nrow(data)),data[,grepl('^X',colnames(data))])
  Z=data[,grepl('^Z',colnames(data))]
  grp=factor(row.names(data))
  if(!'pdMat'%in%names(match.call()[3])) pdMat="pdSym"
  BIC_vec<-BIC_DIFF<-BIC_DIFF_I<-Inf
  i=0
  #pb <- txtProgressBar(min = 0,max = max.iter,initial = 0,style = 3)
  for(i in 1:max.iter){
    #setTxtProgressBar(pb, i)
    capture.output({
      suppressWarnings({object <- lmmlasso::lmmlasso(x=X,y=y,z=Z,grp=grp,lambda=lambda[i],pdMat = pdMat,...)})
    })
    out<-c(object$coefficients,sqrt(diag(object$Psi)), sigma=object$sigma,bic.opt=object$bic,lambda.opt=lambda[i])
    BIC_vec<-c(BIC_vec,object$bic)
    if(i>1){
      BIC_DIFF=c(BIC_DIFF,BIC_vec[i]-BIC_vec[i-1])
      if(abs(BIC_DIFF[i])<1e-4) break 
    }
  }
  #close(pb)
  out
}
