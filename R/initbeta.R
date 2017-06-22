#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param dat PARAM_DESCRIPTION
#' @param method PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #init.beta(sim.data[[1]])
#'  }
#' }
#' @seealso 
#'  \code{\link[glmnet]{cv.glmnet}}
#' @rdname init.beta
#' @export 
#' @importFrom glmnet cv.glmnet
#' @importFrom lme4 lmer
init.beta<-function(dat,method=c('glmnet','lme4')){
  
  dat=as.matrix(dat)
  
  y=matrix(dat[,grepl('^y',colnames(dat))],ncol=1)
  X=dat[,grepl('^X',colnames(dat))]
  
  switch(as.numeric(method=='lme4')+1,
         glmnet={
           cvg=glmnet::cv.glmnet(x = X,y = y,alpha = 0)
           cvg$glmnet.fit$beta[,which(cvg$lambda.min==cvg$glmnet.fit$lambda)]         
         },
         lme4={
           Z=dat[,grepl('^Z',colnames(dat))]
           Z.fit = cbind(rep(1,nrow(Z)), Z)
           subject=as.numeric(rownames(dat))
           init.fit = lme4::lmer(y ~ X -1 + ( Z.fit -1 |subject))
           init.fit@beta
         })
  
}
