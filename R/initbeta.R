#' @title Evaluate fixed effects initial values for lmmen
#' @description Evaluate fixed effects initial values for lmmen via cv.glmnet or lme4. 
#' @param dat data.frame, data to solve initial values
#' @param method character, method to use, c('glmnet','lme4')
#' @return numeric
#' @details cv.glmnet is set to ridge regression. 
#' @examples 
#'  dat <- initialize_example(n.i = 5,n = 30,q=4,seed=1)
#'  init.beta(dat,method='glmnet')
#'  init.beta(dat,method='lme4')
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
         }
         )
  
}
