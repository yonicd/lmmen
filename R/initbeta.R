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
# @importFrom nmle glmmPQL
init.beta<-function(dat,method=c('glmnet','lme4','glmmPQL'),...){
  
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
         #, glmmPQL={
         #   dat<-scale(dat,center=T,scale=T)
         #   PQL<-nmle::glmmPQL(data=dat,...)
         #   Delta.start<-c(as.numeric(PQL$coef$fixed),rep(0,6),as.numeric(t(PQL$coef$random$team)))
         #   Q.start<-as.numeric(VarCorr(PQL)[1,1])
         #   
         # }
         )
  
}
