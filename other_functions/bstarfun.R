unwrap=function(g,coef="fixed"){
  nr=length(g)
  k=which(coef==names(g[[1]]))
  nc=length(g[[1]][[k]])
  out=data.frame(iter=rep(seq(1,nr),times=rep(nc,nr)),type=rep(coef,nr*nc),id=rep(seq(1,nc),nr),
                 ldply(g, function(meas) data.frame(est=meas[[k]],rep.row(meas$frac,length(meas[[k]])),check.names=F))
  )
  names(out)[-c(1:4)]=colnames(g[[1]]$frac)
  return(out)
}

glmmen.cv=function(gin1,seqmat){
  ldply(gin1,function(gi){
    din=list(Y=gi$Y,X=gi$X,Z=gi$Z,subject=gi$subject)
    init.beta=gi$init.beta
    gout=foreach(seqvec=iter(seqmat, by='row'),.packages=c("quadprog","mvtnorm"),
                 .export=c("GLMMEN4"))%dopar% {
                   GLMMEN4(din, init.beta, frac =seqvec, eps = 10^(-2))}
    Coefs=rbind(unwrap(gout,"fixed"),unwrap(gout,"stddev"),unwrap(gout,"BIC"))
    return(Coefs)
  },
  .progress = "text")}

bstar.fun=function(prdata,seqmat,type="glmnet"){
  rmX=c('TgjUserId__1','y','v','TgjUserId','ta','gdMean',
        'value','Ifpid','inst_type','created','activated',
        'closed','UserQuesId','agegroup')
  
  Y=as.matrix(prdata%>%select(value))
  X=as.matrix(prdata%>%
                select_(.dots = -rmX))[,(1:40)]
  
  switch(type,
          glmnet={
            cvg=cv.glmnet(x = X,y = Y,parallel = T)
            out=cvg$glmnet.fit$beta[,which(cvg$lambda.min==cvg$glmnet.fit$lambda)]
            out
            },
          lmmnet={
            Z=as.matrix(prdata%.%select(starts_with("gr"))%.%select(-grit))
            subject=as.matrix(prdata$subject)
            cvg=cv.glmnet(x = X,y = Y,parallel = T,alpha = 0)
            init.beta=cvg$glmnet.fit$beta[,which(cvg$lambda.min==cvg$glmnet.fit$lambda)]
            din=list(Y=Y,X=X,Z=Z,subject=subject)
            gout=foreach(seqvec=iter(seqmat, by='row'),.packages=c("quadprog","mvtnorm"),
                         .export=c("GLMMEN4"))%dopar% {
                           GLMMEN4(din, init.beta, frac =seqvec, eps = 10^(-2))}
            out=rbind(unwrap(gout,"fixed"),unwrap(gout,"stddev"),unwrap(gout,"BIC"))
            out
            },
          penlme={
            Z.0=as.matrix(prdata%.%select(starts_with("gr"))%.%select(-grit))
            subject=as.matrix(prdata$subject)
            n.i = tabulate(prdata$subject)
            n.tot=sum(n.i)
            Z = as.matrix(Z.0,nrow=n.tot)
            Z.fit = cbind(rep(1,n.tot), Z)
            init.fit = lmer(Y ~ X -1 + (Z.fit -1 |subject),
                            control=lmerControl(optCtrl=list(maxfun=30)))
            est = VarCorr(init.fit)
            beta.hat=as.matrix(fixef(init.fit))
            D.lme = as.matrix(est$subject)/((attributes(est)$sc)^2)
            
            init.fit=list(est=est,beta.hat=beta.hat,D.lme=D.lme)
            gout=foreach(tf=iter(seqmat),.packages=c("quadprog","mvtnorm","MASS","lme4"),
                         .export=c("Pen.LME"))%dopar%{Pen.LME(prdata,init.fit,tf,eps = 10^(-3))}
            t1=list(fixed=data.frame(frac=t.fracs,ldply(gout, function(meas){meas$fixed}),check.names=F,row.names=NULL),
                    stddev=data.frame(frac=t.fracs,ldply(gout, function(meas){meas$stddev}),check.names=F,row.names=NULL),
                    bic=data.frame(frac=t.fracs,ldply(gout, function(meas){meas$BIC}),check.names=F,row.names=NULL))
            out=ldply(t1,function(df){reshape(df,direction="long",varying=list(names(df)[-1]),v.names="est",idvar=c("frac"),timevar="id")})
            names(out)[1]="type"
            out}
         )
}