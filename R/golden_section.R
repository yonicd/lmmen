#' @title Golden section grid search on a lmmen penalty
#' @description Solve for local minimum with one dimensional golden 
#' section on a L1 lmmen penalty.
#' @param dat matrix, matrix that includes y (response),X (population covariates),
#' Z (random effects covariates (not incl random intercept))
#' @param init.beta numeric, initial fixed effects estimates
#' @param pen.effect character,which penalty to search on 
#' c('FE.L1','RE.L1',FE.L2'','RE.L2'), Default: 'FE.L1'
#' @param opt.lb numeric, start of search interval, Default: 0
#' @param opt.ub numeric, end of search interval, Default: 1
#' @param opt.maxiter numeric, maximum iterations to search, Default: 100
#' @param opt.tol numeric, accuracy value, Default: 0.1
#' @param opt.tau numeric, golden proportion coefficient (~0.618) Default: (sqrt(5) - 1)/2
#' @return lmmen list object inluding lmmen fit object of min BIC solution
#'  and summary statistics from the grid searc
#' @examples
#' \dontrun{ 
#'  dat <- initialize_example(n.i = 5,n = 30,q=4,seed=1)
#'  init <- init.beta(dat,method='glmnet')
#'  golden_section(dat,init,pen.effect = 'FE.L1')
#'  }
#' @rdname golden_section
#' @export 

golden_section = function(dat, 
                          init.beta, 
                          pen.effect='FE.L1', 
                          opt.lb=0, 
                          opt.ub=1, 
                          opt.maxiter=100, 
                          opt.tol=0.1,
                          opt.tau=(sqrt(5)-1)/2)
{

gs.k=0    # number of iterations
gs.x1=opt.lb+(1-opt.tau)*(opt.ub-opt.lb)             # computing x values
gs.x2=opt.lb+opt.tau*(opt.ub-opt.lb)

fracVec=function(x,pen.effect){
  vout <- rep(1,4) 
  names(vout)<-c('FE.L1','RE.L1','FE.L2','RE.L2')
  vout[pen.effect] <- x
  vout
}

gs.fit1.0=gs.fit1=lmmen(dat, init.beta = init.beta, frac = fracVec(gs.x1,pen.effect))
gs.fit2.0=gs.fit2=lmmen(dat, init.beta = init.beta, frac = fracVec(gs.x2,pen.effect))

gs.f_x1=gs.fit1$BIC
gs.f_x2=gs.fit2$BIC
opt.ubIC.min=min(c(gs.f_x1,gs.f_x2))
gs.x.iter=gs.x12=c(gs.x1,gs.x2)
gs.fx.iter=c(gs.f_x1,gs.f_x2)
gs.x=gs.x12[which.min(c(gs.f_x1,gs.f_x2))]

while ((abs(opt.ub-opt.lb)>opt.tol)*(gs.k<opt.maxiter)){

  if(gs.f_x1<gs.f_x2)
    {
      cat('.')
      opt.ub=gs.x2
      gs.x1=opt.lb+(1-opt.tau)*(opt.ub-opt.lb)
      gs.x2=opt.lb+opt.tau*(opt.ub-opt.lb)
      
      gs.fit2=gs.fit1
      gs.fit1=lmmen(dat, init.beta = init.beta, frac = fracVec(gs.x1,pen.effect))
    }else{
      cat('.')
      opt.lb=gs.x1
      gs.x1=opt.lb+(1-opt.tau)*(opt.ub-opt.lb)
      gs.x2=opt.lb+opt.tau*(opt.ub-opt.lb)
      gs.fit1=gs.fit2
      gs.fit2=lmmen(dat, init.beta = init.beta, frac = fracVec(gs.x2,pen.effect))
    }
  
  gs.x.iter=rbind(gs.x.iter,c(gs.x1,gs.x2))
  gs.f_x1=gs.fit1$BIC  
  gs.f_x2=gs.fit2$BIC
  
  
  gs.fx.iter=rbind(gs.fx.iter,c(gs.f_x1,gs.f_x2))
  
  gs.k=gs.k+1
  opt.ubIC.min=c(opt.ubIC.min,min(c(gs.f_x1,gs.f_x2)))
  gs.fit12=list(gs.fit1,gs.fit2)
  gs.fit=gs.fit12[[which.min(c(gs.f_x1,gs.f_x2))]]
  gs.x12=c(gs.x1,gs.x2)
  gs.x=c(gs.x,gs.x12[which.min(c(gs.f_x1,gs.f_x2))])
}

names(gs.fit$fixed) <- names(init.beta)

l.out <- list(fit=gs.fit,
              BIC.path=opt.ubIC.min,
              frac.path=gs.x,
              iter.star=gs.k)

l.out <- structure(l.out,class=c('lmmen','list'))
cat('\n')
return(l.out)
}