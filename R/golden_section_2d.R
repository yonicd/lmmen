#' @title Golden section two dimensional grid search on L1 lmmen penalties
#' @description Solve for local minimum with two dimensional golden 
#' section on L1 lmmen penalties.
#' @param dat matrix, matrix that includes y (response),X (population covariates),
#' Z (random effects covariates (not incl random intercept))
#' @param init.beta numeric, initial fixed effects estimates
#' @param l2 numeric, L2 penalty levels Default: c(1, 1)
#' @param opt.lb numeric, start of interval for L1 fixed and L1 random effects, Default: c(0, 0)
#' @param opt.ub numeric, end of interval for L1 fixed and L1 random effects Default: c(1, 1)
#' @param opt.maxiter numeric, maximum iterations to search, Default: 100
#' @param opt.tol numeric, accuracy value, Default: 0.1
#' @param opt.tau numeric, golden proportion coefficient (~0.618) Default: (sqrt(5) - 1)/2
#' @return lmmen list object inluding lmmen fit object of min BIC solution
#'  and summary statistics from the grid search
#' @examples 
#' \dontrun{
#'  dat <- initialize_example(n.i = 5,n = 30,q=4,seed=1)
#'  init <- init.beta(dat,method='glmnet')
#'  golden_section_2d(dat,init)}
#' @rdname golden_section_2d
#' @export 

golden_section_2d = function(dat, 
                             init.beta, 
                             l2=c(1,1), 
                             opt.lb=c(0, 0), 
                             opt.ub=c(1, 1), 
                             opt.maxiter=100, 
                             opt.tol=0.1, 
                             opt.tau=(sqrt(5)-1)/2)
{

# opt.tau=(sqrt(5)-1)/2      # golden proportion coefficient, around 0.618  
# opt.lb=c(0,0)                            # start of interval
# opt.ub=c(1,1)
# opt.tol=0.1               # accuracy value
# opt.maxiter= 10                       # maximum number of iterations

gs.k=0                            # number of iterations

gs.y=gs.x=rep(0,2)

gs.x[1]=opt.lb[1]+(1-opt.tau)*(opt.ub[1]-opt.lb[1])
gs.x[2]=opt.lb[1]+opt.tau*(opt.ub[1]-opt.lb[1])
gs.y[1]=opt.lb[2]+(1-opt.tau)*(opt.ub[2]-opt.lb[2])
gs.y[2]=opt.lb[2]+opt.tau*(opt.ub[2]-opt.lb[2])
#gs.y[2]=opt.ub[2]

gs.z=rbind(c(gs.x[1],gs.y[1]),c(gs.x[1],gs.y[2]),c(gs.x[2],gs.y[1]),c(gs.x[2],gs.y[2]))

gs.fit.0=gs.fit=vector("list", 4)

gs.fit.0[[1]]=gs.fit[[1]]=lmmen(dat, init.beta, frac = c(gs.z[1,1],gs.z[1,2],l2))
gs.fit.0[[2]]=gs.fit[[2]]=lmmen(dat, init.beta, frac = c(gs.z[2,1],gs.z[2,2],l2))
gs.fit.0[[3]]=gs.fit[[3]]=lmmen(dat, init.beta, frac = c(gs.z[3,1],gs.z[3,2],l2))
gs.fit.0[[4]]=gs.fit[[4]]=lmmen(dat, init.beta, frac = c(gs.z[4,1],gs.z[4,2],l2))

# gs.fit1=gs.fit1.0
# gs.fit2=gs.fit2.0

gs.f_z=c(gs.fit[[1]]$BIC,gs.fit[[2]]$BIC,gs.fit[[3]]$BIC,gs.fit[[4]]$BIC)

opt.ubIC.min=min(gs.f_z)
gs.z.iter=gs.z
gs.fz.iter=gs.f_z
gs.z.star.iter=gs.z.min=gs.z.star=gs.z[which.min(gs.f_z),]


while ((abs(opt.ub[1]-opt.lb[1])>opt.tol)*(abs(opt.ub[2]-opt.lb[2])>opt.tol)*(gs.k<opt.maxiter)){
  cbind(gs.f_z,gs.z)
  
  if((which.min(gs.f_z))==1)
    {
      cat('.')
      opt.ub[1]=gs.x[2]
      opt.ub[2]=gs.y[2]

      gs.fit[[4]]=gs.fit[[1]]
      
      gs.x[1]=opt.lb[1]+(1-opt.tau)*(opt.ub[1]-opt.lb[1])
      gs.x[2]=opt.lb[1]+opt.tau*(opt.ub[1]-opt.lb[1])
      gs.y[1]=opt.lb[2]+(1-opt.tau)*(opt.ub[2]-opt.lb[2])
      gs.y[2]=opt.lb[2]+opt.tau*(opt.ub[2]-opt.lb[2])
      
      gs.z=rbind(c(gs.x[1],gs.y[1]),c(gs.x[1],gs.y[2]),c(gs.x[2],gs.y[1]),c(gs.x[2],gs.y[2]))
      
      gs.fit[[1]]=lmmen(dat, init.beta, frac = c(gs.z[1,1],gs.z[1,2],l2))
      gs.fit[[2]]=lmmen(dat, init.beta, frac = c(gs.z[2,1],gs.z[2,2],l2))
      gs.fit[[3]]=lmmen(dat, init.beta, frac = c(gs.z[3,1],gs.z[3,2],l2))
    }
  if((which.min(gs.f_z))==2)
    {
    cat('.')
    opt.ub[1]=gs.x[2]
    opt.lb[2]=gs.y[1]

    gs.fit[[3]]=gs.fit[[2]]
    
    gs.x[1]=opt.lb[1]+(1-opt.tau)*(opt.ub[1]-opt.lb[1])
    gs.x[2]=opt.lb[1]+opt.tau*(opt.ub[1]-opt.lb[1])
    gs.y[1]=opt.lb[2]+(1-opt.tau)*(opt.ub[2]-opt.lb[2])
    gs.y[2]=opt.lb[2]+opt.tau*(opt.ub[2]-opt.lb[2])
    
    gs.z=rbind(c(gs.x[1],gs.y[1]),c(gs.x[1],gs.y[2]),c(gs.x[2],gs.y[1]),c(gs.x[2],gs.y[2]))
    
    gs.fit[[1]]=lmmen(dat, init.beta, frac = c(gs.z[1,1],gs.z[1,2],l2))
    gs.fit[[2]]=lmmen(dat, init.beta, frac = c(gs.z[2,1],gs.z[2,2],l2))
    gs.fit[[4]]=lmmen(dat, init.beta, frac = c(gs.z[4,1],gs.z[4,2],l2))
    }
  if((which.min(gs.f_z))==3)
    {
    cat('.')
    opt.lb[1]=gs.x[1]
    opt.ub[2]=gs.y[2]
    
    gs.fit[[2]]=gs.fit[[3]]
    
    gs.x[1]=opt.lb[1]+(1-opt.tau)*(opt.ub[1]-opt.lb[1])
    gs.x[2]=opt.lb[1]+opt.tau*(opt.ub[1]-opt.lb[1])
    gs.y[1]=opt.lb[2]+(1-opt.tau)*(opt.ub[2]-opt.lb[2])
    gs.y[2]=opt.lb[2]+opt.tau*(opt.ub[2]-opt.lb[2])
    
    gs.z=rbind(c(gs.x[1],gs.y[1]),c(gs.x[1],gs.y[2]),c(gs.x[2],gs.y[1]),c(gs.x[2],gs.y[2]))
    
    gs.fit[[1]]=lmmen(dat, init.beta, frac = c(gs.z[1,1],gs.z[1,2],l2))
    gs.fit[[3]]=lmmen(dat, init.beta, frac = c(gs.z[3,1],gs.z[3,2],l2))
    gs.fit[[4]]=lmmen(dat, init.beta, frac = c(gs.z[4,1],gs.z[4,2],l2))
  }
  if((which.min(gs.f_z))==4)
    {
    cat('.')
    opt.lb[1]=gs.x[1]
    opt.lb[2]=gs.y[1]

    gs.fit[[1]]=gs.fit[[4]]
    
    gs.x[1]=opt.lb[1]+(1-opt.tau)*(opt.ub[1]-opt.lb[1])
    gs.x[2]=opt.lb[1]+opt.tau*(opt.ub[1]-opt.lb[1])
    gs.y[1]=opt.lb[2]+(1-opt.tau)*(opt.ub[2]-opt.lb[2])
    gs.y[2]=opt.lb[2]+opt.tau*(opt.ub[2]-opt.lb[2])
    
    gs.z=rbind(c(gs.x[1],gs.y[1]),c(gs.x[1],gs.y[2]),c(gs.x[2],gs.y[1]),c(gs.x[2],gs.y[2]))
    
    gs.fit[[2]]=lmmen(dat, init.beta, frac = c(gs.z[2,1],gs.z[2,2],l2))
    gs.fit[[3]]=lmmen(dat, init.beta, frac = c(gs.z[3,1],gs.z[3,2],l2))
    gs.fit[[4]]=lmmen(dat, init.beta, frac = c(gs.z[4,1],gs.z[4,2],l2))
  }
  
  
  gs.f_z=c(gs.fit[[1]]$BIC,gs.fit[[2]]$BIC,gs.fit[[3]]$BIC,gs.fit[[4]]$BIC)
  
  opt.ubIC.min=c(opt.ubIC.min,min(gs.f_z))
  gs.z.iter=rbind(gs.z.iter,gs.z)
  gs.fz.iter=rbind(gs.fz.iter,gs.f_z)
  gs.fit.star=gs.fit[[which.min(gs.f_z)]]
  gs.z.star.iter=rbind(gs.z.star.iter,gs.z[which.min(gs.f_z),])
  gs.z.star=gs.z[which.min(gs.f_z),]
  gs.z.min=rbind(gs.z.min,gs.z.star)
    
  gs.k=gs.k+1
}

names(gs.fit.star$fixed) <- names(init.beta)

l.out <- list(fit=gs.fit.star,
              BIC.min=opt.ubIC.min,
              xy.star=gs.z.star,
              xy.min=gs.z.min,
              xy.iter=gs.z.iter,
              BIC.iter=gs.fz.iter,
              iter=gs.k)

l.out <- structure(l.out,class=c('lmmen','list'))
cat('\n')
return(l.out)
}