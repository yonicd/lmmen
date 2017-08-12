#' @title Initialize Scenario
#' @description Create a scenario to run the evaluation functions.
#' @param n.i integer, Observations per subject, Default: 5
#' @param n integer, Number of subjects, Default: 30
#' @param q integer, Number of random effects, Default: 4
#' @param total.beta integer, Number of simulated fixed effects, Default: 9
#' @param true.beta numeric, True of fixed effects indicies, Default: c(1,1,1)
#' @param seed integer, set a seed for reproducibility, Default: NULL
#' @return (n.i*n) x (1+total.beta+q) matrix containing where the subjects index are the matrix rownames
#' \tabular{lcc}{
#' \strong{Description} \tab \strong{Parameter} \tab \strong{Dimension}\cr
#' Response \tab y \tab (n.i*n) x 1 \cr
#' Fixed \tab X \tab (n.i*n) x total.beta\cr
#' Random \tab Z \tab (n.i*n) x q \cr
#'}
#' @examples 
#'  initialize_example(n.i = 5,n = 30,q=4,seed=1)
#'  initialize_example(n.i = 10,n = 60,q=4,seed=1)
#'  initialize_example(n.i = 5,n = 60,q=10,seed=1)
#' @seealso 
#'  \code{\link[mvtnorm]{rmvnorm}}
#' @rdname initialize_example
#' @export 
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats runif
initialize_example=function(n.i = 5,n = 30,q=4,total.beta=9,true.beta=c(1,1,1),seed=NULL){
  
  if(!is.null(seed)) set.seed(seed)
  
  y = NULL 
  X = NULL 
  Z = NULL 
  subject = kronecker(1:n, rep(1,n.i))
  
  true.beta = c(true.beta, rep(0,c(total.beta-length(true.beta))))
  
  if(q<10){
    true.D = matrix(c(9,4.8,0.6,0, 4.8,4,1,0, .6,1,1,0, 0,0,0,0), nrow=q, ncol=q)
  }else{
    true.D = matrix(c(9,4.8,0.6,0, 4.8,4,1,0, .6,1,1,0, 0,0,0,0), nrow=4, ncol=4)
    true.D=rbind(cbind(true.D,matrix(0,nrow=4,ncol=q-4)),matrix(0,nrow=q-4,ncol=q))
  }
  for (i in 1:n)  
  { 
    Xi = cbind(matrix(stats::runif(length(true.beta)*n.i,-2,2), nrow=n.i)) 
    Zi=switch((q<10)+1,{
      cbind(1,Xi)
    },
    {
      cbind(1,matrix(stats::runif((nrow(true.D)-1)*n.i,-2,2), nrow=n.i))}
    )
    X = rbind(X, Xi) 
    Z = rbind(Z, Zi[,-1]) 
    S = Zi%*%true.D%*%t(Zi) + diag(n.i) 
    y = rbind(y, t(mvtnorm::rmvnorm(1, Xi%*%true.beta, S)) ) 
  }

  m=cbind(y,X,Z)
  colnames(m)=c('y',sprintf('X%s',1:ncol(X)),sprintf('Z%s',1:ncol(Z)))
  row.names(m)=subject
  m
}