#' @export
print.lmmen <- function(x,...){
  if(inherits(x,'list')){
    x <- x$fit
  }
  print(x[c('fixed','stddev','BIC','sigma.2')])
}