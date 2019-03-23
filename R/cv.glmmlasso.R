#' @title Cross Validation for glmmLasso package
#' @description Cross Validation for glmmLasso package as shown in example xxx
#' @param dat data.frame, containing y,X,Z and subject variables
#' @param form.fixed formula, fixed param formula, Default: NULL
#' @param form.rnd list, named list containing random effect formula
#' @param lambda numeric, vector containing lasso penalty levels
#' @param family family, family function that defines the distribution link of the glmm, Default: gaussian(link = "identity")
#' @param progress logical
#' @return list of a fitted glmmLasso object and the cv BIC path
#' @examples
#' \dontrun{cv.glmmLasso(initialize_example(seed=1))}
#' @seealso 
#'  \code{\link[glmmLasso]{glmmLasso}}
#'  @references 
#'  Variable selection for generalized linear mixed models by ell 1-penalized estimation. 
#'  Statistics and Computing, pages 1–18, 2014.
#' @rdname cv.glmmLasso
#' @export 
#' @importFrom glmmLasso glmmLasso
#' @importFrom stats gaussian as.formula
cv.glmmLasso=function(dat,
                      form.fixed=NULL,
                      form.rnd=NULL,
                      lambda=seq(500,0,by=-5),
                      family=stats::gaussian(link = "identity"),
                      progress = TRUE
                      ) {

    stopifnot(inherits(dat, "data.frame"))
    d.size <- length(unique(dat$subject))*(sum(grepl('^Z',names(dat)))+1)+
        (sum(grepl('^X',names(dat)))+1)
  
    if (is.null(form.fixed)) {
        form.fixed<- reformulate(grep("^X",names(dat), value=TRUE),
                                 response="y")
    }
    if(is.null(form.rnd)) {
        ## form.rnd <- reformulate(grep("^Z",names(dat), value=TRUE),
        ## response="subject")
        form.rnd <- setNames(list(reformulate(grep("^Z",names(dat), value=TRUE))),
                             "subject")
    }
    
    if (progress) pb <- txtProgressBar(style=3, max=length(lambda))

    get.lastpars <- function(f) {
        list(start=f$Deltamatrix[f$conv.step,],
             q.start=f$Q_long[[f$conv.step]])
    }
             
    fit.list <- list()
    fit.list[[1]] <-
        try(glmmLasso::glmmLasso(fix = form.fixed,
                                 rnd = form.rnd,
                                 data = dat,
                                 lambda = lambda[1],
                                 switch.NR = FALSE,
                                 final.re=FALSE),
            silent=TRUE)
    if (inherits(fit.list[[1]], "try-error")) {
        stop("first lambda value failed")
    }        
    if (progress) setTxtProgressBar(pb, 1)
    next.start <- get.lastpars(fit.list[[1]])
        
    for (j in seq_along(lambda)[-1]) {
        if (progress) setTxtProgressBar(pb, j)        
        fn <- suppressMessages(suppressWarnings(
            try(glmmLasso::glmmLasso(fix = form.fixed,
                                     rnd = form.rnd,
                                     data = dat,
                                     lambda = lambda[j],
                                     switch.NR = FALSE,
                                     final.re=FALSE,
                                     control =
                                         glmmLassoControl(start=next.start$start,
                                                         q.start=next.start$q.start)),
                                     silent=TRUE)))
        if (!inherits(fn,"try-error")) {
            next.start <- get.lastpars(fn)
            fit.list[[j]] <- fn
        }
    }

    BIC_vec <- vapply(fit.list,
                      function(x) if (is.null(x)) NA_real_ else x[["bic"]],
                      numeric(1))
    
    fit.opt <- fit.list[[which.min(BIC_vec)]]
    list(fit.opt=fit.opt, BIC_path=BIC_vec)
}  
