
#' @title linear mixed model elastic net
#' @description Solve a linear mixed model elastic net optimization problem.
#' @param data matrix, data
#' @param init.beta numeric, initial values for fixed effects coefficients
#' @param frac PARAM_DESCRIPTION
#' @param eps numeric, tolerance level to pass to solve.QP, Default: 10^(-4)
#' @param verbose boolean, show output during optimization Default: FALSE
#' @return list
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @importFrom quadprog solve.QP
#' @export 

lmmen = function(data, init.beta, frac, eps = 10^(-4),verbose=FALSE)
{
  
  tc <- textConnection(NULL, "w") 
	if (eps <=0) {return(cat("ERROR: Eps must be > 0. \n"))}
  y=matrix(data[,grepl('^y',colnames(data))],ncol=1)
  X=data[,grepl('^X',colnames(data))]
  Z.0=data[,grepl('^Z',colnames(data))]
  subject=as.numeric(rownames(data))
  init.beta=c(as.matrix(init.beta))
	round.set=6
	n.i = tabulate(subject)
	n.tot = sum(n.i)
	n = length(n.i)
	Z = as.matrix(Z.0,nrow=n.tot)
	Z = cbind(rep(1,n.tot), Z)
	p = ncol(X)
	q = ncol(Z)	
	hyper=1

  
	if (qr(X)$rank < p) {return(cat("ERROR: Design matrix for fixed effects is not full rank. \n"))}
	if (qr(Z)$rank < q) {return(cat("ERROR: Design matrix for random effects is not full rank. \n"))}

	eps.tol = 1e-8
  sigma.hat=1
  beta.hat=init.beta
	beta.hatp=(abs(beta.hat))^hyper
	beta.p=t(1/rbind(beta.hatp,beta.hatp))

  D.lme=diag(q)
  stdev.init=rep(1,q)

	junk = t(chol(D.lme+eps.tol*diag(q)))

	lambda.hat = diag(junk)
	lambda.hat = pmax(lambda.hat, eps.tol)
	gamma.init = diag(as.vector(1/lambda.hat))%*%junk
	gamma.hat = gamma.init[lower.tri(gamma.init)]
	lambda.p = t(as.matrix(1/(lambda.hat^hyper)))

	#Block Diagonal of Z
  Z.bd = matrix(0,nrow=(n.tot),ncol=(n*q))
	W.bd = matrix(0,nrow=(n*q),ncol=(n*q))
	start.point = 1
	for (i in 1:n){
		end.point = start.point + (n.i[i] - 1) 
		Z.bd[start.point:end.point,(q*(i-1)+1):(q*i)] = Z[start.point:end.point,]
		start.point = end.point + 1
	}
	
	W.bd = t(Z.bd)%*%Z.bd

	new.beta = beta.hat
	new.lambda = lambda.hat
	new.gamma = gamma.hat
	sigma.2.current = as.numeric(sigma.hat)

	X.star = cbind(X,-X)
	X.star.quad = t(X.star)%*%X.star

	
  ####################################################
  A.trans = rbind(diag(2*p + q), -c(beta.p,rep(0,q)),-c(rep(0,2*p),lambda.p))
	####################################################
    
  cr.full.k = rep(1,n)%x%diag(q)

	K.matrix=rbind(
	                t(diag(rep(1,q-1))%x%rep(1,q-1))[,which(lower.tri(matrix(1,q-1,q-1),diag = TRUE))],
            	    rep(0,q*(q-1)/2)
	                )
	
	ident.tilde=t(rep(1,q-1)%x%diag(rep(1,q-1)))[,which(lower.tri(matrix(1,q-1,q-1),diag = TRUE))]

	beta.est = NULL
	lambda.est = NULL
	gamma.est = NULL
	sigma.est = NULL
	BIC.value = NULL

#LASSO Penalties
 	t.bound.l1.f = frac[1]*(sum(abs(init.beta)))
 	t.bound.l1.r = frac[2]*(sum(abs(lambda.hat)))

#Ridge Penalties
  #alpha=l1.f/l1.f+l2.f
  #l1.f=alpha*l1.f+alpha*l2.f
  #l2.f=(l1.f-alpha*l1.f)/alpha
  alpha=frac[3:4]
  t.bound.l2.f=((t.bound.l1.f-alpha[1]*t.bound.l1.f)/alpha[1])*10^(4)
  t.bound.l2.r=((t.bound.l1.r-alpha[2]*t.bound.l1.r)/alpha[2])*10^(4)
  
	b.0 = c(rep(0,2*p + q), -t.bound.l1.f, -t.bound.l1.r)
	outer.converge = FALSE
	n.iter = 0

 		while ((outer.converge==FALSE) && (n.iter < 200)){
			n.iter = n.iter + 1

			beta.current = beta.iterate = new.beta
			lambda.current = new.lambda
			gamma.current = new.gamma
			gamma.mat.current = diag(q)
			gamma.mat.current[lower.tri(gamma.mat.current)] = gamma.current
			full.gamma.mat = diag(n)%x%gamma.mat.current

			n.iter1 = 0
			inner.converge = FALSE
 			while ((inner.converge==FALSE) && (n.iter1 < 100)){	
				beta.current = new.beta
				lambda.current = new.lambda
				n.iter1 = n.iter1 + 1
				resid.vec.current = y-(X%*%beta.current)
				full.gamma.mat = diag(n)%x%gamma.mat.current

				full.D.mat = diag(n)%x%diag(as.vector(lambda.current))
				Cov.mat.temp = as.matrix(Z.bd%*%full.D.mat%*%full.gamma.mat)
				sigma.2.current = as.numeric(t(resid.vec.current)%*%solve(Cov.mat.temp%*%t(Cov.mat.temp)+diag(n.tot))%*%resid.vec.current/n.tot)

				full.inv.Cov.mat = solve(t(Cov.mat.temp)%*%Cov.mat.temp + diag(n*q))
				exp.bhat = full.inv.Cov.mat%*%t(Cov.mat.temp)%*%resid.vec.current
				exp.Uhat = full.inv.Cov.mat*sigma.2.current
				exp.Ghat = exp.Uhat + exp.bhat%*%t(exp.bhat)
		
				right.side.mat = as.matrix(Z.bd%*%diag(as.vector(full.gamma.mat%*%exp.bhat))%*%cr.full.k)
        
        #Ridge Pen for RE
				  lower.diag.mat = (as.matrix(t(cr.full.k)%*%(W.bd * (full.gamma.mat%*%exp.Ghat%*%t(full.gamma.mat)))
                                     %*%cr.full.k)+diag(rep(t.bound.l2.r,q))/(1+t.bound.l2.r))

				  full.right.side = as.matrix(t(X.star)%*%right.side.mat)

				#Ridge Pen for FE
				  D.quadratic.prog = rbind(cbind(as.matrix((X.star.quad+t.bound.l2.f*diag(rep(1,(2*p)))/(1+t.bound.l2.f))), full.right.side),cbind(t(full.right.side), lower.diag.mat))

				d.linear.prog = as.vector(t(y)%*%cbind(as.matrix(X.star), right.side.mat))*1e-9

				D.quadratic.prog = D.quadratic.prog*1e-9+eps.tol*diag(nrow(D.quadratic.prog))
				
        ##############################
				beta.lambda = quadprog::solve.QP(D.quadratic.prog, d.linear.prog, t(A.trans), bvec=b.0)
        ##############################
        
				new.beta = round(beta.lambda$solution[1:p]-beta.lambda$solution[(p+1):(2*p)],6)
				new.lambda = round(beta.lambda$solution[-(1:(2*p))],6)

				diff = abs(beta.current-new.beta)
			if (max(c(diff))<eps) inner.converge = TRUE
			
 		}

			E.A = NULL
			start.point = 1
			d.d.t = new.lambda%*%t(new.lambda)

			full.A.t.A.matrix = 0*diag(q*(q-1)/2)
			T.vec = rep(0,q*(q-1)/2)

			for (i in 1:n)
			{
				end.point = start.point + (n.i[i] - 1) 
				E.Ai = as.matrix((as.matrix(rep(1,n.i[i]),ncol=1)%*%exp.bhat[(q*(i-1)+1):(q*i)]%*%K.matrix)*
				                   (Z.bd[start.point:end.point,(q*(i-1)+2):(q*i)]%*%diag(new.lambda[-1])%*%ident.tilde))
				E.A = rbind(E.A, E.Ai)
				start.point = end.point + 1

				G.i = exp.Ghat[(q*(i-1)+1):(q*i),(q*(i-1)+1):(q*i)]
				Z.Z.i = W.bd[(q*(i-1)+1):(q*i),(q*(i-1)+1):(q*i)]
				B.matrix = as.matrix(Z.Z.i * d.d.t)
				A.t.A.matrix = NULL
				Cross.matrix = NULL

				for (j in 1:(q-1))
				{
					U.j.matrix = diag(q)[-(1:j),]
					Cross.matrix = rbind(Cross.matrix, U.j.matrix%*%B.matrix)
					A.t.A.row = NULL
					for (k in 1:(q-1))
					{
						V.k.matrix = diag(q)[,-(1:k)]
						A.t.A.row = cbind(A.t.A.row, U.j.matrix%*%B.matrix%*%V.k.matrix*G.i[j,k])  
					}
					A.t.A.matrix = rbind(A.t.A.matrix, A.t.A.row)
				}
				T.vec = T.vec + as.vector(((t(K.matrix)%*%G.i)*Cross.matrix)%*%rep(1,q))
				full.A.t.A.matrix = full.A.t.A.matrix + A.t.A.matrix
			}

			A.eigen = eigen(full.A.t.A.matrix)
			A.eigen.vals = round(A.eigen$values,5)
			A.eigen.vecs = A.eigen$vectors
			eig.A = A.eigen.vals^(-1)
			eig.A[is.infinite(eig.A)] = 0

			A.t.A.inv = (A.eigen.vecs%*%diag(eig.A)%*%t(A.eigen.vecs))
			lin.term = t(E.A)%*%resid.vec.current-T.vec
      
      #Gamma Star##############################
			new.gamma = round(A.t.A.inv%*%lin.term,round.set)
      #########################################
      
			counter.1 = 1
			for (j in 1:(q-1))
			{
				for (k in (j+1):q)
				{
					new.gamma[counter.1] = new.gamma[counter.1]*(new.lambda[j]>0)
					counter.1 = counter.1 + 1
				}
			}

			diff = abs(beta.iterate-new.beta)
			if (max(c(diff))<eps) outer.converge = TRUE
			
		}

		beta.est = as.matrix(cbind(beta.est, new.beta))
		lambda.est = as.matrix(cbind(lambda.est, new.lambda))
		gamma.est = as.matrix(cbind(gamma.est, new.gamma))
		resid.vec = y-(X%*%new.beta)
		gamma.mat = diag(q)
		gamma.mat[lower.tri(gamma.mat)] = new.gamma

		full.gamma.mat = diag(n)%x%gamma.mat

		full.D.mat = diag(n)%x%diag(as.vector(new.lambda))
		Cov.mat.temp = Z.bd%*%full.D.mat%*%full.gamma.mat
		Full.cov.mat = as.matrix(Cov.mat.temp%*%t(Cov.mat.temp)+diag(n.tot))
		new.sigma.2 = as.numeric(t(resid.vec)%*%solve(Full.cov.mat)%*%resid.vec/n.tot)

		sigma.est = c(sigma.est, new.sigma.2)
		Full.Cov.est = new.sigma.2*Full.cov.mat
		Mean.est = X%*%new.beta

		loglikes = -2*(dmvnorm(x=as.vector(y),mean=as.vector(Mean.est),sigma=as.matrix(Full.Cov.est),log=TRUE))
		df.par = sum(new.beta!=0)+sum(new.lambda!=0)*((sum(new.lambda!=0)+1)/2)
    BIC.frac=loglikes + df.par*log(n.tot)
		BIC.value = c(BIC.value,BIC.frac)

    if(verbose==TRUE){
		message("Finished bound of (" ,frac[1:2],",",frac[3:4],"), BIC:",BIC.frac,"\n")
    }

	min.BIC = which.min(BIC.value)
	beta.BIC = beta.est[,min.BIC]
	lambda.BIC = lambda.est[,min.BIC]
	sigma.2.BIC = sigma.est[min.BIC]
	gamma.BIC = gamma.est[,min.BIC]    
	gamma.BIC.mat = diag(q)
	gamma.BIC.mat[lower.tri(gamma.BIC.mat)] = gamma.BIC
	temp.mat = diag(as.vector(lambda.BIC))%*%gamma.BIC.mat
	Cov.Mat.RE = sigma.2.BIC*temp.mat%*%t(temp.mat)
	fit = NULL
	fit$fixed = beta.BIC
	fit$stddev = sqrt(diag(Cov.Mat.RE))
  fit$lambda=lambda.BIC
	fit$BIC = BIC.value
	fit$frac = frac
	fit$sigma.2 = sigma.2.BIC
  fit$gamma.mat = gamma.BIC.mat
	fit$corr = round(diag(1/(fit$stddev+eps.tol))%*%Cov.Mat.RE%*%diag(1/(fit$stddev+eps.tol)),round.set)
  fit$beta.lambda=beta.lambda
  close(tc) 
	return(fit)
}