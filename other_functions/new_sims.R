.libPaths('../lib')
library(qapply)
library(parallel)
library(lmmlasso)
library(lmmen)

system.time(out<-mclapply(sim.data[201:210],cv.lmmlasso,startValue=2,mc.cores = detectCores()-1))

system.time(p<-qapply(X = sim.data[401:600],
                                  FUN = cv.lmmlasso,fargs = list(startValue=0),
                                  workDir = '../qapply',
                                  tag = 'lmmlasso',nCores = 80,internalize = FALSE))

sim.lmmlasso3 <- qinternalize("../qapply/out/lmmlasso")
sim.lmmlasso3<-do.call('rbind',sim.lmmlasso3)
save(sim.lmmlasso3,file='sim.output/lmmlassoScen3.Rdata')

set.seed(1)
c0=rnorm(1,0,0.5)
data(sim.data)
scen4<-lapply(sim.data[201:400],function(dat,a,c0){
  dat[,'X3']=a*dat[,'X1']+(1-a)*dat[,'X2']+c0
  dat
},a=.3,c0=c0)

sim.data<-c(sim.data,scen4)

do.call('rbind',out)

cv.lmmlasso(sim.data[[601]],startValue=0)


#############

load('sim.output/graph_data2.rdata')
load('sim.output/glmmlassoScen1.Rdata')
load('sim.output/glmmlassoScen2.Rdata')

load('sim.output/lmmlassoScen1.Rdata')
load('sim.output/lmmlassoScen2.Rdata')
load('sim.output/lmmlassoScen3.Rdata')

load('sim.output/scadScen1.Rdata')
load('sim.output/scadScen2.Rdata')
load('sim.output/scadScen3.Rdata')

colnames(sim.lmmlasso3)[1:20]<-c('(Intercept)',paste0('X',1:9),paste0('V',1:10))
load('sim.output/lmmlassoScen4.Rdata')
colnames(sim.lmmlasso4)[1:14]<-c('(Intercept)',paste0('X',1:9),paste0('V',1:4))

colnames(sim.scad1)<-gsub('^z','V',colnames(sim.scad1))
colnames(sim.scad2)<-gsub('^z','V',colnames(sim.scad2))
colnames(sim.scad3)<-gsub('^z','V',colnames(sim.scad3))

sim.lmmlasso<-list(sim.lmmlasso1,sim.lmmlasso2,sim.lmmlasso3,sim.lmmlasso4)
sim.glmmlasso<-list(glmmlasso.scen1,glmmlasso.scen2)
sim.scad<-list(sim.scad1,sim.scad2,sim.scad3)

names(sim.lmmlasso)<-c(1:4)
names(sim.glmmlasso)<-c(1:2)
names(sim.scad)<-c(1:3)

sim.lmmlasso<-plyr::ldply(sim.lmmlasso,.fun=function(data){
  x<-data%>%data.frame()%>%
    mutate(iterid=1:nrow(.),iter=lambda.opt/5)%>%select(-lambda.opt)%>%
    reshape2::melt(.,c('iterid','iter'),value.name='est')%>%filter(grepl('[0-9]|bic',variable))
  
  x$type<-as.character(x$variable)
  x$type[grepl('^X',x$type)]<-'FIXED'
  x$type[grepl('^V',x$type)]<-'STDDEV'
  x$type[grepl('^bic',x$type)]<-'BIC'
  
  x=x%>%group_by(iterid,type)%>%mutate(id=1:n())
  x$simtype='lmmlasso'
  
  x<-x[,head(names(fit.min.all),-1)]  
},.id='ex')%>%mutate(ex=as.numeric(as.character(ex)))

sim.scad<-plyr::ldply(sim.scad,.fun=function(data){
  x<-data%>%data.frame()%>%
    mutate(iterid=1:nrow(.),iter=1)%>%
    reshape2::melt(.,c('iterid','iter'),value.name='est')%>%
    filter(grepl('[1-9]|bic',variable))
  
  x$type<-as.character(x$variable)
  x$type[grepl('^X',x$type)]<-'FIXED'
  x$type[grepl('^V',x$type)]<-'STDDEV'
  x$type[grepl('^bic',x$type)]<-'BIC'
  
  x=x%>%group_by(iterid,type)%>%mutate(id=1:n())
  x$simtype='scad'
  
  x<-x[,head(names(fit.min.all),-1)]  
},.id='ex')%>%mutate(ex=as.numeric(as.character(ex)))

sim.glmmlasso<-plyr::ldply(sim.glmmlasso,.fun=function(data){
  x<-data%>%data.frame()%>%
    mutate(iterid=1:nrow(.),iter=lambda.opt/5)%>%select(-lambda.opt)%>%
    reshape2::melt(.,c('iterid','iter'),value.name='est')%>%filter(grepl('[0-9]|subject|bic',variable))
  
  x$type<-as.character(x$variable)
  x$type[grepl('^X',x$type)]<-'FIXED'
  x$type[grepl('^subject',x$type)]<-'STDDEV'
  x$type[grepl('^bic',x$type)]<-'BIC'
  
  x$est[x$type=='STDDEV']=sqrt(x$est[x$type=='STDDEV'])
  
  x=x%>%group_by(iterid,type)%>%mutate(id=1:n())
  x$simtype='glmmlasso'
  
  x<-x[,head(names(fit.min.all),-1)]  
},.id='ex')%>%mutate(ex=as.numeric(as.character(ex)))

fit.min.all<-rbind(fit.min.all,sim.lmmlasso,sim.glmmlasso,sim.scad)

###########

dat.df<-sim.data[[1]]
dat.df$Subject<-as.numeric(rownames(dat.df))
rownames(dat.df)<-NULL
dat.df$iter<-NULL
fit.lme4 <- glmer(y ~ X - 1 + (Z - 1 | cluster), data = dat,family = "poisson")
sat.fit <- build.start.fit(fit.lme4, gamma = 2)

dat<-initialize_example(n.i = 5,n = 30,q=4,total.beta = 100,true.beta = rep(1,20),seed=1)
library(rpql)
cv.scad<-function(dat,pen.type='scad',...){
  dat<-as.matrix(dat)
  if('adl'%in%pen.type){  
    ids <- as.numeric(rownames(dat))
    y <- dat[,'y']
    X <- cbind(X0=1,dat[,grepl('^X',colnames(dat))])
    Z <- cbind(Z0=1,dat[,grepl('^Z',colnames(dat))])
    fit.lmer <- lme4::lmer(y ~ X - 1 + (Z - 1 | ids))
    sat.fit <- rpql::build.start.fit(fit.lmer, gamma = 2)
  }
  
  ids <- list(subject=as.numeric(rownames(dat)))
  rownames(dat)<-NULL
  y <- dat[,'y']
  X <- cbind(X0=1,dat[,grepl('^X',colnames(dat))])
  Z <- list(subject=cbind(Z0=1,dat[,grepl('^Z',colnames(dat))]))
  lambda.seq <- lseq(1e-4,1,length=100)
  
  if('adl'%in%pen.type){
    fit <- rpql::rpqlseq(y = y, X = X, Z = Z, id = ids, lambda = lambda.seq, pen.type = pen.type, 
                         pen.weights = sat.fit$pen.weights, start = sat.fit,...)    
  }else{
    fit <- rpql::rpqlseq(y = y, X = X, Z = Z, id = ids, lambda = lambda.seq, pen.type = pen.type,...) 
  }
  #c(fit$best.fit[[3]]$fixef,sqrt(diag(fit$best.fit[[3]]$ran.cov$subject)),bic.opt=fit$best.fits[[3]]$ics[[3]])
  fit$best.fit
}

x<-cv.scad(sim.data[[1]],pen.type = c('adl'))
out<-sapply(x,function(y) c(y$fixef,sqrt(diag(y$ran.cov$subject))))
View(out)

system.time(sim.scad_adl1<-qapply(X = sim.data[1:200],FUN = cv.scad,tag = 'rpql',fargs = list(pen.type=c('scad','adl')),workDir ='../qapply',nCores = 32 ,clearWd = FALSE))
sim.scad_adl1<-do.call('rbind',sim.scad_adl1)

system.time(sim.scad_adl2<-qapply(X = sim.data[201:400],FUN = cv.scad,tag = 'rpql',fargs = list(pen.type=c('scad','adl')),workDir ='../qapply',nCores = 32 ,clearWd = FALSE))
sim.scad_adl2<-do.call('rbind',sim.scad_adl2)

system.time(sim.scad_adl3<-qapply(X = sim.data[401:432],FUN = cv.scad,tag = 'rpql',fargs = list(pen.type=c('scad','adl')),workDir ='../qapply',nCores = 32 ,clearWd = FALSE))
sim.scad_adl3<-do.call('rbind',sim.scad_adl3)

save(sim.scad_adl1,file='sim.output/scad_adl_Scen1.Rdata')
save(sim.scad_adl2,file='sim.output/scad_adl_Scen2.Rdata')
save(sim.scad_adl3,file='sim.output/scad_adl_Scen3.Rdata')
