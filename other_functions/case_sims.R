library(parallel)
out<-mclapply(sim.data[401:402],cv.lmmlasso,mc.cores = 3)
do.call('rbind',out)

a<-cv.lmmlasso(initialize_example(n.i=5,n=30,q=4,total.beta = 90,true.beta = 20,seed = 1))


case.setup<-function(prdata){
rmX=c('TgjUserId__1','y','v','TgjUserId','ta','gdMean',
      'value','Ifpid','inst_type','created','activated',
      'closed','UserQuesId','agegroup')

zs<-c("gr1B","gr1C","gr2A","gr2B","gr2C","gr4A","gr4B","gr4C")

Y=as.matrix(prdata%>%select(value))
X=as.matrix(scale(prdata[,!names(prdata)%in%rmX]))[,(1:40)]
Z=as.matrix(prdata[,names(prdata)%in%zs])
ids<-as.numeric(prdata$subject)

dat<-cbind(Y,X,Z)
colnames(dat)<-c('y',paste0('X',1:40),paste0('Z',1:8))
row.names(dat)<-ids
dat
}

case.scad<-parallel::mclapply(case.data,function(df,pen.type='scad'){
  cv.scad(case.setup(df),pen.type = pen.type)
},mc.cores=parallel::detectCores()-1)

case.adl<-parallel::mclapply(case.data,function(df,pen.type='adl'){
  cv.scad(case.setup(df),pen.type = pen.type)
},mc.cores=parallel::detectCores()-1)

case.scad.mat<-sapply(1:5,function(idx){
  
  mat<-sapply(case.scad,function(fit,idx){
  c(fit[[idx]]$fixef[-1],sqrt(diag(fit[[idx]]$ran.cov$subject)))  
},idx=idx)

apply(mat,1,function(x) sum(x!=0))
})

colnames(case.scad.mat)<-names(case.scad[[1]])
case.scad.df=as.data.frame(case.scad.mat)
case.scad.df$var<-row.names(case.scad.mat)

lmmen.star<-rbind(glmmnet$bstar,glmmnet$rstar)

plot.data<-case.scad.df%>%mutate(lmmen=c(apply(lmmen.star,1,function(x) sum(x!=0))))%>%
  reshape2::melt(.,'var')%>%
  mutate(type=ifelse(grepl('^X',var),'FIXED','RANDOM'),
         variable=factor(variable,labels=c('AIC','BIC1','BIC2','IC1','IC2','LMMEN')))

plot.data%>%filter(type=='RANDOM')%>%
  ggplot(aes(x=variable,y=value,group=factor(var),fill=variable))+
  geom_bar(stat='identity',position='dodge',colour='black')

summary(x$best.fits[[3]])

out<-sapply(x,function(y) c(y$fixef,sqrt(diag(y$ran.cov$subject))))

library(lmmen)
x1<-cv.lmmlasso(case.data[[1]])
