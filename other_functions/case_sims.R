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

case.scad.mat<-sapply(3,function(idx){
  
  mat<-sapply(case.scad,function(fit,idx){
  c(fit[[idx]]$fixef[-1],sqrt(diag(fit[[idx]]$ran.cov$subject)))  
},idx=idx)

apply(mat,1,function(x) sum(x!=0))
})

colnames(case.scad.mat)<-names(case.scad[[1]][3])
case.scad.df=as.data.frame(case.scad.mat)
case.scad.df$var<-row.names(case.scad.mat)

lmmen.star<-rbind(glmmnet$bstar,glmmnet$rstar)

plot.data<-case.scad.df%>%mutate(lmmen=c(apply(lmmen.star,1,function(x) sum(x!=0))))%>%
  reshape2::melt(.,'var')%>%
  mutate(type=ifelse(grepl('^X',var),'FIXED','RANDOM'),
         variable=factor(variable,labels=c('SCAD','LMMEN')))


library(sqldf)
library(plyr)
library(dplyr)
load('~/Dropbox/GJ_paper/r/paper/casedata.rdata')
load('~/Dropbox/GJ_paper/r/paper/predout_glmmnet.rdata')
load('sim.output/case_lmmen.Rdata')
load('sim.output/case_scad.Rdata')

bstar.lmmen<-glmmnet$bstar
bstar.scad<-sapply(case.scad,function(fit) fit[[3]]$fixef)

lmmen.pred<-plyr::ldply(case.data,.fun=function(prdata){
  Y.temp=prdata%>%select(Ifpid,value)
  X.temp=prdata%>%select(-c(y,v,TgjUserId,TgjUserId__1,ta,gdMean,value,Ifpid,inst_type,created,activated,closed,UserQuesId,agegroup))
  
  w.est=as.matrix(X.temp[1:40])%*%bstar.lmmen
  w.est=exp(3*w.est/sd(w.est))
  a=a=t(apply(w.est,1,quantile,c(0.2,0.8)))
  a1=data.frame(id=seq(1:nrow(a)),lb=a[,1],ub=a[,2])
  west=data.frame(id=rep(seq(1,ncol(w.est)),rep(nrow(w.est),ncol(w.est))),w=c(w.est))
  
  w.est=sqldf::sqldf('select a.id,max(b.lb,min(a.w,b.ub)) as w from west as a left join a1 as b where a.id=b.id')
  estdata=data.frame(id=w.est$id,Ifpid=rep(prdata$Ifpid,length(unique(w.est$id))),y=rep(prdata$value,length(unique(w.est$id))),w=w.est$w)
  pred=sqldf::sqldf("select a.*,b.gm from (select id,Ifpid,sum(y*w)/sum(w) as pr from estdata group by id,Ifpid) as a left join
                    (select Ifpid,avg(Value) as gm from prdata group by Ifpid) as b where a.Ifpid=b.Ifpid")
  pred$bs=(pred$pr-1)^2
  pred$bsgm=(pred$gm-1)^2
  
  pred  
},.progress = 'text')
scad.pred<-plyr::ldply(case.data,.fun=function(prdata){
  Y.temp=prdata%>%select(Ifpid,value)
  X.temp=prdata%>%select(-c(y,v,TgjUserId,TgjUserId__1,ta,gdMean,value,Ifpid,inst_type,created,activated,closed,UserQuesId,agegroup))
  
  w.est=cbind(1,as.matrix(X.temp[1:40]))%*%bstar.scad
  w.est=exp(3*w.est/sd(w.est))
  a=t(apply(w.est,1,quantile,c(0.2,0.8)))
  a1=data.frame(id=seq(1:nrow(a)),lb=a[,1],ub=a[,2])
  west=data.frame(id=rep(seq(1,ncol(w.est)),rep(nrow(w.est),ncol(w.est))),w=c(w.est))
  
  
  w.est=sqldf::sqldf('select a.id,max(b.lb,min(a.w,b.ub)) as w from west as a left join a1 as b where a.id=b.id')
  estdata=data.frame(id=w.est$id,Ifpid=rep(prdata$Ifpid,length(unique(w.est$id))),y=rep(prdata$value,length(unique(w.est$id))),w=w.est$w)
  pred=sqldf::sqldf("select a.*,b.gm from (select id,Ifpid,sum(y*w)/sum(w) as pr from estdata group by id,Ifpid) as a left join
                    (select Ifpid,avg(Value) as gm from prdata group by Ifpid) as b where a.Ifpid=b.Ifpid")
  pred$bs=(pred$pr-1)^2
  pred$bsgm=(pred$gm-1)^2
  
  pred  
},.progress = 'text')

pred<-rbind(lmmen.pred%>%mutate(type='LMMEN'),scad.pred%>%mutate(type='SCAD'))%>%
  group_by(type,Ifpid)%>%summarise_at(.vars=vars('pr','bs','gm','bsgm'),.funs = funs(mean,sd))


p3<-pred%>%reshape2::dcast(Ifpid~type,value.var='bs_mean')%>%
  ggplot(aes(x=LMMEN,y=SCAD))+geom_point()+
  geom_abline(intercept=0,slope=1,linetype=2)+
  theme_bw(base_size = 21)+
  labs(x='LMMEN Brier Scores',y='SCAD Brier Scores')

p2=plot.data%>%filter(type=='RANDOM')%>%
  ggplot(aes(x=var,y=value/100,fill=variable))+
  geom_bar(stat='identity',position='dodge',colour='black',
           show.legend = FALSE)+
  labs(x='Random Covariate',y='% Persistency')+
  theme_bw(base_size = 21)+
  scale_x_discrete(labels=c('(Intercept)',names(case.data[[1]])[55:62]))

p1=plot.data%>%filter(type=='FIXED')%>%
  left_join(plot.data%>%filter(type=='FIXED'&variable=='LMMEN')%>%
             arrange(desc(value))%>%mutate(x=factor(1:n()))%>%select(var,x),by=c('var'))%>%
  ggplot(aes(x=x,y=value/100,group=factor(var),fill=variable))+
  geom_bar(stat='identity',position='dodge',colour='black',
           show.legend = FALSE)+
  facet_grid(variable~.)+
  labs(x='Fixed Covariate',y='% Persistency')+
  theme_bw(base_size = 21)+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

pdf('~/Dropbox/GJ_paper/paper/lmmen_paper-casepanel001.pdf',width=18,height=13)
plot(as.ggedit(list(p1+ggtitle('Fixed Effects Persistence (c)'),
                    p2+ggtitle('Random Effects Persistence (b)'),
                    p3+ggtitle('Prediction Accuracy (a)'))),
     plot.layout = list(list(rows=2,cols=1),
                        list(rows=2,cols=2),
                        list(rows=1,cols=1:2)))

dev.off()

pred%>%
  reshape2::melt(.,c('type','Ifpid'))%>%
  mutate(type=ifelse(grepl('gm',variable),'GM',type),
         variable=as.character(variable),
         variable=gsub('^bsgm','bs',variable),
         variable=gsub('^gm','bs',variable))%>%
  distinct%>%
  group_by(type,variable)%>%
  summarise(stat=mean(value))%>%
  reshape2::dcast(variable~type,value.var='stat')%>%
  select(-GM)%>%slice(1:2)