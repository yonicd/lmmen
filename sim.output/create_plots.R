library(ggplot2);library(scales);library(dplyr);library(grid)
load('./sim.output/graph_data2.rdata')
load('./sim.output/glmmlassoScen1.Rdata')
load('./sim.output/glmmlassoScen2.Rdata')

load('./sim.output/lmmlassoScen1.Rdata')
load('./sim.output/lmmlassoScen2.Rdata')
load('./sim.output/lmmlassoScen3.Rdata')

load('./sim.output/scadScen1.Rdata')
load('./sim.output/scadScen2.Rdata')
load('./sim.output/scadScen3.Rdata')

colnames(sim.lmmlasso3)[1:20]<-c('(Intercept)',paste0('X',1:9),paste0('V',1:10))
load('./sim.output/lmmlassoScen4.Rdata')
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


vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)

fs=21

# real_param=data.frame(ex=c(rep(seq(1,5),c(rep(9,4),200)),rep(seq(1,5),c(4,4,10,4,4))),
#                       id=c(rep(seq(1,9),4),seq(1,200),rep(seq(1,4),2),seq(1,10),rep(seq(1,4),2)),
#                       type=rep(c("FIXED","STDDEV"),rep(c(236,26))),
#                       p=c(rep(c(1,1,rep(0,7)),2),1,0,1,rep(0,6),rep(1,3),rep(0,6),rep(1,20),rep(0,180),
#                           rep(c(rep(1,3),0),2),rep(1,3),rep(0,7),rep(c(rep(1,3),0),2)))

temp=left_join(
  fit.min.all%>%
    filter(type!="BIC")%>%
    mutate(type=factor(type)),
  real_param,
  by=c("ex","type","id"))%>%
  mutate(nz=est!=0)%>%
  group_by(ex,simtype,type,iterid)%>%
  summarise(nz=sum(p))
  
  
  mutate(nz=as.numeric(as.numeric(est!=0)==p))%>%
  group_by(ex,simtype,type,iterid)%>%
  summarise(truez=sum(nz))

temp=left_join(temp,real_param%>%group_by(ex,type)%>%
                 summarise(n=n()),by=c("ex","type"))%>%
  mutate(truez=ifelse((truez/n)>1,truez/2,truez),pct=truez/n)

temp1=temp%>%group_by(simtype,ex,iterid)%>%
  summarise(truez=sum(truez),n=sum(n))%>%mutate(pct=truez/n)

a1<-rbind(
  temp1%>%group_by(simtype,ex)%>%summarise(pct=mean(pct))%>%mutate(type="MODEL"),
  temp%>%group_by(simtype,type,ex)%>%summarise(pct=mean(pct)))%>%
  mutate(pct.type="Mean Percent")

a2<-rbind(temp1%>%
            mutate(type="MODEL",all=as.numeric(pct==1))%>%
            group_by(type,simtype,ex,all)%>%
            summarise(pct=1-(n()/200))%>%
            filter(all==0)%>%
            select(-all),
          temp%>%
            mutate(all=as.numeric(pct==1))%>%
            group_by(type,simtype,ex,all)%>%
            summarise(pct=1-(n()/200))%>%
            filter(all==0)%>%
            select(-all))%>%
  ungroup%>%mutate(pct.type="Oracle")

a3<-expand.grid(type=c("MODEL","FIXED","STDDEV"),pct.type=c("Mean Percent","Oracle"))%>%mutate(simtype="penlme",ex=5,pct=0)

temp.compare=bind_rows(a1,a2,a3)%>%ungroup%>%
  arrange(simtype,type,pct.type,ex)

temp.compare$type=factor(temp.compare$type,levels=unique(temp.compare$type),
                         labels=c("Model","Fixed Effects","Random Effects"))

temp.compare$simtype=factor(temp.compare$simtype,
                            labels=c('GLMMLASSO',"LMMEN",'LMMLASSO',"Pen.LME","SCAD"))

p.compare= temp.compare%>%filter(!is.na(pct))%>%mutate(ex=factor(ex))%>%
  ggplot(aes(x=ex,y=pct,fill=simtype))+geom_bar(stat="identity",position="dodge")+facet_grid(pct.type~type)+theme_bw(base_size = fs)+
  ylab("Percent")+xlab("Scenario")+
  scale_fill_discrete(name="Simulation Type")+theme(legend.position="bottom")+
  scale_y_continuous(labels = percent)+theme(legend.position="bottom")

# fit.min.ex5=fit.min.all%>%filter(type!="BIC",ex==5)%>%mutate(id=ifelse(id>20,21,id))
# fit.min.ex5$id=factor(fit.min.ex5$id,labels=c(seq(1,20),"+20"))
# 
# ggplot(fit.min.ex5%>%filter(type=="FIXED"),aes(x=factor(id),y=est,fill=simtype))+ geom_boxplot()+facet_wrap(ex~type,scales="free",ncol=1)+theme(axis.ticks=element_blank(),axis.text.x=element_blank())

# case.full=data.frame(simtype=rep(c("glmnet","lme4","lmmen"),rep(49,3)),id=rep(c(1:40,1:9),3),type=rep(rep(c("fixed","stddev"),c(40,9)),3))

temp.case=case.fit.min.all%>%filter(type!="BIC")%>%mutate(nz=ifelse(est!=0,1,0))%>%group_by(simtype,type,id,nz)%>%summarise(pct=n()/100,mest=mean(est))%>%filter(nz==1)%>%select(-nz)

temp.case$simtype=factor(temp.case$simtype,levels=c("lme4","glmnet","lmmen"),labels=c("lme4","GLMNET","LMMEN"))
temp.case$type=toupper(temp.case$type)

p.case=ggplot(temp.case%>%filter(id!=36),aes(x=factor(id),y=mest,fill=pct))+geom_bar(stat="identity",position="dodge")+facet_grid(simtype~type,scales="free_x")+theme_bw()+scale_fill_discrete(name="Persistency",breaks=seq(0,1,0.1),labels=percent)+scale_x_discrete(labels=seq(1:39))+ylab("Mean Parameter Estimate")+xlab("Parameter")

gen.plot=fit.min.all%>%
  dplyr::mutate(id=ifelse(id>20,21,id))%>%
  dplyr::filter(type!="BIC")

gen.plot$type=factor(gen.plot$type,labels=c("Fixed Effects","Random Effects"))

gen.plot$simtype.lbl=factor(gen.plot$simtype,labels=c('GLMMLASSO',"LMMEN",'LMMLASSO',"Pen.LME","SCAD"))

gen.plot=gen.plot%>%mutate(id=factor(id,labels=c(1:20,"21+")),flabel=sprintf('%s: %s',ex,type))

gen.plot$simtype.lbl<-factor(gen.plot$simtype.lbl,levels=levels(gen.plot$simtype.lbl)[c(2,4,1,3,5)])

yint<-gen.plot%>%filter(ex<=3)%>%
  mutate(type=as.character(type))%>%
  dplyr::distinct(flabel,type)%>%
  dplyr::left_join(data.frame(yint=c(0,1,3,2,1,0),type=c(rep('Fixed Effects',2),rep('Random Effects',4)),stringsAsFactors = FALSE),by='type')

p4.1=gen.plot%>%filter(ex<=3)%>%
  ggplot(aes(x=id,y=est,fill=simtype.lbl))+
    geom_boxplot()+
    facet_wrap(~flabel,scales="free",ncol=2)+
    xlab("Parameter")+
    ylab("Estimate")+
    scale_fill_discrete(name="Simulation Type")+
    theme_bw(base_size = fs)+
    theme(legend.position="bottom")+
    geom_hline(aes(yintercept=yint),linetype=2,data=yint)

yint<-gen.plot%>%filter(ex>3)%>%
  mutate(type=as.character(type))%>%
  dplyr::distinct(flabel,type)%>%
  dplyr::left_join(data.frame(yint=c(0,1,3,2,1,0),type=c(rep('Fixed Effects',2),rep('Random Effects',4)),stringsAsFactors = FALSE),by='type')

p4.2=gen.plot%>%filter(ex>3)%>%
  ggplot(aes(x=id,y=est,fill=simtype.lbl))+
  geom_boxplot()+
  facet_wrap(~flabel,scales="free",ncol=2)+
  xlab("Parameter")+
  ylab("Estimate")+
  scale_fill_discrete(name="Simulation Type")+
  theme_bw(base_size = fs)+
  theme(legend.position="bottom")+
  geom_hline(aes(yintercept=yint),linetype=2,data=yint)


p_scat=ggplot(case.panel$obs, aes(lmmen,glmnet))+ geom_point()+theme_bw(base_size = fs)+
  #geom_text(vjust=1,size=3)+
  geom_abline(intercept=0,slope=1,linetype="dashed")+xlim(min(case.panel$obs[,1:2]),max(case.panel$obs[,-1]))+ylim(min(case.panel$obs[,1:2]),max(case.panel$obs[,-1]))+ylab("GLMNET")+xlab("LMMEN")

case.panel$mae$type=factor(case.panel$mae$type,labels=c("LMMEN","GLMNET"))

p_mae=ggplot(case.panel$mae,aes(x=type,y=err,fill=type))+
  geom_boxplot()+
  xlab('Selection Type')+ylab('Mean Absolute Error')+theme_bw(base_size = fs)+
  scale_fill_discrete(guide=FALSE)

case.panel$fixed$Type=factor(case.panel$fixed$Type,labels=c("LMMEN","GLMNET"))

p_eps_f=ggplot(case.panel$fixed,aes(x=var,y=stat,fill=Type))+theme_bw(base_size = fs)+
  geom_bar(stat = "identity")+ theme(axis.text.x = element_blank(),axis.ticks=element_blank())+
  facet_wrap(~Type,nrow=2,scales="free_y")+ylab('% Persistence')+xlab('Fixed Effect')+
  scale_fill_discrete(guide=FALSE)+
  scale_y_continuous(labels = percent)

case.panel$ran$Type=factor(case.panel$ran$Type,labels=c("LMMEN"))

p_eps_r=ggplot(case.panel$ran,aes(x=factor(var),y=stat,fill=Type))+geom_bar(stat = "identity")+ylab('% Persistence')+xlab('Random Effect')+theme_bw(base_size = fs)+
  facet_wrap(~Type)+scale_y_continuous(limits=c(0,1),labels=percent)+  scale_fill_discrete(guide=FALSE)

gggsave(filename = 'figs/lmmen_paper-scenario.pdf',plot = list(p4.1,p4.2),onefile=FALSE,width=14,height=11)

ggsave(filename = 'figs/lmmen_paper-oracle001.pdf',plot = p.compare,width=14,height=8)


pdf('figs/lmmen_paper-casepanel001.pdf',width=18,height=13)
grid.newpage()
pushViewport(viewport(layout=grid.layout(2,2)))
print(p_eps_f+ggtitle("Fixed Effects Persistence (a)"),vp=vplayout(1,1))
print(p_eps_r+ggtitle("Random Effects Persistence (b)"),vp=vplayout(1,2))
print(p_scat+ggtitle("Brier Score Comparison (c)"),vp=vplayout(2,1))
print(p_mae+ggtitle("Mean Absolute Prediction Error (d)"),vp=vplayout(2,2))
dev.off()