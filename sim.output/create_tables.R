temp=left_join(
  real_param,
  fit.min.all%>%
    filter(type!="BIC")%>%
    mutate(type=factor(type)),
  by=c("ex","type","id"))

#Selection Table

#Choose at least real params
partial<-
  left_join(
    #Model
      temp%>%filter(p==1)%>%mutate(check=(est!=0))%>%group_by(ex,simtype)%>%summarise(MODEL=sum(check)/n()),
    #Effects
      temp%>%filter(p==1)%>%mutate(check=(est!=0))%>%group_by(ex,simtype,type)%>%summarise(checksum=sum(check)/n())%>%reshape2::dcast(ex+simtype~type,value.var='checksum'),
      by=c('ex','simtype'))%>%ungroup()

#Choose exactly real params
oracle<-left_join(
  #Model
    temp%>%mutate(check=(p==1&est!=0))%>%group_by(ex,simtype)%>%summarise(MODEL=sum(check)/n()),
  #Effects
    temp%>%mutate(check=(p==1&est!=0))%>%group_by(ex,simtype,type)%>%summarise(checksum=sum(check)/n())%>%reshape2::dcast(ex+simtype~type,value.var='checksum'),
  by=c('ex','simtype'))%>%ungroup()

names(partial)<-c('Scenario','Method','C (subset)','CF (subset)','CR (subset)')
names(oracle)<-c('Scenario','Method','C (oracle)','CF (oracle)','CR (oracle)')

tbl<-partial%>%left_join(oracle,by=c('Scenario','Method'))%>%mutate(Scenario=factor(Scenario))

xtbl<-xtable::xtable(tbl)

texPreview::texPreview(xtbl,stem = 'selection',fileDir = 'references/tbls',
                       margin = list(left =10, top = 5, right = 400, bottom = 5),
                       print.xtable.opts = list(include.rownames = FALSE))


#Estimation Table

#FE 1-4
  tbl_est_fe<-temp%>%
    filter(ex<5&type=='FIXED')%>%
    group_by(ex,simtype,id)%>%
    summarise(med=round(median(est,na.rm = TRUE),1),
              ci=sprintf('(%s,%s)',
                           round(quantile(est,probs = 0.25,na.rm = TRUE),2),
                           round(quantile(est,probs = 0.75,na.rm = TRUE),2)
                           ))%>%
    reshape2::melt(.,id=c('ex','simtype','id'))%>%
    reshape2::dcast(ex+simtype+variable~id,value.var='value')%>%
    mutate(simtype=ifelse(variable=='ci',NA,sprintf('\\multirow{2}{*}{%s}',simtype)))%>%
    select(-variable)
  
  names(tbl_est_fe)<-c('Scenario','Method',sprintf('$\\hat{\\beta}_{%s}$',1:9))
  
  tbl_est_fe<-tbl_est_fe%>%mutate(Scenario=factor(Scenario))
  
  xtbl_est_fe<-xtable::xtable(tbl_est_fe)
  
  texPreview::texPreview(xtbl_est_fe,stem = 'fe14',fileDir = 'references/tbls',
                         margin = list(left =10, top = 5, right = 10, bottom = 5),
                         print.xtable.opts = list(include.rownames = FALSE,
                                                  sanitize.colnames.function=identity,
                                                  sanitize.text=identity),
                         imgFormat = 'svg')
  

  #FE 5
  tbl_est_fe5<-temp%>%
    mutate(id5=ifelse(id<=20,id,20.1))%>%
    filter(ex==5&type=='FIXED')%>%
    group_by(ex,simtype,id5)%>%
    summarise(avg=round(mean(est,na.rm = TRUE),1),
              ci=sprintf('(%s,%s)',
                         round(quantile(est,probs = 0.25,na.rm = TRUE),2),
                         round(quantile(est,probs = 0.75,na.rm = TRUE),2)
              ))%>%
    reshape2::melt(.,id=c('ex','simtype','id5'))%>%
    mutate(id2=ifelse(id5<=10,1,2),id5.header=ifelse(id2==2,id5-10,id5),
    lbl=sprintf('$\\hat{\\beta}_{%s}$',id5))%>%
    select(-id5)%>%rename(stat=variable)%>%
    reshape2::melt(.,id=c('ex','simtype','id2','id5.header','stat'))%>%
    do(.,cbind(row=rep(c(2,4,1,3),c(10,11,10,11)),.))%>%
    mutate(row=ifelse(id5.header>10.1,row+2,row),
           value=ifelse(id5.header>10.1&row==5,'$\\hat{\\beta}_{21-200}$',value),
           id5.header=ifelse(id5.header>10.1,1,id5.header))%>%
    reshape2::dcast(ex+simtype+row~id5.header,value.var='value')%>%select(-c(row))
  
  
  tbl_scen5<-function(temp,ids){
  
  tbl_est_fe5b<-temp%>%
    filter(ex==5&type=='FIXED'&id%in%ids)%>%
    group_by(ex,simtype,id)%>%
    summarise(med=round(mean(est,na.rm = TRUE),1),
              ci=sprintf('(%s,%s)',
                         round(quantile(est,probs = 0.25,na.rm = TRUE),2),
                         round(quantile(est,probs = 0.75,na.rm = TRUE),2)
              ))%>%
    mutate(lbl=sprintf('$\\hat{\\beta}_{%s}$',id),id=ifelse(id%/%10>0,id-(10*id%/%10),id),
           id=ifelse(ids%in%c(10,20),10,id))%>%
    reshape2::melt(.,id=c('ex','simtype','id'))%>%
    reshape2::dcast(ex+simtype+variable~id,value.var='value')%>%
    mutate(simtype=ifelse(variable=='lbl'&1%in%ids,sprintf('\\multirow{9}{*}{%s}',simtype),NA))%>%
    select(-variable)
  
  tbl_est_fe5b=tbl_est_fe5b[c(3,1,2),]
  
  tbl_est_fe5b<-tbl_est_fe5b%>%mutate(nrow=1:n(),n=n())%>%
    ungroup()%>%
    mutate(ex=ifelse(nrow==1&1%in%ids,sprintf('\\multirow{%s}{*}{%s}',3*n,ex),NA))%>%
    select(-c(nrow,n))
  
  tbl_est_fe5b
  
  }
  
  tbl_est_fe5<-bind_rows(tbl_scen5(temp,1:10),tbl_scen5(temp,11:20),tbl_scen5(temp,21))
  
  names(tbl_est_fe5)<-c('Scenario','Method',rep(' ',10))
  
  xtbl_est_fe5<-xtable::xtable(tbl_est_fe5)
  
  texPreview::texPreview(xtbl_est_fe5,stem = 'fe5',fileDir = 'references/tbls',
                         margin = list(left =10, top = 5, right = 400, bottom = 5),
                         print.xtable.opts = list(include.colnames=FALSE,
                                                  include.rownames = FALSE,
                                                  sanitize.text.function=identity),
                         imgFormat = 'svg')
    
#RE
  tbl_est_re<-temp%>%
    filter(type=='STDDEV')%>%
    group_by(ex,simtype,id)%>%
    summarise(med=round(median(est,na.rm = TRUE),1),
              ci=sprintf('(%s,%s)',
                         round(quantile(est,probs = 0.25,na.rm = TRUE),2),
                         round(quantile(est,probs = 0.75,na.rm = TRUE),2)
              ))%>%
    reshape2::melt(.,id=c('ex','simtype','id'))%>%
    reshape2::dcast(ex+simtype+variable~id,value.var='value')%>%
    mutate(simtype=ifelse(variable=='ci',NA,sprintf('\\multirow{2}{*}{%s}',simtype)))%>%
    select(-variable)%>%
    group_by(ex)%>%
    mutate(nrow=1:n(),n=n())%>%
    ungroup()%>%
    mutate(ex=ifelse(nrow==1,sprintf('\\multirow{%s}{*}{%s}',n,ex),NA))%>%
    select(-c(nrow,n))
  
  names(tbl_est_re)<-c('Scenario','Method',sprintf('$\\hat{\\sigma}_{%s}$',1:10))
  
  tbl_est_re<-tbl_est_re%>%mutate(Scenario=factor(Scenario))
  
  xtbl_est_re<-xtable::xtable(tbl_est_re)
  
  texPreview::texPreview(xtbl_est_re,stem = 're',fileDir = 'references/tbls',
                         margin = list(left =10, top = 5, right = 600, bottom = 5),
                         print.xtable.opts = list(include.rownames = FALSE,
                                                  sanitize.colnames.function=identity,
                                                  sanitize.text.function=identity),
                         imgFormat = 'svg')