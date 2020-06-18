
setwd('D:/me_impute/code/githubTest')
source('./glmnet/R/glmnet.R')

source('./functions.R')
#source('/home/fyu/me_impute/realData_PBM/realData_v03/s2_realData_v03_noPreprocess_159_tempSource0.R')
library(TSPred)
library(snpar)
library(glmnet)
library(flexmix)
library(data.table)
library(nnls)


# 1.load and format annotated WGBS data #######
df1 = read.table('./df_beta_3_annot_N=19.txt', header=T, sep='\t', stringsAsFactors=F)

#head(df1)
df1$region = df1$tss0
Imin=50
I=40
N=19
K=4
L=2
M=50
regionID = as.integer(names(table(df1$region))[which(table(df1$region)>=Imin)])
regionID = regionID[order(table(df1$region[which(df1$region %in% regionID)]))]
M1=length(regionID)
subjID = names(df1)[2:20]


df1$regionID1='0'
xx = c()
for(m in 1:length(regionID)){
  xx1=which(df1$region==regionID[m])
  xx4 = df1$tss0[xx1][1]
  xx5 = df1$start[xx1]
  #xx2=slidingWindows(xx5,Imin) 
  xx6 = rank(c(xx4, xx5))[1]
  if(xx6<25){
    xx9=1
  }else if(length(xx5)-xx6<25){
    xx9=length(xx5)-49
  }else{
    xx9=xx6-24
  }
  
  xx10 = xx5[xx9:(xx9+49)]
  rank(c(xx4, xx10))
  df1[which(df1$start %in% xx10),'regionID1']=m
  
}
table(df1$regionID1)



dfIn$promoter=NULL
dfIn$tss0=NULL

set.seed(1234567)
regionID = sample(unique(df1$regionID1[-1]), size=M)

df=list()
for(m in 1:M){
  df[[m]]=list()
  for(n in 1:N){
    df[[m]][[n]]=data.frame(POS=dfIn[which(dfIn$regionID1==regionID[m]),'start'], y=dfIn[which(dfIn$regionID1==regionID[m]),subjID[n]])
    df[[m]][[n]]$x=scale01(df[[m]][[n]]$POS)
    df[[m]][[n]]$subj=n
    df[[m]][[n]]$region = m
    df[[m]][[n]]$regionID = regionID[m]
    if(length(which(is.na(df[[m]][[n]]$y)))>0){
      df[[m]][[n]]=df[[m]][[n]][-which(is.na(df[[m]][[n]]$y)),]
    }
    
    
  }
  
}


#getwd()





#2.MRM #######################################
library(doParallel)

#2.1 Regional model##########################

k.region = 2:4  #.region.max #N
k.subj = c(2,5,10) #M
pctMi=0.6
simx=1
set.iter.max = 5
set.makeCluster=30


cl <- makeCluster(set.makeCluster)
registerDoParallel(cl)

#time.start = Sys.time()
df_regionBased<- foreach(m=1:M, .packages = c('data.table', 'flexmix') ) %dopar%{
  print(m)
  subjMi = c(1:N)
  subjRef = setdiff(c(1:N),subjMi)
  df2 =rbindlist(df[[m]])
  df2 = as.data.frame(df2)
  df2$train = 1
  set.seed(12345)
  for(n_temp in 1: N){
    I=nrow(df[[m]][[n_temp]])
    x_test = sample(I, size = round(pctMi*I), replace = F)
    df2[which(df2$subj== subjMi[n_temp]),][x_test,'train'] = 0 
  }
  H = createHmatrix(x=df2$x,gamma=-10, mu=seq(-1,1,length.out = Imin))
  
  df2_train = data.frame(df2, H)[which(df2$train==1),]
  df2_test = data.frame(df2, H)[which(df2$train==0),]
  
  xnam = colnames(H)
  fmla <- as.formula(paste("y ~ ", paste(xnam, collapse= "+"),'|subj'))
  
  Model=FLXMRglmnet(lambda=10^seq(from=-5, to=5, by=1/3),alpha=0)
  
  set.seed(123456789)
  m2 <- initFlexmix(fmla, data = df2_train, k = k.region, 
                      model = Model, 
                      control = list(iter.max = set.iter.max), nrep=1)
  if(length(k.region)>1){
      modelFit = getModel(m2, which = "ICL")  
    }else{
      modelFit = m2
    }
    
  
  y_hat_l = fitted(modelFit)
  y_c= clusters(modelFit)
  y_hat=c()
  for( i_t in 1:length(y_c)){
    y_hat[i_t] = y_hat_l[,y_c[i_t]][i_t]
    
  }
  df2_train$y_hat = y_hat
  df2_train$C.pred = modelFit@cluster
  cor(df2_train$y,y_hat)
  
  
  y_hat_l = predict(modelFit, newdata=df2_test)
  y_c01 = unique(df2_train[df2_train$subj %in% subjMi,c('subj','C.pred')])
  
  y_c = rep(y_c01$C.pred, times=table(df2_test$subj))
  y_hat=c()
  for( i_t in 1:length(y_c)){
    y_hat[i_t] = y_hat_l[[y_c[i_t]]][i_t]
    
  }
  df2_test$y_hat = y_hat
  df2_test$C.pred = y_c
  #cor(df2_test$y, df2_test$y_hat)
  #plot(df2_test$y, df2_test$y_hat)
  #rmse(df2_test$y, df2_test$y_hat)
  df2 = rbind(df2_train, df2_test)
  df2$res = df2$y_hat - df2$y
 
  df2
}
#time.end = Sys.time()
#time.run1 = difftime(time.end, time.start, unit='mins')

df3 = list()
subjMi = c(1:N)
for(n in subjMi){
  for(m in 1:M){
    df3[[m]]=df_regionBased[[m]][which(df_regionBased[[m]]$subj==n),]
    df3[[m]] = df3[[m]][,c('POS','y','x','subj','region','regionID','train','y_hat','C.pred','res')]
  }
  
}

#2.2 Subject model##############
cl <- makeCluster(30)
registerDoParallel(cl)
temp_df = list()

#time.start = Sys.time() 
temp_df<- foreach(n=1:N, .packages = c('data.table', 'flexmix') ) %dopar%{
  
  print(n)
  
  for(m in 1:M){
    df3[[m]]=df_regionBased[[m]][which(df_regionBased[[m]]$subj==n),]
    #df3[[m]]$region = m
    df3[[m]] = df3[[m]][,c('POS', 'y','x','subj','region','regionID','train','y_hat','C.pred','res')]
    
    #df3[[m]] = df3[[m]][,c( 'y','x','subj','region','train','y_hat','C.pred','res')]
  }
  
  
  df4 = as.data.frame(rbindlist(df3))
  #head(df4)
  table(df4$train)
  df4$y2In = ifelse(df4$train==1, df4$y, df4$y_hat)
  
  
  H = createHmatrix(x=df4$x,gamma=-10, mu=seq(-1,1,length.out = Imin))
  
  df4 = data.frame(df4, H)
  
  
  
  #m4 <- flexmix(fmla, data=df4, k=4)
  #xnam = paste0('H',c(1:I))
  xnam = colnames(H)
  fmla <- as.formula(paste("y2In ~ ", paste(xnam, collapse= "+"),'|region'))
  #fmla <- as.formula(paste("y ~ ", paste(xnam, collapse= "+")))
  #df4_train = df4[which(df4$train==1),]
  #df4_test = df4[which(df4$train==0),]
  
  #x01=cbind(df4[,xnam])
  #y01= df4$y
  
  #fo <- sample(rep(seq(10), length = nrow(df4)))
  #data01=data.frame(y01,x01)
  Model=FLXMRglmnet(lambda=10^seq(from=-5, to=5, by=1/3))
  
  
  if(TRUE){
    if(FALSE){
      k_cluster=4
      z <- sample(1:k_cluster, nrow(df4), replace = TRUE)
      m2 <- flexmix(fmla, data = df4, k = 4, 
                    model = Model, cluster=z,
                    control = list(iter.max = 2))
      parameters(m2)  
    }
    
    
    m2 <- initFlexmix(fmla, data = df4, k = k.subj, 
                      model = Model, 
                      control = list(iter.max = set.iter.max), nrep=1)
    
    #plot(m2)
    
	if(length(k.subj)>1){
      modelFit = getModel(m2, which = "ICL")  
    }else{
      modelFit = m2
    }
	
    y_hat_l = fitted(modelFit)
    
    y_c= clusters(modelFit)
    y_hat = c()
    for( i in 1:nrow(df4)){
      y_hat[i] = y_hat_l[,y_c[i]][i]
      
    }
    nLoop = 2
    y_hat_name = paste0('y_hat_', nLoop)
    df4[,y_hat_name] =y_hat
    #df_subjBased[[n]]=df4
    
  }else{
    nLoop = 2
    y_hat_name = paste0('y_hat_', nLoop)
    df4[,y_hat_name] =df4$y_hat
    #df_subjBased[[n]]=df4
  }
  
  #names(df4)
  #cor(df4$y[which(df4$train==0)],df4$y_hat[which(df4$train==0)])
  #cor(df4$y[which(df4$train==0)],df4$y_hat_2[which(df4$train==0)])
  df4
}
#time.end = Sys.time()
#time.run2  =  difftime(time.end, time.start, unit='mins')
df_subjBased = temp_df



#2.3 stacked model#######################

iter=100
registerDoParallel(cl)
#time.start = Sys.time()
temp_df<- foreach(m=1:M, .packages = c('nnls','data.table', 'flexmix') ) %dopar%{
  print(m)
  #dfws_out[[m]] = list()
  dfws_out.1 = list()
  for(n in 1:N){
    
    dfws = df_subjBased[[n]]
    dfws1  = dfws[which(dfws$region==m),]
    
    xx=matrix(,nrow=iter,ncol=2)
    for(i_iter in 1:iter){
      dfws2 = dfws1[which(dfws1$train==1),]
      #dfws2 = dfws2[sample(nrow(dfws2),size=round(0.8*nrow(dfws2))),]
      dfws2 = dfws2[sample(nrow(dfws2),size=round(0.9*nrow(dfws2))),]
      A=as.matrix(dfws2[,c('y_hat','y_hat_2')])
      b=dfws2[,'y']
      xx1=nnls(A,b)
      xx[i_iter,]=xx1$x
      
    }
    w1=mean(xx[,1])
    w2=mean(xx[,2])
    dfws1$y_hat_3 = w1*dfws1$y_hat + w2*dfws1$y_hat_2
    dfws1 = dfws1[order(dfws1$subj, dfws1$x),]
    #print(c(cor(dfws1$y,dfws1$y_hat), cor(dfws1$y,dfws1$y_hat_3)))
    #print(c(cor(dfws1$y[which(dfws1$train==0)],dfws1$y_hat[which(dfws1$train==0)]), cor(dfws1$y[which(dfws1$train==0)],dfws1$y_hat_3[which(dfws1$train==0)])))
    
    #dfws_out[[m]][[n]] = dfws1 
    dfws_out.1[[n]] = dfws1
  }
  dfws_out.1
}
#time.end = Sys.time()
#time.run3 =  difftime(time.end, time.start, unit='mins')
dfws_out = temp_df




dft = list()
for(n in 1:N){
  dft[[n]] = list()
  for(m in 1:M){
    dft[[n]][[m]]=dfws_out[[m]][[n]]
    
  }
  
}



MRM_out = list()
for(m in 1:M){
  print(m)
  df1 = rbindlist(dfws_out[[m]])
  df1 = df1[order(df1$subj, df1$x),]
  df1[,c(10:61)]<-NULL
  MRM_out[[m]] = df1
}




