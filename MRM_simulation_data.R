

set.wd = 'D:/me_impute/code/githubTest'#-------work directory
setwd(set.wd)


# 1. simulate methylation data by MRM ######################
library(flexmix)
library(data.table)
    
    N = 20 #----number of subject
    noise0.var = 0.2 #----variance of noise level
    pctMi = 0.2 #----missing rate
    I=50 #----number of CpGs per region
    M=100 #----number of regions
    K=4 # number of clusters of subjects
    #L=2 
    
    N_mu = 10 #---- # of rfb centers 
    w =list()
    #a = 1
    #b =10
    #beta_var = a*b/(a+b)^2/(a+b+1)
    mu = seq(-1,1,length.out=N_mu) #---- rfb centers 
    gamma = -10
    
    #noise0.var=0.2
    
    set.seed(1234567)
    df = list()
    for(m in 1: M){
      
      for(k in 1:K){
        #w[[k]]= sample(1:N_mu, N_mu )
        w[[k]]=runif(N_mu, -1, 1)
        #w[[k]]=rbeta(N_mu, a,b,0)
        #w[[k]]=c(1,1,1,1,100)
      }
      C = sample(1:K,size = N, replace = T, prob =c(0.4, 0.3, 0.2, 0.1) )
      
      O3 = list()
      
      for(n in 1:N){
        h = matrix(NA,nrow=I,ncol=N_mu)
        x = sort(runif(I,-1,1))
        p=c()
        for(i in 1:I){
          h[i,] = exp(gamma*(x[i]-mu)^2) 
          noise0 = rnorm(1,0,noise0.var)
          #noise1 = rbeta(1,a,b) 
          p[i] =  w[[C[n]]] %*% h[i,] +   noise0 
        }
        p = (p-min(p))/(max(p)-min(p))
        #plot(x,p)
        #O1 =  data.frame(y = p, x=x)
        #O2 = list(h=h, C=C,w=w)
        O3[[n]] =  data.frame(y = p, x=x,POS=c(1:I), h=h, subj=n, region=m, regionID=m)
      }
      
      df[[m]] = list(data = O3, par = list(C=C, w=w))
      names(df)[m] = paste0('region_', m)
    }
    
    #df----simulated data
    
    


#2. MRM- region ##########################################

source(paste0(set.wd, '/glmnet/R/glmnet.R'))
source(paste0(set.wd,'/functions.R'))
library(glmnet)
library(flexmix)
library(data.table)
library(nnls)
library(doParallel)



simx=1
Imin=50
subjMi = c(1:N)
set.iter.max = 10
set.makeCluster=30
k.region=20
k.subj=5



cl <- makeCluster(set.makeCluster)
registerDoParallel(cl)



time.start = Sys.time() 
temp_df<- foreach(m=1:M, .packages = c('data.table', 'flexmix') ) %dopar%{
  #print(m)
  subjRef = setdiff(c(1:N),subjMi)
  df2 =rbindlist(df[[m]]$data)
  df2 = as.data.frame(df2)
  df2$train = 1
  set.seed(12345)
  for(i_temp in 1: N){
    I=nrow(df[[m]]$data[[i_temp]])
    x_test = sample(I, size = round(pctMi*I), replace = F)
    print(x_test)
    df2[which(df2$subj== subjMi[i_temp]),][x_test,'train'] = 0 
  }
  
  #table(df2$train)
  H = createHmatrix(x=df2$x,gamma=-10, mu=seq(-1,1,length.out = Imin))
  df2_train = data.frame(df2, H)[which(df2$train==1),]
  df2_test = data.frame(df2, H)[which(df2$train==0),]
  xnam = colnames(H)
  fmla <- as.formula(paste("y ~ ", paste(xnam, collapse= "+"),'|subj'))
  Model=FLXMRglmnet(lambda=10^seq(from=-5, to=5, by=1/3),alpha=0)
  m2 <- initFlexmix(fmla, data = df2_train, k = k.region, 
                      model = Model, 
                      control = list(iter.max =set.iter.max), nrep=1)
    #plot(m2)
    if(length(k.region)>1){
      modelFit = getModel(m2, which = "ICL")  
    }else{
      modelFit = m2
    }
    #modelFit = getModel(m2, which = "ICL") 
    #modelFit =m2@models$`3`
    #parameters(modelFit)
    
  y_hat_l = fitted(modelFit)
  y_c= clusters(modelFit)
  y_hat=c()
  for( i_t in 1:length(y_c)){
    y_hat[i_t] = y_hat_l[,y_c[i_t]][i_t]
    
  }
  df2_train$y_hat = y_hat
  df2_train$C.pred = modelFit@cluster
  #cor(df2_train$y,y_hat)
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
  df2 = df2[order(df2$subj, df2$x),]
  df2$res = df2$y_hat - df2$y
  #rmse(df2$y,df2$y_hat)
  
  df2
  
}
time.end = Sys.time()
time.run1  =  difftime(time.end, time.start, unit='mins')

df_regionBased = temp_df


df3 = list()
for(n in subjMi){
  for(m in 1:M){
    df3[[m]]=df_regionBased[[m]][which(df_regionBased[[m]]$subj==n),]
    df3[[m]] = df3[[m]][,c('POS','y','x','subj','region','regionID','train','y_hat','C.pred','res')]
  }
  
}


#3. MRM-subject##########
registerDoParallel(cl)
temp_df = list()
time.start = Sys.time() 

temp_df<- foreach(n=1:N, .packages = c('data.table', 'flexmix') ) %dopar%{
  #print(n)
  for(m in 1:M){
    df3[[m]]=df_regionBased[[m]][which(df_regionBased[[m]]$subj==n),]
    df3[[m]] = df3[[m]][,c('POS', 'y','x','subj','region','regionID','train','y_hat','C.pred','res')]
  }
  df4 = as.data.frame(rbindlist(df3))
  df4$y2In = ifelse(df4$train==1, df4$y, df4$y_hat)
  H = createHmatrix(x=df4$x,gamma=-10, mu=seq(-1,1,length.out = Imin))
  df4 = data.frame(df4, H)
  xnam = colnames(H)
  fmla <- as.formula(paste("y2In ~ ", paste(xnam, collapse= "+"),'|region'))
  Model=FLXMRglmnet(lambda=10^seq(from=-5, to=5, by=1/3))
  
  m2 <- initFlexmix(fmla, data = df4, k = k.subj, 
                      model = Model, 
                      control = list(iter.max = set.iter.max), nrep=1)
    
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
  
  df4
}
time.end = Sys.time()
time.run2  =  difftime(time.end, time.start, unit='mins')
df_subjBased = temp_df



# 4. MRM-stack (region specific weights)######

iter=100
registerDoParallel(cl)
time.start = Sys.time()
temp_df<- foreach(m=1:M, .packages = c('nnls','data.table', 'flexmix') ) %dopar%{
  print(m)
  dfws_out.1 = list()
  for(n in 1:N){
    dfws = df_subjBased[[n]]
    dfws1  = dfws[which(dfws$region==m),]
    xx=matrix(,nrow=iter,ncol=2)
    for(i_iter in 1:iter){
      dfws2 = dfws1[which(dfws1$train==1),]
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
    dfws_out.1[[n]] = dfws1
  }
  dfws_out.1
}
time.end = Sys.time()
time.run3  =  difftime(time.end, time.start, unit='mins')
dfws_out = temp_df



dft = list()
for(n in 1:N){
  dft[[n]] = list()
  for(m in 1:M){
    dft[[n]][[m]]=dfws_out[[m]][[n]]
  }
}



#5. comparison with other methods ##############################################


library(randomForest)
library(impute)


# .1 knn ####

knn_out = list()

for(m in 1:M){
  print(m)
  dfknn1 = rbindlist(dfws_out[[m]])
  dfknn1 = dfknn1[order(dfknn1$subj, dfknn1$x),]
  dim(dfknn1)
  which(dfknn1$train == 0)
  xx = reshape(dfknn1, direction = 'wide', idvar ='POS',timevar ='subj', v.names='y', drop = names(dfknn1)[c(7:63)]  )
  xx_tr = reshape(dfknn1, direction = 'wide', idvar ='POS',timevar ='subj', v.names ='train', drop = names(dfknn1)[c(2,8:63)]  )
  xx1 = xx
  x.missing = matrix(FALSE, nrow = nrow(xx), ncol = ncol(xx))
  x.missing[,5:(5+N-1)] = xx_tr[,5:(5+N-1)] == 0
  xx1[x.missing] = NA
  x.missing2 = x.missing[,5:(5+N-1)]
  xx_imp = impute.knn(data = as.matrix(xx1[,5:(5+N-1)]),  rowmax = 1, colmax = 1)
  xx_imp1 = cbind(xx[,1],xx_imp$data)
  xx_imp2 = reshape(xx_imp1,direction = 'long', varying = list(names(xx_imp1)[2:(2+N-1)]), timevar = 'subj')[,1:3]
  names(xx_imp2)[3] = 'y_hat_knn'
  xx_imp2$subj_pos = paste0(xx_imp2$subj,'_', xx_imp2$POS)
  dfknn1$subj_pos = paste0(dfknn1$subj,'_', dfknn1$POS)
  dfknn2 = merge(dfknn1, xx_imp2, by = 'subj_pos', all.x = T, all.y = F)
  dfknn2 = dfknn2[,c(2:12,63,64,67)]
  names(dfknn2)[c(1,4)] = c('POS','subj')
  knn_out[[m]] = dfknn2
}




#.2 weighted average of one upsteam and one downstream ################

dfknn1 = rbindlist(dfws_out[[m]])
dfupd.l = list()#upstream and downstream
for(n in 1:N){
  dfupd.l[[n]] = list()
  for(m in 1:M){
    dfupd.l[[n]][[m]] = dfws_out[[m]][[n]]
    
  }
}


dfupd.l2 = list()
for(n in 1:N){
  print(n)
  dfupd = as.data.frame(rbindlist(dfupd.l[[n]]))
  dfupd = dfupd[order(dfupd$region, dfupd$x),]
  dfupd = dfupd[,c(1:11,62,63)]
  dfupd$y_hat_upd = NA
  for(t in which(dfupd$train==0)){
    if(t==1)  next
    x1 = t-1
    while(dfupd[x1,'train']==0 &&x1>0){
      x1=x1-1
    }
    if(x1<1) next
    
    x2 = t+1
    while(dfupd[x2,'train']==0 &&x2<nrow(dfupd)){
      x2=x2+1
    }
    if(x2==nrow(dfupd)) next
    
    y1 = dfupd[x1,'y']
    y2 = dfupd[x2,'y']
    d1 = t-x1
    d2 = x2-t
    dfupd$y_hat_upd[t] =d1/(d1+d2)*y1+d2/(d1+d2)*y2 
    
  }
  dfupd.l2[[n]] = dfupd  
}





#.3 RF ###############

dfrf.l = list()#upstream and downstream

for(n in 1:N){
  dfrf.l[[n]] = list()
  for(m in 1:M){
    dfrf.l[[n]][[m]] = dfws_out[[m]][[n]]
    
  }
}


dfrf.l2 = list()
cl <- makeCluster(set.makeCluster)
registerDoParallel(cl)

dfrf.l2<- foreach(n=1:N, .packages = c('data.table', 'randomForest') ) %dopar%{
  print(n)
  dfrf = as.data.frame(rbindlist(dfrf.l[[n]]))
  dfrf = dfrf[order(dfrf$region, dfrf$x),]
  dfrf = dfrf[,c(1:11,62,63)]
  dfrf = dfrf[order(dfrf$region,dfrf$POS),]
  dfrf$y1 = NA
  dfrf$y2 = NA
  dfrf$d1 = NA
  dfrf$d2 = NA
  for(t in 1:nrow(dfrf)){
    if(t==1)  next
    x1 = t-1
    while(dfrf[x1,'train']==0 &&x1>0){
      x1=x1-1
    }
    if(x1<1) next
    
    x2 = t+1
    while(dfrf[x2,'train']==0 &&x2<nrow(dfrf)){
      x2=x2+1
    }
    if(x2==nrow(dfrf)) next
    
    dfrf$y1[t] = dfrf[x1,'y']
    dfrf$y2[t] = dfrf[x2,'y']
    dfrf$d1[t] = t-x1
    dfrf$d2[t] = x2-t
      }
  
  model.1 = y~ y1+y2+d1+d2
  train.comp = dfrf[which(dfrf$train==1),c('region','POS','y', 'y1','y2', 'd1', 'd2')]
  train.comp = train.comp[which(complete.cases(train.comp)),]
  fit.1 <- randomForest(model.1 , data=train.comp, importance=TRUE,ntree=1000)
  y_hat_rf<-predict(fit.1,newdata=dfrf)
  dfrf = cbind(dfrf, y_hat_rf)
  dfrf
}





#.4 Melissa ####################################

library(Melissa)

dft.2 = df_regionBased
dft.train = list()
for(n in 1:N){
  dft.train[[n]] = list()
  for(m in 1:M){
    dft.1 = dft.2[[m]][which(dft.2[[m]]$subj==n),]
    dft.1 = dft.1[which(dft.1$train==1),c('x','y','POS')]
    row.names(dft.1)<- dft.1$POS
    dft.1$POS<-NULL
    dft.1$y = ifelse(dft.1$y>0.5, 1, 0)
    dft.train[[n]][[m]]=as.matrix(dft.1)
  }
  
}


dft.test = list()
for(n in 1:N){
  dft.test[[n]] = list()
  for(m in 1:M){
    dft.1 = dft.2[[m]][which(dft.2[[m]]$subj==n),]
    #dft.1 = dft.1[which(dft.1$train==0),c('x','y')]
    dft.1 = dft.1[which(dft.1$train==0),c('x','y','POS')]
    #row.names(dft.1)<- paste0(n, '_',dft.1$POS)
    row.names(dft.1)<- paste0(n, '_',m, '_',dft.1$POS)
    dft.1$POS<-NULL
    dft.1$y = ifelse(dft.1$y>0.5, 1, 0)
    dft.test[[n]][[m]]=as.matrix(dft.1)
  }
  
}

if(T){
  dft.test.c = list()
  for(n in 1:N){
    dft.test.c[[n]] = list() #-----------continuous
    for(m in 1:M){
      dft.1 = dft.2[[m]][which(dft.2[[m]]$subj==n),]
      dft.1 = dft.1[which(dft.1$train==0),c('x','y','POS')]
      row.names(dft.1)<- dft.1$POS
      dft.1$POS<-NULL
      #dft.1 = dft.1[which(dft.1$train==0),c('x','y')]
      #dft.1$y = ifelse(dft.1$y>0.5, 1, 0)
      dft.test.c[[n]][[m]]=as.matrix(dft.1)
    }
    
  }
  
}
time.start.mel = Sys.time()
basis_obj <- BPRMeth::create_rbf_object(M = 3)
melissa_obj <- melissa(X = dft.train,
                       K=4,
                       basis = basis_obj, 
                       vb_max_iter = 10,
                       vb_init_nstart = 1,
                       is_parallel = FALSE, 
                       is_verbose = FALSE)
imputation_obj <- impute_met_state(obj = melissa_obj, test = dft.test)
time.end.mel = Sys.time()
time.run.mel = difftime(time.end.mel, time.start.mel, units = 'mins')

#6. format data###############################################


dfrf_allsubj=rbindlist(dfrf.l2)
dfknn_allsubj = rbindlist(knn_out)
dfupd_allsubj = rbindlist(dfupd.l2)
  
dfrf_allsubj$POS1 = paste0(dfrf_allsubj$subj, '_', dfrf_allsubj$region, '_', dfrf_allsubj$POS)
dfupd_allsubj$POS1 = paste0(dfupd_allsubj$subj, '_', dfupd_allsubj$region, '_', dfupd_allsubj$POS)
dfknn_allsubj$POS1 = paste0(dfknn_allsubj$subj, '_',dfknn_allsubj$region, '_', dfknn_allsubj$POS)
dfrf_allsubj.1 = dfrf_allsubj[,c('POS1', 'y_hat_rf')]
dfupd_allsubj.1 = dfupd_allsubj[,c('POS1', 'y_hat_upd')]
df_methods = merge(dfknn_allsubj, dfrf_allsubj.1, by.x = 'POS1', by.y='POS1', all=T )
df_methods = merge(df_methods, dfupd_allsubj.1, by.x = 'POS1', by.y='POS1', all=T )
  
df01 = df_methods[which(df_methods$train==0),]
df02 = df01[order(df01$subj,df01$region, df01$POS),]
  
pred_obs = matrix(NA, nrow=length(imputation_obj$pred_obs), ncol=2)
pred_obs[,2] = imputation_obj$pred_obs
pred_obs[,1] = names(imputation_obj$pred_obs)
pred_obs=as.data.frame(pred_obs, stringsAsFactors =F)
names(pred_obs)<-c('POS1', 'y_mel')
pred_obs$y_mel<-as.numeric(pred_obs$y_mel)
  
  
  
  df03 = merge(df02,pred_obs, by.x='POS1', by.y='POS1', all.x=T)
  #which(is.na(df03$y_mel))
  #which(is.na(df03$y_hat))
  if(length(which(is.na(df03$y_mel)))>0){
    df04 = df03[-which(is.na(df03$y_mel)),]
  }else{
    df04= df03
  }
  
head(df04) 
names(df04)

#7 calculate correlation, RMSE, AUC###########################################

library(data.table)
library(randomForest)
library(impute)
library(doParallel)
library(ROCR)
library(pROC)


col1 = names(df04)[c(9,14,15,17,16,18)]
aa = matrix(NA, ncol = length(col1), nrow = N )
aa=data.frame(aa)
colnames(aa) <- col1

#.1 AUC ####
for(n in 1:N){
  print(n)
  dft2 = as.data.frame(df04[which(df04$subj==n),])
  
  for(ic in 1: length(col1)){
    xx1 = dft2[,col1[ic]]
    xx2=ifelse(dft2[,'y']>0.5,1,0)
    xx3 = which(is.na(xx1))
    if(length(xx3)>0){
      xx1=xx1[-xx3]
      xx2=xx2[-xx3]
    }
    pred <- prediction( xx1, xx2)
    perf <- performance( pred, "tpr", "fpr" )
    aa[n,col1[ic]] = performance(pred,"auc")@y.values[[1]]
    
  }
}  
aa_out = aa
names(aa_out) = paste0('auc_', names(aa))
#f_O_a = 'compareAUC.txt'
#write.table(aa_out, f_O_a, col.names=T, row.names=F, sep='\t', quote=F)

#.2 correlation#########

col1 = names(df04)[c(9,14,15,17,16,18)]
aa = matrix(NA, ncol = length(col1), nrow = N )
aa=data.frame(aa)
colnames(aa) <- col1

for(n in 1:N){
  print(n)
  dft2 = as.data.frame(df04[which(df04$subj==n),])
  
  for(ic in 1:length(col1)){
    xx1 = dft2[,col1[ic]]
    xx2 = dft2[,'y']
    xx3 = which(is.na(xx1))
    if(length(xx3)>0){
      xx1=xx1[-xx3]
      xx2=xx2[-xx3]
    }
    aa[n,col1[ic]]<-cor(xx1,xx2)
  }
}

cc_out = aa
names(cc_out) = paste0('cor_', names(aa))
#f_O_c = 'compareCor.txt'
#write.table(cc_out, f_O_c, col.names=T, row.names=F, sep='\t', quote=F)


#.3 rmse######
col1 = names(df04)[c(9,14,15,17,16,18)]
aa = matrix(NA, ncol = length(col1), nrow = N )
aa=data.frame(aa)
colnames(aa) <- col1

for(n in 1:N){
  print(n)
  dft2 = as.data.frame(df04[which(df04$subj==n),])
  for(ic in 1:length(col1)){
    xx1 = dft2[,col1[ic]]
    xx2 = dft2[,'y']
    xx3 = which(is.na(xx1))
    if(length(xx3)>0){
      xx1=xx1[-xx3]
      xx2=xx2[-xx3]
    }
  aa[n,col1[ic]]<-rmse(xx1,xx2)
  }
}

rr_out = aa
names(rr_out) = paste0('rmse_', names(aa))
#f_O_r = 'compareRMSE.txt'
#write.table(rr_out, f_O_r, col.names=T, row.names=F, sep='\t', quote=F)




#.4 t-test ##########################
library(psych)

t1 = apply(aa_out, 2, mean)
t2 = apply(aa_out, 2, sd)
t3 = rep(NA,ncol(aa_out))
t4 = rep(NA,ncol(aa_out))
for(i in 1:length(t3)){
  xx = t.test(aa_out$auc_region_k4, aa_out[,i],mu=0, alternative = 'two.sided')
  t3[i] = xx$statistic
  t4[i] = xx$p.value
}
aa_out2 = rbind( t1,t2,t3,t4)
rownames(aa_out2)<-c('mean','sd','t','p.value')
aa_out2 = t(aa_out2)
#f_O_a2 = 'compareAUC_ttest.txt'
#write.table(aa_out2, f_O_a2, col.names=T, row.names=T, sep='\t', quote=F)



t1 = apply(cc_out, 2, mean)
t2 = apply(cc_out, 2, sd)
t3 = rep(NA,ncol(cc_out))
t4 = rep(NA,ncol(cc_out))
for(i in 1:length(t3)){
  xx = t.test(cc_out$cor_region_k4, cc_out[,i],mu=0, alternative = 'two.sided')
  t3[i] = xx$statistic
  t4[i] = xx$p.value
}
cc_out2 = rbind( t1,t2,t3,t4)
rownames(cc_out2)<-c('mean','sd','t','p.value')
cc_out2 = t(cc_out2)
#f_O_c2 = 'compareCor_ttest.txt'
#write.table(cc_out2, f_O_c2, col.names=T, row.names=T, sep='\t', quote=F)




t1 = apply(rr_out, 2, mean)
t2 = apply(rr_out, 2, sd)
t3 = rep(NA,ncol(rr_out))
t4 = rep(NA,ncol(rr_out))
for(i in 1:length(t3)){
  xx = t.test(rr_out[,3], rr_out[,i],mu=0, alternative = 'two.sided')
  t3[i] = xx$statistic
  t4[i] = xx$p.value
}
rr_out2 = rbind( t1,t2,t3,t4)
rownames(rr_out2)<-c('mean','sd','t','p.value')
rr_out2 = t(rr_out2)
#f_O_c2 = 'compareRmse_ttest.txt'
#write.table(rr_out2, f_O_r2, col.names=T, row.names=T, sep='\t', quote=F)



# 8. specific test for cor/rmse/auc and format ###############################################
library(pROC)

# .1. Delong's test for auc ####################


test.auc =matrix(NA,ncol = 5, nrow=6)
rownames(test.auc) = names(aa_out)
test.auc [1,] = NA


test1 = roc.test(response = ifelse(df04$y>0.5,1,0), predictor1 = df04$y_hat, predictor2 = df04$y_hat)
test.auc['auc_region',1:2] = test1$estimate
test.auc['auc_region',3] = test1$statistic
test.auc['auc_region',4] = test1$p.value
test.auc['auc_region',5] = test1$method


test1 = roc.test(response = ifelse(df04$y>0.5,1,0), predictor1 = df04$y_hat, predictor2 = df04$y_hat_3)
test.auc['auc_ws',1:2] = test1$estimate
test.auc['auc_ws',3] = test1$statistic
test.auc['auc_ws',4] = test1$p.value
test.auc['auc_ws',5] = test1$method

test1 = roc.test(response = ifelse(df04$y>0.5,1,0), predictor1 = df04$y_hat, predictor2 = df04$y_hat_knn)
test.auc['auc_knn',1:2] = test1$estimate
test.auc['auc_knn',3] = test1$statistic
test.auc['auc_knn',4] = test1$p.value
test.auc['auc_knn',5] = test1$method

test1 = roc.test(response = ifelse(df04$y>0.5,1,0), predictor1 = df04$y_hat, predictor2 = df04$y_hat_upd)
test.auc['auc_upd',1:2] = test1$estimate
test.auc['auc_upd',3] = test1$statistic
test.auc['auc_upd',4] = test1$p.value
test.auc['auc_upd',5] = test1$method

test1 = roc.test(response = ifelse(df04$y>0.5,1,0), predictor1 = df04$y_hat, predictor2 = df04$y_hat_rf)
test.auc['auc_rf',1:2] = test1$estimate
test.auc['auc_rf',3] = test1$statistic
test.auc['auc_rf',4] = test1$p.value
test.auc['auc_rf',5] = test1$method

test1 = roc.test(response = ifelse(df04$y>0.5,1,0), predictor1 = df04$y_hat, predictor2 = df04$y_mel)
test.auc['auc_mel',1:2] = test1$estimate
test.auc['auc_mel',3] = test1$statistic
test.auc['auc_mel',4] = test1$p.value
test.auc['auc_mel',5] = test1$method


colnames(test.auc) = c('estimate1','estimate2','statistic','p.value','method')

test.auc = test.auc[c(1,2,4,3,5,6),]
rownames(test.auc) <- c('MRM_region', 'MRM_stacked','average', 'KNN', 'RF', 'Melissa' )
test.auc

t1 = gregexpr('_', lsName)[[1]][1]
t2 = gregexpr('RD', lsName)[[1]][1]
f_O_test.auc = paste0('', substr(lsName,t1,t2-1 ), 'txt')
f_O_test.auc0 = 'test.auc.txt'

#cat(paste0(f_O_test.auc, '\n'), file = f_O_test.auc0, append=T)
#write.table(test.auc, f_O_test.auc0, col.names = NA, row.names = T, sep='\t', quote=F, append=T)



# .2 Fisher's test Correlations  ###############  

#install.packages('psych', dependencies = T)
#install.packages('D:/me_impute/code/VB/mnormt_1.5-7.zip', type='binary', repos = NULL)
#install.packages('tmvnsim')


library(psych)

#r.test(n, r12, r34 = NULL, r23 = NULL, r13 = NULL, r14 = NULL, r24 = NULL, n2 = NULL,pooled=TRUE, twotailed = TRUE)

test.cor =matrix(NA,ncol = 5, nrow=6)
rownames(test.cor) = names(cc_out)


r12=cor(df04$y, df04$y_hat)
r13 = cor(df04$y, df04$y_hat)
r23 = cor(df04$y_hat, df04$y_hat)
test1 = r.test(n=nrow(df04), r12=r12, r34 = NULL, r23 = r23, r13 = r13 , r14 = NULL, r24 = NULL,  n2 = NULL, pooled=TRUE, twotailed = TRUE)
test.cor['cor_region',1:2] = c(r12, r13)
test.cor['cor_region',3] = test1$t
test.cor['cor_region',4] = test1$p
test.cor['cor_region',5] = test1$Test

r12=cor(df04$y, df04$y_hat)
r13 = cor(df04$y, df04$y_hat_3)
r23 = cor(df04$y_hat, df04$y_hat_3)
test1 = r.test(n=nrow(df04), r12=r12, r34 = NULL, r23 = r23, r13 = r13 , r14 = NULL, r24 = NULL,  n2 = NULL, pooled=TRUE, twotailed = TRUE)
test.cor['cor_ws',1:2] = c(r12, r13)
test.cor['cor_ws',3] = test1$t
test.cor['cor_ws',4] = test1$p
test.cor['cor_ws',5] = test1$Test

r12=cor(df04$y, df04$y_hat)
r13 = cor(df04$y, df04$y_hat_knn)
r23 = cor(df04$y_hat, df04$y_hat_knn)
test1 = r.test(n=nrow(df04), r12=r12, r34 = NULL, r23 = r23, r13 = r13 , r14 = NULL, r24 = NULL,  n2 = NULL, pooled=TRUE, twotailed = TRUE)
test.cor['cor_knn',1:2] = c(r12, r13)
test.cor['cor_knn',3] = test1$t
test.cor['cor_knn',4] = test1$p
test.cor['cor_knn',5] = test1$Test

r12=cor(df04$y, df04$y_hat)
r13 = cor(df04$y, df04$y_hat_upd, use='complete.obs')
r23 = cor(df04$y_hat, df04$y_hat_upd, use='complete.obs')
test1 = r.test(n=nrow(df04), r12=r12, r34 = NULL, r23 = r23, r13 = r13 , r14 = NULL, r24 = NULL,  n2 = NULL, pooled=TRUE, twotailed = TRUE)
test.cor['cor_upd',1:2] = c(r12, r13)
test.cor['cor_upd',3] = test1$t
test.cor['cor_upd',4] = test1$p
test.cor['cor_upd',5] = test1$Test

r12=cor(df04$y, df04$y_hat)
r13 = cor(df04$y, df04$y_hat_rf,use='complete.obs')
r23 = cor(df04$y_hat, df04$y_hat_rf,use='complete.obs')
test1 = r.test(n=nrow(df04), r12=r12, r34 = NULL, r23 = r23, r13 = r13 , r14 = NULL, r24 = NULL,  n2 = NULL, pooled=TRUE, twotailed = TRUE)
test.cor['cor_rf',1:2] = c(r12, r13)
test.cor['cor_rf',3] = test1$t
test.cor['cor_rf',4] = test1$p
test.cor['cor_rf',5] = test1$Test

r12=cor(df04$y, df04$y_hat)
r23 = cor(df04$y_hat, df04$y_mel)
r13 = cor(df04$y, df04$y_mel)
test1 = r.test(n=nrow(df04), r12=r12, r34 = NULL, r23 = r23, r13 = r13 , r14 = NULL, r24 = NULL,  n2 = NULL, pooled=TRUE, twotailed = TRUE)
test.cor['cor_mel',1:2] = c(r12, r13)
test.cor['cor_mel',3] = test1$t
test.cor['cor_mel',4] = test1$p
test.cor['cor_mel',5] = test1$Test


colnames(test.cor) = c('estimate1','estimate2','statistic','p.value','method')

test.cor = test.cor[c(1,2,4,3,5,6),]
rownames(test.cor) <- c('MRM_region', 'MRM_stacked','average', 'KNN', 'RF', 'Melissa' )
test.cor

t1 = gregexpr('_', lsName)[[1]][1]
t2 = gregexpr('RD', lsName)[[1]][1]
f_O_test.cor = paste0('', substr(lsName,t1,t2-1 ), 'txt')
f_O_test.cor0 = 'test.cor.txt'

#cat(paste0(f_O_test.cor, '\n'), file = f_O_test.cor0, append=T)
#write.table(test.cor, f_O_test.cor0, col.names = NA, row.names = T, sep='\t', quote=F, append=T)



# .3 Diebold-Mariano testrmse ############################################
#Diebold-Mariano test compares the forecast accuracy of two forecast methods
#install.packages('forecast')
library(forecast)


test.rmse =matrix(NA,ncol = 5, nrow=6)
rownames(test.rmse) = names(rr_out)


r12=rmse(df04$y, df04$y_hat)
r13 = rmse(df04$y, df04$y_hat_3)
test1 = dm.test(  e1= df04$y_hat- df04$y,  e2= df04$y_hat_3 - df04$y, alternative = "two.sided", h = 1 )
test.rmse['rmse_ws',1:2] = c(r12, r13)
test.rmse['rmse_ws',3] = test1$statistic
test.rmse['rmse_ws',4] = test1$p.value
test.rmse['rmse_ws',5] = test1$method



r12=rmse(df04$y, df04$y_hat)
r13 = rmse(df04$y, df04$y_hat_3)
test1 = dm.test(  e1= df04$y_hat- df04$y,  e2= df04$y_hat_3 - df04$y, alternative = "two.sided", h = 1 )
test.rmse['rmse_ws',1:2] = c(r12, r13)
test.rmse['rmse_ws',3] = test1$statistic
test.rmse['rmse_ws',4] = test1$p.value
test.rmse['rmse_ws',5] = test1$method


r12=rmse(df04$y, df04$y_hat)
r13 = rmse(df04$y, df04$y_hat_knn)
test1 = dm.test(  e1= df04$y_hat- df04$y,  e2= df04$y_hat_knn - df04$y, alternative = "two.sided", h = 1 )
test.rmse['rmse_knn',1:2] = c(r12, r13)
test.rmse['rmse_knn',3] = test1$statistic
test.rmse['rmse_knn',4] = test1$p.value
test.rmse['rmse_knn',5] = test1$method



if(length(is.na(df04$y_hat_upd>0))){
  df04.1 = df04[-which(is.na(df04$y_hat_upd>0)),]
}else{
  df04.1 = df04
}
r12=rmse(df04$y, df04$y_hat)
r13 = rmse(df04.1$y, df04.1$y_hat_upd)
test1 = dm.test(  e1= df04.1$y_hat- df04.1$y,  e2= df04.1$y_hat_upd - df04.1$y, alternative = "two.sided", h = 1 )
test.rmse['rmse_upd',1:2] = c(r12, r13)
test.rmse['rmse_upd',3] = test1$statistic
test.rmse['rmse_upd',4] = test1$p.value
test.rmse['rmse_upd',5] = test1$method




if(length(is.na(df04$y_hat_rf>0))){
  df04.1 = df04[-which(is.na(df04$y_hat_rf>0)),]
}else{
  df04.1 = df04
}
r12=rmse(df04$y, df04$y_hat)
r13 = rmse(df04.1$y, df04.1$y_hat_rf)
test1 = dm.test(  e1= df04.1$y_hat- df04.1$y,  e2= df04.1$y_hat_rf - df04.1$y, alternative = "two.sided", h = 1 )
test.rmse['rmse_rf',1:2] = c(r12, r13)
test.rmse['rmse_rf',3] = test1$statistic
test.rmse['rmse_rf',4] = test1$p.value
test.rmse['rmse_rf',5] = test1$method


r12=rmse(df04$y, df04$y_hat)
r13 = rmse(df04$y, df04$y_mel)
test1 = dm.test(  e1= df04.1$y_hat- df04.1$y,  e2= df04.1$y_mel - df04.1$y, alternative = "two.sided", h = 1 )
test.rmse['rmse_mel',1:2] = c(r12, r13)
test.rmse['rmse_mel',3] = test1$statistic
test.rmse['rmse_mel',4] = test1$p.value
test.rmse['rmse_mel',5] = test1$method



colnames(test.rmse) = c('estimate1','estimate2','statistic','p.value','method')

test.rmse = test.rmse[c(1,2,4,3,5,6),]
rownames(test.rmse) <- c('MRM_region', 'MRM_stacked','average', 'KNN', 'RF', 'Melissa' )
test.rmse

t1 = gregexpr('_', lsName)[[1]][1]
t2 = gregexpr('RD', lsName)[[1]][1]
f_O_test.rmse = paste0('', substr(lsName,t1,t2-1 ), 'txt')
f_O_test.rmse0 = 'test.rmse.txt'

#cat(paste0(f_O_test.rmse, '\n'), file = f_O_test.rmse0, append=T)
#write.table(test.rmse, f_O_test.rmse0, col.names = NA, row.names = T, sep='\t', quote=F, append=T)







