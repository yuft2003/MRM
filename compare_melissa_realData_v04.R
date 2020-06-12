

#setwd('D:/me_impute/code/sim_v03')
#setwd('D:/me_impute/code/realData_v03')
#setwd('D:/me_impute/code/realData_v03_2')
#setwd('D:/me_impute/code/realData_v03_2/condition2')

setwd('D:/me_impute/code/sim_v04')
ls1 = list.files('./', pattern='compare_condition2')
load('./compare_0.2.RData')

#run patch1 first
#source('..')


rmse <- function(x,y){
  return(sqrt(sum((x-y)^2)/length(x)))
}

# 1.melissa to compute predicted value##############

#setwd('D:/me_impute/code/sim_v03_fixSeed')
#setwd('D:/me_impute/code/sim_v03')
library(data.table)
library(Melissa)
library(ROCR)

#ls1 = read.table('./ls.txt', header=F, sep='\t', stringsAsFactors = F)



#ls1 = list.files('./', pattern='compare_condition2')
#ix=1

#aucMel.l = list()

#for(ix in 1:length(ls1)){
  ##print(ix)
  #lsName = ls1[ix]#'compare_1_20_0.1_0.2.RData'
  #load(lsName)
  
  pctMi 
  
  
  dft.2 = df_regionBased
  
  n=1
  m=1
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
      row.names(dft.1)<- paste0(n, '_',dft.1$POS)
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
                         #K=4,
                         basis = basis_obj, 
                         vb_max_iter = 10,
                         vb_init_nstart = 1,
                         is_parallel = FALSE, 
                         is_verbose = FALSE)
  
  imputation_obj <- impute_met_state(obj = melissa_obj, test = dft.test)
  
  time.end.mel = Sys.time()
  
  cor(imputation_obj$act_obs, imputation_obj$pred_obs)
  
  dft.3 = list()
  dft.3.nrow=c()
  for(n in 1:N){
    #dft.3[[n]] = rbind(dft.test[[n]])
    dft.3[[n]]=as.data.frame(do.call(rbind, dft.test.c[[n]]))
    dft.3.nrow[n] = nrow(dft.3[[n]])
  }
  
  
  
  dft.4 = list()
  dft.4.nrow=c()
  dft.4.name = list()
  for(n in 1:N){
    #dft.3[[n]] = rbind(dft.test[[n]])
    dft.4[[n]]=as.data.frame(do.call(rbind, dft.test[[n]]))
    dft.4.nrow[n] = nrow(dft.4[[n]])
    dft.4.name[[n]] = row.names(as.data.frame(do.call(rbind, dft.test[[n]])))
  }
  
  
  
  
  
  imputation_obj$subj = rep(c(1:N), times = dft.3.nrow)
  imputation_obj$POS=names(imputation_obj$act_obs)
  imputation_obj$act_obs_c = as.data.frame(do.call(rbind, dft.3))$y
  
  
  imputation_obj.1 = as.data.frame(do.call(cbind,imputation_obj))
  
  imputation_obj.1$act_obs = as.numeric(as.character(imputation_obj.1$act_obs))
  imputation_obj.1$pred_obs = as.numeric(as.character(imputation_obj.1$pred_obs))
  imputation_obj.1$subj = as.numeric(as.character(imputation_obj.1$subj))
  imputation_obj.1$act_obs_c = as.numeric(as.character(imputation_obj.1$act_obs_c))
  imputation_obj.1$POS = as.character(imputation_obj.1$POS)
  
  
  
  #.1combine cor rmse of mel with previous output#####################
  
  
  
  
  
  cc=data.frame( cor_mel=NA)
  for(n in 1:N){
    dft1= imputation_obj.1[which(imputation_obj.1$subj==n),]
    #dft2 = dft1[which(dft1$train==0),]
    #cc[n,'cor_region'] = cor(dft2$y, dft2$y_hat)
    
    #cc[n,'cor_ws'] = cor(dft2$y, dft2$y_hat_3)
    cc[n,'cor_mel'] = cor(dft1$act_obs, dft1$pred_obs)
    
  }  
  cc
  cc_out = cbind(cc_out, cc)
  cc_out
  
  
  
  rr=data.frame( rmse_mel=NA)
  for(n in 1:N){
    dft1= imputation_obj.1[which(imputation_obj.1$subj==n),]
    #dft2 = dft1[which(dft1$train==0),]
    #cc[n,'cor_region'] = cor(dft2$y, dft2$y_hat)
    
    #cc[n,'cor_ws'] = cor(dft2$y, dft2$y_hat_3)
    rr[n,'rmse_mel'] = rmse(dft1$act_obs_c, dft1$pred_obs)
    
  }  
  rr
  rr_out = cbind(rr_out,rr)
  
  getwd()
  
  t1 = gregexpr('_', lsName)[[1]][1]
  t2 = gregexpr('RD', lsName)[[1]][1]
  f_O_c = paste0('compareCor6', substr(lsName,t1,t2-1 ), 'txt')
  f_O_r = paste0('compareRmse6', substr(lsName,t1,t2-1 ), 'txt')
  
  #write.table(cc_out, f_O_c, col.names=T, row.names=F, sep='\t', quote=F)
  
  #write.table(rr_out, f_O_r, col.names=T, row.names=F, sep='\t', quote=F)
  
  
  
  time.end.mel = Sys.time()
  time.run.mel = difftime(time.end.mel, time.start.mel, units='mins')
  time.run.mel
  
  
  
  
  t1 = apply(cc_out, 2, mean)
  t2 = apply(cc_out, 2, sd)
  t3 = c(1:ncol(cc_out))
  t4 = c(1:ncol(cc_out))
  for(i in 1:length(t3)){
    xx = t.test(cc_out[,1], cc_out[,i],mu=0, alternative = 'two.sided')
    t3[i] = xx$statistic
    t4[i] = xx$p.value
  }
  cc_out2 = rbind(cc_out, t1,t2,t3,t4)
  rownames(cc_out2)[1:4+nrow(cc_out)]<-c('mean','sd','t','p.value')
  cc_out2 = cc_out2[,c(1,2,4,3,5,6)]
  colnames(cc_out2) <- paste0('cor',c('MRM_region', 'MRM_stacked','average', 'KNN', 'RF', 'Melissa' ))
  
  
  
  t1 = apply(rr_out, 2, mean)
  t2 = apply(rr_out, 2, sd)
  t3 = c(1:ncol(rr_out))
  t4 = c(1:ncol(rr_out))
  for(i in 1:length(t3)){
    xx = t.test(rr_out[,1], rr_out[,i],mu=0, alternative = 'two.sided')
    t3[i] = xx$statistic
    t4[i] = xx$p.value
  }
  rr_out2 = rbind(rr_out, t1,t2,t3,t4)
  rownames(rr_out2)[1:4+nrow(rr_out)]<-c('mean','sd','t','p.value')
  rr_out2 = rr_out2[,c(1,2,4,3,5,6)]
  colnames(rr_out2) <- paste0('rmse_',c('MRM_region', 'MRM_stacked','average', 'KNN', 'RF', 'Melissa' ))
  
  t1 = gregexpr('_', lsName)[[1]][1]
  t2 = gregexpr('RD', lsName)[[1]][1]
  f_O_c2 = paste0('compareCor2', substr(lsName,t1,t2-1 ), 'txt')
  f_O_r2 = paste0('compareRmse2', substr(lsName,t1,t2-1 ), 'txt')
  
  write.table(cc_out2, f_O_c2, col.names=T, row.names=T, sep='\t', quote=F)
  write.table(rr_out2, f_O_r2, col.names=T, row.names=T, sep='\t', quote=F)
  
  
  
  
  #2 combine predictied value of methods #############################
  
  #ls1 = list.files('./', pattern='compare_condition2')
  #ix=1
  #print(ix)
  #lsName = ls1[ix]#'compare_1_20_0.1_0.2.RData'
  #lsName = 'compare_condition2_0.2_3.RData' 
  #load(lsName)
  
  
  aa=data.frame(auc_region=NA,  auc_ws=NA, auc_knn=NA)
  
  dfrf_allsubj=rbindlist(dfrf.l2)
  dfknn_allsubj = xxdf
  dfupd_allsubj = rbindlist(dfupd.l2)
  
  
  
  dfrf_allsubj$POS1 = paste0(dfrf_allsubj$subj, '_', dfrf_allsubj$POS)
  dfupd_allsubj$POS1 = paste0(dfupd_allsubj$subj, '_', dfupd_allsubj$POS)
  dfknn_allsubj$POS1 = paste0(dfknn_allsubj$subj, '_', dfknn_allsubj$POS)
  dfrf_allsubj.1 = dfrf_allsubj[,c('POS1', 'y_hat_rf')]
  dfupd_allsubj.1 = dfupd_allsubj[,c('POS1', 'y_hat_upd')]
  
  
  
  df_methods = merge(dfknn_allsubj, dfrf_allsubj.1, by.x = 'POS1', by.y='POS1', all=T )
  df_methods = merge(df_methods, dfupd_allsubj.1, by.x = 'POS1', by.y='POS1', all=T )
  
  
  
  df01 = df_methods[which(df_methods$train==0),]
  df02 = df01[order(df01$subj,df01$POS),]
  df02$POS1 = paste0(df02$subj, '_',df02$POS)
  
  
  pred_obs = matrix(NA, nrow=length(imputation_obj$pred_obs), ncol=2)
  pred_obs[,2] = imputation_obj$pred_obs
  pred_obs[,1] = names(imputation_obj$pred_obs)
  pred_obs=as.data.frame(pred_obs, stringsAsFactors =F)
  names(pred_obs)<-c('POS1', 'y_mel')
  pred_obs$y_mel<-as.numeric(pred_obs$y_mel)
  
  
  
  df03 = merge(df02,pred_obs, by.x='POS1', by.y='POS1', all.x=T)
  which(is.na(df03$y_mel))
  which(is.na(df03$y_hat))
  
  if(length(which(is.na(df03$y_mel)))>0){
    df04 = df03[-which(is.na(df03$y_mel)),]
  }else{
    df04= df03
  }
  
  
if(FALSE){
  xx=which(df04$subj==1)
  cor(df04[xx,'y'], df04[xx,'y_hat'])
  cor(df04[xx,'y'], df04[xx,'y_hat_knn'])
  cor(df04[xx,'y'], df04[xx,'y_hat_upd'], use='complete.obs')
  cor(df04[xx,'y'], df04[xx,'y_hat_rf'],use='complete.obs')
  cor(df04[,'y'], df04[,'y_mel'],use='complete.obs')
  
}
  

#3.compare auc ##############################################

##############################################,


library(data.table)
library(randomForest)
library(impute)
library(doParallel)
library(ROCR)
library(pROC)


#install.packages('impute')
#install.packages('pROC')

#if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

#BiocManager::install("impute")





#ls1 = list.files('./', pattern='compare_condition2')
#ix=1

  #print(ix)
  #lsName = ls1[ix]#'compare_1_20_0.1_0.2.RData'
  #lsName = 'compare_condition2_0.2_3.RData' 
  #load(lsName)
  
  
  
 
  
  aa=data.frame(auc_region=NA,  auc_ws=NA, auc_knn=NA)
  
  #.1 knn regional--no change from v03####
  for(n in 1:N){
    #dft1= xxdf[which(xxdf$subj==n),]
    dft2 = df04[which(df04$subj==n),]
    
    
    head(dft2)
    
    xx1=dft2[,'y_hat']
    xx2=ifelse(dft2[,'y']>0.5,1,0)
    pred <- prediction( xx1, xx2)
    perf <- performance( pred, "tpr", "fpr" )
    #plot(perf, col=colors()[33], lwd=1.5)
    aa[n,'auc_region'] = performance(pred,"auc")@y.values[[1]]
    
    
    xx1=dft2[,'y_hat_3']
    xx2=ifelse(dft2[,'y']>0.5,1,0)
    pred <- prediction( xx1, xx2)
    perf <- performance( pred, "tpr", "fpr" )
    #plot(perf, col=colors()[33], lwd=1.5)
    aa[n,'auc_ws'] = performance(pred,"auc")@y.values[[1]]
    
    xx1=dft2[,'y_hat_knn']
    xx2=ifelse(dft2[,'y']>0.5,1,0)
    pred <- prediction( xx1, xx2)
    perf <- performance( pred, "tpr", "fpr" )
    #plot(perf, col=colors()[33], lwd=1.5)
    aa[n,'auc_knn'] = performance(pred,"auc")@y.values[[1]]
    
  }  
  aa
  mean(aa$auc_region)
  mean(aa$auc_knn)
  var(aa$auc_region)
  var(aa$auc_knn)
  aa_knn=aa
  
  
  
  
  #.2 one upsteam and one downstream: non-regional version#####
  
  
  if(FALSE){ #------------------deleted from v03
    dfknn1 = rbindlist(dfws_out[[m]])
    dfupd.l = list()#upstream and downstream
    
    for(n in 1:N){
      dfupd.l[[n]] = list()
      for(m in 1:M){
        dfupd.l[[n]][[m]] = dfws_out[[m]][[n]]
        
      }
    }
    
    
    
    n=1
    dfupd.l2 = list()
    for(n in 1:N){
      print(n)
      dfupd = as.data.frame(rbindlist(dfupd.l[[n]]))
      
      dfupd = dfupd[order(dfupd$region, dfupd$x),]
      #dfupd$POS=rep(c(1:I),M)
      #dfupd$regionID = dfupd$region
      #names(dfknn1)
      #dfupd = dfupd[,c(62,1:4,63,5:61)]
      
      dfupd = dfupd[,c(1:11,62,63)]
      #dfupd = dfupd[order(dfupd$POS),]
      #names(dfupd)
      #table(dfupd$subj)
      #dim(dfupd)
      #length(unique(dfupd$POS))
      #dfupd = dfupd[-which(duplicated(dfupd$POS)),]
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
  }
  
  
  
  #aa=data.frame()
  for(n in 1:N){
    print(n)
    #dft1= dfupd.l2[[n]]
    dft2 = df04[which(df04$subj==n),]
    head(dft2)
    
    
    xx1=dft2[,'y_hat_upd']
    xx2=ifelse(dft2[,'y']>0.5,1,0)
    xx3 = which(is.na(xx1))
    if(length(xx3)>0){
      xx1=xx1[-xx3]
      xx2=xx2[-xx3]
    }
    pred <- prediction( xx1, xx2)
    perf <- performance( pred, "tpr", "fpr" )
    #plot(perf, col=colors()[33], lwd=1.5)
    aa[n,'auc_upd'] = performance(pred,"auc")@y.values[[1]]
  }  
  

  
  
  mean(aa$auc_region)
  mean(aa$auc_ws)
  mean(aa$auc_upd)
  
  
  
  
  
  
  #.3  RF ####
  
  
  
  #cl=makeCluster(6)
  # library(doParallel)
  
  if(FALSE){ #-------------------deleted
    dfknn1 = rbindlist(dfws_out[[m]])
    dfupd.l = list()#upstream and downstream
    
    for(n in 1:N){
      dfupd.l[[n]] = list()
      for(m in 1:M){
        dfupd.l[[n]][[m]] = dfws_out[[m]][[n]]
        
      }
    }
    
    
    
    
    n=1
    dfupd.l2 = list()
    
    registerDoParallel(cl)
    
    dfupd.l2<- foreach(n=1:N, .packages = c('data.table', 'randomForest') ) %dopar%{
      print(n)
      dfupd = as.data.frame(rbindlist(dfupd.l[[n]]))
      
      dfupd = dfupd[order(dfupd$region, dfupd$x),]
      #dfupd$POS=rep(c(1:I),M)
      #dfupd$regionID = dfupd$region
      #names(dfknn1)
      #dfupd = dfupd[,c(62,1:4,63,5:61)]
      
      dfupd = dfupd[,c(1:11,62,63)]
      dfupd = dfupd[order(dfupd$region,dfupd$POS),]
      #names(dfupd)
      #table(dfupd$subj)
      #dim(dfupd)
      #length(unique(dfupd$POS))
      #dfupd = dfupd[-which(duplicated(dfupd$POS)),]
      #dfupd$y_hat_upd = NA
      dfupd$y1 = NA
      dfupd$y2 = NA
      dfupd$d1 = NA
      dfupd$d2 = NA
      for(t in 1:nrow(dfupd)){
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
        
        dfupd$y1[t] = dfupd[x1,'y']
        dfupd$y2[t] = dfupd[x2,'y']
        dfupd$d1[t] = t-x1
        dfupd$d2[t] = x2-t
        #dfupd$y_hat_upd[t] =d1/(d1+d2)*y1+d2/(d1+d2)*y2 
        
        
        
        
      }
      
      model.1 = y~ y1+y2+d1+d2
      train.comp = dfupd[which(dfupd$train==1),c('region','POS','y', 'y1','y2', 'd1', 'd2')]
      train.comp = train.comp[which(complete.cases(train.comp)),]
      #test.comp = dfupd[which(dfupd$train==0),c('region','POS','y', 'y1','y2', 'd1', 'd2')]
      #which(!complete.cases(test.comp))
      fit.1 <- randomForest(model.1 , data=train.comp, importance=TRUE,ntree=1000)
      
      y_hat_rf<-predict(fit.1,newdata=dfupd)
      dfupd = cbind(dfupd, y_hat_rf)
      #dfupd.l2[[n]] = dfupd  
      dfupd
    }
    
    
    
  }
  
  
  
  for(n in 1:N){
    
    dft2 =  df04[which(df04$subj==n),]
    
    xx1=dft2[,'y_hat_rf']
    xx2=ifelse(dft2[,'y']>0.5,1,0)
    xx3 = which(is.na(xx1))
    if(length(xx3)>0){
      xx1=xx1[-xx3]
      xx2=xx2[-xx3]
    }
    
    pred <- prediction( xx1, xx2)
    perf <- performance( pred, "tpr", "fpr" )
    #plot(perf, col=colors()[33], lwd=1.5)
    aa[n,'auc_rf'] = performance(pred,"auc")@y.values[[1]]  
  }  
  
  mean(aa$auc_region)
  mean(aa$auc_rf)
  
 
  

  
  
  # .4 auc of mel ############
  
  for(n in 1:N){
    dft1= imputation_obj.1[which(imputation_obj.1$subj==n),]
    
    head(dft1)
    #colnames(xx)<-c("start", "sim", "Y", "pred1", "pred2")
    xx1=as.numeric(as.character(dft1[,'pred_obs']))
    xx2=as.numeric(as.character(dft1[, 'act_obs']))
    pred <- prediction( xx1, xx2)
    perf <- performance( pred, "tpr", "fpr" )
    plot(perf, col=colors()[33], lwd=1.5)
    
    aa[n,'auc_mel'] = performance(pred,"auc")@y.values[[1]]
    
  }  
  
  
    aa_out = aa
  
  
  
  #write.table(aa_out, f_O_a, col.names=T, row.names=F, sep='\t', quote=F)
  
  
  t1 = apply(aa_out, 2, mean)
  t2 = apply(aa_out, 2, sd)
  t3 = c(1:ncol(aa_out))
  t4 = c(1:ncol(aa_out))
  for(i in 1:length(t3)){
   xx = t.test(aa_out[,1], aa_out[,i],mu=0, alternative = 'two.sided')
   t3[i] = xx$statistic
   t4[i] = xx$p.value
  }
  aa_out2 = rbind(aa_out, t1,t2,t3,t4)
  rownames(aa_out2)[1:4+nrow(aa_out)]<-c('mean','sd','t','p.value')
  aa_out2 = aa_out2[,c(1,2,4,3,5,6)]
  colnames(aa_out2) <- paste0('auc_',c('MRM_region', 'MRM_stacked','average', 'KNN', 'RF', 'Melissa' ))
  
  
  t1 = gregexpr('_', lsName)[[1]][1]
  t2 = gregexpr('RD', lsName)[[1]][1]
  f_O_a = paste0('compareAUC', substr(lsName,t1,t2-1 ), 'txt')
  f_O_a2 = paste0('compareAUC2', substr(lsName,t1,t2-1 ), '.txt')
  write.table(aa_out2, f_O_a2, col.names=T, row.names=T, sep='\t', quote=F)
    

  




# 4. test and format ###############################################
library(pROC)

# .1. auc ####################
#roc.test(response = ifelse(df04$y>0.5,1,0), predictor1 = ifelse(df04$y_hat>0.5,1,0), predictor2 = ifelse(df04$V2>0.5,1,0), method=c('delong'))

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

cat(paste0(f_O_test.auc, '\n'), file = f_O_test.auc0, append=T)
write.table(test.auc, f_O_test.auc0, col.names = NA, row.names = T, sep='\t', quote=F, append=T)



# .2 Correlations  ###############  

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

cat(paste0(f_O_test.cor, '\n'), file = f_O_test.cor0, append=T)
write.table(test.cor, f_O_test.cor0, col.names = NA, row.names = T, sep='\t', quote=F, append=T)



#.3 rmse ############################################
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

cat(paste0(f_O_test.rmse, '\n'), file = f_O_test.rmse0, append=T)
write.table(test.rmse, f_O_test.rmse0, col.names = NA, row.names = T, sep='\t', quote=F, append=T)



###################

des_auc = psych::describe(aa_out)
des_cor = psych::describe(cc_out)
des_rmse = psych::describe(rr_out)
#print(des, digits = 6)

cat(paste0(pctMi, '\n'), file = './des_auc.txt', append=T)
write.table(des_auc, './des_auc.txt', row.names=T, col.names = NA, sep='\t', quote=F, append =T)



cat(paste0(pctMi, '\n'), file = './des_cor.txt', append=T)
write.table(des_cor, './des_cor.txt', row.names=T, col.names = NA, sep='\t', quote=F, append =T)



cat(paste0(pctMi, '\n'), file = './des_rmse.txt', append=T)
write.table(des_rmse, './des_rmse.txt', row.names=T, col.names = NA, sep='\t', quote=F, append =T)












