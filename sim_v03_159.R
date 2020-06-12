#!/bin/bash
#SBATCH --job-name sim08
#SBATCH --qos normal
#SBATCH -N 1
#SBATCH -n 30
#SBATCH -o aa_log08.txt
#SBATCH -e aa_err08.txt
#SBATCH --time=24:00:00
#SBATCH --array=1





R  --quiet --no-save > /dev/null << EOF

setwd('/home/fyu/me_impute/realData_PBM')

simPar = read.table('./sim_v03_simNo.txt', header=T, sep='\t', stringsAsFactors=F)

simNo = 8
N = simPar$N[simNo]
noise0.var = simPar$noise[simNo]
pctMi = simPar$missingRate[simNo]
k.region.max = simPar$k.region.max[simNo]
k.region = 2:k.region.max #N
k.subj = 2:8 #M=100

if(FALSE){
  for(simNo in simPar$simNo){
    print(simNo)
    N = simPar$N[simNo]
    noise0.var = simPar$noise[simNo]
    pctMi = simPar$missingRate[simNo]	
    
    
    library(flexmix)
    library(data.table)
    
    I=50
    M=100
    
    
    K=4
    L=2
    
    N_mu = 10 
    w =list()
    a = 1
    b =10
    #beta_var = a*b/(a+b)^2/(a+b+1)
    mu = seq(-1,1,length.out=N_mu)
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
    
    #df[["region_1"]]["par"]
    
    
    f_O_1 = paste0('./sim_v03/data_', simNo,'_', N,'_',noise0.var,'.RData' )
    save.image( f_O_1, version=2)
  }
  
  
}


f_O_1 = paste0('./sim_v03/data_', simNo,'_', N,'_',noise0.var,'.RData' )
load(f_O_1)



###########################################

## ensemble mixReg + mixReg (v02)

## v02 ##########################################
#install.packages('glmnet')
source('./source/glmnet/R/glmnet.R')
library(glmnet)
library(flexmix)
library(data.table)
library(nnls)
library(doParallel)
createHmatrix <- function(x,   gamma,mu) {
  H = matrix(,ncol = length(mu), nrow = length(x) )
  for(i in 1: length(x)){
    H[i,] = exp(gamma*(x[i]-mu)^2) 
  }
  colnames(H)<-paste0("H", 1:length(mu))
  return(H)
}

rmse <- function(x,y){
  return(sqrt(sum((x-y)^2)/length(x)))
}

FLXMRglmnet <-
  function(formula = .~., family = c("gaussian", "binomial", "poisson"), adaptive = TRUE, select = TRUE, offset = NULL, ...) {
    family <- match.arg(family)
    z <- FLXMRglm(formula = formula, family = family)
    z@preproc.x <- function(x) {
      if (!isTRUE(all.equal(x[, 1], rep(1, nrow(x)), check.attributes = FALSE)))
        stop("The model needs to include an intercept in the first column.")
      x
    }
    z@fit <- function(x, y, w) {
      if (all(!select)) {
        coef <- if (family == "gaussian")
          lm.wfit(x, y, w = w)$coef
        else if (family == "binomial")
          glm.fit(x, y, family = binomial(), weights = w)$coef
        else if (family == "poisson")
          glm.fit(x, y, family=poisson(), weights = w)$coef
      } else {
        if (adaptive) {
          coef <- if (family == "gaussian")
            lm.wfit(x, y, w = w)$coef[-1]
          else if(family == "binomial")
            glm.fit(x, y, family = binomial(), weights = w)$coef[-1]
          else if (family == "poisson")
            glm.fit(x, y, family = poisson(), weights = w)$coef[-1]
          penalty <- mean(w) / abs(coef)
        } else
          penalty <- rep(1, ncol(x) - 1)
        if (any(!select)){
          select <- which(!select)
          penalty[select] <- 0
        }      
        m <-  glmnet::cv.glmnet(x[, -1, drop = FALSE], y, family = family, weights = w, ...)
        coef <- as.vector(coef(m, s = "lambda.min"))
      }
      df <- sum(coef != 0)
      sigma <- if (family == "gaussian") sqrt(sum(w * (y - x %*% coef)^2/mean(w))/(nrow(x) - df)) else NULL
      z@defineComponent(
        list(coef = coef, sigma = sigma, df = df + ifelse(family == "gaussian", 1, 0)))
    }
    z
  }

predFromList=function(y_hat_list, y_c){
  if(length(y_c) != length(y_hat_l[[1]]) ){
    print('error')
  }else{
    
  }
  for( i_t in 1:length(y_c)){
    y_hat[i_t] = y_hat_l[y_c[i_t]][i_t]
  }
  return(y_hat)
}

scale01 <- function(x){(x-0.5*((max(x)-min(x))))/(max(x)-min(x))}



### region based mixture regression
### The last subj has n_missing missing values
#pctMi=0.2
Imin=50
subjMi = c(1:N)

#df_regionBased = list()
cl <- makeCluster(30)
registerDoParallel(cl)

temp_df<- foreach(m=1:M, .packages = c('data.table', 'flexmix') ) %dopar%{
  print(m)
  
  
  subjRef = setdiff(c(1:N),subjMi)
  
  df2 =rbindlist(df[[m]]$data)
  #C = df$region_1$par$C
  #w = df$region_1$par$w
  
  #df2$subj = rep(c(1:N), each=I)
  #df2$C = rep(C, each=I)
  df2 = as.data.frame(df2)
  
  
  #df2 = df2[-which(is.na(df2$y)),]
  
  #indMi = which(df2$subj %in% subjMi)
  df2$train = 1
  set.seed(12345)
  
  for(i_temp in 1: N){
    I=nrow(df[[m]]$data[[i_temp]])
    x_test = sample(I, size = round(pctMi*I), replace = F)
    print(x_test)
    df2[which(df2$subj== subjMi[i_temp]),][x_test,'train'] = 0 
  }
  #which(df2$train==0)
  table(df2$train)
  
  H = createHmatrix(x=df2$x,gamma=-10, mu=seq(-1,1,length.out = Imin))
  
  df2_train = data.frame(df2, H)[which(df2$train==1),]
  df2_test = data.frame(df2, H)[which(df2$train==0),]
  #m4 <- flexmix(y ~ h.1+h.2+h.3+h.4+h.5 | subj, data = df2, k = 4)
  #m4 <- flexmix(fmla, data=df2, k=4)
  xnam = colnames(H)
  fmla <- as.formula(paste("y ~ ", paste(xnam, collapse= "+"),'|subj'))
  #fmla <- as.formula(paste("y ~ ", paste(xnam, collapse= "+")))
  
  #fo <- sample(rep(seq(10), length = nrow(df2)))
  #x01=cbind(df2[,xnam])
  #y01= df2$y
  
  #fo <- sample(rep(seq(10), length = nrow(df2)))
  #data01=data.frame(y01,x01)
  Model=FLXMRglmnet(lambda=10^seq(from=-5, to=5, by=1/3),alpha=0)
  
  if(FALSE){
    k_cluster=4
    z <- sample(1:k_cluster, nrow(df2_train), replace = TRUE)
    m2 <- flexmix(fmla, data = df2_train, k = k_cluster,model = Model, cluster=z, control = list(iter.max = 100))
    modelFit = m2
  }  
  
  
  if(TRUE){
    m2 <- initFlexmix(fmla, data = df2_train, k = k.region, 
                      model = Model, 
                      control = list(iter.max = 100), nrep=1)
    #plot(m2)
    
    modelFit = getModel(m2, which = "ICL") 
    #modelFit =m2@models$`3`
    #parameters(modelFit)
    
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
  cor(df2_test$y, df2_test$y_hat)
  #plot(df2_test$y, df2_test$y_hat)
  rmse(df2_test$y, df2_test$y_hat)
  
  df2 = rbind(df2_train, df2_test)
  
  df2 = df2[order(df2$subj, df2$x),]
  df2$res = df2$y_hat - df2$y
  
  rmse(df2$y,df2$y_hat)
  
  
  
  #summary(modelFit)
  
  #parameters(modelFit)
  #m4@cluster
  #check cluster accuracy
  
  #  table(df2$C, df2$C.pred)
  
  # names(df2)
  
  
  
  df2
  
}
df_regionBased = temp_df
#stopCluster(cl)
#temp=df_regionBased[[2]]
#which(temp$train==0)


#xx = rbindlist(df_regionBased)
#xx1 = xx[which(xx$train==0),]
#dim(xx1)

#cor(xx1$y, xx1$y_hat)


# subject based mixture regression
#######################################

#method 1 checked: using result from first step

df3 = list()
for(n in subjMi){
  for(m in 1:M){
    df3[[m]]=df_regionBased[[m]][which(df_regionBased[[m]]$subj==n),]
    #df3[[m]]$region = m
    df3[[m]] = df3[[m]][,c('POS','y','x','subj','region','regionID','train','y_hat','C.pred','res')]
  }
  
}







#cl <- makeCluster(20)
registerDoParallel(cl)
temp_df = list()

#df_subjBased[[n]]=df4

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
                      control = list(iter.max = 100), nrep=1)
    
    #plot(m2)
    m3 = getModel(m2, which = "ICL")
    parameters(m3)
    
    
    modelFit = m3
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

df_subjBased = temp_df




if(FALSE){
  xx=df4
  xx=df_subjBased[[1]]
  names(xx)
  xx1 = xx[which(xx$train==0),]
  cor(xx$y_hat_2, xx$y)
  cor(xx$y_hat, xx$y)
  
  cor(xx1$y_hat_2, xx1$y)
  cor(xx1$y_hat, xx1$y)
  
  plot(xx$x[which(xx$region==4)],xx$y[which(xx$region==4)])
  
  
}

#train the weight using non-negative least square

## method 1. region specific weights
## sample from subjMi






#dfws = df4
iter=100
#dfws_out=list()

#dopar

registerDoParallel(cl)
#temp_df = list()

#df_subjBased[[n]]=df4

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

dfws_out = temp_df




dft = list()
for(n in 1:N){
  dft[[n]] = list()
  for(m in 1:M){
    dft[[n]][[m]]=dfws_out[[m]][[n]]
    
  }
  
}

if(FALSE){
  cc=data.frame(cor_region=NA, cor_subj=NA, cor_ws=NA)
  for(n in 1:N){
    dft1 = rbindlist(dft[[n]])
    dft2 = dft1[which(dft1$train==0),]
    cc[n,1] = cor(dft2$y, dft2$y_hat)
    cc[n,2] = cor(dft2$y, dft2$y_hat_2)
    cc[n,3] = cor(dft2$y, dft2$y_hat_3)
    
  }  
  
  cc
  
  mean(cc[,1])
  mean(cc[,2])
  mean(cc[,3])
  names(dft1)
  
  getwd()
}


f_O_2 = paste0('./sim_v03/compute_', simNo,'_', N,'_',noise0.var,'_', pctMi,'.RData' )

save.image(f_O_2, version=2)

#compare ##############################################

##############################################

library(data.table)
library(randomForest)
library(impute)

#install.packages('impute')

#if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

#BiocManager::install("impute")


#knn regional####

knn_out = list()

M
m=1

for(m in 1:M){
  print(m)
  
  
  dfknn1 = rbindlist(dfws_out[[m]])
  dfknn1 = dfknn1[order(dfknn1$subj, dfknn1$x),]
  #dfknn1$POS=rep(c(1:I),N)
  #dfknn1$regionID = dfknn1$region
  #names(dfknn1)
  #dfknn1 = dfknn1[,c(62,1:4,63,5:61)]
  
  dim(dfknn1)
  
  which(dfknn1$train==0)
  
  
  #names(dfknn1)
  #dfknn2 = dfknn1[which(dfknn1$region==m),]
  
  
  #table(dfknn2$subj)
  
  
  xx=reshape(dfknn1, direction = 'wide', idvar='POS',timevar='subj', v.names='y', drop=names(dfknn1)[c(7:63)]  )
  
  xx_tr=reshape(dfknn1, direction = 'wide', idvar='POS',timevar='subj', v.names='train', drop=names(dfknn1)[c(2,8:63)]  )
  
  xx1=xx
  x.missing=matrix(FALSE, nrow=nrow(xx), ncol=ncol(xx))
  x.missing[,5:(5+N-1)] = xx_tr[,5:(5+N-1)] == 0
  xx1[x.missing]=NA
  
  x.missing2 = x.missing[,5:(5+N-1)]
  
  
  xx_imp= impute.knn(data=as.matrix(xx1[,5:(5+N-1)]),  rowmax=1, colmax=1)
  
  
  
  xx_imp1 = cbind(xx[,1],xx_imp$data)
  
  xx_imp2 = reshape(xx_imp1,direction='long', varying=list(names(xx_imp1)[2:(2+N-1)]), timevar='subj')[,1:3]
  
  names(xx_imp2)[3]='y_hat_knn'
  xx_imp2$subj_pos = paste0(xx_imp2$subj,'_', xx_imp2$POS)
  
  dfknn1$subj_pos = paste0(dfknn1$subj,'_', dfknn1$POS)
  
  dfknn2 = merge(dfknn1, xx_imp2, by='subj_pos', all.x=T, all.y=F)
  
  dfknn2=dfknn2[,c(2:12,63,64,67)]
  
  #xx2=as.matrix(xx[,5:23])
  #xx3=xx2[x.missing2]
  #knn_out[[m]] = data.frame(xx3, xx_imp1,region=m)
  names(dfknn2)[c(1,4)] = c('POS','subj')
  
  knn_out[[m]] = dfknn2
}

xxdf=rbindlist(knn_out)
#names(xxdf)
#xxdf1 = xxdf[which(xxdf$subj==1),]

cc=data.frame(cor_region=NA,  cor_ws=NA, cor_knn=NA)
for(n in 1:N){
  dft1= xxdf[which(xxdf$subj==n),]
  dft2 = dft1[which(dft1$train==0),]
  cc[n,'cor_region'] = cor(dft2$y, dft2$y_hat)
  
  cc[n,'cor_ws'] = cor(dft2$y, dft2$y_hat_3)
  cc[n,'cor_knn'] = cor(dft2$y, dft2$y_hat_knn)
  
}  
cc
mean(cc$cor_region)
mean(cc$cor_knn)
var(cc$cor_region)
var(cc$cor_knn)
cc_knn=cc


rr=data.frame(rmse_region=NA,  rmse_ws=NA, rmse_knn=NA)
for(n in 1:N){
  dft1= xxdf[which(xxdf$subj==n),]
  dft2 = dft1[which(dft1$train==0),]
  rr[n,'rmse_region'] = rmse(dft2$y, dft2$y_hat)
  
  rr[n,'rmse_ws'] = rmse(dft2$y, dft2$y_hat_3)
  rr[n,'rmse_knn'] = rmse(dft2$y, dft2$y_hat_knn)
  
}  
rr_knn = rr




library(e1071)

if(FALSE){
  cc.region=data.frame(cor_region=NA,  cor_ws=NA, cor_knn=NA, skew=NA,var=NA)
  for(m in 1:M){
    dft1= xxdf[which(xxdf$region==m),]
    dft2 = dft1[which(dft1$train==0),]
    cc.region[m,'cor_region'] = cor(dft2$y, dft2$y_hat)
    
    cc.region[m,'cor_ws'] = cor(dft2$y, dft2$y_hat_3)
    cc.region[m,'cor_knn'] = cor(dft2$y, dft2$y_hat_knn)
    cc.region[m,'skew'] = skewness(dft1$y)
    cc.region[m,'var'] = var(dft1$y)
    cc.region[m,'range'] = max(dft1$y)-min(dft1$y)
    
  }  
  cc.region
  
  mean(cc.region$cor_region)
  mean(cc.region$cor_knn)
  
  cc.region$dif = cc.region$cor_region -cc.region$cor_knn
  cor(cc.region$skew+10,cc.region$dif,use='complete.obs')
  cor(cc.region$var,cc.region$dif,use='complete.obs')
  
  plot(cc.region$var,cc.region$dif,use='complete.obs')
  
  
  length(which(cc.region$cor_region>cc.region$cor_knn))
  length(which(cc.region$cor_region<cc.region$cor_knn))
  
  cc.region[which(cc.region$cor_region>cc.region$cor_knn),]
  cc.region[which(cc.region$cor_region<cc.region$cor_knn),]
  
  cc.region.1 = cc.region[which(cc.region$cor_knn>0.7),]
  
  cc.region.1 = cc.region[which(cc.region$cor_knn>0.7& cc.region$cor_region>0.7),]
  
  cc.region.1 = cc.region[which(cc.region$skew>-1),]
  
  which(cc.region$cor_knn>0.7)
  
  #cc.region.1 = cc.region[which(cc.region$cor_region>0.7),]
  length(which(cc.region$cor_region>0.7))
  length(which(cc.region.1$cor_region>cc.region.1$cor_knn))
  mean(cc.region.1$cor_region)
  mean(cc.region.1$cor_knn)
  
  
  cor(xx$xx3,xx$xx_imp1,use='complete.obs')
  
}



#######KNN non- regional ######################################
if(FALSE){
  xx1.l = list()
  xx.l = list()
  dfknn1.l = list()
  
  for(m in 1:M){
    print(m)
    
    
    dfknn1 = rbindlist(dfws_out[[m]])
    #names(dfknn1)
    #dfknn2 = dfknn1[which(dfknn1$region==m),]
    
    
    #table(dfknn2$subj)
    
    
    xx=reshape(dfknn1, direction = 'wide', idvar='x',timevar='subj', v.names='y', drop=names(dfknn1)[c(7:63)]  )
    
    xx_tr=reshape(dfknn1, direction = 'wide', idvar='x',timevar='subj', v.names='train', drop=names(dfknn1)[c(2,8:63)]  )
    
    xx1=xx
    x.missing=matrix(FALSE, nrow=nrow(xx), ncol=ncol(xx))
    x.missing[,5:23] = xx_tr[,5:23] == 0
    xx1[x.missing]=NA
    
    x.missing2 = x.missing[,5:23]
    
    xx1.l[[m]]=xx1
    xx.l[[m]] = xx
    dfknn1.l[[m]] = dfknn1
  }
  
  xx1.rb = as.data.frame(rbindlist(xx1.l)  )
  xx1.rb = xx1.rb[order(xx1.rb$POS),]
  #xx1.rb = xx1.rb[-which(duplicated(xx1.rb$POS)),]
  xx.rb =  as.data.frame(rbindlist(xx.l)  )
  dfknn1 = as.data.frame(rbindlist(dfknn1.l)  )
  
  
  
  xx_imp= impute.knn(data=as.matrix(xx1.rb[,5:23]), k=10, rowmax=1, colmax=1)
  
  
  
  xx_imp1 = as.data.frame(cbind(xx1.rb[,1],xx_imp$data))
  
  head(xx_imp1)
  names(xx_imp1)[1]='POS'  
  typeof(xx_imp1)  
  dim(xx_imp1)
  
  xx_imp2 = reshape(xx_imp1,direction='long', varying=list(names(xx_imp1)[2:20]), timevar='subj')[,1:3]
  
  names(xx_imp2)[3]='y_hat_knn'
  xx_imp2$subj_pos = paste0(xx_imp2$subj,'_', xx_imp2$POS)
  
  dfknn1$subj_pos = paste0(dfknn1$subj,'_', dfknn1$POS)
  #dfknn1 = dfknn1[-which(duplicated(dfknn1$subj_pos)),]
  
  dfknn2 = merge(dfknn1, xx_imp2, by='subj_pos', all.x=T, all.y=F)
  
  dfknn2=dfknn2[,c(2:12,63,64,67)]
  
  #xx2=as.matrix(xx[,5:23])
  #xx3=xx2[x.missing2]
  #knn_out[[m]] = data.frame(xx3, xx_imp1,region=m)
  names(dfknn2)[c(1,4)] = c('POS','subj')
  
  
  
  
  xxdf=dfknn2
  
  names(xxdf)
  xxdf1 = xxdf[which(xxdf$subj==1),]
  
  cc=data.frame(cor_region=NA,  cor_ws=NA, cor_knn=NA)
  for(n in 1:N){
    dft1= xxdf[which(xxdf$subj==n),]
    dft2 = dft1[which(dft1$train==0),]
    cc[n,'cor_region'] = cor(dft2$y, dft2$y_hat)
    
    cc[n,'cor_ws'] = cor(dft2$y, dft2$y_hat_3)
    cc[n,'cor_knn'] = cor(dft2$y, dft2$y_hat_knn)
    
  }  
  
  cc
  
  
  
  
  
}







#one upsteam and one downstream: non-regional version#####




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


cc=data.frame(cor_region=NA,  cor_ws=NA, cor_upd=NA)
for(n in 1:N){
  dft1= dfupd.l2[[n]]
  dft2 = dft1[which(dft1$train==0),]
  cc[n,'cor_region'] = cor(dft2$y, dft2$y_hat)
  
  cc[n,'cor_ws'] = cor(dft2$y, dft2$y_hat_3)
  cc[n,'cor_upd'] = cor(dft2$y, dft2$y_hat_upd, use='complete.obs')
  
}  

cc
cc_upd = cc

mean(cc$cor_region)
mean(cc$cor_ws)
mean(cc$cor_upd)

rr=data.frame(rmse_region=NA,  rmse_ws=NA, rmse_upd=NA)
for(n in 1:N){
  dft1= dfupd.l2[[n]]
  dft2 = dft1[which(dft1$train==0),]
  rr[n,'rmse_region'] = rmse(dft2$y, dft2$y_hat)
  
  rr[n,'rmse_ws'] = rmse(dft2$y, dft2$y_hat_3)
  
  if(length(which(is.na(dft2$y)|(is.na(dft2$y_hat_upd))))>0){
    rr[n,'rmse_upd'] = rmse(dft2$y[-which(is.na(dft2$y)|(is.na(dft2$y_hat_upd)))], dft2$y_hat_upd[-which(is.na(dft2$y)|(is.na(dft2$y_hat_upd)))])
  }else{
    rr[n,'rmse_upd'] = rmse(dft2$y, dft2$y_hat_upd)
  }
  
  
}  
rr_upd = rr


##################################################

#one upsteam and one downstream: RF (unfinished)####

######################################################


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



cc=data.frame(cor_region=NA,  cor_ws=NA, cor_rf=NA)
for(n in 1:N){
  dft1= dfupd.l2[[n]]
  dft2 = dft1[which(dft1$train==0),]
  cc[n,'cor_region'] = cor(dft2$y, dft2$y_hat)
  
  cc[n,'cor_ws'] = cor(dft2$y, dft2$y_hat_3)
  cc[n,'cor_rf'] = cor(dft2$y, dft2$y_hat_rf, use='complete.obs')
  
}  

mean(cc$cor_region)
mean(cc$cor_rf)

cc_rf = cc


rr=data.frame(rmse_region=NA,  rmse_ws=NA, rmse_rf=NA)
for(n in 1:N){
  dft1= dfupd.l2[[n]]
  dft2 = dft1[which(dft1$train==0),]
  rr[n,'rmse_region'] = rmse(dft2$y, dft2$y_hat)
  
  rr[n,'rmse_ws'] = rmse(dft2$y, dft2$y_hat_3)
  
  if(length(which(is.na(dft2$y)|(is.na(dft2$y_hat_rf))))>0){
    rr[n,'rmse_rf'] = rmse(dft2$y[-which(is.na(dft2$y)|(is.na(dft2$y_hat_rf)))], dft2$y_hat_rf[-which(is.na(dft2$y)|(is.na(dft2$y_hat_rf)))])
  }else{
    rr[n,'rmse_rf'] = rmse(dft2$y, dft2$y_hat_rf)
  }
  
  
}  

mean(rr$rmse_region)
mean(rr$rmse_rf)

rr_rf = rr

stopCluster(cl)

cc_out = data.frame(cc_knn, cor_upd=cc_upd$cor_upd, cor_rf=cc_rf$cor_rf)
rr_out = data.frame(rr_knn, rmse_upd=rr_upd$rmse_upd, rmse_rf=rr_rf$rmse_rf)


f_O_3 = paste0('./sim_v03/compare_', simNo,'_', N,'_',noise0.var,'_', pctMi,'.RData' )
save.image(f_O_3, version=2)

f_O_cor = paste0('./sim_v03/compareCor_', simNo,'_', N,'_',noise0.var,'_', pctMi,'.txt' )

f_O_rmse = paste0('./sim_v03/compareRmse_', simNo,'_', N,'_',noise0.var,'_', pctMi,'.txt' )

write.table(cc_out, f_O_cor, col.names=T, row.names=F, quote=FALSE, sep='\t')

write.table(rr_out, f_O_rmse, col.names=T, row.names=F, quote=FALSE, sep='\t')

EOF



