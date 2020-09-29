# 1. simulate methylation data by MRM ######################


#N = 20           #----number of subject
#noise0.var = 0.2 #----variance of noise level
#pctMi = 0.2      #----missing rate
#I=50             #----number of CpGs per region
#M=100            #----number of regions
#K=4              #----number of clusters of subjects
#N_mu = 10        #----number of RBF centers 
#mu = seq(-1,1,length.out=N_mu) #---- RBF centers 
#gamma = -10                    #---- RBF scale 
#noise0.var=0.2


#The following code create 'df', a 2-level list storing simulated data.
# Level 1: M lists represents M regions, each is a list of N dataframe 
# Level 2: List of N dataframe, each storing methyation data of subject n region m, 
# columns of N: y-methylation value, x-scaled position, POS-order of CpGs, h1 to h10-RBFs
# df_par: a 2-level list storing simulation parameters.
# level 1: list of M 
# level 2: list of 2 objects. C-cluster information of subjects, w-regression coeffecient of 4 clusters



set.seed(1234567)
df = list()
w =list()
df_par = list()
for(m in 1: M){
  
  for(k in 1:K){
    #w[[k]]= sample(1:N_mu, N_mu )
    w[[k]]=runif(N_mu, -1, 1)
    #w[[k]]=rbeta(N_mu, a,b,0)
    #w[[k]]=c(1,1,1,1,100)
  }
  C = sample(1:K,size = N, replace = T, prob =c(0.4, 0.3, 0.2, 0.1) )
  
  df[[m]] = list()
  
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
    df[[m]][[n]] =  data.frame(y = p, x=x,POS=c(1:I), h=h, subj=n, region=m, regionID=m)
  }
  
  df_par[[m]] = list(C=C, w=w)
  names(df)[m] = paste0('region_', m)
}


