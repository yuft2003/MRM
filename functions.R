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

#scale01 <- function(x){(x-0.5*((max(x)-min(x))))/(max(x)-min(x))}

scale00 <- function(x){(x-min(x))/(max(x)-min(x))}

scale01 <- function(x){2*(x-(0.5*((max(x)-min(x)))+min(x)))/(max(x)-min(x)+2)}

runs.test.1 <-  function(x, exact = FALSE, alternative = c("two.sided", "less", "greater"))  {
  DNAME <- deparse(substitute(x))
  METHOD <- ifelse(exact,"Exact runs test", "Approximate runs rest")
  alternative <- match.arg(alternative)
  x <- x[is.finite(x)]
  N <- as.integer(length(x))
  if (N < 1L) 
    stop("not enough 'x' data")
  if (!is.numeric(x)) 
    stop("'x' must be numeric")
  x.s <- if (length(unique(x)) == 2) (x == unique(x)[1])*1 
  else (x > mean(as.matrix(x)))*1
  m <- sum(1 - x.s)
  n <- N - m
  R <- 1
  for (i in 1:(N - 1))
  {
    if (x.s[i] != x.s[i+1])  R <- R + 1
  }
  STATISTIC <- setNames(R, "Runs")
  P.even <- function(m,n,r) 2*choose(m-1,r-1)*choose(n-1,r-1)/choose(m + n,n)
  P.odd <- function(m,n,r)
    (choose(m-1,r-1)*choose(n-1,r) + choose(n-1,r-1)*choose(m-1,r))/choose(m + n,n) 
  if (exact)
  {
    if (any(is.na(P.even(m,n,1:floor(R/2)))) || any(is.na(P.odd(m,n,1:floor(R/2)))))
      stop("can't calculate exact p-value; please use approximate method")
    if (R%%2 == 0)
    {
      p.val <- sum(P.even(m,n,1:(R/2))) + sum(P.odd(m,n,1:(R/2 - 1)))
      p.val1 <- 1 - p.val + P.even(m,n,R/2)
    }
    else 
    {
      p.val <- sum(P.even(m,n,1:floor(R/2))) + sum(P.odd(m,n,1:floor(R/2))) 
      p.val1 <- 1 - p.val + P.odd(m,n,floor(R/2))
    }        
  }
  else
    Z <- (R - 2*m*n/N - 1)/sqrt(2*m*n*(2*m*n - N)/(N^2*(N - 1)))
  P.VAL <- switch(alternative, two.sided = ifelse(exact, 2*min(p.val,1 - p.val), 
                                                  2*min(pnorm(Z),1- pnorm(Z))), less = ifelse(exact, p.val, pnorm(Z)),
                  greater = ifelse(exact, p.val1, 1 - pnorm(Z)))
  RVAL <- list(data.name = DNAME, method = METHOD, alternative = alternative,
               statistic = STATISTIC, p.value = P.VAL)
  class(RVAL) <- "htest"
  return(RVAL)
}

