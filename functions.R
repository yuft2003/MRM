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