#' Stability estimator
#' @param x A feature selection matrix of dimension m and p
#' @param  K Number of training subsampling sets
#' @param  Verbose A logical and if TRUE,more information of the package gStab is returned.
#' @return Stability score over m subsampling experiments
#' @export
##  selection function
library(MASS);
selection <- function(x){
  j <- 1:length(x)
  ifelse(x[j]!= 0, 1, 0)
}
##
gStab <- function(x,K, Verbose=FALSE){
  ## m <- dim(x)[1];
  p <- dim(x)[2];
  ### selection
  M <- t(apply(x,1,selection));
  ## frequency vec function
  ### "K" represents the number of subsampling experiments
  Sn <- colSums(M);
  an <- rep(0,K);
  for(j in 1:K){
    In <- ifelse(j==Sn[1:length(Sn)],1,0);
    an[j] <- sum(In);
  }
  ## selected set sizes in "m" subsampling experiments
  r <- dim(M)[1];
  N <- rep(0,r);
  for(k in 1:r){
    N[k] <- length(M[k,][M[k,]!= 0]);
  }
  ## stability in selected feature set sizes
  h <- length(N);
  stab_size <- rep(0,h-1);
  for(i in 1:h-1){
    stab_size[i] <- ifelse(N[i+1] > N[i],N[i]/N[i+1],ifelse(N[i+1] < N[i],N[i+1]/N[i],1));
  }
  alp <-  c(N[1]/p,stab_size);
  ## proposed stability estimator equation
  m <- length(an);
  lambda <- rowSums(M);
  score <- rep(0,m);
  for(j in 1:m){
    score[j] <- (j^2)*an[j]*alp[j]/lambda[j];
  }
  stab <- (1/(K^2))*sum(score);
  if(Verbose){
    return(list(Sn=Sn,an=an,N=N,stab_size=stab_size,alpha=alp,stab=stab))
  }
  else{return(stab)};
}
