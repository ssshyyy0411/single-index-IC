gen.data<-function(n,gt,tau,alpha_true,beta_true,fun){
  
  X1 <- W1 <- runif(n,-1,1)
  X2 <- W2 <- runif(n,-1,1)
  X3 <- runif(n,-1,1)
  X4 <- runif(n,-1,1)
  
  X <- cbind(X1,X2,X3,X4)
  W <- cbind(W1,W2)
  
  A <- rbinom(n, 1, 0.5) #treatment assignment
  
  xb <- X%*%alpha_true + A*phi0(as.numeric(W%*%beta_true),fun=fun)   
  exb <- exp(xb)
  U <- runif(n,0,1)
  Lambt <- -log(U)*exb^(-1)  
  Ti <- lambda*Lambt^(1/r)
  

  d1 = d2 = d3 = rep(0,n)
  L = R = c()

  U1 <- runif(n,0,2*tau/5)
  U2 <- pmin((0.1 + U1 + rexp(n,1)*tau/2),tau)
  for(i in 1:n){
    if(Ti[i]< U1[i]) {
      d1[i] = 1           #d2[i] = d3[i] = 0
      L[i] = 0
      R[i] = U1[i]
    }
    if(Ti[i]> U2[i]) {
      d3[i] = 1           #d2[i] = d1[i] = 0
      L[i] = U2[i]
      R[i] = Inf
    }
    if(U1[i]< Ti[i] && Ti[i]<= U2[i]){
      d2[i] = 1           #d3[i] = d1[i] = 0
      L[i] = U1[i]
      R[i] = U2[i]
    }
  }

  L.prob = mean(d1)
  print(L.prob)
  I.prob = mean(d2)
  R.prob = mean(d3)
  print(R.prob)
  
  L[d1==1] <- NA
  R[d3==1] <- NA
  
  sdata <- data.frame(Li=L,Ri=R,d1=d1,d2=d2,d3=d3,X1=X1,X2=X2,X3=X3,X4=X4,Ti=Ti,W1=W1,W2=W2,A=A)
  
  return(sdata)
}