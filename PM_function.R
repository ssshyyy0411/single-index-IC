#precision medicine

library(Rsolnp)
library(nloptr)
library(survival)
library(ICsurv)
library(splines)
library(splines2)
library(quadprog)
library(Matrix)


##################### necessary functions###############

#baseline hazard function ~ weibull
lambdat <- function(lambda,r,t)
{
  return(r*t^(r-1)/lambda^r)
}

#cumulative hazard function
Lambdat <- function(lambda,r,t)   
{
  return(t^r/lambda^r)
}


#link function
phi0<-function(u,fun)
{
  
  if(fun == 1){return(2*u + 0.1)}           #linear 
  if(fun == 2){return(exp(u) - 1.8)}        #exponential  
  if(fun == 3){return(3*u^2 - 2.5)}         #misspecified
}


###########################################

gen.initial<-function(sdata,order,n.lam,n.phi,t.seq){
  
  n<-nrow(sdata)
  d1<-sdata$d1
  d2<-sdata$d2
  d3<-sdata$d3
  Li<-sdata$Li
  Ri<-sdata$Ri
  A <-sdata$A
  
  X<-cbind(sdata$X1,sdata$X2,sdata$X3,sdata$X4)
  W<-cbind(sdata$W1,sdata$W2)
  AW<-A*W
  Xp<-cbind(sdata$X1,sdata$X2,sdata$X3,sdata$X4,AW,sdata$A)
  d<-dim(Xp)[2]
  
  #################### initialization ##################
  
  tol<-0.005
  e0<-rep(1,order+n.lam)
  b0<-rep(0,d)
  fit.cox<-fast.PH.ICsurv.EM(d1, d2, d3, Li, Ri, Xp, n.lam, order, e0, b0, tol, t.seq, equal = TRUE)
  alpha0<-fit.cox$b[1:ncol(X)]
  beta0<-fit.cox$b[(ncol(X)+1):(d-1)]
  c<-fit.cox$b[d]                     #coefficient of A
  eta0<-fit.cox$g
  

  x<-W%*%beta0
  x.max <- max(x) + 1e-05
  x.min <- min(x) - 1e-05
  y = W%*%beta0+c
  id <- seq(0, 1, length.out = (n.phi + 2))
  id <- id[-c(1, (n.phi + 2))]
  knots <- quantile(x, id)
  bs_x<-bs(x,knots=knots,degree=degree,Boundary.knots=c(x.min,x.max))
  bspline.fit = lm(y ~ bs_x)

  gamma0<-as.vector(coef(bspline.fit))    #initial of gamma 
 
  ini.par<-list(alpha0=alpha0, beta0=beta0, eta0=eta0, gamma0=gamma0)
  
  return(ini.par)
}

##########################################################

loglik <-function(beta,alpha,eta,gamma,sdata,Xp,W,Ml.Li,Ml.Ri){
  GRi<-Ml.Ri%*%eta
  GLi<-Ml.Li%*%eta
  u<-W%*%beta
  id <- seq(0, 1, length.out = (n.phi + 2))
  id <- id[-c(1, (n.phi + 2))]
  u.max <- max(u) + 1e-05
  u.min <- min(u) - 1e-05
  knots <- quantile(u, id)
  Bw<-bs(x=u, knots = knots, degree = degree, intercept =T,Boundary.knots=c(u.min,u.max))
  xb<-Xp%*%alpha+sdata$A*(Bw%*%gamma)
  l1<-(1-exp(-(GRi*exp(xb))))^sdata$d1
  l1[which(l1<=0)] <-10^-3
  l2<-((1-exp(-(GRi*exp(xb))))-(1-exp(-(GLi*exp(xb)))))^sdata$d2
  l2[which(l2<=0)] <-10^-3
  l3<-(1-(1-exp(-(GLi*exp(xb)))))^sdata$d3
  l3[which(l3<=0)] <-10^-3
  lpart1<-log(l1)
  lpart2<-log(l2)
  lpart3<-log(l3)
  
  return(sum(lpart1+lpart2+lpart3))
}

#########################################################

Q1<-function(bb,aa,gg,X,W,A,EWil,EYil,Ml.Li,Ml.Ri,d1,d2,d3){  
  N<-length(d1)
  L<-ncol(Ml.Li)
  r.mat<-function(x) matrix(x,ncol=N,nrow=L,byrow=TRUE)
  c.mat<-function(x) matrix(x,ncol=N,nrow=L)
  unew<-W%*%bb
  id <- seq(0, 1, length.out = (n.phi + 2))
  id <- id[-c(1, (n.phi + 2))]
  unew.max <- max(unew) + 1e-05
  unew.min <- min(unew) - 1e-05
  knots <- quantile(unew, id)
  Bwnew<-bs(x=unew, knots = knots, degree = degree, intercept =T,Boundary.knots=c(unew.min,unew.max))
  exb<-exp(X%*%aa+A*(Bwnew%*%gg))
  numer<-EWil+r.mat(1-d1)*EYil  
  denom<-(r.mat(1-d3)*t(Ml.Ri)+r.mat(d3)*t(Ml.Li))*r.mat(exb)    
  e1<- apply(numer,1,sum)/apply(denom,1,sum)
  return(e1)
}

#############################################
#iterative EM procedure

surv.EM.opt<-function(sdata, ini.par,degree,order,n.lam,n.phi,t.seq,max.iter, cov.rate, equal ,Cn) 
{
  n<-nrow(sdata)
  d1<-sdata$d1;d2<-sdata$d2; d3<-sdata$d3
  Li<-sdata$Li;Ri<-sdata$Ri
  A <-sdata$A
  
  X<-cbind(sdata$X1,sdata$X2,sdata$X3,sdata$X4)
  W<-cbind(sdata$W1,sdata$W2)
  
  G<-order+n.lam
  
  ##################################################
  #splines for cumulative baseline hazard function
  
  Li[ d1 == 1] <-  Ri[ d1 == 1]  
  Ri[ d3 == 1] <-  Li[ d3 == 1]  
  ti <- c(Li[d1 == 0], Ri[d3 == 0])
  
  if (equal == TRUE) {         
    ti.max <- max(ti) + .00001
    ti.min <- min(ti) - .00001
    knots <- seq(ti.min, ti.max, length.out = (n.lam + 2))
  }
  if (equal == FALSE) {       
    id <- seq(0, 1, length.out = (n.lam + 2))
    id <- id[-c(1, (n.lam + 2))]
    ti.max <- max(ti) + .00001
    ti.min <- min(ti) - .00001
    knots <- c(ti.min, quantile(ti, id), ti.max)
  }
  
  Ml.Li<-Ispline(Li,order=order,knots=knots) 
  Ml.Ri<-Ispline(Ri,order=order,knots=knots) 
  Mt <- t(Ispline(t.seq,order=order,knots=knots) )
  Ml.Li=t(Ml.Li)
  Ml.Ri=t(Ml.Ri)
  
  ###########################################
  alpha=ini.par$alpha0
  beta=ini.par$beta0
  eta=ini.par$eta0
  gamma=ini.par$gamma0
  
  dd=1
  iter<-0
  ll<-numeric()
  
  # loop
  while(dd>cov.rate & iter < max.iter) {
    
    # update iter
    iter<-iter + 1
    
    # link function B-splines
    u<-W%*%beta
    id <- seq(0, 1, length.out = (n.phi + 2))
    id <- id[-c(1, (n.phi + 2))]
    u.max <- max(u) + 1e-05
    u.min <- min(u) - 1e-05
    knots <- quantile(u, id)
    Bw<-bs(x=u, knots = knots, degree = degree, intercept =T,Boundary.knots=c(u.min,u.max))
    diff_Bw<-bSpline(x=u, knots = knots, degree = degree, intercept =T,deriv=1,Boundary.knots=c(u.min,u.max))  #derivative of B-splines
    
    #conditional expextations
    exb0<-exp(X%*%alpha+A*(Bw%*%gamma))
    Lamd.Li<-Ml.Li%*%eta
    Lamd.Ri<-Ml.Ri%*%eta
    dw<- 1-exp(-Lamd.Ri*exb0)
    dw[which(dw==0)] <- 1
    dy<-1-exp((Lamd.Li-Lamd.Ri)*exb0)
    dy[which(dy==0)] <- 1
    
    r.mat<-function(x) matrix(x,ncol=n,nrow=G,byrow=TRUE)
    c.mat<-function(x) matrix(x,ncol=n,nrow=G)
    
    EWi<-d1*exb0*Lamd.Ri/dw  
    EWil<-(r.mat(d1)*t(Ml.Ri))*c.mat(eta)*r.mat(exb0)/r.mat(dw)
    
    EYi<-d2*exb0*(Lamd.Ri-Lamd.Li)/dy
    EYil<-(r.mat(d2)*t(Ml.Ri-Ml.Li))*c.mat(eta)*r.mat(exb0)/r.mat(dy)
    

    #update beta
    fn1=function(beta01)
    {
      exb <-exp( X%*%alpha+A*(Bw%*%gamma+diff_Bw%*%gamma*W%*%(beta01-beta)))
      p1<-sum((EWil+r.mat(d2+d3)*EYil)*r.mat(X%*%alpha+A*(Bw%*%gamma+diff_Bw%*%gamma*W%*%(beta01-beta))))
      denom<-apply(((r.mat(1-d3)*t(Ml.Ri)+r.mat(d3)*t(Ml.Li))*r.mat(exb) ),1,sum)
      denom[denom==0]=0.01
      p2 = sum((EWil+r.mat(d2+d3)*EYil) * c.mat(log(denom)))
      if (is.nan(p2)) cat("NaN in p2\n")
      return(-(p1-p2))
    }
    
    eqn1=function(x){
      mod_beta= sqrt(sum(x^2))
      return(mod_beta)
    }
    # ineqfun = function(x) {
    #   return(x[1])
    # }
    control <- list(trace = FALSE, tol = 1e-3)
    powell=solnp(beta+0.0001, fun = fn1, 
                 eqfun = eqn1, eqB = c(1),
                 # ineqfun = ineqfun, ineqLB = c(1e-6), ineqUB = c(Inf),
                 control=control)
    
    beta.est<-powell$pars
    
    
    #update eta,alpha,gamma
    #s.t.-Cn<=gamm_1<=...<=gamma_jn<=Cn
    Theta<-c(gamma,alpha)
    
    eval_f<-function(Theta){
      X<-cbind(sdata$X1,sdata$X2,sdata$X3,sdata$X4)
      W<-cbind(sdata$W1,sdata$W2)
      
      gamma<-Theta[1:length(gamma)]
      alpha<-Theta[(length(gamma)+1):(length(gamma)+length(alpha))]
      e1<- Q1(bb=beta.est,aa=alpha,gg=gamma, X,W,A,EWil=EWil,EYil=EYil,Ml.Li=Ml.Li,Ml.Ri=Ml.Ri,d1=d1,d2=d2,d3=d3)
      e1[e1==0]=0.01
      l_eta<-log(e1)
      unew<-W%*%beta.est
      id <- seq(0, 1, length.out = (n.phi + 2))
      id <- id[-c(1, (n.phi + 2))]
      unew.max <- max(unew) + 1e-05
      unew.min <- min(unew) - 1e-05
      knots <- quantile(unew, id)
      Bwnew<-bs(x=unew, knots = knots, degree = degree, intercept =T,Boundary.knots=c(unew.min,unew.max))
      BX <- cbind(A * Bwnew, X)
      exbnew<-exp(BX %*% Theta)
      p1<-sum((EWil+r.mat(d2+d3)*EYil)*r.mat(BX %*% Theta))
      p2<-sum(EWil*c.mat(l_eta)+r.mat(d2+d3)*EYil*c.mat(l_eta))
      p3<-sum(c.mat(e1)*r.mat(exbnew)*(r.mat(d2+d1)*t(Ml.Ri)+r.mat(d3)*t(Ml.Li)))
      Q<-p1+p2-p3
      return(-Q)
    }
    
    #constraint functions
    eval_g_ineq<-function(Theta){
      gamma<-Theta[2:(length(gamma))]
      Jn=length(gamma)
      A2<--diag(Jn)
      A3=diag(Jn)*0
      diag(A3[-1,-Jn]) = 1
      A4<-rep(0,Jn)
      A4[length(A4)]=1
      A5<-rbind(A2+A3,A4)
      con<-A5%*%gamma
      con[1]=con[1]-Cn
      con[length(con)]=con[length(con)]-Cn
      constr<-c(con)
      return(constr)
    }
    

    res <- nloptr( x0=Theta,
                   eval_f=eval_f,
                   eval_g_ineq=eval_g_ineq,
                   opts=list( "algorithm" = "NLOPT_LN_COBYLA",xtol_rel = 1e-3)
    )
    gamma.est<-res$solution[1:length(gamma)]
    alpha.est<-res$solution[(length(gamma)+1):(length(gamma)+length(alpha))]
    
    eta.est<-Q1(bb=beta.est,aa=alpha.est,gg=gamma.est, X,W,A,EWil=EWil,EYil=EYil,Ml.Li=Ml.Li,Ml.Ri=Ml.Ri,d1=d1,d2=d2,d3=d3)
    eta.est[eta.est==0]=0.01
    
    ll<-c(ll,loglik(beta=beta.est,alpha=alpha.est,eta=eta.est,gamma=gamma.est,sdata,Xp=X,W=W,Ml.Li,Ml.Ri))
    if(iter>2) dd<-abs(ll[iter]-ll[iter-1])
    loglikli<-ll[iter]
    
    beta = beta.est
    eta = eta.est
    alpha = alpha.est
    gamma = gamma.est
    
    print(iter)
    print(dd)
    print(alpha)
    print(beta)
  }
  Lam<- Mt%*%eta
  AIC<-2*(length(beta)+length(gamma)+length(eta)+length(alpha))-2*loglikli
  BIC<-(length(beta)+length(gamma)+length(eta)+length(alpha))*log(n)-2*loglikli
  
  u.seq=seq(-1.2,1.2,length.out =100)
  id <- seq(0, 1, length.out = (n.phi + 2))
  id <- id[-c(1, (n.phi + 2))]
  u.max <- max(u.seq) + 1e-05
  u.min <- min(u.seq) - 1e-05
  knots <- quantile(u.seq, id)
  Bwf<-bs(x=u.seq, knots = knots, degree = degree, intercept =T,Boundary.knots=c(u.min,u.max))
  Phi<-Bwf%*%gamma
  
  
  #average treatment effect beta_av
  uf<-W%*%beta
  uf.max <- max(uf) + 1e-05
  uf.min <- min(uf) - 1e-05
  write.table(t(c( max(uf), min(uf))),paste(filename,"uf.txt",sep=""),row.names=FALSE,col.names=FALSE,append=TRUE)
  knots <- quantile(uf, id)
  Bw_hat<-bs(x=uf, knots = knots, degree = degree, intercept =T,Boundary.knots=c(uf.min,uf.max))
  beta_av=mean(Bw_hat%*%gamma)
  
  # return
  ans<-list(a=alpha, b=beta, e=eta, g=gamma, Lam=Lam, Phi=Phi,beta_av=beta_av,AIC=AIC,BIC=BIC,iter=iter)
  
  return(ans)
}


#############################################
