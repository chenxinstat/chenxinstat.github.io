########################################################################
# Code to do auto-adaptive M-estimation
#
#  must-have inputs: y = the data
#
#  optional inputs: mr = the mesh ratio (default value 10)
#                   alpha = type I error rate for confidence interval (defualt value 0.05)
#                   penvals = vector of candidate values of penalty (o/w will get one automatically)
#                   pen = the penalty constant (o/w will use cross-validation to find one)
#                   mu0 = the starting point of mean (o/w will use median of the data)
#
#  returns: bhat = the estiamted spline density coefficients
#           confidence.intervals = the confidence interval of mean of the data
#           fhat = estimated error density function
#           knots = planced knots for spline density estiamtion
#           muhat = estiamted mean for the data
#           pen = selected penalty constant by cross-validation
#           penvals = automatically generated candidate values of penalty
#           risk = estimated risks for cross-validation
#
########################################################################
## need library(quadprog)
########################################################################

onesamp=function(y, mr=10, mu0=NA, alpha=0.05, pen=NA, penvals=NA){
  ans=new.env()
  n=length(y)
# use median as the starting value, if not given by user
  if (is.na(mu0)) mu0=median(y)    
# place knots
  xs=abs(y-mu0)
  S=1.5*quantile(xs,0.98)  
  J=round(n^(1/7)*10)
  kn=S*(mr^(seq(0,J-1)/(J-2))-1)/(mr^((J-1)/(J-2))-1)
  k=J-2
  x=sort(xs)
# make candiate values of penalty, if not given by user
  if ((is.na(penvals) && is.na(pen) )) {
    penbasis=c(n^(-12/7)*exp(1)^seq(-25,10,length=30))   
    penvals=penbasis*quantile(x,0.9)^5    
  }
# Create the penalty D matrix
  etamat=matrix(0,nrow=k+1,ncol=2)
  for(j in 1:k){
    etamat[j,1]=-2/(kn[j+2]-kn[j])/(kn[j+1]-kn[j])
    etamat[j,2]=2/(kn[j+2]-kn[j])/(kn[j+2]-kn[j+1])
  }
  etamat[k+1,1]=-2/(kn[k+2]-kn[k+1])^2
  dmat=matrix(0,nrow=J-2,ncol=J)
  dmat[1,1]=etamat[1,2]-etamat[1,1]
  dmat[1,2]=etamat[2,1]
  for(i in 2:(J-2)){
    dmat[i,i-1]=-etamat[i-1,2]
    dmat[i,i]=etamat[i,2]-etamat[i,1]
    dmat[i,i+1]=etamat[i+1,1]
  }
# create the H matrix
  hmat=matrix(nrow=k+2,ncol=k+2)
  for(j in 1:k){
    hmat[j,j] = kn[j+1]-2*(kn[j+1]-kn[j])^2/(kn[j+2]-kn[j])/3 +(kn[j+1]-kn[j])**3/(kn[j+2]-kn[j])^2/5+(kn[j+2]-kn[j+1])^3/(kn[j+2]-kn[j])**2/5
  }
  hmat[k+1,k+1] = (8*kn[k+2]+7*kn[k+1])/15
  hmat[k+2,k+2]=kn[k+2]
  if(k>1){for(j in 1:(k-1)){
    hmat[j,j+1]=kn[j+1] + ((kn[j+2]-kn[j+1])^2 - (kn[j+1]-kn[j])^2) / (kn[j+2]-kn[j])/3 -(kn[j+2]-kn[j+1])^3 / (kn[j+3]-kn[j+1])/(kn[j+2]-kn[j]) /30
  }}
  hmat[k,k+1]=kn[k+1]+ ((kn[k+2]-kn[k+1])^2 - (kn[k+1]-kn[k])^2) / (kn[k+2]-kn[k])/3 -(kn[k+2]-kn[k+1])^2 /(kn[k+2]-kn[k]) /30
  if(k>1){for(j in 1:(k-1)){
    for(l in (j+2):(k+1)){
      hmat[j,l] = (kn[j]+kn[j+1]+kn[j+2])/3
    }
  }}
  for(j in 1:k){
    hmat[j,k+2]=(kn[j]+kn[j+1]+kn[j+2])/3
  }	
  hmat[k+1,k+2]=(2*kn[k+2]+kn[k+1])/3
  for(i in 2:(k+2)){
    for(j in 1:(i-1)){
      hmat[i,j]=hmat[j,i]
    }
  }
# create delta matrix
  delta=matrix(0,nrow=J,ncol=n)	
  delta[J,]=1:n*0+1
  for(j in 1:k){
    ind=x<=kn[j]
    delta[j,ind] = 1
    ind=x>kn[j]&x<=kn[j+1]
    delta[j,ind] = 1 - (x[ind]-kn[j])^2 / (kn[j+2]-kn[j]) / (kn[j+1]-kn[j])
    ind=x>kn[j+1]&x<=kn[j+2]
    delta[j,ind] = (x[ind]-kn[j+2])^2/(kn[j+2]-kn[j+1])/(kn[j+2]-kn[j])
    ind=x>kn[j+2]
    delta[j,ind]=0
  }
  ind=x<=kn[k+1]
  delta[k+1,ind]=1
  ind=x>kn[k+1]
  delta[k+1,ind]=1-(x[ind]-kn[k+1])^2/(kn[k+2]-kn[k+1])^2
# create the z vector
  one=1:n*0+1/n
  zvec=delta%*%one
# find penalty using 10-fold cross-validation
  if (is.na(pen)) {
    nv=trunc(n/10)
    fold=rep(1:10,nv)
    if(nv<n/10){fold=c(fold,1:(n-10*nv))}
    nf=1:10;for(i in 1:10){nf[i]=sum(fold==i)}
# scramble
    fold=sample(fold,n)
    risk=1:length(penvals)*0
    for(ell in 1:length(penvals)){
      hpen=(hmat+penvals[ell]*t(dmat)%*%dmat)
      hj=hpen[,J]
      b0=rep(1,J)/sum(hj)
      wmat=rbind(hj[-1],diag(-hj[1],J-1))
      wmat=t(t(wmat)/sqrt(apply((wmat)^2,2,sum)))
      qmat=t(wmat)%*%hpen%*%wmat
      cvec=t(wmat)%*%(zvec-hpen%*%b0)
      betavec=solve.QP(qmat,cvec,t(wmat),-b0)$solution
      bhat=(wmat)%*%betavec+b0
      for(i in 1:10){
        del=delta[,fold!=i]
        nm=n-nf[i]
        onef=1:nm*0+1/nm
        zvecf=del%*%onef
        cvecf=t(wmat)%*%(zvecf-hpen%*%b0)
        betavec=solve.QP(qmat,cvecf,t(wmat),-b0)$solution
        bhatf=(wmat)%*%betavec+b0
        gminus=(t(delta)%*%bhatf)[fold==i]
        risk[ell]=risk[ell]+sum(gminus)
      }
      risk[ell]=t(bhat)%*%hmat%*%bhat-2/n*risk[ell]
    }
# find the largest local minimizer of candiate values
    pen=penvals[which.min(risk)]
    lp=length(penvals)
    for(i in 0:(lp-2)){
      if ((risk[lp-1-i]-risk[lp-i])/abs(risk[lp-i]) > -1e-4) {
        target=lp-i
        pen=penvals[target]
        break
      }
    }
# output risk for diagnostics of cross-validation
    ans$risk=risk                                                                   
  }
# estimate error density
  hpen=(hmat+pen*t(dmat)%*%dmat)
  hj=hpen[,J]
  b0=rep(1,J)/sum(hj)
  wmat=rbind(hj[-1],diag(-hj[1],J-1))
  wmat=t(t(wmat)/sqrt(apply((wmat)^2,2,sum)))
  qmat=t(wmat)%*%hpen%*%wmat
  cvec=t(wmat)%*%(zvec-hpen%*%b0)
  betavec=solve.QP(qmat, cvec, t(wmat), -b0)$solution
  bhat=(wmat)%*%betavec+b0
# estimate the mean
  sy=sort(y)
  lwr=sy[round(n*.2)];upr=sy[round(n*.8)]
  fans=optimize(findmu,bhat=bhat,y=y,kn=kn,lower=lwr,upper=upr) 
  muhat=fans$minimum
# compute confidence interval
  b=bhat/2
  de=integrate(fd,0,max(kn),kn=kn,k=k,b=b,stop.on.error=FALSE)$value^2*4*n
  nu=integrate(fn,0,max(kn),kn=kn,k=k,b=b,stop.on.error=FALSE)$value*2
  p.var=nu/de      
  CI=c(muhat-qnorm(1-alpha/2)* sqrt(p.var),muhat+qnorm(1-alpha/2)*sqrt(p.var))
# outputs
  ans$confidence.intervals=CI
  ans$knots=kn
  ans$bhat=bhat
  ans$muhat=muhat
  ans$pen=pen                                
  ans$penvals=penvals
# output estimated density function g
  density.out=function(dp){
    dp=abs(dp)
    if(any(dp>max(kn))) return("Input is out of support")
    k=length(kn)-2
    b=bhat/2
    n=length(dp)
    delta=matrix(0,nrow=k+2,ncol=n)
    for(j in 1:k){
      ind=dp<=kn[j]
      delta[j,ind] = 1
      ind=dp>kn[j]&dp<=kn[j+1]
      delta[j,ind]=1-(dp[ind]-kn[j])^2/(kn[j+2]-kn[j])/(kn[j+1]-kn[j])
      ind=dp>kn[j+1]&dp<=kn[j+2]
      delta[j,ind]=(kn[j+2]-dp[ind])^2/(kn[j+2]-kn[j+1])/(kn[j+2]-kn[j])
    }
    ind=dp<=kn[k+1]
    delta[k+1,ind]=1
    ind=dp>kn[k+1]&dp<=kn[k+2]
    delta[k+1,ind]=1-(dp[ind]-kn[k+1])^2/(kn[k+2]-kn[k+1])^2
    g=(t(b)%*%delta)
    as.numeric(g)
  }
  ans$fhat=density.out
  ans
}


#################################
# objective function: for finding the mean
#################################
findmu=function(mu,bhat,y,kn){
  n=length(y)
  x=sort(abs(y-mu))
  J=length(kn);k=J-2
  delta=matrix(0,nrow=J,ncol=n)	
  delta[J,]=1:n*0+1
  for(j in 1:k){
    ind=x<=kn[j]
    delta[j,ind] = 1
    ind=x>kn[j]&x<=kn[j+1]
    delta[j,ind] = 1 - (x[ind]-kn[j])^2 / (kn[j+2]-kn[j]) / (kn[j+1]-kn[j])
    ind=x>kn[j+1]&x<=kn[j+2]
    delta[j,ind] = (x[ind]-kn[j+2])^2/(kn[j+2]-kn[j+1])/(kn[j+2]-kn[j])
    ind=x>kn[j+2]
    delta[j,ind]=0
  }
  ind=x<=kn[k+1]
  delta[k+1,ind]=1
  ind=x>kn[k+1]
  delta[k+1,ind]=1-(x[ind]-kn[k+1])^2/(kn[k+2]-kn[k+1])^2
  -sum(t(bhat)%*%delta)
}


###########################
# integral functions: for confidence intervals
###########################
fn=function(x,kn,k,b){
  n=length(x)
  delta=matrix(0,nrow=k+2,ncol=n)
  delta_1=matrix(0,nrow=k+2,ncol=n)
  for(j in 1:k){
    ind=x<=kn[j]
    delta[j,ind] = 1
    ind=x>kn[j]&x<=kn[j+1]
    delta[j,ind]=1-(x[ind]-kn[j])^2/(kn[j+2]-kn[j])/(kn[j+1]-kn[j])
    delta_1[j,ind] = 2*(x[ind]-kn[j]) / (kn[j+2]-kn[j]) / (kn[j+1]-kn[j])
    ind=x>kn[j+1]&x<=kn[j+2]
    delta[j,ind]=(kn[j+2]-x[ind])^2/(kn[j+2]-kn[j+1])/(kn[j+2]-kn[j])
    delta_1[j,ind] = -2*(x[ind]-kn[j+2]) /(kn[j+2]-kn[j+1])/(kn[j+2]-kn[j])
  }
  ind=x<=kn[k+1]
  delta[k+1,ind]=1
  ind=x>kn[k+1]&x<=kn[k+2]
  delta[k+1,ind]=1-(x[ind]-kn[k+1])^2/(kn[k+2]-kn[k+1])^2
  delta_1[k+1,ind]= 2*(x[ind]-kn[k+1]) /(kn[k+2]-kn[k+1])^2
  g=t(b)%*%delta
  g1=t(b)%*%delta_1
  result=g1^2*g
  result
}

fd=function(x,kn,k,b){
  n=length(x)
  delta=matrix(0,nrow=k+2,ncol=n)
  delta_2=matrix(0,nrow=k+2,ncol=n)
  for(j in 1:k){
    ind=x<=kn[j]
    delta[j,ind] = 1
    ind=x>kn[j]&x<=kn[j+1]
    delta[j,ind]=1-(x[ind]-kn[j])^2/(kn[j+2]-kn[j])/(kn[j+1]-kn[j])
    delta_2[j,ind] = 2 / (kn[j+2]-kn[j]) / (kn[j+1]-kn[j])
    ind=x>kn[j+1]&x<=kn[j+2]
    delta[j,ind]=(kn[j+2]-x[ind])^2/(kn[j+2]-kn[j+1])/(kn[j+2]-kn[j])
    delta_2[j,ind] = -2 /(kn[j+2]-kn[j+1])/(kn[j+2]-kn[j])
  }
  ind=x<=kn[k+1]
  delta[k+1,ind]=1
  ind=x>kn[k+1]&x<=kn[k+2]
  delta[k+1,ind]=1-(x[ind]-kn[k+1])^2/(kn[k+2]-kn[k+1])^2
  delta_2[k+1,ind]= 2 /(kn[k+2]-kn[k+1])^2
  g=t(b)%*%delta
  g2=t(b)%*%delta_2
  result=g2*g
  result
}
