########################################################################
# Code to do auto-adaptive M-estimation (osar)
#
#  must-have inputs: y = the data
#
#  optional inputs: p = the order of autoregressive coefficients
#                   mr = the mesh ratio (default value 10)
#                   pen = the penalty constant (o/w will use cross-validation to find one)
#                   penvals = vector of candidate values of penalty (o/w will get one automatically)
#                   nloops = maximum number of Cochrane-Orcutt flavored loops 
#
#  returns: bhat = the estiamted spline density coefficients
#           fhat = estimated error density function
#           knots = planced knots for spline density estiamtion
#           muhat = estiamted mean for the data
#           phihat = estiamted autoregressive phi for the data
#           mvec = record of every muhat in loops
#           pmat = record of every phihat in loops
#           pen = selected penalty constant by cross-validation
#           penvals = automatically generated candidate values of penalty
#           risk = estimated risks for cross-validation
#
#
########################################################################
## need library(quadprog)
########################################################################

osar=function(y, p=1, mr=10, pen=NA, penvals=NA, nloops=20, figures=FALSE){ # mr:10->15
  ans=new.env()
  n=length(y)
  sy=sort(y)
  lwr=sy[round(n*.2)];upr=sy[round(n*.8)] # may need to update again later!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  # burn-in stage 0
  mu0=optimize(find0,lower=lwr,upper=upr,phi=rep(0,p),y=y,p=p)$minimum
  phi0=optim(par=rep(0,p),fn=find0,lower=rep(-1,p),upper=rep(1,p),method=ifelse(p>1,"L-BFGS-B","Brent"),y=y,p=p,mu=mu0)$par
  mvec=mu0
  pmat=matrix(phi0,p)
  for (i in 1:nloops) {
    m.old=mvec[length(mvec)]
    p.old=pmat[,dim(pmat)[2]]
    m.new=optimize(find0,lower=lwr,upper=upr,phi=p.old,y=y,p=p)$minimum
    mvec=c(mvec,m.new)
    p.new=optim(par=p.old,fn=find0,lower=rep(-1,p),upper=rep(1,p),method=ifelse(p>1,"L-BFGS-B","Brent"),y=y,p=p,mu=m.new)$par
    pmat=cbind(pmat,p.new)
    if(abs((m.new-m.old)/m.old)<1e-3 & all(abs((p.new-p.old)/p.old)<1e-3) ) break
  } 
  # place knots
  mu0=mvec[length(mvec)]
  phi0=pmat[,dim(pmat)[2]]
  temp.y=y[(p+1):n]
  for (i in 1:p) {
    temp.y=temp.y-phi0[i]*y[(p+1-i):(n-i)]
  }
  xs=abs(temp.y-(1-sum(phi0))*mu0)
  n=length(xs)
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
  # burn-in stage 1
  for (i in 1:nloops) {
    m.old=mvec[length(mvec)]
    p.old=pmat[,dim(pmat)[2]]
    m.new=optimize(find1,lower=lwr,upper=upr,phi=p.old,y=y,p=p,kn=kn,bhat=bhat)$minimum
    mvec=c(mvec,m.new)
    p.new=optim(par=p.old,fn=find1,lower=rep(-1,p),upper=rep(1,p),method=ifelse(p>1,"L-BFGS-B","Brent"),y=y,p=p,mu=m.new,kn=kn,bhat=bhat)$par
    pmat=cbind(pmat,p.new)
    if(abs((m.new-m.old)/m.old)<1e-3 & all(abs((p.new-p.old)/p.old)<1e-3) ) break
  } 
  mu0=mvec[length(mvec)]
  phi0=pmat[,dim(pmat)[2]]
  # final joint-estimation after all burn-in preparations
  fans=optim(par=c(mu0,phi0),fn=find2,lower=c(lwr,rep(-1,p)),upper=c(upr,rep(1,p)),method="L-BFGS-B",y=y,p=p,kn=kn,bhat=bhat)
  muhat=fans$par[1]
  phihat=fans$par[-1]
  mvec=c(mvec,muhat)
  pmat=cbind(pmat,phihat)
  #---------------- plots
  if (figures==TRUE & p==1){ # marginal surface
    par(mfrow=c(1,3),mar=c(3,2.5,2,1),cex=1,cex.main=1,las=0,mgp=c(1.5,0.5,0))
    mu.grid=seq(min(y),max(y),length=200)
    phi.grid=seq(-1,1,length=200)
    findmu.plot=findphi.plot=rep(NA,200)
    for (ii in 1:200) {
      findmu.plot[ii]=find1(mu=mu.grid[ii],phi=phihat,y=y,p=p,bhat=bhat,kn=kn)
      findphi.plot[ii]=find1(mu=muhat,phi=phi.grid[ii],y=y,p=p,bhat=bhat,kn=kn)
    }
    plot(findmu.plot~mu.grid,type="l",yaxt='n',ylab = "",xlab=expression(mu),main="(a)",lwd=2)
    plot(findphi.plot~phi.grid,type="l",yaxt='n',ylab = "",xlab=expression(phi),main="(b)",lwd=2)
    mu.grid=seq(lwr,upr,length=100)
    phi.grid=seq(-1,1,length=100)
    crit=matrix(NA,length(mu.grid),length(phi.grid))
    for (i in 1:length(mu.grid)) {
      for (j in 1:length(phi.grid)) {
        crit[i,j]=find2(c(mu.grid[i],phi.grid[j]),y=y,p=p,bhat=bhat,kn=kn) 
      }
    }
    contour(mu.grid, phi.grid, crit,nlevels = 50, xlab=expression(mu),ylab=expression(phi),main="(c)")
  }
  
  # outputs
  ans$knots=kn
  ans$bhat=bhat
  ans$muhat=muhat
  ans$phihat=phihat
  ans$mvec=mvec
  ans$pmat=pmat
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
# objective functions
#################################
find0=function(mu,phi,p,y){
  n=length(y)
  temp.y=y[(p+1):n]
  for (i in 1:p) {
    temp.y=temp.y-phi[i]*y[(p+1-i):(n-i)]
  }
  x=sort(abs(temp.y-(1-sum(phi))*mu))
  sum(abs(x))
}

find1=function(mu,phi,p,y,bhat,kn){
  n=length(y)
  temp.y=y[(p+1):n]
  for (i in 1:p) {
    temp.y=temp.y-phi[i]*y[(p+1-i):(n-i)]
  }
  x=sort(abs(temp.y-(1-sum(phi))*mu))
  n=length(x)
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

find2=function(est,p,y,bhat,kn){
  mu=est[1]
  phi=est[-1]
  n=length(y)
  temp.y=y[(p+1):n]
  for (i in 1:p) {
    temp.y=temp.y-phi[i]*y[(p+1-i):(n-i)]
  }
  x=sort(abs(temp.y-(1-sum(phi))*mu))
  n=length(x)
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
