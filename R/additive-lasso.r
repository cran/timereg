lasso.boostR<-function(D,d,lambda,max.it=10,beta=0,prnt=0)
{ nkol<-nrow(D); if (sum(beta)==0) beta<-matrix(0,nkol,1); 
  i<-1; outbreak<-9; 

  while ( i<=max.it)
    {
      Db<-D%*%beta; bd<-sum(beta*d); bDb<-t(beta)%*%Db; 

      g<--(Db-d); mg<-apply(cbind(g,-g),1,min); 
      index<-order(mg,decreasing=0)[1]; gamma<-sign(g[index]); 
      if (abs(g[index])<0.00000001) {outbreak<-0; break;}
      lg<-gamma*lambda
      if (prnt==1) print(index); 

      k<-(Db[index,1]*lg-bDb+bd-lg*d[index])/
        (-bDb-lg^2*D[index,index]+2*lg*Db[index,1]); k<-c(k)
      if (prnt==1) print(c(k,lg)); 
      if (prnt==1) print(c(bDb,bd)); 
      if (prnt==1) print(c(D[index,index],d[index])); 

      L1<-(1/2)*lg^2*D[index,index]-lg*d[index]
      L0<-(1/2)*bDb-bd
      Lk<-(1/2)*((1-k)^2*bDb+k^2*lg^2*D[index,index]+
                 2*k*(1-k)*lg*Db[index,1])-(1-k)*bd-k*lg*d[index]

      if (k<0 | k >1) Lk<-L1+1; 
      mi<-order(c(L0,Lk,L1),decreasing=FALSE)[1]
      k<-c(0,k,1)[mi]
      if (prnt==1) print(round(c(L0,L1,Lk,k),3)); 
                                        #if (k==0) {outbreak<-1; break;}

      beta<-(1-k)*beta; beta[index]<-beta[index]+k*lg
      i<-i+1
    }
  out<-list(beta=c(beta),outbreak=outbreak,loops=i,Lbeta=L0,
            L1=L1,Lk=Lk,k=k,index=index,absdL=mg[index],l1=sum(abs(beta)))
  return(out)
}

L<-function(D,d,beta) {
  return(1/2*beta %*% D %*% beta-sum(beta*d))}
dL<-function(D,d,beta){return( D %*% beta -d)}

lasso.boost<-function(D,d,lambda,max.it=10000,beta=0,detail=0)
{
  nkol<-nrow(D); if (sum(beta)==0) beta<-matrix(0,nkol,1); 
                                        #dyn.load("kim.so"); 
  out<-.C("l1boost", as.double(D),as.integer(nkol),as.double(d),
          as.double(lambda),as.integer(max.it),as.double(beta),
          as.integer(detail),PACKAGE="timereg")
  beta<-out[[6]]; 
  return( list(beta=beta,L=L(D,d,beta),l1=sum(abs(beta))))
}

###############################################
#################### LASSO ####################
###############################################
lasso.add.hazard<-function(D,d,lambda,ridge,start.sign,max.it=10)
{
  #require(quadprog);
  nkol<-nrow(D); Amat<-matrix(start.sign,nkol,1)
  b0<-c(); l1<-lambda+1; i<-1

  while (l1 > lambda & i<max.it)
    {
      print(c(i,l1))
      b0<-c(b0,lambda);
      i<-i+1
      if (i>2) Amat<-cbind(Amat,sign(ud1$solution))

      ud1<-solve.QP(D+diag(ridge,nkol),d,-Amat,-b0)
      l1<-sum(abs(ud1$solution))
    }
  return(list(sol=ud1,D=D,d=d,ridge=ridge,lambda=lambda,l1=l1,Amat=Amat)
         )
}

surv.lars<-function(time, status, z, l1.weights=NULL, ...) 
{
    antpers=length(time); p<-ncol(z); 
    if (is.null(l1.weights)==FALSE) adap<-1 else adap<-0; 

    out<-des.aalen(Surv(time,status)~const(z))
    y<-out$Y; x<-as.matrix(out[,1:p]);
    if (adap==1) x<-t(t(x)/l1.weights)
    fit <- mylars(x, y , trace = trace, ...)
    if (adap==1) fit$coef<-t(t(fit$coef)/l1.weights)
return(fit)
}

surv.lars.cv<-function(time, status, z, l1.weights=NULL, K = 10, 
fraction = seq(from = 0, to = 1, length = 100), trace = FALSE, 
plot.it = TRUE, se = TRUE, ...) 
{
    antpers=length(time); p<-ncol(z); all.folds <- cv.folds(antpers, K)
    residmat <- matrix(0, length(fraction), K); 
    if (is.null(l1.weights)==FALSE) adap<-1 else adap<-0; 
    for (i in seq(K)) {
      omit <- all.folds[[i]]
      fit<-surv.lars(time[-omit], status[-omit], z[-omit], 
      l1.weights=l1.weights) 
 #     out<-des.aalen(Surv(time[-omit],status[-omit])~const(z[-omit,]))
 #     y<-out$Y; x<-as.matrix(out[,1:p]);
 #     if (adap==1) x<-t(t(x)/l1.weights)
 #     fit <- mylars(x, y , trace = trace, ...)
 #     if (adap==1) fit$beta<-t(t(fit$beta)/l1.weights)
      out<-des.aalen(Surv(time[omit],status[omit])~const(z[omit,]))
      y<-out$Y; x<-as.matrix(out[,1:p]);
      fit <- mypredict.lars(fit, x,mode="fraction", s = fraction)$fit
      if (length(omit) == 1) fit <- matrix(fit, nrow = 1)
      residmat[, i] <- apply((y - fit)^2, 2, mean)

      if (trace) cat("\n CV Fold", i, "\n\n")
    }
    cv <- apply(residmat, 1, mean); 
    cv.error <- sqrt(apply(residmat, 1, var)/K)
    object <- list(fraction = fraction, cv = cv, cv.error = cv.error,
                   cv.frac=fraction[order(cv)][1])
    if (plot.it) plotCVLars(object, se = se)
    invisible(object)
}
#beta0<-c(0.1,0.3,0.5); 
#ud<-makeadd(100,9,beta=beta0,covs=3,cens=1.51,rho=0.3,sig=0)
#p<-ncol(ud$dat); des<-as.matrix(ud$dat[,3:p]);
#outa<-aalen(Surv(time,status)~const(des),ud$dat,n.sim=0,robust=0)
#cvala<-surv.lars.cv(ud$dat$time,ud$dat$status,des,c(outa$gamma))
#par(mfrow=c(1,2)); plot(scale(cvala$cv)); plot(scale(cvala$cv2))

surv.lars.gcv<-function(time, status, z, K = 10, 
fraction = seq(from = 0, to = 1, length = 100), 
trace = FALSE, plot.it = TRUE, se = TRUE, ...) 
{   
    cat(" does not work !!! \n"); 
    cvg <-pl<- rep(0, length(fraction))
    out<-des.aalen(Surv(time,status)~const(z))
    p<-ncol(z); y<-out$Y; x<-as.matrix(out[,1:p]);
    fitl <- lars(x, y , trace = trace); # , ...)
    totsum=sum(abs(fitl$beta[p,]))
    fit <- predict(fitl, x, mode = "fraction", s = fraction)$fit
    if (length(fit) == 1) fit <- matrix(fit, nrow = 1)
    lambda<-fraction*totsum; 
    betas<-coef.lars(fitl,mode = "fraction", s = fraction)
    k=0; 
    for (ss in fraction) {
       k<-k+1; W<-diag(abs(betas[k,])); 
       pl[k]<-sum(diag(
       z %*% solve( t(z) %*% z + lambda[k]* ginv(W) )%*% t(z) )) }
    rss<-apply((y - fit)^2, 2, mean);
    cvg <- rss/ (1 - pl/length(time))^2; 
    if (plot.it) plot(fraction,cvg, type="l",...)
    object <- list(fraction = fraction, rss=rss, cvg = cvg,
                   cvg.frac=fraction[order(cvg)][1])
}

ginv = function(X, tol = sqrt(.Machine$double.eps)){
  s = svd(X)
  nz = s$d > tol * s$d[1]
  if(any(nz)) s$v[,nz] %*% (t(s$u[,nz])/s$d[nz]) else X*0
}
