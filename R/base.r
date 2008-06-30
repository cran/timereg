coefBase<- function(object, digits=3, d2logl=0) {
  res <- cbind(object$gamma,
               diag(object$var.gamma)^0.5,
               diag(object$robvar.gamma)^0.5)
  if (d2logl==1) res<-cbind(res,diag(object$D2linv)^.5)
  wald <- object$gamma/diag(object$var.gamma)^0.5
  waldp <- (1 - pnorm(abs(wald))) * 2
  res <- as.matrix(cbind(res, wald, waldp))
  if (d2logl==1) colnames(res) <- c("Coef.", "SE", "Robust SE","D2log(L)^-1","z","P-val") else colnames(res) <- c("Coef.", "SE", "Robust SE", "z", "P-val")
  prmatrix(signif(res, digits))
}

timetest<-function(object,digits=3)
{ 
  cat("Test for nonparametric terms \n")
  if (is.null(object$conf.band)==TRUE)  mtest<-FALSE else mtest<-TRUE;
  if (mtest==FALSE) cat("Test not computed, sim=0 \n\n")
  if (mtest==TRUE) {
  test0<-cbind(object$obs.testBeq0,object$pval.testBeq0)
  testC<-cbind(object$obs.testBeqC,object$pval.testBeqC)
  colnames(test0)<-c("sup|  hat B(t)/SD(t) |","p-value H_0: B(t)=0")
  colnames(testC)<-c("sup| B(t) - (t/tau)B(tau)|","p-value H_0: B(t)=b t")
  if (is.null(object$obs.testBeqC.is)!=TRUE)  {
  testCis<-cbind(object$obs.testBeqC.is,object$pval.testBeqC.is)
  colnames(testCis) <-c("int  (B(t)-(t/tau)B(tau))^2dt","p-value H_0: B(t)=b t")
  }
  cat("\n")
  cat("Test for non-significant effects \n")
  prmatrix(signif(test0,digits))
  cat("\n")
  cat("Test for time invariant effects \n")
  prmatrix(signif(testC,digits))
  if (is.null(object$obs.testBeqC.is)!=TRUE)  prmatrix(signif(testCis,digits))
  cat("\n")
}
}
