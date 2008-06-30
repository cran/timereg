cum.residuals<-function(object,data=sys.parent(),modelmatrix=0,cum.resid=0,
n.sim=500,weighted.test=1,start.design=1)
{
  if (!(class(object)!="aalen" | class(object)!="timecox" |
        class(object)!="cox.aalen" ))
    stop ("Must be output from aalen() timecox() or cox.aalen() functions\n") 
  if (class(object)=="timecox") if (object$method!="basic") 
    stop("Residuals available only for method=basic\n")
  if (class(object)=="timecox") if (is.null(object$gamma)==FALSE) 
    stop("Residuals available only for timecox model with no const terms\n")
  if (class(object)=="aalen") if (is.null(object$gamma)==FALSE) 
    stop("Residuals available only for Aalen model with no const terms\n")
  if (is.null(object$residuals$dM)==TRUE) 
    stop("Residuals not computed, add option residuals=1\n");

  if (sum(modelmatrix)==0 && cum.resid==0) 
    stop("No modelmatrix or continous covariates given to cumulate residuals\n"); 

  if (class(object)=="cox.aalen") {
    dcum<-apply(as.matrix(object$cum[,-1]),2,diff); 
    beta<-object$gamma; coxaalen<-1; } else {
      dcum<-0; beta<-0; coxaalen<-0; pg<-0;designG<-0; }

  id<-attr(object,"id"); cluster<-attr(object,"cluster"); 
  formula<-attr(object,"Formula"); 
  start.time<-attr(object,"start"); 
  pers<-unique(id); antpers<-length(pers); 
  clust<-unique(cluster); antclust<-length(clust); 

if (class(object)=="cox.aalen") 
  ldata<-aalen.des(formula,data,model="cox.aalen") else ldata<-aalen.des(formula,data) 

  X<-ldata$X; 
  time<-ldata$time; time2<-ldata$time2; 
  covar<-X; 

  status<-ldata$status; 
  if (coxaalen==1) {
    designG<-ldata$Z;  covnamesZ<-ldata$covnamesZ; 
    pg<-ncol(designG); }
  Ntimes <- sum(status); 

  if (sum(duplicated(time2[status==1]))>0) {
    cat("Ties may cause difficulties, break them ! \n"); 
  }

  if (ldata$type == "right") {
    cat("Cumulative martingale residuals for Right censored survival times\n");
  } else {cat("Counting process style data\n"); }
  times<-c(start.time,time2[status==1]); times<-sort(times);

  antpers=length(time); ntot<-nrow(X); px<-ldata$px

  lmgresids<-length(object$residuals$time); 

  if (coxaalen==1) gamma.iid<-object$residuals$gamma.iid else gamma.iid<-0;

  if (coxaalen==1) {covar<-cbind(X,designG);
                    ptot<-px+pg;} else {covar<-covar; ptot<-px;}

  covnames0<-colnames(covar); 
  antal<-0;  maxval<-0; intercept<-0; antal2<-0;  maxval2<-0; 
  xvals<-list();  ant<-rep(0,ptot);  
  for (i in 1:ptot) xvals[[i]]<-c(); 
  rani<-round(runif(1)*10000);
  test<-matrix(0,n.sim,ptot); uni.test<-matrix(0,n.sim,ptot); 
  testOBS<-rep(0,4*ptot); 
                                        # testOBS size and location of supremum, test simulated sup's
  unitime.test<-time.test<-mult.test<-multtime.test<-test
  unitime.testOBS<-uni.testOBS<-time.testOBS<-mult.testOBS<-
    multtime.testOBS<-testOBS

  keepcumz<-c()
  for (j in 1:ptot) { 
    z<-unique(covar[,j]); z<-sort(z); 
    xvals[[j]]<-z; antal<-antal+1; ant[j]<-length(z); 
    if (ant[j]>2) keepcumz<-c(keepcumz,j); 
    maxval<-max(maxval,length(z)); }
  xval<-matrix(0,maxval,ptot); k<-1; 
  for (i in 1:ptot) {xval[1:ant[i],k]<-xvals[[i]]; k<-k+1;}
  if (sum(keepcumz)==0 && cum.resid==1) 
    stop(" No continous covariates given to cumulate residuals \n"); 

  univar.proc<-matrix(0,maxval,ptot); 
  sim.univar.proc<-matrix(0,maxval,50*ptot); 
  time.proc<-matrix(0,lmgresids,2); sim.time.proc<-matrix(0,lmgresids,50); 

  if (sum(modelmatrix)==0) {modelmatrix<-0;model<-0;pm<-1;}  else 
  {model<-1; modelmatrix<-as.matrix(modelmatrix); pm<-ncol(modelmatrix); 
   test<-matrix(0,n.sim,3*pm); testOBS<-rep(0,2*pm); 
   covnames<-colnames(modelmatrix); 
 }

  Ut<-cummgt<-matrix(0,lmgresids,pm+1); 
  robvarcum<-matrix(0,lmgresids,pm+1); 
  simUt<-matrix(0,lmgresids,pm*50); 

                                        #dyn.load("linmgresid.so");

  mgout<- .C("mgresid",
             as.double(X),as.integer(ntot),as.integer(px), 
             as.integer(antpers),as.double(time),as.double(time2),
             as.integer(status),as.integer(id),as.double(object$residuals$time),
             as.integer(lmgresids),as.double(object$residuals$dM),as.integer(n.sim),
             as.integer(rani),as.double(xval), as.integer(ant), 
             as.double(univar.proc),as.double(time.proc),as.double(sim.univar.proc),
             as.double(sim.time.proc),as.double(uni.test),as.double(uni.testOBS),
             as.double(time.test),as.double(time.testOBS),as.double(unitime.test),
             as.double(unitime.testOBS),as.double(modelmatrix),as.integer(model),
             as.integer(pm),as.double(cummgt),as.double(object$residuals$dM.iid),
             as.double(robvarcum),as.double(testOBS),as.double(test),
             as.double(simUt),as.double(Ut),as.integer(cum.resid),
             as.integer(maxval),as.integer(start.design),as.integer(coxaalen),
             as.double(dcum),as.double(beta),as.double(designG),
             as.integer(pg),as.double(gamma.iid),as.integer(cluster),
             as.integer(antclust),PACKAGE="timereg")


  if (model==1) {
    cum<-matrix(mgout[[29]],lmgresids,pm+1); 
    robvar.cum<-matrix(mgout[[31]],lmgresids,pm+1);
    var.cum<-robvar.cum; Ut<-matrix(mgout[[35]],lmgresids,pm+1);

    colnames(Ut)<-colnames(cum)<-colnames(var.cum)<-
      colnames(robvar.cum)<- c("time",covnames)
    test.procBeq0<-Ut; 

    sim<-1; 
    if (sim>=1) {
      Uit<-matrix(mgout[[34]],lmgresids,50*pm);
      UIt<-list(); for (i in (0:49)*pm) UIt[[i/pm+1]]<-as.matrix(Uit[,i+(1:pm)]);
      testOBS<-mgout[[32]];
      test<-matrix(mgout[[33]],n.sim,3*pm);
      testval<-c(); unifCI<-c(); 
      for (i in 1:(2*pm)) testval<-c(testval,pval(test[,i],testOBS[i]))
      for (i in 1:pm) unifCI<-as.vector(c(unifCI,percen(test[,2*pm+i],0.95)));
      obs.testBeq0<-as.vector(testOBS[1:pm]);
      obs.testBeq0.is<-as.vector(testOBS[(pm+1):(2*pm)]);
      pval.testBeq0<-as.vector(testval[1:pm]);
      pval.testBeq0.is<-as.vector(testval[(pm+1):(2*pm)]);
      sim.testBeq0<-test[,(2*pm+1):(3*pm)]; 

      sim.test.procBeq0<-UIt; 
      names(unifCI)<- names(pval.testBeq0)<- names(obs.testBeq0)<- 
        names(pval.testBeq0.is)<- names(obs.testBeq0.is)<- covnames
    } else {
      test<-unifCI<-Ut<-UIt<-pval.testBeq0<-pval.testBeq0.is<-obs.testBeq0<-
        obs.testBeq0.is<-sim.testBeq0<-NULL; }
  } else {
    cum<-robvar.cum<-test<-unifCI<-Ut<-UIt<-pval.testBeq0<-
      pval.testBeq0.is<-obs.testBeq0<-obs.testBeq0.is<-sim.testBeq0<-NULL; 
  } 

  if (cum.resid>=1) {
    univar.p<-matrix(mgout[[16]],maxval,ptot)[,keepcumz];
    univar.p<-as.matrix(univar.p); 
    univar.proc<-list(); 
    k<-1
    for (i in keepcumz) {
      univar.proc[[k]]<-cbind(xvals[[i]],univar.p[1:ant[i],k]); 
      colnames(univar.proc[[k]])<-c(covnames0[i],"cum. martingale residual"); 
      k<-k+1; 
    }
    Uiz<-matrix(mgout[[18]],maxval,50*ptot);
    UIz<-list(); 
    k<-1; 
    for (i in keepcumz) {UIz[[k]]<-matrix(Uiz[1:ant[i],i+(0:49)*ptot],ncol=50);
                         k<-k+1;}
    uni.test<-matrix(mgout[[20]],n.sim,ptot)[,keepcumz]; 
    uni.test<-as.matrix(uni.test); 
    uni.testOBS<-mgout[[25]][1:ptot][keepcumz]; 
    testval<-c(); 
    lkeep<-length(keepcumz); 

    for (i in 1:lkeep)  testval<-c(testval,pval(uni.test[,i],uni.testOBS[i])) 
    uni.pval<-testval
    names(uni.testOBS)<- names(uni.pval)<- colnames(uni.test)<-covnames0[keepcumz]; 
                                        #uni.testOBS<-uni.testOBS[keepcumz]; uni.pval<-uni.pval[keepcumz]
                                        #uni.test<-uni.test[,keepcumz];
  } else {
    uni.testOBS<-uni.pval<-proc.cumz<-UIz<-unitime.pval<-unitime.testOBS<-NULL;}

  ud<-list(cum=cum,robvar.cum=robvar.cum,
           pval.testBeq0=pval.testBeq0, obs.testBeq0=obs.testBeq0,
           pval.testBeq0.is=pval.testBeq0.is, obs.testBeq0.is=obs.testBeq0.is,
           sim.testBeq0=sim.testBeq0,
           procBeq0=Ut,sim.test.procBeq0=UIt,
           conf.band=unifCI,
           obs.test=uni.testOBS,pval.test=uni.pval,
           sim.test=uni.test,
                                        #pval.testtime=unitime.pval,
           proc.cumz=univar.proc,sim.test.proccumz=UIz)

  attr(ud,"Call")<-sys.call(); 
  class(ud)<-"cum.residuals"
  return(ud); 
}

"print.cum.residuals"<- function (x,...)
{
    object <- x; rm(x);
	if (!inherits(object, 'cum.residuals')) stop ("Must be an MG resid object")

		cat("  Call: \n")
			dput(attr(object, "Call"))
			cat("\n")
}

"plot.cum.residuals" <- function (x,pointwise.ci=1,hw.ci=0,sim.ci=0,
		robust=1, specific.comps=FALSE,level=0.05, start.time = 0, 
	stop.time = 0, add.to.plot=FALSE, mains=TRUE, 
	xlab="Time",ylab ="Cumulative Residuals",ylim=NULL,
	score=0,...) 
{
  object <- x; rm(x);
  if (!inherits(object,'cum.residuals') ) stop ("Must be output from cum.residuals()") 

  if (score <2) {
    B<-object$cum;
    if (sum(B)==0)  {
      cat(" To compute cumulative residuals provide model matrix \n");
    } }
  if (score==2) {
    if (sum(object$obs.test)==0)
      stop("To plot cumulative residuals vs. covariates, cum.resid=1"); 
  }

  if (score==0) {
    B<-object$cum; 
    if (robust>=1) V<-object$robvar.cum else V<-object$var.cum
    p<-dim(B)[[2]]; 

    if (specific.comps==FALSE) comp<-(2:p) else comp<-specific.comps+1
    if (stop.time==0) stop.time<-max(B[,1]);

    med<-B[,1]<=stop.time & B[,1]>=start.time
    B<-B[med,]; Bs<-B[1,];  B<-t(t(B)-Bs); B[,1]<-B[,1]+Bs[1];
    V<-V[med,]; Vs<-V[1,]; V<-t( t(V)-Vs); 
    Vrob<-object$robvar.cum; 
    Vrob<-Vrob[med,]; Vrobs<-Vrob[1,]; Vrob<-t( t(Vrob)-Vrobs); 

    c.alpha<- qnorm(1-level/2)
    for (v in comp) { 
      c.alpha<- qnorm(1-level/2)
      est<-B[,v];ul<-B[,v]+c.alpha*V[,v]^.5;nl<-B[,v]-c.alpha*V[,v]^.5;
      if (add.to.plot==FALSE) 
        {
          plot(B[,1],est,ylim=1.05*range(ul,nl),type="s",xlab=xlab,ylab=ylab) 
          if (mains==TRUE) title(main=colnames(B)[v]); }
      else lines(B[,1],est,type="s"); 
      if (pointwise.ci>=1) {
        lines(B[,1],ul,lty=pointwise.ci,type="s");
        lines(B[,1],nl,lty=pointwise.ci,type="s"); }
      if (robust>=1) {
        lines(B[,1],ul,lty=robust,type="s"); 
        lines(B[,1],nl,lty=robust,type="s"); }
      if (hw.ci>=1) {
        if (level!=0.05) cat("Hall-Wellner bands only 95 % \n");
        tau<-length(B[,1])
        nl<-B[,v]-1.13*V[tau,v]^.5*(1+V[,v]/V[tau,v])
        ul<-B[,v]+1.13*V[tau,v]^.5*(1+V[,v]/V[tau,v])
        lines(B[,1],ul,lty=hw.ci,type="s"); 
        lines(B[,1],nl,lty=hw.ci,type="s"); }
      if (sim.ci>=1) {
        if (level!=0.05) c.alpha<-percen(object$sim.testBeq0[,v-1],1-level)
        else c.alpha<-object$conf.band[v-1];
        nl<-B[,v]-c.alpha*Vrob[,v]^.5; ul<-B[,v]+c.alpha*Vrob[,v]^.5;
        lines(B[,1],ul,lty=sim.ci,type="s"); 
        lines(B[,1],nl,lty=sim.ci,type="s"); }
      abline(h=0)
    }
  } else if (score==1) 
    {
                                        # plot score proces
      dim1<-ncol(object$procBeq0)
      if (sum(specific.comps)==FALSE) comp<-2:dim1 else comp<-specific.comps+1

      for (i in comp)
        {
          ranyl<-range(object$procBeq0[,i]);
          for (j in 1:50) ranyl<-range(c(ranyl,
                                         (object$sim.test.procBeq0[[j]])[,i-1]));
          mr<-max(abs(ranyl));

          if (add.to.plot==FALSE)
            plot(object$procBeq0[,1],object$procBeq0[,i],type="l",ylim=c(-mr,mr),
                 lwd=2,xlab=xlab,ylab=ylab)
          else
            lines(object$procBeq0[,1],object$procBeq0[,i])
          if (mains==TRUE) title(main=colnames(object$procBeq0)[i]);
          for (j in 1:50)
            lines(object$procBeq0[,1],as.matrix(object$sim.test.procBeq0[[j]])[,i-1],
                  col="grey",lwd=1,lty=1)
          lines(object$procBeq0[,1],object$procBeq0[,i],lwd=2)
        } 
    } else if (score==2) {
                                        # plot score proces
      dim1<-length(object$obs.test)
      if (sum(specific.comps)==FALSE) comp<-1:dim1 else comp<-specific.comps

      for (i in comp)
        {
          if (nrow(object$proc.cumz[[i]])==1) TYPE<-"p" else TYPE<-"l"; 

          if (TYPE=="l") 
            {
              ranyl<-range(object$proc.cumz[[i]][,2]);
              for (j in 1:50) ranyl<-range(c(ranyl,(object$sim.test.proccumz[[i]])[,j])); 
              mr<-max(abs(ranyl));

              if (add.to.plot==FALSE)
                plot(object$proc.cumz[[i]][,1],object$proc.cumz[[i]][,2],type=TYPE,
                     ylim=c(-mr,mr),lwd=2,xlab=colnames(object$proc.cumz[[i]])[1],
                     ylab="Cumulative residuals")
              else
                lines(object$proc.cumz[[i]][,1],object$proc.cumz[[i]][,2],type="l")
              if (mains==TRUE) title(main=colnames(object$proc.cumz[[i]])[1]); 
              if (TYPE=="l") for (j in 1:50)
                lines(object$proc.cumz[[i]][,1],object$sim.test.proccumz[[i]][,j],
                      col="grey",lwd=1,lty=1,type="l")
              if (TYPE=="p") for (j in 1:50)
                points(object$proc.cumz[[i]][,1],object$sim.test.proccumz[[i]][,j],pch=".")
              lines(object$proc.cumz[[i]][,1],object$proc.cumz[[i]][,2],lwd=2); 
            } }
    }

}

"summary.cum.residuals" <-
function (object,digits=3,...) 
{
  if (!inherits(object, 'cum.residuals')) stop ("Must be an cum.residuals object")

                                        # We print information about object:  
  cat("Test for cumulative MG-residuals \n\n")

  mtest<-(sum(object$conf.band)>0)

  if (mtest==FALSE) { 
    cat("Grouped cumulative residuals not computed, you must provide\n")
    cat("modelmatrix to get these (see help) \n\n")
  }
  if (mtest==TRUE) { 
    test0<-cbind(object$obs.testBeq0,object$pval.testBeq0)
    test0.is<-cbind(object$obs.testBeq0.is,object$pval.testBeq0.is) 
    colnames(test0)<- c("sup|  hat B(t) |","p-value H_0: B(t)=0")
    colnames(test0.is)<- c("int ( B(t) )^2 dt","p-value H_0: B(t)=0")  

    cat("Grouped Residuals consistent with model \n\n")
    prmatrix(round(test0,digits))
    cat("\n")
    prmatrix(round(test0.is,digits))
    cat("\n")
  }

  cumtest<-!is.null(object$obs.test)
  if (cumtest==FALSE) { 
    cat("Cumulative tests versus covariates not computed \n\n")
    cat("cum.resid=1 to compute these \n\n"); 
  }
  if (cumtest==TRUE) { 
    test0<-cbind(object$obs.test,object$pval.test)
    colnames(test0)<- c("sup|  hat B(t) |","p-value H_0: B(t)=0")

    cat("Residual versus covariates consistent with model \n\n")
    prmatrix(round(test0,digits))
  }
  cat("  \n");cat("  Call: \n");dput(attr(object, "Call"));
  cat("\n");
}
