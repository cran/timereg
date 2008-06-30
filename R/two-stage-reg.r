prop<-function(x) x

two.stage.reg<-function(formula=formula(data),data=sys.parent(),
beta=0,Nit=10,detail=0,start.time=0,max.time=NULL,id=NULL, 
clusters=NULL, robust=1,
rate.sim=1,beta.fixed=0,theta=NULL,theta.des=NULL)
{
  ratesim<-rate.sim; inverse<-0; 
  call <- match.call()
  m <- match.call(expand=FALSE)
  m$robust<-m$start.time<-m$beta<-m$Nit<-m$detail<-m$max.time<-m$clusters<-m$rate.sim<-m$beta.fixed<-m$theta<-m$theta.des<-NULL

  if (robust==0) cat("When robust=0 no variance estimate\n"); 

  special <- c("prop","cluster")
  Terms <- if(missing(data)) terms(formula, special)
  else          terms(formula, special, data=data)
  m$formula <- Terms
  m[[1]] <- as.name("model.frame")
  m <- eval(m, sys.parent())
  mt <- attr(m, "terms")
  intercept<-attr(mt, "intercept")
  Y <- model.extract(m, "response")
  if (!inherits(Y, "Surv")) stop("Response must be a survival object")

  des<-read.design(m,Terms,model="cox.aalen")
  X<-des$X; Z<-des$Z; npar<-des$npar; px<-des$px; pz<-des$pz;
  covnamesX<-des$covnamesX; covnamesZ<-des$covnamesZ

  if(is.null(clusters)) clusters <- des$clusters  
  
  pxz <- px + pz;

  survs<-read.surv(m,id,npar,clusters,start.time,max.time)
  times<-survs$times;id<-id.call<-survs$id.cal;
  clusters<-cluster.call<-survs$clusters; 
  time<-survs$start; time2<-survs$stop; status<-survs$status;
  ldata<-list(start=survs$start,stop=survs$stop,
              antpers=survs$antpers,antclust=survs$antclust);

  if (npar==FALSE) covar<-data.matrix(cbind(X,Z)) else 
  stop("Both multiplicative and additive model needed");

  Ntimes <- sum(status); 
  times<-c(start.time,time2[status==1]); times<-sort(times);
  if (is.null(max.time)==TRUE) maxtimes<-max(times)+0.1 else maxtimes<-max.time; 
  times<-times[times<maxtimes]

  if ((sum(beta)==0) & (beta.fixed==0)) beta<-coxph(Surv(time,time2,status)~Z)$coef; 

  if (px==0) stop("No nonparametric terms (needs one!)");
  ud<-two.stageBase.reg(times,ldata,X,Z,
                        status,id,clusters,Nit=Nit,detail=detail,beta=beta,
                        robust=robust,ratesim=ratesim,namesX=covnamesX,namesZ=covnamesZ,
                        beta.fixed=beta.fixed,theta=theta,theta.des=theta.des);

  if (px>0) {
    colnames(ud$cum)<-colnames(ud$var.cum)<- c("time",covnamesX)
    if (robust==1) colnames(ud$robvar.cum)<- c("time",covnamesX) }

  rownames(ud$gamma)<-c(covnamesZ); colnames(ud$gamma)<-"estimate"; 
  rownames(ud$score)<-c(covnamesZ); colnames(ud$score)<-"score"; 
                                        #namematrix(ud$var.gamma,covnamesZ); 
                                        #namematrix(ud$robvar.gamma,covnamesZ); 
                                        #namematrix(ud$D2linv,covnamesZ); 

                                        #colnames(ud$var.gamma)<-c(covnamesZ); 
                                        #rownames(ud$var.gamma)<-c(covnamesZ); 
                                        #colnames(ud$robvar.gamma)<-c(covnamesZ); 
                                        #rownames(ud$robvar.gamma)<-c(covnamesZ); 
                                        #colnames(ud$D2linv)<-c(covnamesZ); 
                                        #rownames(ud$D2linv)<-c(covnamesZ); 

  ptheta<-length(ud$theta); 
  if (ptheta>1) rownames(ud$theta)<-colnames(theta.des) else rownames(ud$theta)<-"intercept"

  attr(ud,"Call")<-sys.call(); 
  class(ud)<-"two.stage.reg"
  attr(ud,"Formula")<-formula;
  attr(ud,"id")<-id.call;
  attr(ud,"cluster")<-cluster.call;
  attr(ud,"start")<-start.time; 
  attr(ud,"time2")<-time2; 
  attr(ud,"inverse")<-inverse

  return(ud); 
}

two.stageBase.reg<-function (times, fdata, designX, designG, status,
id, clusters, Nit = 5, beta = 0, detail = 0, robust = 1, 
ratesim = 1, namesZ=NULL,namesX=NULL,beta.fixed=0,theta=NULL,
theta.des=NULL,inverse=0) 
{
    additive.resamp <-0; ridge <- 0; XligZ <- 0;
    Ntimes <- length(times)
    designX <- as.matrix(designX); designG <- as.matrix(designG)
    if (is.matrix(designX) == TRUE) px <- as.integer(dim(designX)[2])
    if (is.matrix(designX) == TRUE) nx <- as.integer(dim(designX)[1])
    if (is.matrix(designG) == TRUE) pg <- as.integer(dim(designG)[2])
    if (is.matrix(designG) == TRUE) ng <- as.integer(dim(designG)[1])
    if (nx != ng) print(" A design og B designs er ikke ens\n")

    cumint <- matrix(0, Ntimes, px + 1)
    vcum <- matrix(0, Ntimes, px + 1)
    Rvcu <- matrix(0, Ntimes, px + 1)
    if (sum(abs(beta)) == 0) betaS <- rep(0, pg) else betaS <- beta
    score <- betaS
    Varbeta <- matrix(0, pg, pg)
    Iinv <- matrix(0, pg, pg)
    RVarbeta <- matrix(0, pg, pg)
    if (is.null(theta.des)==TRUE) ptheta<-1; 
    if (is.null(theta.des)==TRUE) theta.des<-matrix(1,ng,ptheta) else
    theta.des<-as.matrix(theta.des); 
    ptheta<-ncol(theta.des); 
    if (is.null(theta)==TRUE) theta<-rep(0.1,ptheta); 
    if (length(theta)!=ptheta) theta<-rep(theta[1],ptheta); 
    theta.score<-rep(0,ptheta);Stheta<-var.theta<-matrix(0,ptheta,ptheta); 

    cluster.size<-as.vector(table(clusters));

    #dyn.load("two-stage-reg.so"); 

    nparout <- .C("twostagereg", 
        as.double(times), as.integer(Ntimes), 
        as.double(designX), as.integer(nx), as.integer(px), 
	as.double(designG), as.integer(ng), as.integer(pg), 
	as.integer(fdata$antpers),as.double(fdata$start),as.double(fdata$stop),
	as.double(betaS), as.integer(Nit), as.double(cumint), 
	as.double(vcum),  as.double(Iinv), as.double(Varbeta), 
	as.integer(detail), as.double(Rvcu), as.double(RVarbeta), 
         as.integer(id), as.integer(status), as.integer(ratesim), 
	as.double(score), as.integer(robust), as.integer(clusters),
        as.integer(fdata$antclust), as.integer(beta.fixed),
        as.double(theta),as.double(var.theta),as.double(theta.score),
        as.integer(inverse), as.integer(cluster.size), as.double(theta.des),
        as.integer(ptheta), as.double(Stheta),PACKAGE = "timereg")

    gamma <- matrix(nparout[[12]], pg, 1)
    cumint <- matrix(nparout[[14]], Ntimes, px + 1)
    vcum <- matrix(nparout[[15]], Ntimes, px + 1)
    Iinv <- matrix(nparout[[16]], pg, pg)
    Varbeta <- -matrix(nparout[[17]], pg, pg)
    Rvcu <- matrix(nparout[[19]], Ntimes, px + 1)
    RVarbeta <- -matrix(nparout[[20]], pg, pg)
    score <- matrix(nparout[[24]], pg, 1)


   theta<-matrix(nparout[[29]],ptheta,1);  
   var.theta<-matrix(nparout[[30]],ptheta,ptheta); 
   theta.score<-nparout[[31]]; 
   Stheta<-matrix(nparout[[32]],ptheta,ptheta); 

   ud <- list(cum = cumint, var.cum = vcum, robvar.cum = Rvcu, 
       gamma = gamma, var.gamma = Varbeta, robvar.gamma = RVarbeta, 
       D2linv = Iinv, score = score,  theta=theta,var.theta=var.theta,
       S.theta=Stheta)
   return(ud)
}

summary.two.stage.reg <- function(object,...){
  summary.two.stage(object,...)
}

print.two.stage.reg <- function(x,...){
  print.two.stage(x,...)
}

coef.two.stage.reg <- function(object,...){
  coef.two.stage(object,...)
}

plot.two.stage.reg <- function(x,...){
  plot.two.stage(x,...)
}

