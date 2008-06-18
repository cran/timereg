additive.plsR<-function (formula = formula(data),
                         data = sys.parent(), start.time = 0, max.time = NULL, 
                         robust=0, id=NULL, clusters=NULL, residuals = 0, n.sim = 0,  
                         weighted.test=0,
                         covariance=0,resample.iid=0,deltaweight=1,pls.dim=1,
                         scale=FALSE,semi.pls=1){
###deltaweight<-1; # always default
  if (n.sim == 0) sim <- 0 else sim <- 1
  if (resample.iid==1 & robust==0) {
    cat("When robust=0 no iid representaion computed\n"); 
    resample.iid<-0;}
  if (covariance==1 & robust==0) {
    cat("When robust=0 no covariance computed \n"); 
    cat("Covariance based on robust iid representation\n")
    covariance<-0;}
  if (sim==1 & robust==0) {
    cat("When robust=0, No simulations \n"); 
    cat("n.sim set to 0\n"); n.sim<-0;}
  if (residuals==1 & robust==0) {
    cat("When robust=0, no martingale residuals \n"); 
    residuals<-0;}
  if (n.sim>0 & n.sim<50) {n.sim<-50 ; cat("Minimum 50 simulations\n");}
  call <- match.call()
  m <- match.call(expand = FALSE)
  m$start.time <- m$weighted.test <- m$max.time <- m$robust <- 
    m$sim <- m$residuals <- m$n.sim <- m$id <- m$covariance <- 
      m$resample.iid <- m$clusters <- m$deltaweight<-
        m$pls.dim<- m$semipls<- NULL
  special <- c("const","cluster")
  Terms <- if (missing(data)) terms(formula, special)
  else terms(formula, special, data = data)
  m$formula <- Terms
  m[[1]] <- as.name("model.frame")
  m <- eval(m, sys.parent())
  mt <- attr(m, "terms")
  intercept <- attr(mt, "intercept")
  Y <- model.extract(m, "response")
  if (!inherits(Y, "Surv")) stop("Response must be a survival object")

  des<-read.design(m,Terms)
  X<-des$X; Z<-des$Z; npar<-des$npar; px<-des$px; pz<-des$pz;
  covnamesX<-des$covnamesX; covnamesZ<-des$covnamesZ

  if(is.null(clusters)) clusters <- des$clusters  
  
  pxz <- px + pz; 
  npar.pls<-npar; 

  survs<-read.surv(m,id,FALSE,clusters,start.time,max.time)
  times<-survs$times;id<-id.call<-survs$id.cal;
  clusters<-cluster.call<-survs$clusters; time2<-survs$stop
  status<-survs$status; 
  ldata<-list(start=survs$start,stop=survs$stop,
              antpers=survs$antpers,antclust=survs$antclust);

  px<-px-1; 
  beta<-matrix(0,pls.dim,px); pls.comp<-matrix(0,nrow(X),pls.dim); 

  Yorig<-as.matrix(X[,1]); 
  if (npar==FALSE) Z<-scale(Z,scale=scale); X<-scale(X,scale=scale); 
  pls.cov<-as.matrix(X[,-1]); Zorig<-Z; 
  Xorig<-as.matrix(X[,-1]); 


  pvals<-var.gamma<-c()
  for (j in 1:pls.dim)
    {
      cat(paste("PLS component ",j,"\n")) 
      for (i in 1:px)
        {
          if (npar==TRUE ) {
            if (j==1) Z<-as.matrix(Xorig[,i]) else Z<-
              cbind(Xorig[,i],pls.comp[,1:(j-1)]); }
          if (npar==FALSE) {
            if (j==1) Z<-cbind(Xorig[,i],Zorig) else 
            Z<-cbind(Xorig[,i],pls.comp[,1:(j-1)],Zorig); }
                                        #print(Z[1:10,]); 

          print(cbind(Yorig,Z)[1:10,])

                                        #cat("Semiparametric Additive Risk Model\n")
          ud <- semiaalen(times, ldata, Yorig , Z, 
                          status, id , clusters, robust = robust, sim = sim, antsim = n.sim, 
                          weighted.test = weighted.test, retur =
                          residuals,covariance=covariance,
                          resample.iid=resample.iid,namesX=covnamesX,namesZ=covnamesZ,
                          deltaweight=deltaweight)

          pvals<-c(pvals,1-pchisq(ud$gamma[1]^2/ud$var.gamma,1))
          var.gamma<-c(var.gamma,ud$var.gamma)

          beta[j,i]<-c(ud$gamma)[1];
          pls.comp[,j]<-pls.comp[,j]+beta[j,i]*Xorig[,i]; 
        } 
    }
  cnames<-c()
  for (j in 1:pls.dim) cnames<-c(cnames,paste("PLS-",j)); 
  colnames(pls.comp)<-cnames

  if (npar==TRUE) Z<-pls.comp else Z<-cbind(pls.comp,Zorig)
  udf <- semiaalen(times, ldata, Yorig , Z, 
                   status, id , clusters, robust = robust, sim = sim, antsim = n.sim, 
                   weighted.test = weighted.test, retur =
                   residuals,covariance=covariance,
                   resample.iid=resample.iid,namesX=covnamesX,namesZ=covnamesZ,
                   deltaweight=deltaweight)

                                        #print(udf$gamma); 

  ud<-list(pls.comp=pls.comp,beta.pls=beta,
           beta=udf$gamma,baseline=udf$cum,pvals=pvals,
           var.gamma=var.gamma)
}

additive.pls<-function (formula = formula(data), data = sys.parent(), 
start.time=0,max.time=NULL,id=NULL,
pls.dim=1, scale=FALSE,weighted.pls=0,constant=1) 
{ semipls=1; deltaweight<-1; # always !  
call <- match.call()
    m <- match.call(expand = FALSE)
    m$start.time <- m$max.time <-  m$id <- m$pls.dim<- m$semipls<- m$scale<-m$weighted.pls<-m$constant<-NULL
    special <- c("const")
    Terms <- if (missing(data)) terms(formula, special)
    else terms(formula, special, data = data)
    m$formula <- Terms
    m[[1]] <- as.name("model.frame")
    m <- eval(m, sys.parent())
    mt <- attr(m, "terms")
    intercept <- attr(mt, "intercept")
    Y <- model.extract(m, "response")
    if (!inherits(Y, "Surv")) stop("Response must be a survival object")


des<-read.design(m,Terms)
X<-des$X; Z<-des$Z; npar<-des$npar; px<-des$px; pz<-des$pz;
covnamesX<-des$covnamesX; covnamesZ<-des$covnamesZ
pxz <- px + pz; 

clusters=NULL;
survs<-read.surv(m,id,FALSE,clusters,start.time,max.time)
times<-survs$times;id<-id.call<-survs$id.cal;
clusters<-cluster.call<-survs$clusters; time2<-survs$stop
status<-survs$status; 
ldata<-list(start=survs$start,stop=survs$stop,
            antpers=survs$antpers,antclust=survs$antclust);

px<-px-1; 
beta<-matrix(0,pls.dim,px); weights<-beta 
pls.comp<-matrix(0,nrow(X),pls.dim); 

Yorig<-as.matrix(X[,1]); 
if (npar==FALSE) Z<-scale(Z,scale=scale); 
X<-scale(X,scale=scale); pls.cov<-as.matrix(X[,-1]); Zorig<-Z;  


print(npar); 

for (j in 1:pls.dim)
{
if (j<=pls.dim) cat(paste("PLS component ",j,"\n"))
 
Xmat<-as.matrix(Yorig); 
if (npar==TRUE ) {
if (j==1) Z<-as.matrix(pls.cov[,1]) else  {
if (constant==1)  {Z<-cbind(pls.comp[,1:(j-1)],pls.cov[,1]);} else
                  {Z<-as.matrix(pls.cov[,1]); 
                   Xmat<-as.matrix(cbind(Yorig,pls.comp[,1:(j-1)])); 
		   } } }
if (npar==FALSE) {
if (j==1) Z<-cbind(Zorig,pls.cov[,1]) else  {

if (constant==1)  {Z<-cbind(Zorig,pls.comp[,1:(j-1)],pls.cov[,1]);} else
                  {Z<-as.matrix(cbind(Zorig,pls.cov[,1])); 
                   Xmat<-as.matrix(cbind(Yorig,pls.comp[,1:(j-1)]));}}}

Z<-as.matrix(Z); Xmat<-as.matrix(Xmat); 
pg<-ncol(as.matrix(Z)); px<-ncol(as.matrix(Xmat)); nx<-nrow(Yorig); 

dimplscov<-ncol(pls.cov);
betapls<-rep(0,dimplscov); plscomp<-Yorig[,1]; 

Nalltimes <- length(times);
Ntimes<-sum(status[(ldata$stop>times[1]) & 
            (ldata$stop<=times[Nalltimes])])+1;
cum<-matrix(0,Ntimes,px);  

#print(Z[1:10,]); print(Yorig[1:10,]); print(pls.cov[1:10,]); 

semiout<-.C("plssemiadd",
as.double(times),as.integer(Nalltimes),as.integer(Ntimes),
as.double(Xmat),as.integer(nx),as.integer(px),
as.double(Z),as.integer(nx),as.integer(pg),
as.integer(ldata$antpers),as.double(ldata$start),as.double(ldata$stop),
as.integer(id),as.double(cum),as.integer(status),
as.integer(deltaweight),as.double(pls.cov),
as.integer(dimplscov),as.double(betapls),as.double(plscomp),
as.integer(semipls),as.integer(weighted.pls)
,PACKAGE="timereg")

plscomp<-semiout[[20]]; betapls<-semiout[[19]];

beta[j,]<-betapls; 
#if (weighted.pls==1) weights[j,]<-intZHZ;
pls.comp[,j]<-plscomp;
}
cnames<-c()
for (j in 1:pls.dim) cnames<-c(cnames,paste("PLS-",j)); 
#colnames(pls.comp)<-cnames

colnames(beta)<-covnamesX[-1]; 

if (npar==TRUE) Z<-pls.comp else Z<-cbind(pls.comp,Zorig)

udf<-semiaalen(times, ldata, Yorig , Z, 
    status, id , clusters, robust = 0, sim = 0, 
    antsim = 0, weighted.test = 0 , retur = 0 ,
    covariance= 0 , resample.iid=0 ,
    namesX=covnamesX,namesZ=covnamesZ, 
    deltaweight=deltaweight);

betaf<-c(udf$gamma)
tbeta.pls<-apply(beta*betaf[1:pls.dim],2,sum)

ud<-list(pls.comp=pls.comp,beta.pls=beta,beta=betaf,
baseline=udf$cum,tbeta.pls=tbeta.pls)
}
