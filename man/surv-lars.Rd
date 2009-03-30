\name{surv.lars}
\alias{mylars}
\alias{surv.lars}
\alias{surv.lars.cv}
\alias{mycoef.lars}
\alias{mypredict.lars}
\title{Fits LASSO model for  additive hazards model by Lars
algorithm}
\description{
Fits the LASSO estimator for the additive risk
model based on the least squares fitting criterion

\deqn{
L(\beta,D,d) = \beta^T D \beta - 2 \beta^T d
}
where \eqn{D=\int Z H Z dt} and  \eqn{d=\int Z H dN}.

This is equivalent to an appropriate normal least squares
problem on least squares data. 
}
\usage{
surv.lars(time,status,z,l1.weights=NULL,...)
}
\arguments{
\item{time}{covariates for artificial least squares data.}
\item{status}{responses for artificial least squares data.}
\item{z}{covariatesspecifies the weights for the L1 penalty.}
\item{l1.weights}{specifies the weights for the L1 penalty.}
\item{...}{unused arguments - for S3 compatibility}
}
\details{
Modified version of standard lars program in the LARS package.
Essentially all scalings are removed from the lars program. 

The matrices D and d can be computed using the aalen.test() function, but
we here rely on an artificial least squares representation of the data. 
This approach takes up a lot of memory and building the algorithm directly
on D and d should improve the speed considerably.
}
\value{
see lars. 
}
\references{
Martinussen and Scheike, Model selection for the the additive
risks hazards model, submitted. 

Martinussen and Scheike, Dynamic Regression Models for
Survival Data, Springer (2006).

LARS, Efron et al. 
}
\author{Thomas Scheike}
\examples{
## makes data for pbc complete case
data(mypbc)
pbc<-mypbc
pbc$time<-pbc$time+runif(418)*0.1; pbc$time<-pbc$time/365
pbc<-subset(pbc,complete.cases(pbc)); 
covs<-as.matrix(pbc[,-c(1:3,6)])
covs<-cbind(covs[,c(1:6,16)],log(covs[,7:15]))

## computes the matrices needed for the least squares 
## criterion 
#out<-des.aalen(Surv(time,status>=1)~const(covs),pbc)
#x<-as.matrix(out[,1:17]); y<-out$Y

## lasso for survival data 
# must load lars library for it to work 
#fit<-surv.lars(pbc$time,pbc$status,covs); 
#fit<-mylars(x,y)
#cv<-surv.lars.cv(time,status,covs)
#myplot.lars(fit)

#cv.frac<-cv$cv.frac
#beta<-mycoef.lars(fit,s=cv.frac,mode="fraction");

#full<-fit$beta[18,] # least squares solution

## adaptive lasso for survival data 
#fit<-surv.lars(time,status,covs,l1.weights=1/abs(full)); 
#cva<-surv.lars.cv(time,status,covs,l1.weights=1/abs(full))

#cva.frac<-cva$cv.frac
#abeta<-mycoef.lars(afit,s=cv.frac,mode="fraction");
#myplot.lars(afit)
}
\keyword{survival}
