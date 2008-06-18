\name{lasso.boost}
\alias{lasso.boost}
\alias{lasso.boostR}
\alias{L}
\alias{dL}
\alias{lasso.add.hazard}
\title{Fits LASSO model for  additive hazards model}
\description{
Fits the LASSO estimator for the additive risk
model based on the least squares fitting criterion

\deqn{
L(\beta,D,d) = \beta^T D \beta - 2 \beta^T d
}
where \eqn{D=\int Z H Z dt} and  \eqn{d=\int Z H dN}.
}
\usage{
lasso.boost(D,d,lambda,max.it=10000,beta=0,detail=0)
}
\arguments{
\item{D}{defined above}
\item{d}{defined above}
\item{lambda}{l1 regularization}
\item{max.it}{number of steps in l1 boosting algorithm}
\item{detail}{prints details}
\item{beta}{starting value for algorithm.}
}
\details{
lasso.boost       is the boosting algorithm 

lasso.add.hazard  computes the exact solution using the 
                  quadprog package. 
}
\value{
returns a list with the following arguments:
\item{beta}{regression coefficients}
\item{L}{value of the fitting criterion}
\item{l1}{sum of the absolute value of the coeffcients}
}
\references{
Martinussen and Scheike, The Aalen additive hazards model 
with high-dimensional regressors, submitted. 

Martinussen and Scheike, Dynamic Regression Models for
Survival Data, Springer (2006).

Kim and Kim (2004), Gradient LASSO for feature selection,
In  Proceedings of the 21st International Conference on Machine
Learning.
}
\author{Thomas Scheike}
\examples{
## makes data for pbc complete case
data(pbc)
pbc$time<-pbc$time+runif(418); pbc$time<-pbc$time/365
pbc[pbc==-9]<-NA; rsum<-apply(pbc,1,sum); 
pbc<-pbc[!is.na(rsum),]
covs<-as.matrix(pbc[,-c(10,17)])
lcovs<-covs;
lcovs[,c(3,5,6,10,11,13,17,18)]<-
log(lcovs[,c(3,5,6,10,11,13,17,18)])

## computes the matrices needed for the least squares 
## criterion 
out2<-aalen.test(Surv(time,status)~const(lcovs),pbc,robust=0,n.sim=0)
S=out2$intZHZ; s=out2$intZHdN;

## lambda=0.39
out.pbc<-lasso.boost(S,s,0.39,max.it=20000)

## exact solution, slow sometimes !!
# library(quadprog)
# out.ex<-lasso.add.hazard(S,s,0.39,0,1,max.it=20000)
# print(round(cbind(out.pbc$beta,out.ex$sol$solution),4))
}
\keyword{survival}
