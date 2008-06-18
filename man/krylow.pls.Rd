\name{krylow.pls}
\alias{krylow.pls}
\title{Fits Krylow based PLS for  additive hazards model}
\description{
Fits the PLS estimator for the additive risk
model based on the least squares fitting criterion

\deqn{
L(\beta,D,d) = \beta^T D \beta - 2 \beta^T d
}
where \eqn{D=\int Z H Z dt} and  \eqn{d=\int Z H dN}.
}
\usage{
krylow.pls(D,d,dim)
}
\arguments{
\item{D}{defined above}
\item{d}{defined above}
\item{dim}{number of pls dimensions}
}
\value{
returns a list with the following arguments:
\item{beta}{PLS regression coefficients}
}
\references{
Martinussen and Scheike, The Aalen additive hazards model 
with high-dimensional regressors, submitted. 

Martinussen and Scheike, Dynamic Regression Models for
Survival Data, Springer (2006).
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
out<-aalen.test(Surv(time,status)~const(lcovs),pbc,robust=0,n.sim=0)
S=out$intZHZ; s=out$intZHdN;

out<-krylow.pls(S,s,2)
}
\keyword{survival}
