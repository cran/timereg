\name{predict.comprisk}
\alias{predict.aalen}
\alias{predict.comprisk}
\alias{predict.cox.aalen}
\title{Predictions for Survival and Competings Risks Regression} 
\description{
Make predictions based on the survival models 
(Aalen and Cox-Aalen) and the competing risks
models for the cumulative incidence function
(comp.risk). Computes confidence intervals and
confidence bands based on resampling. 
}
\usage{
\method{predict}{comprisk}(object,newdata=NULL,X=NULL,Z=NULL,
n.sim=500, uniform=TRUE,se=TRUE,alpha=0.05,...)
}
\arguments{
  \item{object}{an object belonging to one of the following
    classes: comprisk, aalen or cox.aalen}
  \item{newdata}{specifies the data at which the
    predictions are wanted.}
  \item{X}{alternative to newdata, specifies the
    nonparametric components for predictions.}
  \item{Z}{alternative to newdata, specifies the
    parametric components of the model for predictions.}
  \item{n.sim}{number of simulations in resampling.}
  \item{uniform}{computes resampling based uniform
    confidence bands.}
  \item{se}{computes pointwise standard errors}
  \item{alpha}{specificies the significance levelwhich cause we consider.}
  \item{...}{unused arguments - for S3 compatability}
}
\references{
Scheike, Zhang and Gerds (2007), Predicting cumulative
incidence probability by direct binomial regression,
Biometrika, to appear. 

Scheike and Zhang (2007), Flexible competing risks regression modelling and goodness of fit, work in progress. 

Martinussen and Scheike (2006), Dynamic regression
models for survival data, Springer. 
}
\author{Thomas Scheike, Jeremy Silver}
\examples{
data(bmt); 
times<-bmt$time[bmt$cause==1]; 

add<-comp.risk(Surv(time,cause>0)~platelet+age+tcell,bmt,
bmt$cause,times[-1],causeS=1,resample.iid=1)
summary(add)

par(mfrow=c(2,4))
plot(add); plot(add,score=1)

ndata<-data.frame(platelet=c(1,0,0),age=c(0,1,0),tcell=c(0,0,1))
par(mfrow=c(2,3))
out<-predict(add,ndata,uniform=1,n.sim=1000)
par(mfrow=c(2,2))
plot(out,multiple=0,uniform=1,col=1:3,lty=1,se=1)
# see comp.risk for further examples. 

## SURVIVAL predictions aalen function
data(sTRACE)
out<-aalen(Surv(time,status==9)~const(age)+const(sex)+
const(diabetes)+chf+vf,
sTRACE,max.time=7,n.sim=0,resample.iid=1)

out<-predict(out,X=rbind(c(1,0,0),c(1,1,0)),
Z=rbind(c(55,0,1),c(60,1,1)))
par(mfrow=c(2,2))
plot(out,multiple=1,se=0,uniform=0,col=1:2,lty=1:2)
plot(out,multiple=0,se=1,uniform=1,col=1:2)

data(sTRACE)
out<-cox.aalen(Surv(time,status==9)~prop(age)+prop(sex)+
prop(diabetes)+chf+vf,
sTRACE,max.time=7,n.sim=0,resample.iid=1)

out<-predict(out,X=rbind(c(1,0,0),c(1,1,0)),
Z=rbind(c(55,0,1),c(60,1,1)))
par(mfrow=c(2,2))
plot(out,multiple=1,se=0,uniform=0,col=1:2,lty=1:2)
plot(out,multiple=0,se=1,uniform=1,col=1:2)
}
\keyword{survival}
