\name{comp.risk}
\alias{comp.risk}
\title{Competings Risks Regression} 
\description{
Fits a semiparametric model for the cause-specific quantities :
\deqn{
P(T \leq t, cause=1 | x) = P_1(t,x) = 1 - \exp(- g(t,x) )
}
for the probability of dying from cause 1 in a situation with 
competing causes of death. 

The model is considered in two situations : 
\deqn{
g(t,x,z) = x^T A(t) + (diag(t^p) z)^T \beta
}
the additive situation (with additive subdistribution hazard) and 
the proportional setting that includes the Fine & Gray (FG) model and
some extensions
\deqn{
g(t,x,z) = \exp(x^T A(t) + (diag(t^p) z)^T \beta)
}
The FG model is obtained when \eqn{x=1}. Where p is 1 for the additive
model and 0 for the proportional model. In general p may be
powers of the same length as z. 
}
\usage{
comp.risk(formula,data=sys.parent(),cause,times,Nit=50,
clusters=NULL,gamma=0,n.sim=500,weighted=0,model="additive",
causeS=1,cens.code=0,detail=0,interval=0.01,resample.iid=1,
cens.model="KM",time.pow=0)
}
\arguments{
\item{formula}{a formula object, with the response on the left of a '~'
  operator, and the terms on the right. The response must be a survival
  object as returned by the `Surv' function. The status indicator is not
  important here. Time-invariant regressors are specified by the wrapper
  const(), and cluster variables (for computing robust variances) by the
  wrapper cluster().
}
\item{data}{a data.frame with the variables.}
\item{cause}{specifies the causes  related to the death
	times, the value 0 is the censoring value.}
\item{times}{specifies the times at which the estimator is
	considered. This is typically all cause "1" jump times.}
\item{Nit}{number of iterations for Newton-Raphson algorithm.}
\item{clusters}{specifies cluster structure, for backwards compability.}
\item{gamma}{starting value for constant effects.}
\item{n.sim}{number of simulations in resampling.}
\item{weighted}{Not implemented. To compute a variance weighted version of the test-processes used for testing time-varying effects.}
\item{model}{"additive" or "prop"ortional.}
\item{causeS}{specificies which cause we consider.}
\item{cens.code}{specificies the code for the censoring.}
\item{detail}{if 0 no details are printed during iterations, if 1 details are given.}
\item{interval}{specifies that we only consider timepoints where the Kaplan-Meier of the censoring distribution  is larger than this value.}
\item{resample.iid}{to return the iid decomposition, that can be used to construct confidence bands for predictions}
\item{cens.model}{specified which model to use for the ICPW, KM is Kaplan-Meier alternatively it may be "cox"}
\item{time.pow}{specifies that the power at which the time-arguments  is transformed, for each of the arguments of the const() terms, default is 1 for the additive model and 0 for the proportional model.}
}
\value{returns an object of type 'comprisk'. With the following arguments:
\item{cum}{cumulative timevarying regression coefficient estimates are 
computed within the estimation interval.}
\item{var.cum}{pointwise variances estimates.  }
\item{gamma}{estimate of proportional odds parameters of model.}
\item{var.gamma}{variance for gamma.  }
\item{score}{sum of absolute value of scores.}
\item{gamma2}{estimate of constant effects based on the
non-parametric estimate. Used for testing of constant effects.}
\item{obs.testBeq0}{observed absolute value of supremum of 
cumulative components scaled with the variance.}
\item{pval.testBeq0}{p-value for covariate effects based on supremum test.}
\item{obs.testBeqC}{observed absolute value of supremum of difference between observed cumulative process and estimate under null of constant effect.}
\item{pval.testBeqC}{p-value based on resampling.}
\item{obs.testBeqC.is}{observed integrated squared differences between 
observed cumulative and estimate under null of constant effect.}
\item{pval.testBeqC.is}{p-value based on resampling.}
\item{conf.band}{resampling based constant to construct 95\% uniform confidence bands.}
\item{B.iid}{list of iid decomposition of non-parametric effects.}
\item{gamma.iid}{matrix of iid decomposition of parametric effects.}
\item{test.procBeqC}{observed test process for testing of time-varying effects}
\item{sim.test.procBeqC}{50 resample processes for for testing of time-varying effects}
}
\references{
Scheike, Zhang and Gerds (2007), Predicting cumulative
incidence probability by direct binomial regression,
Biometrika, to appear. 

Scheike and Zhang (2007), Flexible competing risks regression
modelling and goodness of fit, work in progress. 
}
\author{Thomas Scheike}
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
out<-predict(add,ndata,uniform=1,n.sim=100)
par(mfrow=c(2,2))
plot(out,multiple=0,uniform=1,col=1:3,lty=1,se=1)

## fits additive model with some constant effects 
add.sem<-comp.risk(Surv(time,cause>0)~
const(platelet)+const(age)+const(tcell),bmt,
bmt$cause,times[-1],causeS=1,resample.iid=1)
summary(add.sem)

out<-predict(add.sem,ndata,uniform=1,n.sim=100)
par(mfrow=c(2,2))
plot(out,multiple=0,uniform=1,col=1:3,lty=1,se=0)

## Fine & Gray model 
fg<-comp.risk(Surv(time,cause>0)~
const(platelet)+const(age)+const(tcell),bmt,
bmt$cause,times[-1],causeS=1,resample.iid=1,model="prop")
summary(fg)

out<-predict(fg,ndata,uniform=1,n.sim=100)

par(mfrow=c(2,2))
plot(out,multiple=1,uniform=0,col=1:3,lty=1,se=0)

## extended model with time-varying effects
fg.npar<-comp.risk(Surv(time,cause>0)~platelet+age+const(tcell),
bmt,bmt$cause,times[-1],causeS=1,resample.iid=1,model="prop")
summary(fg.npar); 

out<-predict(fg.npar,ndata,uniform=1,n.sim=100)

par(mfrow=c(2,2))
plot(out,multiple=1,uniform=0,col=1:3,lty=1,se=0)
}
\keyword{survival}