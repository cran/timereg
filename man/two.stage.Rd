\name{two.stage}
\alias{two.stage}
\title{Fit Clayton-Oakes-Glidden Two-Stage model}
\description{
Fit Clayton-Oakes-Glidden Two-Stage model with 
Cox-Aalen marginals and Gamma distributed dependence frailties

Fits an Cox-Aalen survival model as the marginal model. 
Time dependent variables and counting process data
(multiple events per subject) are not possible !

The marginal baselines are on the Cox-Aalen form
\deqn{ \lambda_{i}(t) = Y_i(t) 
( X_{i}^T(t) \alpha(t) )  \exp(Z_{i}^T \beta ) 
}

The model thus contains the Cox's regression model and the 
additive hazards model as special cases. (see cox.aalen function
for more on this).

The modelling formula uses the standard survival modelling given in the
\bold{survival} package.
}
\usage{
two.stage(formula=formula(data),data=sys.parent(),beta=0,Nit=10,
detail=0, start.time=0,max.time=NULL,id=NULL,clusters=NULL,
robust=1,rate.sim=1,beta.fixed=0,theta=NULL,inverse=0)
}
\arguments{
\item{formula}{
  a formula object with the response on the left of a '~' operator, and
  the independent terms on the right as regressors. The response must be
  a survival object as returned by the `Surv' function. Terms with a
  proportional effect are specified by the wrapper prop(), and cluster
  variables (for computing robust variances) by the wrapper cluster().
}
\item{data}{a data.frame with the variables.}
\item{start.time}{start of observation period where estimates are computed.}
\item{max.time}{end of observation period where estimates are computed. 
Estimates thus computed from [start.time, max.time]. Default is max of data.}
\item{robust}{to compute robust variances and construct processes for
resampling. May be set to 0 to save memory.}
\item{id}{For timevarying covariates the variable must 
associate each record with the id of a subject.}
\item{clusters}{cluster variable for computation of robust variances.}
\item{beta}{starting value for relative risk estimates}
\item{Nit}{number of iterations for Newton-Raphson algorithm.}
\item{detail}{if 0 no details is printed during iterations, if 1 details are given.}
\item{rate.sim}{rate.sim=1 such that resampling of residuals is based on estimated martingales and thus valid in rate case, rate.sim=0 means that resampling is based on counting processes and thus only valid in intensity case.}
\item{beta.fixed}{option for keeping beta in the Cox-Aalen model fixed.}
\item{theta}{starting values for the frailty variance (default=0.1).}
\item{inverse}{inverse=0, indicates that the parameter is the variance, in contrast to inverse=1 where the inverse variance is used as the parameter. }
}
\details{
The data for a subject is presented as multiple rows or 'observations',
each of which applies to an interval of observation (start, stop]. 
For counting process data with the )start,stop] notation is used
the 'id' variable is needed to identify the records for each subject.
Only one record per subject is allowed in the current
implementation for the estimation of theta. 
The program assumes that there are no ties, and if such are present 
random noise is added to break the ties.
}
\value{
returns an object of type "two.stage". With the following arguments:
\item{cum}{cumulative timevarying regression coefficient estimates are 
computed within the estimation interval.
}
\item{var.cum}{the martingale based pointwise variance estimates.  }
\item{robvar.cum}{robust pointwise variances estimates.  }
\item{gamma}{estimate of parametric components of model.  }
\item{var.gamma}{variance for gamma.  }
\item{robvar.gamma}{robust variance for gamma.  }
\item{D2linv}{inverse of the derivative of the score function.}
\item{score}{value of score for final estimates.}
\item{theta}{estimate of Gamma variance for frailty.}
\item{var.theta}{estimate of variance of theta.}
}
\references{
Glidden (2000), A Two-Stage estimator of the dependence parameter
for the Clayton Oakes model.

Martinussen and Scheike, Dynamic Regression Models for Survival Data,
Springer (2006).
}
\author{Thomas Scheike}
\examples{
library(survival)
data(diabetes)
# Marginal Cox model  with treat as covariate

fit<-two.stage(Surv(time,status) ~ prop(treat) + cluster(id),
diabetes,theta=0.95,detail=0,inverse=0)
summary(fit)

fit<-two.stage(Surv(time,status) ~ prop(treat) + cluster(id),
diabetes,theta=0.95,detail=0,inverse=1)
summary(fit)

# to fit model without covariates, fix beta=0, but still need prop term !
fit<-two.stage(Surv(time,status) ~ prop(treat) + cluster(id),
diabetes,theta=0.95,detail=0,beta.fixed=1,inverse=0)
summary(fit)
}
\keyword{survival}
