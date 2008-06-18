\name{additive.pls}
\alias{additive.pls}
\alias{additive.plsR}
\title{Fits PLS for additive hazards model}
\description{
Fits the partial least squares estimator for the additive risk
model.
Time dependent variables and counting process data 
(multiple events per subject) are possible.

The modelling formula uses the standard survival modelling given in the 
\bold{survival} package.

Covariates \eqn{Z_1,...,Z_p}, fixed covariates \eqn{F(t)}. 
Algorithm : 
1) 
For pls components \eqn{X_1,...,X_K},  fits 
\deqn{
\lambda_0(t) + \alpha^T(t) F(t) + 
\sum_{j=1}^K X_j(t) \gamma_j(t) + Z_i(t) \beta_i^{K+1}
}
for \eqn{i=1,...,p}
2) compute new pls components \eqn{X_{K+1} = \sum \beta_i^{K+1} Z_i(t)}
   and iterate. 
}
\usage{
additive.pls(formula = formula(data), data =
sys.parent(), start.time=0,max.time=NULL,id=NULL, pls.dim=1, 
scale=FALSE, weighted.pls=0, constant=1)
}
\arguments{
\item{formula}{
a formula object with the response on the left of a '~' 
operator, and the independent terms on the right as regressors. 
The response must be a 
survival object as returned by the `Surv' function.

The const terms are kept as fixed covariates that are not involved in the 
pls variable reduction. This may covariates that are know to be of clinical 
importance. 
}
\item{data}{a data.frame with the variables.}
\item{start.time}{start of observation period where estimates are computed.}
\item{max.time}{end of observation period where estimates are computed. 
Estimates thus computed from [start.time, max.time]. Default is max of data.}
\item{id}{For timevarying covariates the variable must 
associate each record with the id of a subject.}
\item{pls.dim}{number of pls components}
\item{scale}{to center and scale the covariates}
\item{weighted.pls}{Gartwaith weights if 1.}
\item{constant}{When 0 the coefficient \eqn{\gamma_j(t)} does not 
depend on time.}
}
\details{
const() specifies the \eqn{F(t)} covariates in the above model. 
}
\value{
returns an object with the following arguments:
\item{baseline}{baseline of the semparametric additive risk
model}
\item{pls.comp}{the pls components, i.e. the covariates multipled on the beta
coefficients.}
\item{beta}{risk regression coefficients related to the pls components.}
\item{beta.pls}{regression coefficients that defines the pls componets.}
\item{tbeta.pls}{the combined regression coefficients from the pls components and the regression coefficients, these leads to risk predictions when applied to
new covariates.}
}
\references{
Martinussen and Scheike, Dynamic Regression Models for
Survival Data, Springer (2006).

Martinussen and Scheike, The Aalen additive hazards model
with high-dimensional regressors, submitted.
}
\author{Thomas Scheike}
\examples{
library(survival)
data(sTRACE)

### 2 PLS components for the 5 covariates 
### based on Gui-Li approach with constant effects
out<-additive.pls(Surv(time,status==9)~age+sex+
vf+diabetes+chf,sTRACE,max.time=7,pls.dim=2,constant=1)

### 2 PLS components for the 5 covariates 
### based on Gui-Li approach with timevarying effects
out<-additive.pls(Surv(time,status==9)~age+sex+
vf+diabetes+chf,sTRACE,max.time=7,pls.dim=2,constant=0)

### 2 PLS components for the 3 covariates 
### age and sex are fixed covariates 
out<-additive.pls(Surv(time,status==9)~const(age)+const(sex)+
vf+diabetes+chf,sTRACE,max.time=7,pls.dim=2)
}
\keyword{survival}
