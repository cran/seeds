% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/modelClass.R
\name{nominalSol}
\alias{nominalSol}
\alias{nominalSol,odeModel-method}
\title{Calculate the nominal solution of the model}
\usage{
nominalSol(odeModel)

\S4method{nominalSol}{odeModel}(odeModel)
}
\arguments{
\item{odeModel}{a object of the class ode model describing the experiment}
}
\value{
a matrix with the numeric solution to the nominal ode equation
}
\description{
After an model is defined it can be evaluated. This returns the numerical solution
for the state equation before hidden inputs are calculated.
}
\examples{

lotka_voltera <- function (t, x, parameters) {
with(as.list(c(x,parameters)), {
 dx1 = x[1]*(alpha - beta*x[2])
  dx2 = -x[2]*(gamma - delta*x[1])
 return(list(c(dx1, dx2)))
})
}

pars <- c(alpha = 2, beta = .5, gamma = .2, delta = .6)
init_state <- c(x1 = 10, x2 = 10)
time <- seq(0, 100, by = 1)
lotVolModel = odeModel(func = lotka_voltera, parms = pars, times = time, y = init_state)
nominalSol(lotVolModel) 

}
