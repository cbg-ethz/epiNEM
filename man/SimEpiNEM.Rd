% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/epiN.R
\name{SimEpiNEM}
\alias{SimEpiNEM}
\title{Compare algorithms.}
\usage{
SimEpiNEM(
  runs = 10,
  do = c("n", "e"),
  random = list(FPrate = 0.1, FNrate = c(0.1, 0.5), single = 3, double = 1, reporters =
    10, replicates = 2),
  maxTime = FALSE,
  forcelogic = TRUE,
  epinemsearch = "greedy",
  bnemsearch = "genetic",
  ...
)
}
\arguments{
\item{runs}{number simulation runs}

\item{do}{string vector of algorithms to
compare: e (epiNEM), n (Nested Effects Models),
b (B-NEM), p (PC algorithm), a (Aracne), e.g. c("e", "n", "p")}

\item{random}{list of false poitive rate FPrate, false negative rates
FNrate, number of single knock-downs single, number of double
knock-downs double, number of effect reporters reporters and number
of replicates replicates}

\item{maxTime}{TRUE if the algorithms are bound to a maximum running
time in respect to epiNEM}

\item{forcelogic}{if TRUE the randomly sampled ground truth network
includes a complex logic with probability 1}

\item{epinemsearch}{greedy or exhaustive search for epiNEM}

\item{bnemsearch}{genetic or greedy search for B-NEM}

\item{...}{additional parameters}
}
\value{
returns list of specificity and sensitivity of inferred edges
(spec, sens) and inferred expected data (spec2, sens2) and accuracy
of logics (logics) and running time (time)
}
\description{
Compares different network reconstruction algorithm on
simulated data.
}
\examples{
res <- SimEpiNEM(runs = 1)
}
\author{
Martin Pirkl
}
