% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/update_beta_MM_ls.R
\name{update_beta_MM_ls}
\alias{update_beta_MM_ls}
\title{Beta update in L2E multivariate regression - MM}
\usage{
update_beta_MM_ls(y, X, beta, tau, max_iter = 100, tol = 1e-04)
}
\arguments{
\item{y}{Response}

\item{X}{Design matrix}

\item{beta}{Initial vector of regression coefficients}

\item{tau}{Precision estimate}

\item{max_iter}{Maximum number of iterations}

\item{tol}{Relative tolerance}
}
\value{
Returns a list object containing the new estimate for beta (vector) and the number of iterations (scalar) the update step utilized
}
\description{
\code{update_beta_MM_ls} updates beta for L2E multivariate regression using MM
}
