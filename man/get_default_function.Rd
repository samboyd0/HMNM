% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_integrated_network.R
\name{get_default_function}
\alias{get_default_function}
\title{Get a built-in function from keyword}
\usage{
get_default_function(f)
}
\arguments{
\item{f}{character string: 'binary', 'shift_scale', 'p_value', 'exp', or NULL}
}
\value{
a function
}
\description{
Return a default function used to compute seed values for random walk with restart.
}
\details{
For the input \emph{f}, NULL means no transformation is done.

'binary' coerces to numeric and then sets all positive values to 1, and to 0 otherwise.

'shift_scale' is for data types with a range centered about zero and takes two parameters: \emph{DOI} (direction of interest: 1 for positive, -1 for negative, and 0 for both) and \emph{w} (numeric value between 0 and 1). It takes the absolute values of the data and then down-weights nodes that weren't in DOI by \emph{w} (if \emph{DOI}=0, \emph{w} is coerced to 1).

'p_value' assumes the data are p-values and calculates \emph{-log10(x)} .

'exp' exponentiates the values and has a \emph{DOI} argument. For \emph{DOI}=-1 or 1, it is exp(\emph{DOI} * \emph{x}). For \emph{DOI}=0, it is exp(abs( \emph{x})).
}
\examples{
get_default_function("p_value")
get_default_function(NULL)

}
