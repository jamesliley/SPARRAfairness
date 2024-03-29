% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SPARRAfairness_functions.R
\name{plot_decomp}
\alias{plot_decomp}
\title{plot_decomp}
\usage{
plot_decomp(decomp1, decomp2, threshold, labels, inc_died = TRUE)
}
\arguments{
\item{decomp1}{matrix for first group; see specification in description}

\item{decomp2}{matrix for second group; see specification in description}

\item{threshold}{score threshold to plot (between 0 and 1)}

\item{labels}{labels for group 1 and group 2}

\item{inc_died}{set to TRUE to include a second panel showing 'death' type admissions}
}
\value{
Silently return ggplot object
}
\description{
Plots a bar graph of decomposition of FORP by cause of admission
}
\details{
Takes two matrices as input with the following specifications:
Each matrix corresponds to one group
Columns are named with the admission types to be plotted. Any admission
types including the string 'Died' are counted as deaths
If the matrix has N rows, these are interpreted as corresponding to N
score quantiles in increasing order.
The (i,j)th entry of the matrix is the number of people admitted for
reason i with a score greater than or equal to (j-1)/N and less than (j/N)
who are in that group
}
\examples{

# See vignette
}
