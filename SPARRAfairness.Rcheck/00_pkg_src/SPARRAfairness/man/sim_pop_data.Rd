% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SPARRAfairness_simulation.R
\name{sim_pop_data}
\alias{sim_pop_data}
\title{sim_pop_data}
\usage{
sim_pop_data(
  npop,
  coef_adjust = 4,
  offset = 1,
  vcor = NULL,
  coefs = c(2, 1, 0, 5, 3, 0, 0),
  seed = 12345,
  incl_id = TRUE,
  incl_reason = TRUE
)
}
\arguments{
\item{npop}{population size}

\item{coef_adjust}{inverse scale for all (true) coefficients (default 4): lower means that hospital admissions are more predictable from covariates.}

\item{offset}{offset for logistic model (default 1): higher means a lower overall prevalence of admission}

\item{vcor}{a valid 5x5 correlation matrix (default NULL), giving correlation between variables. If 'NULL', values roughly represents realistic data.}

\item{coefs}{coefficients of age, male sex, non-white ethnicity, number of previous admissions, and deprivation decile on hospital admissions, Default (2,1,0,5,3). Divided through by coef_adjust.}

\item{seed}{random seed (default 12345)}

\item{incl_id}{include an ID column (default TRUE)}

\item{incl_reason}{include a column indicating reason for admission.}
}
\value{
data frame with realistic values.
}
\description{
Simulates population data with a reasonably realistic joint distribution
}
\details{
Simulates data for a range of people for the variables
\itemize{
\item Age (\code{age})
\item Sex (\code{sexM}; 1 if male)
\item Race/ethnicity (\code{raceNW}: 1 if non-white ethnicity)
\item Number of previous hospital admissions (\code{PrevAdm})
\item Deprivation decile (\code{SIMD}: 1 most deprived, 10 least deprived. NOTE - opposite to English IMD)
\item Urban-rural residence status (\code{urban_rural}: 1 for rural)
\item Mainland-island residence status (\code{mainland_island}: 1 for island)
\item Hospital admission (\code{target}: 1/TRUE if admitted to hospital in year following prediction date)
}

Can optionally add an ID column.

Optionally includes an admission reason for samples with \code{target=1}. These admission reasons
roughly correspond to the first letters of ICD10 categories, and can either correspond to an
admission or death. Admission reasons are simulated with a non-constant multinomial distribution
which varies across age/sex/ethnicity/urban-rural/mainland-island/PrevAdm values in a randomly-
chosen way. The distributions of admission reasons are \emph{not} however chosen to reflect real
distributions, nor are systematic changes in commonality of admission types across categories
intended to appear realistic.
}
\examples{

# Simulate data
dat=sim_pop_data(10000)
pairs(dat[,1:7])
cor(dat[,1:7])

# See vignette
}
