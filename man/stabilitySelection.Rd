% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stabilitySelection.R
\name{stabilitySelection}
\alias{stabilitySelection}
\title{stabilitySelection}
\usage{
stabilitySelection(
  modelSpecification,
  subsampleSize,
  numberOfSubsamples = 100,
  threshold = 70,
  maxTries = 10 * numberOfSubsamples
)
}
\arguments{
\item{modelSpecification}{a call to one of the penalty functions in lessSEM. See
examples for details}

\item{subsampleSize}{number of subjects in each subsample. Must be smaller than
the number of subjects in the original data set}

\item{numberOfSubsamples}{number of times the procedure should subsample and
recompute the model. According to Meinshausen & Bühlmann (2010), 100 seems to
work quite well and is also the default in regsem}

\item{threshold}{percentage of models, where the parameter should be contained in order
to be in the final model}

\item{maxTries}{fitting models in a subset may fail. maxTries sets the maximal
number of subsets to try.}
}
\value{
estimates for each subsample and aggregated percentages for each parameter
}
\description{
Provides rudimentary stability selection for regularized SEM. Stability
selection has been proposed by Meinshausen & Bühlmann (2010) and was
extended to SEM by Li & Jacobucci (2021). The problem that stabiltiy selection
tries to solve is the instability of regularization procedures: Small changes in
the data set may result in different parameters being selected. To address
this issue, stability selection uses random subsamples from the initial data
set and fits models in these subsamples. For each parameter, we can now check
how often it is included in the model for a given set of tuning parameters.
Plotting these probabilities can provide an overview over which of the parameters
are often removed and which remain in the model most of the time. To get
a final selection, a threshold t can be defined: If a parameter is in the model
t\% of the time, it is retained.
}
\section{References}{
\itemize{
\item Li, X., & Jacobucci, R. (2021). Regularized structural equation modeling with
stability selection. Psychological Methods, 27(4), 497–518. https://doi.org/10.1037/met0000389
\item Meinshausen, N., & Bühlmann, P. (2010). Stability selection. Journal of the
Royal Statistical Society: Series B (Statistical Methodology), 72(4), 417–473.
https://doi.org/10.1111/j.1467-9868.2010.00740.x
}
}

\examples{
library(lessSEM)

# Identical to regsem, lessSEM builds on the lavaan
# package for model specification. The first step
# therefore is to implement the model in lavaan.

dataset <- simulateExampleData()

lavaanSyntax <- "
f =~ l1*y1 + l2*y2 + l3*y3 + l4*y4 + l5*y5 + 
     l6*y6 + l7*y7 + l8*y8 + l9*y9 + l10*y10 + 
     l11*y11 + l12*y12 + l13*y13 + l14*y14 + l15*y15
f ~~ 1*f
"

lavaanModel <- lavaan::sem(lavaanSyntax,
                           data = dataset,
                           meanstructure = TRUE,
                           std.lv = TRUE)

# Stability selection
stabSel <- stabilitySelection(
  # IMPORTANT: Wrap your call to the penalty function in an rlang::expr-Block:
  modelSpecification = 
    rlang::expr(
      lasso(
        # pass the fitted lavaan model
        lavaanModel = lavaanModel,
        # names of the regularized parameters:
        regularized = paste0("l", 6:15),
        # in case of lasso and adaptive lasso, we can specify the number of lambda
        # values to use. lessSEM will automatically find lambda_max and fit
        # models for nLambda values between 0 and lambda_max. For the other
        # penalty functions, lambdas must be specified explicitly
        nLambdas = 50)
    ),
  subsampleSize = 80,
  numberOfSubsamples = 5, # should be set to a much higher number (e.g., 100)
  threshold = 70
)
stabSel
plot(stabSel)
}
