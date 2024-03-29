% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mixedPenalty.R
\name{addScad}
\alias{addScad}
\title{addScad}
\usage{
addScad(mixedPenalty, regularized, lambdas, thetas)
}
\arguments{
\item{mixedPenalty}{model of class mixedPenalty created with the mixedPenalty function (see ?mixedPenalty)}

\item{regularized}{vector with names of parameters which are to be regularized.
If you are unsure what these parameters are called, use
getLavaanParameters(model) with your lavaan model object}

\item{lambdas}{numeric vector: values for the tuning parameter lambda}

\item{thetas}{parameters whose absolute value is above this threshold will be penalized with
a constant (theta)}
}
\value{
Model of class mixedPenalty. Use the fit() - function to fit the model
}
\description{
Implements scad regularization for structural equation models.
The penalty function is given by:
\ifelse{html}{
\deqn{p( x_j) = \begin{cases}
\lambda |x_j| & \text{if } |x_j| \leq \theta\\
\frac{-x_j^2 + 2\theta\lambda |x_j| - \lambda^2}{2(\theta -1)} & 
\text{if } \lambda < |x_j| \leq \lambda\theta \\
(\theta + 1) \lambda^2/2 & \text{if } |x_j| \geq \theta\lambda\\
\end{cases}}
where \eqn{\theta > 2}.}{
Equation Omitted in Pdf Documentation.
}
}
\details{
Identical to \pkg{regsem}, models are specified using \pkg{lavaan}. Currently,
most standard SEM are supported. \pkg{lessSEM} also provides full information
maximum likelihood for missing data. To use this functionality,
fit your \pkg{lavaan} model with the argument \code{sem(..., missing = 'ml')}.
\pkg{lessSEM} will then automatically switch to full information maximum likelihood
as well.

scad regularization:
\itemize{
\item Fan, J., & Li, R. (2001). Variable selection via nonconcave penalized
likelihood and its oracle properties. Journal of the American Statistical Association,
96(456), 1348–1360. https://doi.org/10.1198/016214501753382273
}

Regularized SEM
\itemize{
\item Huang, P.-H., Chen, H., & Weng, L.-J. (2017). A Penalized Likelihood Method for Structural Equation Modeling. Psychometrika, 82(2), 329–354. https://doi.org/10.1007/s11336-017-9566-9
\item Jacobucci, R., Grimm, K. J., & McArdle, J. J. (2016). Regularized Structural Equation Modeling. Structural
Equation Modeling: A Multidisciplinary Journal, 23(4), 555–566. https://doi.org/10.1080/10705511.2016.1154793
}

For more details on ISTA, see:
\itemize{
\item Beck, A., & Teboulle, M. (2009). A Fast Iterative Shrinkage-Thresholding
Algorithm for Linear Inverse Problems. SIAM Journal on Imaging Sciences, 2(1),
183–202. https://doi.org/10.1137/080716542
\item Gong, P., Zhang, C., Lu, Z., Huang, J., & Ye, J. (2013).
A General Iterative Shrinkage and Thresholding Algorithm for Non-convex
Regularized Optimization Problems. Proceedings of the 30th International
Conference on Machine Learning, 28(2)(2), 37–45.
\item Parikh, N., & Boyd, S. (2013). Proximal Algorithms. Foundations and
Trends in Optimization, 1(3), 123–231.
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

# We can add mixed penalties as follows:

regularized <- lavaanModel |>
  # create template for regularized model with mixed penalty:
  mixedPenalty() |>
  # add penalty on loadings l6 - l10:
  addScad(regularized = paste0("l", 11:15), 
          lambdas = seq(0,1,.1),
          thetas = 3.1) |>
  # fit the model:
  fit()
}
