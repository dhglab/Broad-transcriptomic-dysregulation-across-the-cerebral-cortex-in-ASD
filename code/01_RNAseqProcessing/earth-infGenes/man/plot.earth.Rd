\name{plot.earth}
\alias{plot.earth}
\alias{earth_plotmodsel}
\title{Plot an earth object}
\description{
Plot an \code{\link{earth}} object.
By default the plot shows model selection, cumulative distribution
of the residuals, residuals versus fitted values, and the residual QQ plot.

This function calls \code{\link[plotmo]{plotres}} internally.
The first arguments are identical to \code{plotres}.
}
\usage{
\method{plot}{earth}(x = stop("no 'x' argument"),

  # the following are identical to plotres arguments

  which = 1:4, info = FALSE, versus = 1, standardize = FALSE, delever = FALSE,
  level = 0, id.n = 3, labels.id = NULL, smooth.col = 2, grid.col = 0,
  jitter = 0, do.par = NULL, caption = NULL,
  trace = 0, npoints = 3000, center = TRUE, type = NULL, nresponse = NA,

  # the following are earth specific

  col.cv = "lightblue", col.grsq = 1, col.rsq = 2, col.infold.rsq = 0,
  col.mean.infold.rsq = 0, col.mean.oof.rsq = "palevioletred",
  col.npreds = if(is.null(object$cv.oof.rsq.tab)) 1 else 0, col.oof.labs = 0,
  col.oof.rsq = "mistyrose2", col.oof.vline = col.mean.oof.rsq,
  col.pch.cv.rsq = 0, col.pch.max.oof.rsq = 0, col.vline = col.grsq,
  col.vseg = 0, lty.grsq = 1, lty.npreds = 2, lty.rsq = 5, lty.vline = "12",
  legend.pos = NULL, \dots)

earth_plotmodsel( # for internal use by plotres
  x, col.rsq = 2, col.grsq = 1, col.infold.rsq = 0,
  col.mean.infold.rsq = 0, col.mean.oof.rsq = "palevioletred",
  col.npreds = NULL, col.oof.labs = 0, col.oof.rsq = "mistyrose2",
  col.oof.vline = col.mean.oof.rsq, col.pch.cv.rsq = 0,
  col.pch.max.oof.rsq = 0, col.vline = col.grsq, col.vseg = 0,
  lty.grsq = 1, lty.npreds = 2, lty.rsq = 5, lty.vline = "12",
  legend.pos=NULL, add = FALSE, jitter = 0,
  max.nterms = length(object$rss.per.subset),
  max.npreds=max(1,get.nused.preds.per.subset(object$dirs,object$prune.terms)),
  ...)
}
\arguments{
\item{x}{
An \code{\link{earth}} object.
This is the only required argument.
(It is called "x" for consistency with the generic \code{\link{plot}}.)
\cr
\cr
}
\item{which,info,versus}{These arguments are identical to \code{\link[plotmo]{plotres}}.
Please see the help page for \code{\link[plotmo]{plotres}}.
}
\item{standardize,delever,level}{.}
\item{id.n,labels.id,smooth.col}{.}
\item{grid.col,jitter}{.}
\item{do.par,caption,trace}{.}
\item{npoints,center}{.}
\item{type,nresponse}{.\cr
\cr
}
\item{col.cv}{Default \code{"lightblue"}.
Color of cross validation line in the residuals plot.
This is the residual of the mean out-fold-predicted value.
\cr
\cr
\bold{The following arguments are for the model selection plot.}
\cr
\cr
}
\item{col.grsq}{Default \code{1}.
Color of GRSq line in the Model Selection plot.
Use \code{0} for no GRSq line.
}
\item{col.rsq}{Default \code{2}.
Color of the RSq line in the Model Selection plot.
Use \code{0} for no RSq line.
}
\item{col.infold.rsq}{
Color of in-fold RSq lines for each fold in the Model Selection plot.
Applies only if \code{nfold} and \code{keepxy} were used
in the original call to \code{earth}.
Default is \code{0}, lines not plotted.
}
\item{col.mean.infold.rsq}{
Color of mean in-fold RSq for each number of terms in the Model Selection plot.
Default is \code{0}, line not plotted.
Applies only if \code{nfold} and \code{keepxy} were used
in the original call to \code{earth}.
}
\item{col.mean.oof.rsq}{Default \code{"palevioletred"}.
Color of mean out-of-fold RSq for each number of terms in the Model Selection plot.
Applies only if \code{nfold} and \code{keepxy} were used
in the original call to \code{earth}.
Use \code{0} to not plot this line.
}
\item{col.npreds}{Color
 of the "number of predictors" plot in the Model Selection plot.
 The default displays the number of predictors unless the \code{oof.rsq}'s are displayed.
 Use \code{0} for no "number of predictors" plot.
}
\item{col.oof.labs}{
Color of fold number labels on the \code{oof.rsq} lines.
Default is \code{0}, no labels.
}
\item{col.oof.rsq}{
Color of out-of-fold RSq lines for each fold in the Model Selection plot.
Applies only if \code{nfold} and \code{keepxy} were used
in the original call to \code{earth}.
Default is \code{"mistyrose2"}, a pale pink.
Use \code{0} to not plot these lines.
May be a vector of colors, which will be recycled if necessary.
}
\item{col.oof.vline}{
Color of vertical line at the maximum \code{oof.rsq} in the Model Selection plot.
Default is \code{col.mean.oof.rsq}.
}
\item{col.pch.cv.rsq}{
Color of point plotted on the \code{oof.rsq} line to indicate the \code{cv.rsq}.
for that fold
(i.e., it is plotted at the number of terms selected by the in-fold GCV).
Default is \code{0}, point not plotted.
}
\item{col.pch.max.oof.rsq}{
Color of point plotted on the \code{oof.rsq} line to indicate the
maximum \code{oof.rsq} for that fold.
Default is \code{0}, point not plotted.
}
\item{col.vline}{
Color of the vertical line at selected model in the Model Selection plot.
Default is \code{col.grsq}.
This will be at the maximum GRSq unless \code{pmethod="none"}.
Use \code{0} for no vertical line.
}
\item{col.vseg}{
Default is \code{0}.
Color of triangular marker at top of vertical line for best GRSq.
}
\item{lty.grsq}{
Line type of GRSq line in the Model Selection plot.
Default is \code{1}
}
\item{lty.npreds}{
Line type of the "number of predictors" plot in the Model Selection plot.
Default is \code{2}.
}
\item{lty.rsq}{
Line type of RSq line in the Model Selection plot.
Default is \code{5}.
}
\item{lty.vline}{
Line type of vertical line at selected model in the Model Selection plot.
Default is \code{"12"}.
}
\item{legend.pos}{
Position of the legend in the Model Selection plot.
Default is \code{NULL} meaning automatic.
Use \code{legend.pos=NA} or \code{0} for no legend.
}
\item{add, max.nterms, max.npreds}{
\code{earth_plotmodsel} arguments for internal use by
\code{plotres}.
\cr
\cr
\cr
}
\item{\dots}{Please see \code{\link[plotmo]{plotres}} for the details
on the dots arguments.

The \code{ylim} argument is treated specially in the model selection plot:
\code{ymin} equal to \code{-1} means use
the smallest GRSq or RSq value, excluding the intercept,
and \code{ymax} equal \code{-1} means use
the largest GRSq or RSq value.
}
}
\note{
For details on interpreting the graphs,
please see the \code{earth} package vignettes
\dQuote{\href{../doc/earth-notes.pdf}{Notes on the earth package}}
and
\dQuote{\href{../doc/earth-varmod.pdf}{Variance models in earth}}.

Note that cross-validation data will not be displayed unless
both \code{nfold} and \code{keepxy} were used in the original call to
\code{earth}.

\code{earth_plotmodsel} is provided for use by \code{\link[plotmo]{plotres}}.
}
\seealso{
  \code{\link{earth}},
  \code{\link{plot.earth.models}},
  \code{\link{plotd}},
  \code{\link[plotmo]{plotmo}}
}
\examples{
data(ozone1)
earth.mod <- earth(O3 ~ ., data = ozone1, degree = 2)
plot(earth.mod)
}
\keyword{models}
