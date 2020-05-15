\name{plotd}
\alias{plotd}
\title{Plot the distribution of predictions for each class}
\description{
Draw a plot of the distribution of the predicted values for each class.
Can be used for \code{\link[earth]{earth}} models, but also for models built by
\code{\link{lm}},
\code{\link{glm}},
\code{\link[MASS]{lda}},
etc.
}
\usage{
plotd(object, hist = FALSE, type = NULL, nresponse = NULL, dichot = FALSE,
      trace = FALSE, xlim = NULL, ylim = NULL, jitter = FALSE, main=NULL,
      xlab = "Predicted Value", ylab = if(hist) "Count" else "Density",
      lty = 1, col = c("gray70", 1, "lightblue", "brown", "pink", 2, 3, 4),
      fill = if(hist) col[1] else 0,
      breaks = "Sturges", labels = FALSE,
      kernel = "gaussian", adjust = 1, zero.line = FALSE,
      legend = TRUE, legend.names = NULL, legend.pos = NULL,
      cex.legend = .8, legend.bg = "white", legend.extra = FALSE,
      vline.col = 0, vline.thresh = .5, vline.lty = 1, vline.lwd = 1,
      err.thresh = vline.thresh, err.col = 0, err.border = 0, err.lwd = 1,
      xaxt = "s", yaxt = "s", xaxis.cex = 1, sd.thresh = 0.01, ...)
}
\arguments{
  To start off, look at the arguments \code{object}, \code{hist}, \code{type}.\cr
  For predict methods with multiple column responses, see the \code{nresponse} argument.\cr
  For factor responses with more than two levels, see the \code{dichot} argument.

  \item{object}{
    Model object.
    Typically a model which predicts a class or a class discriminant.
  }
  \item{hist}{
    \code{FALSE} (default) to call \code{\link{density}} internally.\cr
    \code{TRUE} to call \code{\link{hist}} internally.
  }
  \item{type}{
    Type parameter passed to \code{\link{predict}}.
    For allowed values see the \code{predict} method for
    your \code{object} (such as \code{\link[earth]{predict.earth}}).
    By default, \code{plotd} tries to automatically select a suitable
    value for the model in question.
    (This is \code{"response"} for all objects except \code{rpart} models, where \code{"vector"}
    is used. The choices will often be inappropriate.)
    Typically you would set \code{hist=TRUE} when \code{type="class"}.
  }
  \item{nresponse}{
    Which column to use when \code{predict} returns multiple columns.
    This can be a column index or column name
    (which may be abbreviated, partial matching is used).
    The default is \code{NULL}, meaning use all columns of the predicted response.
  }
  \item{dichot}{
     Dichotimise the predicted response.
     This argument is ignored except for models where the observed response
     is a factor with more than two levels
     and the predicted response is a numeric vector.
     The default \code{FALSE} separates the response into a group for each factor.
     With \code{dichot=TRUE} the response is separated into just two groups:
     the first level of the factor versus the remaining levels.
  }
  \item{trace}{
     Default \code{FALSE}.
     Use \code{TRUE} or \code{1} to trace \code{plotd} ---
     useful to see how \code{plotd}
     partitions the predicted response into classes.
     Use \code{2} for more details.
  }
  \item{xlim}{
    Limits of the x axis.
    The default \code{NULL} means determine these limits automatically,
    else specify \code{c(xmin,xmax)}.
  }
  \item{ylim}{
    Limits of the y axis.
    The default \code{NULL} means determine these limits automatically,
    else specify \code{c(ymin,ymax)}.
  }
  \item{jitter}{
     Jitter the histograms or densities horizontally to minimize overplotting.
     Default \code{FALSE}.
     Specify \code{TRUE} to automatically calculate the jitter,
     else specify a numeric jitter value.
  }
  \item{main}{
    Main title. Values:\cr
       \code{"string"} string\cr
       \code{""} no title\cr
       \code{NULL} (default) generate a title from the call.
  }
  \item{xlab}{
    x axis label.
    Default is \code{"Predicted Value"}.
  }
  \item{ylab}{
    y axis label.
    Default is \code{if(hist) "Count" else "Density"}.
  }
  \item{lty}{
    Per class line types for the plotted lines.
    Default is 1 (which gets recycled for all lines).
  }
  \item{col}{
    Per class line colors.
    The first few colors of the default are intended to be
    easily distinguishable on both color displays and monochrome printers.
  }
  \item{fill}{
    Fill color for the plot for the first class.
    For \code{hist=FALSE}, the default is 0, i.e., no fill.
    For \code{hist=TRUE}, the default is the first element in the \code{col} argument.
  }
  \item{breaks}{
    Passed to \code{\link{hist}}.
    Only used if \code{hist=TRUE}.
    Default is \code{"Sturges"}.
    When \code{type="class"}, setting \code{breaks} to a low number
    can be used to widen the histogram bars
  }
  \item{labels}{
    \code{TRUE} to draw counts on the \code{\link{hist}} plot.
    Only used if \code{hist=TRUE}.
    Default is \code{FALSE}.
  }
  \item{kernel}{
    Passed to \code{\link{density}}.
    Only used if \code{hist=FALSE}.
    Default is \code{"gaussian"}.
  }
  \item{adjust}{
    Passed to \code{\link{density}}.
    Only used if \code{hist=FALSE}.
    Default is \code{1}.
  }
  \item{zero.line}{
    Passed to \code{\link{plot.density}}.
    Only used if \code{hist=FALSE}.
    Default is \code{FALSE}.
  }
  \item{legend}{
    \code{TRUE} (default) to draw a legend, else \code{FALSE}.
  }
  \item{legend.names}{
    Class names in legend.
    The default \code{NULL} means determine these automatically.
  }
  \item{legend.pos}{
    Position of the legend.
    The default \code{NULL} means position the legend automatically,
    else specify \code{c(x,y)}.
  }
  \item{cex.legend}{
    \code{cex} for \code{\link[graphics]{legend}}.
    Default is \code{.8}.
  }
  \item{legend.bg}{
   \code{bg} color for \code{\link[graphics]{legend}}.
   Default is \code{"white"}.
  }
  \item{legend.extra}{
    Show (in the legend) the number of occurrences of each class.
    Default is \code{FALSE}.
  }
  \item{vline.thresh}{
    Horizontal position of optional vertical line.
    Default is \code{0.5}.
    The vertical line is intended to indicate class separation.
    If you use this, don't forget to set \code{vline.col}.
  }
  \item{vline.col}{
    Color of vertical line.
    Default is 0, meaning no vertical line.
  }
  \item{vline.lty}{
    Line type of vertical line.
    Default is \code{1}.
  }
  \item{vline.lwd}{
    Line width of vertical line.
    Default is \code{1}.
  }
  \item{err.thresh}{
    x axis value specifying the error shading threshold.
    See \code{err.col}.
    Default is \code{vline.thresh}.
  }
  \item{err.col}{
    Specify up to three colors to shade the "error areas" of the density plot.
    The default is \code{0}, meaning no error shading.
    This argument is ignored unless \code{hist=FALSE}.
    If there are more than two classes, \code{err.col} uses only the first two.
    This argument is best explained by running an example:\preformatted{
data(etitanic)
earth.mod <- earth(survived ~ ., data=etitanic)
plotd(earth.mod, vline.col=1, err.col=c(2,3,4))
      }
    The three areas are (i) the error area to the left of the threshold,
    (ii) the error area to the right of the threshold, and,
    (iii) the reducible error area.
    If less than three values are specified, \code{plotd} re-uses values in a sensible manner.
    Use values of \code{0} to skip areas.
    Disjoint regions are not handled well by the current implementation.
  }
  \item{err.border}{
    Borders around the error shading.
    Default is \code{0}, meaning no borders, else specify up to three colors.
  }
  \item{err.lwd}{
    Line widths of borders of the error shading.
    Default is \code{1}, else specify up to three line widths.
  }
  \item{xaxt}{
    Default is \code{"s"}.
    Use \code{xaxt="n"} for no x axis.
  }
  \item{yaxt}{
    Default is \code{"s"}.
    Use \code{yaxt="n"} for no y axis.
  }
  \item{xaxis.cex}{
    Only used if \code{hist=TRUE} and \code{type="class"}.
    Specify size of class labels drawn on the x axis.
    Default is 1.
  }
  \item{sd.thresh}{
    Minimum acceptable standard deviation for a density.
    Default is \code{0.01}.
    Densities with a standard deviation less than \code{sd.thresh}
    will not be plotted (a warning will be issued and the legend
    will say \code{"not plotted"}).
  }
  \item{\dots}{
    Extra arguments passed to the predict method for the object.
  }
}
\note{
  This function calls \code{\link{predict}}
  with the data originally used to build
  the model, and with the \code{type} specified above.
  It then separates the predicted values into classes,
  where the class for each predicted value
  is determined by the class of the observed response.
  Finally, it calls \code{\link{density}} (or \code{\link{hist}}
  if \code{hist=TRUE}) for each class-specific set of values,
  and plots the results.

  This function estimates distributions with the
  \code{\link{density}} and \code{\link{hist}} functions,
  and also calls \code{\link{plot.density}} and \code{\link{plot.histogram}}.
  For an overview see Venables and Ripley MASS section 5.6.

  \bold{Partitioning the response into classes}

  Considerable effort is made to partition the predicted response
  into classes in a sensible way.
  This is not always possible for multiple column responses and the \code{nresponse} argument
  should be used where necessary.
  The partitioning details depend on the types and numbers of columns in the observed
  and predicted responses.
  These in turn depend on the model object and the \code{type} argument.

  Use the \code{trace} argument to see how \code{plotd} partitions the
  response for your model.

  \bold{Degenerate densities}

  A message such as\cr
      \code{Warning: standard deviation of "male" density is 0,  density is degenerate?}\cr
  means that the density for that class will not be plotted
  (the legend will say \code{"not plotted"}).

  Set \code{sd.thresh=0} to get rid of this check,
  but be aware that histograms (and sometimes x axis labels)
  for degenerate densities will be misleading.

  \bold{Using plotd for various models}

  This function is included in the \code{\link[earth]{earth}} package
  but can also be used with other models.

  Example with \code{\link{glm}}: \preformatted{
      library(earth); data(etitanic)
      glm.model <- glm(sex ~ ., data=etitanic, family=binomial)
      plotd(glm.model)
  }
  Example with \code{\link{lm}}:\preformatted{
      library(earth); data(etitanic)
      lm.model <- lm(as.numeric(sex) ~ ., data=etitanic)
      plotd(lm.model)
  }
%   Example with \code{\link[rpart]{rpart}}:\preformatted{
%        library(rpart); library(earth); data(etitanic)
%        rpart.model <- rpart(sex ~ ., data = etitanic, method="class")
%        plotd(rpart.model, type="prob", nresponse=1)
%        plotd(rpart.model, type="prob", nresponse=2)
%        plotd(rpart.model, type="class", hist=TRUE, labels=TRUE)
%   }
  \bold{Using plotd with lda or qda}

    The \code{plotd} function has special handling for \code{\link[MASS]{lda}}
    (and \code{\link[MASS]{qda}}) objects.
    For such objects, the \code{type} argument can take one of the
    following values:

    \code{"response"} (default) linear discriminant\cr
    \code{"ld"} same as \code{"response"}\cr
    \code{"class"} predicted classes\cr
    \code{"posterior"} posterior probabilities

    Example:\preformatted{
    library(MASS); library(earth); data(etitanic)
    lda.model <- lda(sex ~ ., data=etitanic)
    plotd(lda.model) # linear discriminant by default
    plotd(lda.model, type="class", hist=TRUE, labels=TRUE)
    }
    This handling of \code{type} is handled internally by \code{plotd}
    and \code{type} is not passed to \code{predict.lda}
    (\code{type} is used merely to select fields in the list
    returned by \code{predict.lda}).
    The type names can be abbreviated down to a single character.

    For objects created with \code{lda.matrix}
    (as opposed to \code{lda.formula}),
    \code{plotd} blindly assumes that the \code{grouping} argument was the second argument.

    \code{plotd} does not yet support objects created with \code{lda.data.frame}.

    For \code{lda} responses with more than two factor levels,
    use the \code{nresponse} argument to
    select a column in the predicted response.
    Thus with the default \code{type=NULL},
    (which gets automatically converted by \code{plotd} to \code{type="response"}),
    use \code{nresponse=1} to select just the first linear discriminant.
    The default \code{nresponse=NULL} selects all columns,
    which is typically not what you want for \code{lda} models.
    Example:\preformatted{
    library(MASS); library(earth);
    set.seed(1)      # optional, for reproducibility
    example(lda)     # creates a model called "z"
    plot(z, dimen=1) # invokes plot.lda from the MASS package
    plotd(z, nresponse=1, hist=1) # equivalent using plotd
                                 # nresponse=1 selects first linear discr.
    }
    The \code{dichot=TRUE} argument is also useful for \code{lda}
    responses with more than two factor levels.\cr

  \bold{TODO}

  Handle degenerate densities in a more useful way.\cr
  Add \code{freq} argument for \code{\link{hist}}.
}
\seealso{
  \code{\link{density}}, \code{\link{plot.density}}\cr
  \code{\link{hist}}, \code{\link{plot.histogram}}\cr
  \code{\link[earth]{earth}}, \code{\link[earth]{plot.earth}}
}
\examples{
if (require(earth)) {
    old.par <- par(no.readonly=TRUE);
    par(mfrow=c(2,2), mar=c(4, 3, 1.7, 0.5), mgp=c(1.6, 0.6, 0), par(cex = 0.8))
    data(etitanic)
    mod <- earth(survived ~ ., data=etitanic, degree=2, glm=list(family=binomial))

    plotd(mod)

    plotd(mod, hist=TRUE, legend.pos=c(.25,220))

    plotd(mod, hist=TRUE, type="class", labels=TRUE, xlab="", xaxis.cex=.8)

    par(old.par)
}
}
\keyword{models}
