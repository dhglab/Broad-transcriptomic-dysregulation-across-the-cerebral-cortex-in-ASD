\name{earth}
\alias{earth}
\alias{earth.default}
\alias{earth.formula}
\alias{earth.fit}
\concept{regression}
\concept{mars}
\concept{Friedman}
\title{Multivariate Adaptive Regression Splines}
\description{
Build a regression model
using the techniques in Friedman's papers "Multivariate Adaptive Regression Splines"
and "Fast MARS".

See the package vignette
\dQuote{\href{../doc/earth-notes.pdf}{Notes on the earth package}}.
}
\usage{
\method{earth}{formula}(formula = stop("no 'formula' argument"), data = NULL,
   weights = NULL, wp = NULL, subset = NULL,
   na.action = na.fail,
   pmethod = c("backward", "none", "exhaustive", "forward", "seqrep", "cv"),
   keepxy = FALSE, trace = 0, glm = NULL, degree = 1, nprune = NULL,
   nfold=0, ncross=1, stratify=TRUE,
   varmod.method = "none", varmod.exponent = 1,
   varmod.conv = 1, varmod.clamp = .1, varmod.minspan = -3,
   Scale.y = (NCOL(y)==1), \dots)

\method{earth}{default}(x = stop("no 'x' argument"), y = stop("no 'y' argument"),
    weights = NULL, wp = NULL, subset = NULL,
    na.action = na.fail,
    pmethod = c("backward", "none", "exhaustive", "forward", "seqrep", "cv"),
    keepxy = FALSE, trace = 0, glm = NULL, degree = 1, nprune = NULL,
    nfold=0, ncross=1, stratify=TRUE,
    varmod.method = "none", varmod.exponent = 1,
    varmod.conv = 1, varmod.clamp = .1, varmod.minspan = -3,
    Scale.y = (NCOL(y)==1), \dots)

\method{earth}{fit}(x = stop("no 'x' argument"), y = stop("no 'y' argument"),
    weights = NULL, wp = NULL, subset = NULL,
    na.action = na.fail,
    pmethod = c("backward", "none", "exhaustive", "forward", "seqrep", "cv"),
    keepxy = FALSE, trace = 0, glm = NULL, degree = 1,
    penalty = if(degree > 1) 3 else 2,
    nk = min(200, max(20, 2 * ncol(x))) + 1,
    thresh = 0.001, minspan = 0, endspan = 0,
    newvar.penalty = 0, fast.k = 20, fast.beta = 1,
    linpreds = FALSE, allowed = NULL,
    nprune = NULL, Object = NULL,
    Scale.y = (NCOL(y)==1), Adjust.endspan = 2, Auto.linpreds = TRUE,
    Force.weights = FALSE, Use.beta.cache = TRUE, Force.xtx.prune = FALSE,
    Get.leverages = NROW(x) < 1e5, Exhaustive.tol = 1e-10, \dots)
}
\arguments{
To start off, look at the arguments
\code{formula},
\code{data},
\code{x},
\code{y},
\code{nk},
\code{degree}, and
\code{trace}.
\cr
If the response is binary or a factor, consider using the \code{glm}
argument.\cr
For cross validation, use the \code{nfold} argument.\cr
For prediction intervals, use the \code{varmod.method} argument.\cr
\cr
Most users will find that the above arguments are all they need,
plus in some cases \code{keepxy} and \code{nprune}.
Unless you are a knowledgeable user, it's best not subvert the
standard algorithm by toying with tuning parameters such as \code{thresh},
\code{penalty}, and \code{endspan}.
\cr

\item{formula}{
     Model formula.
}
\item{data}{
    Data frame for \code{formula}.
}
\item{x}{
     Matrix or dataframe containing the independent variables.
}
\item{y}{
     Vector containing the response variable, or, in the case of multiple responses,
     a matrix or dataframe whose columns are the values for each response.
}
\item{subset}{
     Index vector specifying which cases to use, i.e., which rows in \code{x} to use.
     Default is NULL, meaning all.
}
\item{weights}{
     Case weights.
     Default is NULL, meaning no case weights.
     If specified, \code{weights} must have length equal to \code{nrow(x)}
     before applying \code{subset}.
     Zero weights are converted to a very small nonzero value.
     % Following the precedent set by \code{lm}, \code{earth} first searches for the
     % weights in \code{data}, then in the environment in which \code{earth} is called.
}
\item{wp}{
     Response weights.
     Default is NULL, meaning no response weights.
     If specified, \code{wp} must have an element for each column of
     \code{y} (after \code{\link[=factor]{factors}} in
     \code{y}, if any, have been expanded).
     Zero weights are converted to a very small nonzero value.
%  Earth uses \code{wp} as follows:
%  It normalizes \code{wp} as follows \code{wp <- sqrt(wp / mean(wp))}
%  It multiplies each column of \code{y} by the corresponding
%  element of the normalized \code{wp}.
%  It runs its internal MARS machine with the modified \code{y}.
%  Finally, after running the machine it divides each column of
%  \code{fitted.values}, \code{residuals}, and \code{coefficients} by the
%  corresponding element of the normalized \code{wp}.
}
\item{na.action}{
     NA action. Default is \code{na.fail}, and only \code{na.fail} is supported.
}
\item{keepxy}{
     Default is \code{FALSE}.
     Set to \code{TRUE} to retain the following in the returned value: \code{x} and \code{y} (or \code{data}),
     \code{subset}, and \code{weights}.
     The function \code{\link{update.earth}} and friends will use these
     if present instead of searching for them
     in the environment at the time \code{update.earth} is invoked.\cr
     % The \code{subset} if specified will have been applied to the saved \code{x}, \code{y}, and \code{data}).\cr
     When the \code{nfold} argument is used with \code{keepxy=TRUE},
     \code{earth} keeps more data and calls \code{predict.earth} multiple
     times to generate \code{cv.oof.rsq.tab} and \code{cv.infold.rsq.tab}
     (see the \code{cv.} arguments in the \dQuote{\emph{Value}} section
     below).
     It therefore makes cross-validation significantly slower.
}
\item{trace}{
     Trace \code{earth}'s execution.  Values:\cr
        \code{0} (default) no tracing\cr
        \code{.3} variance model (the \code{varmod.method} arg)\cr
        \code{.5} cross validation (the \code{nfold} arg)\cr
        \code{1} overview\cr
        \code{2} forward pass\cr
        \code{3} pruning\cr
        \code{4} model mats summary, pruning details\cr
        \code{5} full model mats, internal details of operation\cr
}
\item{glm}{
     NULL (default) or a list of arguments to pass on to \code{\link{glm}}.
     See the documentation of \code{\link{glm}} for a description of these arguments
     See \dQuote{\emph{Generalized linear models}} in the vignette.
     Example:\cr
     \code{earth(survived~., data=etitanic, degree=2, glm=list(family=binomial))}
\cr\cr
     \bold{The following arguments are for the forward pass.}
}
\item{degree}{
     Maximum degree of interaction (Friedman's \eqn{mi}).
     Default is \code{1}, meaning build an additive model (i.e., no interaction terms).
}
\item{penalty}{
    Generalized Cross Validation (GCV) penalty per knot.
    Default is \code{if(degree>1) 3 else 2}.
    Simulation studies suggest values in the range of about \code{2} to \code{4}.
    The FAQ section in the vignette has some information on GCVs.\cr
    Special values (for use by knowledgeable users):
The value \code{0} penalizes only terms, not knots.
The value \code{-1} means no penalty, so GCV = RSS/n.
}
\item{nk}{
     Maximum number of model terms before pruning, i.e., the
     maximum number of terms created by the forward pass.
     Includes the intercept.\cr
     The actual number of terms created by the forward pass will often be
     less than \code{nk} because of other stopping conditions.
     See \dQuote{\emph{Termination conditions for the forward pass}}
     in the vignette.\cr
     The default is semi-automatically calculated from the number of predictors
     but may need adjusting.
}
\item{thresh}{
     Forward stepping threshold.
     Default is \code{0.001}.
     This is one of the arguments used to decide when forward stepping
     should terminate:
     the forward pass terminates if adding a term changes RSq by less than \code{thresh}.
     See \dQuote{\emph{Termination conditions for the forward pass}} in the vignette.
}
\item{minspan}{
    Minimum number of observations between knots.
    (This increases resistance to runs of correlated noise in the input data.)\cr
    The default \code{minspan=0} is treated specially and
    means calculate the \code{minspan} internally, as per
    Friedman's MARS paper section 3.8 with \eqn{alpha} = 0.05.
    Set \code{trace>=2} to see the calculated value.\cr
    Use \code{minspan=1} and \code{endspan=1} to consider all x values.\cr
    Negative values of \code{minspan} specify the maximum number of knots
    per predictor.  These will be equally spaced.
    For example, \code{minspan=-3} allows three evenly spaced knots for each predictor.
    As always, knots that fall in the endzones specified by \code{endspan} will be ignored.
}
\item{endspan}{
    Minimum number of observations before the first and after the final knot.\cr
    The default \code{endspan=0} is treated specially and
    means calculate the \code{endspan} internally, as per
    the MARS paper equation 45 with \eqn{alpha} = 0.05.
    Set \code{trace>=2} to see the calculated value.
\cr
Be wary of reducing \code{endspan}, especially if you plan to make
predictions beyond or near the limits of the training data.
Overfitting near the edges of training data is much more
likely with a small \code{endspan}.
The model's \code{RSq} and \code{GRSq} won't indicate when this
overfitting is occurring.
(A \code{\link[plotmo]{plotmo}} plot can help: look for sharp hinges at the
edges of the data).  See also the \code{Adjust.endspan} argumen.
}
\item{newvar.penalty}{
     Penalty for adding a new variable in the forward pass
     (Friedman's \eqn{gamma}, equation 74 in the MARS paper).
     Default is \code{0}, meaning no penalty for adding a new variable.
     Useful non-zero values typically range from about \code{0.01} to \code{0.2}
     and sometimes higher ---
     you will need to experiment.\cr
     A word of explanation.  With the default \code{newvar.penalty=0},
     if two variables have nearly the same effect (e.g. they are
     collinear), at any step in the forward pass \code{earth} will
     arbitrarily select one or the other (depending on noise in the sample).
     Both variables can appear in the
     final model, complicating model interpretation.  On the other hand
     with a non-zero \code{newvar.penalty}, the forward pass will be
     reluctant to add a new variable --- it will rather try to use a
     variable already in the model, if that does not affect RSq too much.
     The resulting final model may be easier to interpret, if you are lucky.
     There will often be a small performance hit (a worse GCV).
}
\item{fast.k}{
     Maximum number of parent terms considered at each step of the forward pass.
     (This speeds up the forward pass. See the Fast MARS paper section 3.0.)\cr
     Default is \code{20}.
     A value of \code{0} is treated specially
     (as being equivalent to infinity), meaning no Fast MARS.
     Typical values, apart from \code{0}, are \code{20}, \code{10}, or \code{5}.\cr
     In general, with a lower \code{fast.k} (say \code{5}), \code{earth} is faster;
     with a higher \code{fast.k}, or with \code{fast.k} disabled (set to \code{0}),
     \code{earth} builds a better model.
     However, because of random variation this general rule often doesn't apply.
}
\item{fast.beta}{
     Fast MARS ageing coefficient, as described in the
     Fast MARS paper section 3.1.
     Default is \code{1}.
     A value of \code{0} sometimes gives better results.
}
\item{linpreds}{
     Index vector specifying which predictors should enter linearly, as in \code{\link{lm}}.
     The default is \code{FALSE}, meaning predictors enter
     in the standard MARS fashion, i.e., in hinge functions.
\cr
\cr
     The linpreds argument does not specifiy that a predictor
     \emph{must} enter the model; only that if it enters, it enters
     linearly.  See \dQuote{\emph{The linpreds argument}} in the
     \href{../doc/earth-notes.pdf}{vignette}.
\cr
     See also the \code{Auto.linpreds} argument below (which describes how
     \code{earth} will \emph{automatically} treat a predictor as linear
     under certain conditions).
\cr
\cr
     Details:
     A predictor's index in \code{linpreds} is the column number in the input matrix \code{x}
     (after factors have been expanded).
\cr
     \code{linpreds=TRUE} makes all predictors enter linearly (the \code{TRUE} gets recycled).
\cr
     \code{linpreds} may be a character vector e.g.
     \code{linpreds=c("wind", "vis")}.  Note: \code{\link{grep}} is used
     for matching. Thus \code{"wind"} will match all variables that have
     \code{"wind"} in their names. Use \code{"^wind$"} to match only the
     variable named \code{"wind"}.
}
\item{allowed}{
     Function specifying which predictors can interact and how.
     Default is NULL, meaning all standard MARS terms are allowed.\cr
     During the forward pass, \code{earth} calls the \code{allowed} function
     before considering a term for inclusion; the term can go into the
     model only if the \code{allowed} function returns \code{TRUE}.
     See \dQuote{\emph{The allowed argument}} in the vignette.
\cr\cr
     \bold{The following arguments are for the pruning pass.}
}
\item{pmethod}{
     Pruning method.
     One of: \code{backward none exhaustive forward seqrep cv}.\cr
     Default is \code{"backward"}.\cr
     Specify \code{pmethod="cv"} to use cross-validation to select the number of terms.
     This selects the number of terms that gives the maximum
     mean out-of-fold RSq on the fold models.
     Requires the \code{nfold} argument.\cr
     Use \code{"none"} to retain all the terms created by the forward pass.\cr
     If \code{y} has multiple columns, then only \code{"backward"} or \code{"none"}
     is allowed.\cr
     Pruning can take a while if \code{"exhaustive"} is chosen and
     the model is big (more than about 30 terms).
     The current version of the \code{\link[leaps]{leaps}} package
     used during pruning does not allow user interrupts
     (i.e., you have to kill your R session to interrupt;
     in Windows use the Task Manager or from the command line use \code{taskkill}).
}
\item{nprune}{
     Maximum number of terms (including intercept) in the pruned model.
     Default is NULL, meaning all terms created by the forward pass
     (but typically not all terms will remain after pruning).
     Use this to enforce an upper bound on the model size (that is less than \code{nk}),
     or to reduce exhaustive search time with \code{pmethod="exhaustive"}.
\cr\cr
     \bold{The following arguments are for cross validation.}
}
\item{nfold}{
     Number of cross-validation folds.
     Default is \code{0}, no cross validation.
     If greater than \code{1}, \code{earth} first builds a standard model as usual with all the data.
     It then builds \code{nfold} cross-validated models,
     measuring R-Squared on the out-of-fold (left out) data each time.
     The final cross validation R-Squared (\code{CVRSq}) is the mean of these
     out-of-fold R-Squareds.\cr
     The above process of building \code{nfold} models is repeated
     \code{ncross} times (by default, once).
     Use \code{trace=.5} to trace cross-validation.\cr
     Further statistics are calculated if \code{keepxy=TRUE} or
     if a binomial or poisson model (specified with the \code{glm} argument).
     See \dQuote{\emph{Cross validation}} in the vignette.
}
\item{ncross}{
    Only applies if \code{nfold>1}.
    Number of cross-validations.  Each cross-validation has \code{nfold} folds.
    Default \code{1}.
}
\item{stratify}{
     Only applies if \code{nfold>1}.
     Default is \code{TRUE}.
     Stratify the cross-validation samples so that
     an approximately equal number of cases with a non-zero response
     occur in each cross validation subset.
     So if the response \code{y} is logical, the \code{TRUE}s will be spread
     evenly across folds.
     And if the response is a multilevel factor, there will be an
     approximately equal number of each factor level in each fold
     (because a multilevel factor response gets expanded to columns of zeros and ones,
     see \dQuote{\emph{Factors}} in  the vignette).
     We say \dQuote{approximately equal} because the number of occurrences of a factor
     level may not be exactly divisible by the number of folds.
\cr\cr
     \bold{The following arguments are for variance models.}
}
\item{varmod.method}{
Construct a variance model.
For details, see \code{\link{varmod}} and the vignette
\dQuote{\href{../doc/earth-varmod.pdf}{Variance models in earth}}.
Use \code{trace=.3} to trace construction of the variance model.
\cr
This argument requires \code{nfold} and \code{ncross}. (We suggest at least \code{ncross=30}
here to properly calculate the variance of the errors --- although
you can use a smaller value, say \code{3}, for debugging.)\cr
The \code{varmod.method} argument should be one of\cr
\bold{\code{"none"}} Default.  Don't build a variance model.
\cr
\bold{\code{"const"}} Assume homoscedastic errors.
\cr
\bold{\code{"lm"}} Use \code{\link{lm}} to estimate standard deviation as a
function of the predicted response.
\cr
\bold{\code{"rlm"}} Use \code{\link{rlm}}.
\cr
\bold{\code{"earth"}} Use \code{\link{earth}}.
\cr
\bold{\code{"gam"}} Use \code{gam}.
This will use either \code{\link[gam]{gam}}
or the \code{mgcv} package, whichever is loaded.
\cr
\bold{\code{"power"}} Estimate standard deviation as
\code{intercept + coef * predicted.response^exponent},
where
\code{intercept}, \code{coef}, and \code{exponent} will be estimated by \code{\link{nls}}.
This is equivalent to \code{varmod.method="lm"} except that \code{exponent} is
automatically estimated instead of being held at the value
set by the \code{varmod.exponent} argument.
\cr
\bold{\code{"power0"}} Same as \code{"power"} but no intercept (offset) term.
\cr
\bold{\code{"x.lm"}},
\bold{\code{"x.rlm"}},
\bold{\code{"x.earth"}},
\bold{\code{"x.gam"}}
Like the similarly named options above,
but estimate standard deviation by regressing on the predictors \code{x}
(instead of the predicted response).
A current implementation restriction is that \code{"x.gam"}
allows only models with one predictor (\code{x} must have only one column).
}
\item{varmod.exponent}{
Power transform applied to the rhs before regressing the
absolute residuals with the specified \code{varmod.method}.
Default is \code{1}.\cr
For example, with \code{varmod.method="lm"}, if you expect the
standard deviance to increase linearly with the mean response, use
\code{varmod.exponent=1}.
If you expect the standard deviance to increase with the square root
of the mean response, use
\code{varmod.exponent=.5}
(where negative response values will be treated as \code{0},
and you will get an error message if more than 20\% of them are negative).
}
%\item{varmod.lambda}{
% The variance model regresses on \code{squared.resids^(varmod.lambda/2)}.
% The default \code{varmod.lambda} is \code{1},
% that is, regress on the absolute value of the residuals.
% }
\item{varmod.conv}{
Convergence criterion for the Iteratively Reweighted Least Squares used
when creating the variance model.\cr
Iterations stop when the mean value of the coefficients of the
residual model change by less than \code{varmod.conv}
percent.
Default is \code{1} percent.\cr
Negative values force the specified number of iterations,
e.g. \code{varmod.conv=-2} means iterate twice.\cr
Positive values are ignored for \code{varmod="const"}
and also currently ignored for \code{varmod="earth"}
(these are iterated just once, the same as using \code{varmod.conv=-1}).
}
\item{varmod.clamp}{
The estimated standard deviation of the main model errors
is forced to be at least a small positive value,
which we call \code{min.sd}.
This prevents negative or absurdly small estimated standard deviations.
Clamping takes place in \code{predict.varmod}, which is called
by \code{predict.earth} when estimating prediction intervals.
The value of \code{min.sd} is determined when building the variance
model as \code{min.sd = varmod.clamp * mean(sd(training.residuals))}.
The default \code{varmod.clamp} is \code{0.1}.
}
\item{varmod.minspan}{
Only applies when \code{varmod.method="earth"} or \code{"x.earth"}.
This is the \code{minspan} used in the internal call to \code{earth}
when creating the variance model (not the main \code{earth} model).\cr
Default is \code{-3}, i.e., three evenly spaced knots per predictor.
Residuals tend to be very noisy, and allowing only this small
number of knots helps prevent overfitting.
\cr\cr
    \bold{The following arguments are for internal or advanced use.}
}
\item{Object}{
Earth object to be updated, for use by \code{\link{update.earth}}.
}
\item{Scale.y}{
\code{\link[=scale]{Scale}} \code{y} internally in the forward pass
for better numeric stability.
This is invisible to the user, up to numerical differences.
Scaling here means subtract the mean and divide by the standard
deviation.  Default is \code{NCOL(y)==1},
i.e., scale \code{y} unless \code{y} has multiple columns.

}
\item{Adjust.endspan}{
In interaction terms, \code{endspan} gets multiplied by this value.
This reduces the possibility of an overfitted interaction term
supported by just a few cases on the boundary of the predictor space
(as sometimes seen in our simulation studies).\cr
The default is \code{2}.
Use \code{Adjust.endspan=1} for compatibility with old
versions of \code{earth}.
}
\item{Auto.linpreds}{
Default is \code{TRUE}, which works as follows
(see \href{../doc/Auto-linpreds-example.pdf}{example}):\cr
At any step in the forward pass, if earth discovers that the best knot
for the best predictor is at the predictor minimum (in the
training data),
then earth adds the predictor to the model as a linear \dQuote{basis
function} (with no hinge).
Compare the following basis functions (printed in bold)
for an example such predictor \code{x}:\cr
\code{Auto.linpreds=TRUE} (default): \bold{\code{x}}\cr
\code{Auto.linpreds=FALSE}: \bold{\code{max(x-99, 0)}} where
\code{99} is the minimum \code{x} in the training data.\cr
Using \code{Auto.linpreds=FALSE} always forces a knot, even when the
knot is at the minimum value of the variable.
This ensures that the basis functions are always expressed as hinge functions
(and will always be non-negative).\cr
Note that \code{Auto.linpreds} affects only how the model behaves \emph{outside}
the training data.
Thus \code{predict.earth} will
make the same predictions from the training data, regardless
of whether the earth model was built with \code{Auto.linpreds} set
\code{TRUE} or \code{FALSE}
(up to possible differences in the size of the model caused by
different GCVs because of the different forms of the terms).
}
\item{Force.weights}{
Default is \code{FALSE}.
For testing the \code{weights} argument.
Force use of the code for handling weights in the \code{earth} code,
even if \code{weights=NULL} or all the weights are the same.
This will not necessarily generate an identical model,
primarily because the non-weighted code requires some tests for
numerical stability that can sometimes affect knot selection.
}
\item{Use.beta.cache}{
Default is \code{TRUE}.
Using the \dQuote{beta cache} takes a little more memory but is faster
(by 20\% and often much more for large models).
The beta cache uses \code{nk * nk * ncol(x) * sizeof(double)} bytes.
(The beta cache is an innovation in this implementation of MARS
and does not appear in Friedman's papers. It is not related to
the \code{fast.beta} argument.  Certain regression coefficients
in the forward pass can be saved and re-used, thus
saving recalculation time.)
}
\item{Force.xtx.prune}{
Default is \code{FALSE}.
This argument pertains to subset evaluation in the pruning pass.
By default,
if \code{y} has a single column then \code{earth} calls the \code{\link[leaps]{leaps}} routines;
if \code{y} has multiple columns then \code{earth} calls \code{EvalSubsetsUsingXtx}.
The \code{leaps} routines are numerically more stable
but do not support multiple responses
(\code{leaps} is based on the QR decomposition and
\code{EvalSubsetsUsingXtx} is based on the inverse of X'X).
Setting \code{Force.xtx.prune=TRUE} forces use of \code{EvalSubsetsUsingXtx}, even
if \code{y} has a single column.
}
\item{Get.leverages}{
Default is \code{TRUE} unless the model has more than 100 thousand cases.
The leverages are the diagonal hat values for the linear regression of
\code{y} on \code{bx}.
(The leverages are needed only for certain model checks, for example
when \code{plotres} is called with \code{versus=4}).
\cr
Details:
This argument was introduced to reduce peak memory usage.
When \code{n >> p}, memory use peaks when \code{earth} is
calculating the leverages.
% If memory thrashing occurs,  model-building will be slow.
% The hard drive light will flicker continuously as the operating
% system pages memory.
}
\item{Exhaustive.tol}{
Default \code{1e-10}.
Applies only when \code{pmethod="exhaustive"}.
If the reciprocal of the condition number of \code{bx}
is less than \code{Exhaustive.tol}, \code{earth} forces \code{pmethod="backward"}.
See \dQuote{\emph{XHAUST returned error code -999}} in the vignette.
}
\item{\dots}{
     Dots are passed on to \code{earth.fit}.
 }
}
\value{
  An S3 model of class \code{"earth"}.
  See \code{\link{earth.object}} for a complete description.
}
\author{
  Stephen Milborrow, derived from \code{mda::\link[mda]{mars}}
  by Trevor Hastie and Robert Tibshirani.

  The approach used for GLMs was motivated by work done by
  Jane Elith and John Leathwick
  (a representative paper is given below).

  The \code{\link{evimp}} function uses ideas from Max Kuhn's \code{caret} package
  \url{https://CRAN.R-project.org/package=caret}.

  Parts of Thomas Lumley's \code{\link[leaps]{leaps}} package have been
  incorporated into \code{earth}, so \code{earth} can directly access
  Alan Miller's Fortran functions without going through hidden functions
  in the \code{leaps} package.
}
\references{
  The primary references are the Friedman papers.
  Readers may find the MARS section in Hastie, Tibshirani,
  and Friedman a more accessible introduction.
  The Wikipedia article is recommended for an elementary introduction.
  Faraway takes a hands-on approach,
  using the \code{\link[=ozone1]{ozone}} data to compare \code{mda::mars} with other techniques.
  (If you use Faraway's examples with \code{earth} instead of \code{mars}, use \code{$bx}
  instead of \code{$x}, and check out the book's errata.)
  Friedman and Silverman is recommended background reading for the MARS paper.
  Earth's pruning pass uses code from the \code{\link[leaps]{leaps}} package
  which is based on techniques in Miller.

  Faraway (2005) \emph{Extending the Linear Model with R}
  \url{http://www.maths.bath.ac.uk/~jjf23}

  Friedman (1991) \emph{Multivariate Adaptive Regression Splines (with discussion)}
  Annals of Statistics 19/1, 1--141
  % \code{https://statistics.stanford.edu/research/multivariate-adaptive-regression-splines}\cr
  \url{http://projecteuclid.org/euclid.aos/1176347963}\cr
  \doi{10.1214/aos/1176347963}

  Friedman (1993) \emph{Fast MARS}
  Stanford University Department of Statistics, Technical Report 110
  % use \code and not \url below else get libcurl error 35 from cran check
  % note also that there is no DOI for this paper as far as I can tell
  \code{https://statistics.stanford.edu/research/fast-mars}

  Friedman and Silverman (1989)
  \emph{Flexible Parsimonious Smoothing and Additive Modeling}
  Technometrics, Vol. 31, No. 1.
  \url{http://links.jstor.org/sici?sici=0040-1706\%28198902\%2931\%3A1\%3C3\%3AFPSAAM\%3E2.0.CO\%3B2-Z}

  Hastie, Tibshirani, and Friedman (2009) \emph{The Elements of Statistical Learning (2nd ed.)}
  \url{http://web.stanford.edu/~hastie/pub.htm}

  Leathwick, J.R., Rowe, D., Richardson, J., Elith, J., & Hastie, T. (2005)
  \emph{Using multivariate adaptive regression splines to predict the distributions
  of New Zealand's freshwater diadromous fish} Freshwater Biology, 50, 2034-2052
  \url{http://web.stanford.edu/~hastie/pub.htm},
  \url{http://www.botany.unimelb.edu.au/envisci/about/staff/elith.html}

  Miller, Alan (1990, 2nd ed. 2002) \emph{Subset Selection in Regression}
  \url{http://wp.csiro.au/alanmiller/index.html}

  Wikipedia article on MARS
  \url{http://en.wikipedia.org/wiki/Multivariate_adaptive_regression_splines}
}
\seealso{
  Start with \code{\link{summary.earth}}, \code{\link{plot.earth}},
  \code{\link{evimp}}, and \code{\link[plotmo]{plotmo}}.

Please see the main package vignette
\dQuote{\href{../doc/earth-notes.pdf}{Notes on the earth package}}.
The vignette can also be downloaded from
\url{http://www.milbo.org/doc/earth-notes.pdf}.

The vignette
\dQuote{\href{../doc/earth-varmod.pdf}{Variance models in earth}}
is also included with the package.
It describes how to build variance models and
generate prediction intervals for \code{earth} models.
}
\examples{
earth.mod <- earth(Volume ~ ., data = trees)
plotmo(earth.mod)
summary(earth.mod, digits = 2, style = "pmax")
}
\keyword{smooth}
\keyword{models}
\keyword{regression}
