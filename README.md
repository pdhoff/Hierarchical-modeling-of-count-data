Hierarchical negative binomial regression
================
Peter Hoff
28 June, 2020



### Summary

Many datasets of current interest include count data from multiple
groups. For example, data available at
<https://github.com/pdhoff/US-counties-C19-data> include time series of
C19 counts for different regions in the US. The file `negbinHGLM.r`
contains a few functions that provide regression modeling of count data
from multiple groups. Specifically, the main function `negbinHGLM` fits
a hierarchical negative binomial regression model using a Markov chain
Monte Carlo approximation. An example use of this code includes fitting
an autoregressive model of state-level time-series of C19 counts.

### Why hierarchical modeling?

It is often the case that the data patterns in one group are similar to,
but not exactly equal to, the data patterns of other groups. Because of
the potential differences, one data analysis approach is to fit a
separate statistical model to the data from each group. However, if the
amount of data from each group is small, the group-specific parameter
estimates, inferences and forecasts will be very noisy. The effects of
noise can be reduced by leveraging the potential similarities among the
groups. This is what hierarchical modeling does - it provides
group-specific estimates and inferences that make use of potential
across-group similarities.

### A hierarchical negative binomial regression model

#### Within-group model

Let \(y_{i,j}\) be the \(i\)th count in group \(j\), and let \(x_{i,j}\)
be a vector of predictors for this count. A group-specific negative
binomial regression model is that \[
  y_{i,j} \sim \text{negative binomial}( \lambda_j e^{\beta_j^\top x_{i,j}},\lambda_j), 
\] independently, within and across groups. Also, the code as currently
written forces the inclusion of an intercept as the first element of the
vector \(x_{i,j}\), even if the user doesn’t include one directly in the
specification of \(x_{i,j}\). So \[
  \beta_j^\top x_{i,j} = \beta_{j,1}  + \beta_{j,2} x_{i,j,2} + \cdots + \beta_{j,p} x_{i,j,p}, 
\] and the first element of each \(\beta_j\)-vector is an intercept
parameter.

The negative binomial distribution is parameterized here as

  - \(E[ y_{i,j}] \equiv \mu_{i,j} = \lambda_j e^{\beta_j^\top x_{i,j}}\)

  - \(V[ y_{i,j} ] = \mu_{i,j}( 1+ \mu_{i,j}/\lambda_j)\).

The \(\lambda_j\) parameter is typically called the *size* parameter -
see the R helpfile for `dnbinom` for more details. The distribution
approaches a Poisson distribution as \(\lambda_j\) increases.

#### Across-group model

A hierarchical model specifies a distribution for the group-specific
parameters \(\{ (\beta_j,\lambda_j),j=1,\ldots,n\}\). The model used
here is that

independently across groups, where the \(W_j\)’s are matrices of known
group-level predictors, and \((\alpha,\Sigma,a,b)\) are unknown
parameters. Simple hierarchical models typically just have \(W_j=I\),
where \(I\) is the identity matrix of dimension equal to the length of
\(\beta_j\). In this case, the \(\alpha\) parameter represents the
across-group average of the regression coefficient vectors
\(\beta_1,\ldots, \beta_n\). Allowing \(W_j\) to vary across groups
permits inclusion of known group-level information the model. In
particular, when the group represent counts in a population, it may be
useful to let the model for the intercept \(\beta_{j,1}\) depend on the
log population in group \(j\).

### Examples

  - [Simulation study](SimStudy.md)

  - [State-level C19 analyses](C19ARModelFittingWithnegbinHGLM.md)

  - [State-level C19 forecasts](C19ARForecast.md)

### Technical comments

  - The MCMC routine in the function `negbinHGLM` makes use of the
    Polya-gamma latent variable representation described in Polson,
    Scott and Windle \[2013\]. To implement this representation in R,
    the package `BayesLogit` needs to have been installed.

  - You need to be careful with the interpretation of the \(\lambda_j\)
    parameter. The way the model is parameterized, \(\lambda_j\)
    controls both the mean and the variance. This creates a bit of
    confounding between the intercept term \(\beta_{j,1}\) and
    \(\lambda_j\), since the \(\lambda_j\) parameters are not very well
    estimated unless the within-group sample size is quite large.

  - I’m using a grid approximation to the full conditional distribution
    of \(\lambda_j\) to update it. However, I also include a
    Metropolis-Hastings updater that you can use instead if you prefer.
    It is implemented with a symmetric proposal distribution for
    \(-\log \lambda_j\). I have found that the mixing of the MCMC in
    terms of the \(\lambda_j\) parameters is generally not very good,
    presumably due to the partial confounding of this parameter with
    \(\beta_{j,1}\). Probably it would be much better to update
    \(\lambda_j\) and \(\beta_{j,1}\) simultaneously. Maybe I’ll get to
    that some weekend.
