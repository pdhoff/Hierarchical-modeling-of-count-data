Hierarchical negative binomial regression
================
Peter Hoff
28 June, 2020

### Summary

Several datasets of current interest include count data from multiple
groups. For example, data available at
<https://github.com/pdhoff/US-counties-C19-data> include time series of
C19 counts for different regions in the US.  
The file `negbinHGLM.r` contains a few functions that provide a Bayesian
fit to a hierarchical negative binomial regression model. An example of
such a model would be an autoregressive model of state-level time-series
of C19 counts.

### Why hierarchical modeling?

It is often the case that the data patterns in one group are similar to,
but not exactly the same as, the data patterns of other groups. Because
of the potential differences, one data analysis approach is to fit a
separate statistical model to the data from each group. However, if the
amount of data from each group is small, the group-specific parameter
estimates, inferences and forecasts will be very noisy. The effects of
noise can be reduced by leveraging the potential similatities among the
groups. This is what hierarchical modeling provides - group-specific
estimates and inferences that make use of potential across-group
similarities.

### The hierarchical negative binomial regression model

Let \(y_{i,j}\) be the \(i\)th count in group \(j\), and let \(x_{i,j}\)
be a vector of predictors for this count. A group-specific negative
binomial regression model is that \[
  y_{i,j} \sim \text{negative binomial}( e^{ x_{i,j}^\top \beta_ j},\lambda_j).
\]

\(\theta\)
