---
title: "Perplexity and Intrinsic Dimensionality"
subtitle: "1. Theory"
date: "March 20, 2018"
output:
  html_document:
    theme: cosmo
    toc: true
    toc_float:
      collapsed: false
---

Up: [Documentation Home](https://jlmelville.github.io/smallvis/).

Much of the following discussion is lifted nearly verbatim from the documentation to
the fore-runner to `smallvis`,
[sneer](http://jlmelville.github.io/sneer/dimensionality.html), except I've tried
to be a lot more rigorous and consistent with the nomenclature and symbols. For
a practical application of all of this, see 
[part 2](https://jlmelville.github.io/smallvis/idp.html).

## Introduction

Choosing the correct perplexity value for an embedding is an open research
problem. The usual approach is to plump for a value around 30. The
[How to Use t-SNE Effectively](http://distill.pub/2016/misread-tsne/) web site
demonstrates that multiple perplexity values can give quite different results
even for simple synthetic datasets, although to be fair 
my own work on 
[perplexity](https://jlmelville.github.io/smallvis/perplexity.html) suggests
that a lot of the sensitivity at low perplexity is due to random initialization.

Lee and co-workers' work on
[multi-scale neighbor embedding](https://dx.doi.org/10.1016/j.neucom.2014.12.095)
attempts to get around this problem by considering multiple perplexities during
optimization. Part of their discussion involves making use of the concept of 
intrinsic dimensionality of the input data and relating that to perplexity.

## Some Definitions

This discussion of dimensionality involves the input probabilities, so we will 
assume the use of gaussian similarities, $v_{ij}$, which represents
the similarity between two observations in the dataset, $i$ and $j$.

$$
v_{ij} = \exp(-\beta_{i} r_{ij}^{2})
$$

where $r_{ij}^{2}$ is the squared distance between point $i$ and $j$ in the 
input space and $\beta_i$ is the precision (inverse of the bandwidth) 
associated with point $i$. $\beta_i$ is the parameter that is directly optimized
to reproduce a specific perplexity during the t-SNE initialization.

We use a point-wise normalization to define a probability, $p_{j|i}$:

$$
p_{j|i} = \frac{v_{ij}}{\sum_{k} v_{ik}}
$$

The Shannon Entropy, $H_U$, of the probability distribution with perplexity $U$ 
is:

$$
H_U = -\sum_{j} p_{j|i} \log p_{j|i}
$$

and the perplexity, $U$, in units of nats is:

$$
U = \exp(H_U)
$$

Lee and co-workers show that, assuming the data is uniformly distributed, the 
relationship between the precision, perplexity and the correlation dimension
around point $i$, $\delta_{i, U}$, is:

$$
U \propto \beta_i^{-\left(\delta_{i, U} / 2 \right)}
$$

This suggests that if you did a log-log plot of the perplexity against
the precision for multiple perplexities, the graph should be linear and the
gradient would be $-\delta_{i, U}/2$.

where $\beta_{i, U}$ is the precision parameter for the gaussian function which
generates the $i$th similarity/weight for some perplexity, $U$.

## Estimating Intrinsic Dimensionality

For estimating the intrinsic dimensionality of the input data, Lee and
co-workers give the following scheme:

* Calculate a "soft" correlation dimension, $\delta_{i,U}$ at each point $i$ for a 
given perplexity, $U$. The term "soft" here refers to the fact that the
calculation uses the Gaussian input weights as usually calculated in t-SNE and
related methods, rather than the hard sphere usually used in correlation
dimension calculations. The calculation itself is based on one-sided finite
differences, which is convenient for multi-scaled neighbor embedding because you
calculate the input probabilities for multiple perplexities:

$$
\delta_{i,U} = \frac{2}{\log_{2}\beta_{i,U} - \log_{2}\beta_{i,U^{+}}}
$$

$\beta_{i,U}$ is the Gaussian function precision for input function
associated with point $i$ and perplexity $U$. The precision is the parameter
that is directly calibrated when trying to create probabilities at the target
perplexity during the t-SNE perplexity calibration procedure. $U^{+}$ just
represents some other (larger) perplexity. In the multi-scaling procedure,
a standard set of perplexities (increasing in powers of 2) is used, so $U^{+}$
is whatever the next-highest perplexity in the list is.

* The correlation dimension for the entire dataset and a specific perplexity,
$\hat{\delta}_{U}$, is then calculated as the mean average of all the point-wise 
estimates.

$$
\hat{\delta}_{U} = \frac{1}{N} \sum_{i}^{N} \delta_{i,U}
$$

* The maximum value that the mean correlation dimension achieves is the estimate
of the intrinsic dimensionality, $D$:

$$D = \max_{U} \hat{\delta}_{U}$$

## An Analytical Expression for Correlation Dimension

The purpose of the rest of this document is to derive an analytical expression
for the per-point correlation dimension, $\delta_{i,U}$ rather than a finite
difference, so you only need the current perplexity calibration to generate the
correlation dimension. The calculation of the of $\hat{\delta}_{U}$ and $D$
remains the same.

I should say that the following is something I just came up with myself. It's
*not* in the paper by Lee and co-workers, so take it for what it's worth
(probably not much), but when we put this into practice in the 
[next part](https://jlmelville.github.io/smallvis/idp.html) of this series,
it does produce results quite close to the finite difference values from the
multi-scaling JSE paper.

The gradient of the log-log plot of the perplexity against the precision is:

$$
\frac{\partial \log U}{\partial \log \beta_i} = 
\beta_i \frac{\partial H_U}{\partial \beta_i} = 
-\frac{\delta_{i, U}}{2} 
$$

so

$$
\delta_{i, U} = -2 \beta_{i, U} \frac{\partial H_U}{\partial \beta_{i, U}} 
$$

For the rest of this discussion, we're doing all these calculations at a fixed
perplexity $U$, so for clarity I will drop the $U$ subscript on $\beta_{i, U}$
and $H_U$.

We need to find an expression for $\partial H / \partial \beta_{i}$.
This calls for use of the chain rule. It's going to be tedious, but not that
difficult.

To start, we can write the derivative as:

$$
 \frac{\partial H}{\partial \beta_i} = 
 \sum_{jklm} 
 \frac{\partial H}{\partial p_{k|j}}
 \frac{\partial p_{k|j}}{\partial v_{lm}}
 \frac{\partial v_{lm}}{\partial \beta_{i}}
$$

$\partial v_{lm} / \partial \beta_{i} = 0$ unless $l = i$. Additionally, due to 
the point-wise normalization, $\partial p_{k|j} / \partial v_{lm} = 0$ unless 
$j = l$, allowing us to simplify the summation to:

$$
 \frac{\partial H}{\partial \beta_i} = 
 \sum_{km} 
 \frac{\partial H}{\partial p_{k|i}}
 \frac{\partial p_{k|i}}{\partial v_{im}}
 \frac{\partial v_{im}}{\partial \beta_{i}}
$$

Now let us regroup that double summation into two single summations, and also
rename $m$ to $j$:

$$
 \frac{\partial H}{\partial \beta_i} = 
 \sum_{j}
 \left[
   \sum_{k} 
   \frac{\partial H}{\partial p_{k|i}}
   \frac{\partial p_{k|i}}{\partial v_{ij}}
 \right]
 \frac{\partial v_{ij}}{\partial \beta_{i}}
$$

The full details of how to derive the expression of the 
[gradient of the weight normalization](http://jlmelville.github.io/sneer/gradients.html) 
are available elsewhere, but we shall jump straight to inserting the result for 
$\partial p_{k|i} / \partial v_{ij}$
to get:

$$
 \frac{\partial H}{\partial \beta_i} = 
 \sum_{j}
 \frac{1}{V_{i}}
 \left[
   \frac{\partial H}{\partial p_{j|i}}
   -\sum_{k} 
   \frac{\partial H}{\partial p_{k|i}}
   p_{k|i}
 \right]
 \frac{\partial v_{ij}}{\partial \beta_{i}}
$$
where $V_{i}$ is:

$$
V_{i} = \sum_{j} v_{ij}
$$

The gradient of the Shannon Entropy with respect to the probability is:

$$
\frac{\partial H}{\partial p_{j|i}} =
- \log \left( p_{j|i} \right) - 1
$$
substituting into the expression in square brackets:

$$
 \left[
   \frac{\partial H}{\partial p_{j|i}}
   -\sum_{k} 
   \frac{\partial H}{\partial p_{k|i}}
   p_{k|i}
 \right]
 = 
 \left[
- \log \left( p_{j|i} \right) - 1
   -\sum_{k} 
   \left\{
   - \log \left( p_{k|i} \right) - 1
   \right\} 
   p_{k|i}
 \right] \\
  = 
 \left[
   - \log \left( p_{j|i} \right) - 1
   +\sum_{k} 
   p_{k|i} \log \left( p_{k|i} \right)
   +\sum_{k} 
   p_{k|i}
 \right]
$$
and because $\sum_k p_{k|i} = 1$, we eventually get to:
$$
 \left[
   \frac{\partial H}{\partial p_{j|i}}
   -\sum_{k} 
   \frac{\partial H}{\partial p_{k|i}}
   p_{k|i}
 \right]
 = 
 \left[
- \log \left( p_{j|i} \right) - H
 \right] \\
 =
 -\left[
  \log \left( p_{j|i} \right) + H
 \right]
$$
At this point:

$$
 \frac{\partial H}{\partial \beta_i} = 
 \sum_{j}
 -\frac{1}{V_{i}}
 \frac{\partial v_{ij}}{\partial \beta_{i}}
  \left[
  \log \left( p_{j|i} \right) + H
 \right]
$$

which leads to:

$$
\delta_{i, U} = \frac{2 \beta_i}{V_i}
 \sum_{j}
 \frac{\partial v_{ij}}{\partial \beta_{i}}
  \left[
  \log \left( p_{j|i} \right) + H
 \right]
$$

The gradient of the exponential weight with respect to the precision parameter, $\beta_i$,
is:

$$
\frac{\partial v_{ij}}{\partial \beta_{i}}
=
-r_{ij}^2 v_{ij}
$$

Substituting these two sub expressions into the total gradient, we are left with:
$$
 \frac{\partial H}{\partial \beta_i} = 
 \sum_{j}
 \frac{r_{ij}^2 v_{ij}}{V_{i}}
 \left[
 \log \left( p_{j|i} \right) + H
 \right]
$$

Additionally, $p_{j|i} = v_{ij} / V_{i}$, so the final expression for the 
gradient is:

$$
 \frac{\partial H}{\partial \beta_i} = 
 \sum_{j}
 r_{ij}^2 p_{j|i}
 \left[
 \log \left( p_{j|i} \right) + H
 \right]
$$

The only thing to left to do is to multiply this expression by $-2 \beta_{i}$
to get this expression for the correlation dimension:

$$
\delta_{i,U} = -2 \beta_i \sum_j r^2_{ij} p_{j|i} \left[\log\left(p_{j|i}\right) + H\right]
$$

which is useful if you've carried out the perplexity-based calibration
on the input weights, as you have already calculated $p_{j|i}$, $H$ and 
$\beta_i$. 


### Un-normalized weights

If you'd rather think in terms of the un-normalized weights and the distances 
only, the usual expression for Shannon entropy can be rewritten in terms of
weights as:

$$
H = \log V_{i} -\frac{1}{V_i} \left( \sum_j v_{ij} \log v_{ij} \right) 
$$

and that can be combined with:

$$
\log \left( p_{j|i} \right) = \log \left( v_{ij} \right) - \log \left(V_i \right) 
$$

to give:

$$
\delta_{i,U} = \frac{2 \beta_i}{V_i^2}
 \sum_{j}
 \frac{\partial v_{ij}}{\partial \beta_{i}}
  \left(
   V_i \log v_{ij} - \sum_k v_{ik} \log v_{ik}
 \right)
$$

*10 October 2025*: pausing to scold myself from seven years in the future, as
I appear to have pulled out a factor of $1/V_i$ in the hope of it simplifying
later, but it doesn't. You'll see I just multiply it back in again later.

We can also express the relation between the weight and the squared distance as:

$$
\log \left( v_{ij} \right) = -\beta_i r_{ij}^2
$$
and the gradient with respect to the precision as:

$$
\frac{\partial v_{ij}}{\partial \beta_{i}}
=
-r_{ij}^2 v_{ij} = \frac{v_{ij} \log v_{ij}}{\beta_i}
$$

and the Shannon entropy expression as:

$$
H = \log V_i + \frac{\beta_i}{V_i} \sum_j r_{ij}^2 v_{ij}
$$

With all that, you can eventually get to two equivalent expressions for $\delta_{i, U}$:

$$
\delta_{i,U} = \frac{2}{V_i}
\left\{
\sum_j v_{ij} \left[ \log \left( v_{ij} \right) \right] ^ 2
-\frac{1}{V_i} \left[ \sum_j v_{ij} \log \left( v_{ij} \right) \right]^2
\right\}
$$
$$
\delta_{i, U} = \frac{2 \beta_i^2}{V_i}
\left[
\sum_j r_{ij}^4 v_{ij}
-\frac{1}{V_i} \left( \sum_j r_{ij}^2 v_{ij} \right)^2
\right]
$$

The first one only requires the $\mathbf{V}$ matrix, although as you've been 
carrying out lots of exponential operations to generate the weights, it seems a
pity to have to then carry out lots of expensive log calculations, in which case
the second expression might be better, but which requires the squared distance
matrix and $\beta_i$ also.

With this expression for the correlation dimension for point $i$, we proceed as
before: the estimate for the entire dataset, $\hat{\delta}_U$ is the average
over all $\delta_{i, U}$. The intrinsic dimensionality, $D$ is then estimated as
the maximum value that $\hat{\delta}_U$ reaches when calculated over a range of
perplexities. In the multi-scale JSE paper, perplexities in increasing powers of
2 are used.

For a practical application, see [part 2](https://jlmelville.github.io/smallvis/idp.html).

Up: [Documentation Home](https://jlmelville.github.io/smallvis/).
