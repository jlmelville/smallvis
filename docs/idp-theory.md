---
title: "Perplexity and Intrinsic Dimensionality"
subtitle: "1. Theory"
date: "March 20, 2017"
output:
  html_document:
    theme: cosmo
---

Much of the following discussion is lifted verbatim from the documentation to
the fore-runner to `smallvis`,
[sneer](http://jlmelville.github.io/sneer/dimensionality.html). For a practical
application of all of this, I've applied the 
[Intrinsic Dimensionality Perplexity to t-SNE](https://jlmelville.github.io/smallvis/idp.html).


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

For estimating the intrinsic dimensionality of the input data, they suggest the
following scheme:

* Calculate a "soft" correlation dimension at each point for a given perplexity.
The term "soft" here refers to the fact that the calculation uses the Gaussian
input weights as usually calculated in t-SNE and related methods, rather than
the hard sphere usually used in correlation dimension calculations. The
calculation itself is based on one-sided finite differences, which is convenient
for multi-scaled neighbor embedding because you calculate the input
probabilities for multiple perplexities.

* The correlation dimension for the entire dataset is then calculated as the
mean average of all the point-wise estimates.

* The maximum value that the mean correlation dimension achieves is the estimate
of the intrinsic dimensionality.

Here I'm going to derive an analytical expression for the correlation dimension
with the use of the exponential kernel (exponential with regard to the squared
distances, or gaussian in the raw distances, if you prefer) that only requires
one perplexity calculation.

First, some definitions. This discussion of dimensionality involves the input
probabilities, so we will assume the use of gaussian similarities, not anything
exotic like the t-distribution or other heavy-tailed functions:

$$
w_{ij} = \exp(-\beta_i d_{ij}^{2})
$$
where $d_{ij}^{2}$ is the squared distance between point $i$ and $j$ in the 
input space and $\beta_i$ is the precision (inverse of the bandwidth) 
associated with point $i$. We use a point-wise normalization to define a 
probability, $p_{j|i}$:

$$
p_{j|i} = \frac{w_{ij}}{\sum_{k} w_{ik}}
$$
The Shannon Entropy, $H$, of the probability distribution is:
$$
H = -\sum_{j} p_{j|i} \log p_{j|i}
$$
and the perplexity, $U$, in units of nats is:

$$
U = \exp(H)
$$

Lee and co-workers show that, assuming the data is uniformly distributed, the 
relationship between the precision, perplexity and correlation dimension
around point $i$, $D_i$, is:

$$
U \propto \beta_i^{-\left({D_i}/{2}\right)}
$$
This suggests that if you did a log-log plot of the perplexity against
the precision for multiple perplexities, the graph should be linear and the
gradient would be $-D_{i}/2$.

With a multi-scaling approach, a finite difference expression using base 2 
logarithms is given as:
$$
D_i = \frac{2}{\log_2\left(\beta_{U,i}\right)-\log_2\left(\beta_{V,i}\right)}
$$

where $\beta_{U,i}$ is the precision parameter for the gaussian function which
generates the $i$th similarity/weight for some perplexity, $U$.

### An Analytical Expression for Correlation Dimension

The following is something I just came up with myself. It's *not* in the paper 
by Lee and co-workers, so take it for what it's worth (probably not much), but
it does produce 
[result quite close to the finite difference values](https://jlmelville.github.io/smallvis/idp.html).

The gradient of the log-log plot of the perplexity against the precision is:

$$
\frac{\partial \log U}{\partial \log \beta_i} = 
\beta_i \frac{\partial H}{\partial \beta_i} = 
-\frac{D_{i}}{2} 
$$
so

$$
D_i = -2 \beta_i \frac{\partial H}{\partial \beta_i} 
$$

We therefore need to find an expression for $\partial H / \partial \beta_{i}$.
This calls for use of the chain rule. It's going to be tedious, but not that
difficult.

To start, we can write the derivative as:

$$
 \frac{\partial H}{\partial \beta_i} = 
 \sum_{jklm} 
 \frac{\partial H}{\partial p_{k|j}}
 \frac{\partial p_{k|j}}{\partial w_{lm}}
 \frac{\partial w_{lm}}{\partial \beta_{i}}
$$

$\partial w_{lm} / \partial \beta_{i} = 0$ unless $l = i$. Additionally, due to 
the point-wise normalization, $\partial p_{k|j} / \partial w_{lm} = 0$ unless 
$j = l$, allowing us to simplify the summation to:

$$
 \frac{\partial H}{\partial \beta_i} = 
 \sum_{km} 
 \frac{\partial H}{\partial p_{k|i}}
 \frac{\partial p_{k|i}}{\partial w_{im}}
 \frac{\partial w_{im}}{\partial \beta_{i}}
$$
Now let us regroup that double summation into two single summations, and also
rename $m$ to $j$:

$$
 \frac{\partial H}{\partial \beta_i} = 
 \sum_{j}
 \left[
   \sum_{k} 
   \frac{\partial H}{\partial p_{k|i}}
   \frac{\partial p_{k|i}}{\partial w_{ij}}
 \right]
 \frac{\partial w_{ij}}{\partial \beta_{i}}
$$

The full details of how to derive the expression of the 
[gradient of the weight normalization](http://jlmelville.github.io/sneer/gradients.html) 
are available elsewhere, but we shall jump straight to inserting the result for 
$\partial p_{k|i} / \partial w_{ij}$
to get:

$$
 \frac{\partial H}{\partial \beta_i} = 
 \sum_{j}
 \frac{1}{S_{i}}
 \left[
   \frac{\partial H}{\partial p_{j|i}}
   -\sum_{k} 
   \frac{\partial H}{\partial p_{k|i}}
   p_{k|i}
 \right]
 \frac{\partial w_{ij}}{\partial \beta_{i}}
$$
where $S_{i}$ is:

$$
S_{i} = \sum_{j} w_{ij}
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
 \right]
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
 \right]
 =
 -\left[
  \log \left( p_{j|i} \right) + H
 \right]
$$
At this point:

$$
 \frac{\partial H}{\partial \beta_i} = 
 \sum_{j}
 -\frac{1}{S_{i}}
 \frac{\partial w_{ij}}{\partial \beta_{i}}
  \left[
  \log \left( p_{j|i} \right) + H
 \right]
$$

which leads to:

$$
D_{i} = \frac{2 \beta_i}{S_i}
 \sum_{j}
 \frac{\partial w_{ij}}{\partial \beta_{i}}
  \left[
  \log \left( p_{j|i} \right) + H
 \right]
$$

The gradient of the exponential weight with respect to the precision parameter, $\beta_i$,
is:

$$
\frac{\partial w_{ij}}{\partial \beta_{i}}
=
-d_{ij}^2 w_{ij}
$$

Substituting these two sub expressions into the total gradient, we are left with:
$$
 \frac{\partial H}{\partial \beta_i} = 
 \sum_{j}
 \frac{d_{ij}^2 w_{ij}}{S_{i}}
 \left[
 \log \left( p_{j|i} \right) + H
 \right]
$$

Additionally, $p_{j|i} = w_{ij} / S_{i}$, so the final expression for the 
gradient is:

$$
 \frac{\partial H}{\partial \beta_i} = 
 \sum_{j}
 d_{ij}^2 p_{j|i}
 \left[
 \log \left( p_{j|i} \right) + H
 \right]
$$

The only thing to left to do is to multiply this expression by $-2 \beta_{i}$
to get this expression for the correlation dimension:

$$
D_i = -2 \beta_i \sum_j d^2_{ij} p_{j|i} \left[\log\left(p_{j|i}\right) + H\right]
$$

which is useful if you've carried out the perplexity-based calibration
on the input weights, as you have already calculated $p_{j|i}$, $H$ and 
$\beta_i$. 


### Un-normalized weights

If you'd rather think in terms of the un-normalized weights and the distances 
only, the usual expression for Shannon entropy can be rewritten in terms of
weights as:

$$
H = \log S_{i} -\frac{1}{S_i} \left( \sum_j w_{ij} \log w_{ij} \right) 
$$

and that can be combined with:

$$
\log \left( p_{j|i} \right) = \log \left( w_{ij} \right) + \log \left(S_i \right) 
$$

to give:

$$
D_{i} = \frac{2 \beta_i}{S_i^2}
 \sum_{j}
 \frac{\partial w_{ij}}{\partial \beta_{i}}
  \left(
   S_i \log w_{ij} - \sum_k w_{ik} \log w_{ik}
 \right)
$$

We can also express the relation between the weight and the squared distance as:

$$
\log \left( w_{ij} \right) = -\beta_i d_{ij}^2
$$
and the gradient with respect to the precision as:

$$
\frac{\partial w_{ij}}{\partial \beta_{i}}
=
-d_{ij}^2 w_{ij} = -\frac{w_{ij} \log w_{ij}}{\beta_i}
$$

and the Shannon entropy expression as:

$$
H = \log S_i + \frac{\beta_i}{S_i} \sum_j d_{ij}^2 w_{ij}
$$

With all that, you can eventually get to two equivalent expressions for $D_i$:

$$
D_{i} = \frac{2}{S_i}
\left\{
\sum_j w_{ij} \left[ \log \left( w_{ij} \right) \right] ^ 2
-\frac{1}{S_i} \left[ \sum_j w_{ij} \log \left( w_{ij} \right) \right]^2
\right\}
$$
$$
D_{i} = \frac{2 \beta_i^2}{S_i}
\left[
\sum_j d_{ij}^4 w_{ij}
-\frac{1}{S_i} \left( \sum_j d_{ij}^2 w_{ij} \right)^2
\right]
$$
The first one only requires the $W$ matrix, although as you've been carrying out
lots of exponential operations to generate the weights, it seems a pity to have
to then carry out lots of expensive log calculations, in which case the second
expression might be better, but which requires the squared distance matrix and
$\beta_i$ also.

With this expression for the correlation dimension for point $i$, the estimate
for the entire dataset is the average over all $D_{i}$. The intrinsic
dimensionality is then estimated as the maximum value that the mean
correlation dimension reaches when calculated over a range of perplexities. In
the multi-scale JSE paper, perplexities in increasing powers of 2 are used.

For a practical application, see [Intrinsic Dimensionality Perplexity with t-SNE](https://jlmelville.github.io/smallvis/idp.html).
