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
assume the use of gaussian similarities (or affinities), $v_{j|i}$, which 
represents the similarity between two observations in the dataset, $i$ and $j$.

$$
v_{j|i} = \exp(-\beta_{i} r_{ij}^{2})
$$

where $r_{ij}^{2}$ is the squared distance between point $i$ and $j$ in the 
input space and $\beta_i$ is the precision (inverse of the bandwidth) 
associated with point $i$. $\beta_i$ is the parameter that is directly optimized
to reproduce a specific perplexity during the t-SNE initialization.

I should also stipulate here that when $i = j$, $v_{j|i} = 0$. In addition, I'm 
slightly abusing conditional probability-style notation to indicate that the
affinities are not symmetric: $v_{j|i} \ne v_{i|j}$, because $\beta_i$ and 
$\beta_j$ are different.

We use a point-wise normalization to define a probability, $p_{j|i}$:

$$
p_{j|i} = \frac{v_{j|i}}{\sum_{k} v_{k|i}}
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
 \frac{\partial p_{k|j}}{\partial v_{m|l}}
 \frac{\partial v_{m|l}}{\partial \beta_{i}}
$$

$\partial v_{m|l} / \partial \beta_{i} = 0$ unless $l = i$. Additionally, due to 
the point-wise normalization, $\partial p_{k|j} / \partial v_{m|l} = 0$ unless 
$j = l$, allowing us to simplify the summation to:

$$
 \frac{\partial H}{\partial \beta_i} = 
 \sum_{km} 
 \frac{\partial H}{\partial p_{k|i}}
 \frac{\partial p_{k|i}}{\partial v_{m|i}}
 \frac{\partial v_{m|i}}{\partial \beta_{i}}
$$

Now let us regroup that double summation into two single summations, and also
rename $m$ to $j$:

$$
 \frac{\partial H}{\partial \beta_i} = 
 \sum_{j}
 \left[
   \sum_{k} 
   \frac{\partial H}{\partial p_{k|i}}
   \frac{\partial p_{k|i}}{\partial v_{j|i}}
 \right]
 \frac{\partial v_{j|i}}{\partial \beta_{i}}
$$

The full details of how to derive the expression of the 
[gradient of the weight normalization](http://jlmelville.github.io/sneer/gradients.html) 
are available elsewhere, but we shall jump straight to inserting the result for 
$\partial p_{k|i} / \partial v_{j|i}$
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
 \frac{\partial v_{j|i}}{\partial \beta_{i}}
$$
where $V_{i}$ is:

$$
V_{i} = \sum_{j} v_{j|i}
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
 \frac{\partial v_{j|i}}{\partial \beta_{i}}
  \left[
  \log \left( p_{j|i} \right) + H
 \right]
$$

which leads to:

$$
\delta_{i, U} = \frac{2 \beta_i}{V_i}
 \sum_{j}
 \frac{\partial v_{j|i}}{\partial \beta_{i}}
  \left[
  \log \left( p_{j|i} \right) + H
 \right]
$$

The gradient of the exponential weight with respect to the precision parameter, $\beta_i$,
is:

$$
\frac{\partial v_{j|i}}{\partial \beta_{i}}
=
-r_{ij}^2 v_{j|i}
$$

Substituting these two sub expressions into the total gradient, we are left with:
$$
 \frac{\partial H}{\partial \beta_i} = 
 \sum_{j}
 \frac{r_{ij}^2 v_{j|i}}{V_{i}}
 \left[
 \log \left( p_{j|i} \right) + H
 \right]
$$

Additionally, $p_{j|i} = v_{j|i} / V_{i}$, so the final expression for the 
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
H = \log V_{i} -\frac{1}{V_i} \left( \sum_j v_{j|i} \log v_{j|i} \right) 
$$

and that can be combined with:

$$
\log \left( p_{j|i} \right) = \log \left( v_{j|i} \right) - \log \left(V_i \right) 
$$

to give:

$$
\delta_{i,U} = \frac{2 \beta_i}{V_i^2}
 \sum_{j}
 \frac{\partial v_{j|i}}{\partial \beta_{i}}
  \left(
   V_i \log v_{j|i} - \sum_k v_{k|i} \log v_{k|i}
 \right)
$$

*10 October 2025*: pausing to scold myself from seven years in the future, as
I appear to have pulled out a factor of $1/V_i$ in the hope of it simplifying
later, but it doesn't. You'll see I just multiply it back in again later.

We can also express the relation between the weight and the squared distance as:

$$
\log \left( v_{j|i} \right) = -\beta_i r_{ij}^2
$$
and the gradient with respect to the precision as:

$$
\frac{\partial v_{j|i}}{\partial \beta_{i}}
=
-r_{ij}^2 v_{j|i} = \frac{v_{j|i} \log v_{j|i}}{\beta_i}
$$

and the Shannon entropy expression as:

$$
H = \log V_i + \frac{\beta_i}{V_i} \sum_j r_{ij}^2 v_{j|i}
$$

With all that, you can eventually get to two equivalent expressions for $\delta_{i, U}$:

$$
\delta_{i,U} = \frac{2}{V_i}
\left\{
\sum_j v_{j|i} \left[ \log \left( v_{j|i} \right) \right] ^ 2
-\frac{1}{V_i} \left[ \sum_j v_{j|i} \log \left( v_{j|i} \right) \right]^2
\right\}
$$
$$
\delta_{i, U} = \frac{2 \beta_i^2}{V_i}
\left[
\sum_j r_{ij}^4 v_{j|i}
-\frac{1}{V_i} \left( \sum_j r_{ij}^2 v_{j|i} \right)^2
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

## The Fisher Information Connection

*10 October 2025*: I'm keeping the above expressions here for completeness, but
there are some better and simpler and more interesting ways of writing these
expressions. In the interests of full disclosure, I will say that these
connections were pointed out to me by ChatGPT 5 Pro. Yes I know. The derivations
below were done by hand by your entirely human author, but I would never have
spotted the relation to Fisher Information without the help of AI.

Let's start with the expression that is entirely in terms of substituting in the
different derivatives expressed as un-normalized affinities, $v_{j|i}$, but 
without any further cancellation or re-arrangement:

$$
\delta_{i,U} = \frac{2 \beta_i}{V_i}
 \sum_{j}
 \frac{v_{j|i} \log v_{j|i}}{\beta_i}
  \left(
  \log v_{j|i} - \frac{1}{V_i} \sum_k v_{k|i} \log v_{k|i}
 \right)
$$

Now let's cancel out the $\beta_i$ terms and bring both of the $1/V_i$ terms
inside the summations:

$$
\delta_{i,U} = 2
 \sum_{j}
 \frac{v_{j|i} \log v_{j|i}}{V_i}
  \left(
  \log v_{j|i} - \frac{\sum_k v_{k|i} \log v_{k|i}}{V_i}
 \right)
$$

We can now reintroduce the normalized probabilities, $p_{j|i} = v_{j|i}/V_i$:

$$
\delta_{i,U} = 2
 \sum_{j} p_{j|i} \log v_{j|i}
  \left(
  \log v_{j|i} - \sum_k p_{k|i} \log v_{k|i}
 \right)
$$

and lets expand the parentheses:

$$
\delta_{i,U} = 2 \left[
 \sum_{j} p_{j|i} (\log v_{j|i})^2
 - \left( \sum_k p_{k|i} \log v_{k|i} 
 \right) \left( \sum_{j} p_{j|i} \log v_{j|i} \right)
\right]
$$

and those two summations on the right are now the same, so we don't need $k$
any more:

$$
\delta_{i,U} = 2 \left[
 \sum_{j} p_{j|i} (\log v_{j|i})^2
 - \left( \sum_{j} p_{j|i} \log v_{j|i} 
 \right)^2
\right]
$$

At the cost of mixing probabilities and the un-normalized weights, we now have
a much less cluttered expression. It also looks a lot like a variance. In fact,
it's the weighted variance of $\log v_{j|i}$, where the weights are provided by 
$p_{j|i}$. I admit I had never thought about variances before, but it's
not any more outlandish than a weighted mean.

So a compact way of expressing the correlation dimension is:

$$
\delta_{i,U} = 2 \, \mathrm{Var}_{p_{j|i}}(\log v_{j|i})
$$

where the subscript on the $\mathrm{Var}$ indicates that the variance is
weighted by whatever's in the subscript.

It's also true that subtracting a constant from every value in a variance
doesn't change it, so because $\log p_{j|i} = \log v_{j|i} - \log V_i$, we can
also write:

$$
\delta_{i,U} = 2 \, \mathrm{Var}_{p_{j|i}}(\log p_{j|i})
$$

and you can do the calculation entirely with probabilities.

Finally, because $\log v_{j|i} = -\beta_i r_{ij}^2$, we can *also* write:

$$
\delta_{i,U} = 2 \beta_i^2 \, \mathrm{Var}_{p_{j|i}}(r_{ij}^2)
$$

I'm going to write this out in full again, because hopefully later you are
going to scroll back up to here and go "oh, yeah":

$$
\delta_{i,U} = 2 \beta_i^2 \left[
 \sum_{j} p_{j|i} (r_{ij}^2)^2
 - \left( \sum_{j} p_{j|i} r_{ij}^2 
 \right)^2
\right]
$$

The probability-weighted variance part turns out to be the Fisher Information 
for $\beta_i$. This may not be immediately obvious even if you are familiar with
Fisher Information. To prove it, let's go through the derivation. I am going to
play a bit fast and loose with notation because in all the examples I see, there
are a ton of parentheses and semi-colons which seems unnecessary for the fairly
simple distribution we have here. Apologies to any statisticians reading this.

The likelihood, $L$ of picking point $j$ given point $i$ is just the probability:

$$
L = p_{j|i}
$$

The log-likelihood, $l$ is:

$$
l = \log p_{j|i} = \log v_{j|i} - \log V_i = -\beta_i r_{ij}^2 - \log V_i
$$

The *score*, $s$, is the gradient of the log-likelihood with respect to the
parameter in question, which is $\beta_i$ in this case:

$$
\frac{\partial l}{\partial \beta_i} 
= -\frac{\partial (\beta_i r_{ij}^2)}{\partial \beta_i} - \frac{\partial \log V_i}{\partial \beta_i} \\
=-r_{ij}^2 - \frac{1}{V_i} \sum_k \frac{\partial v_{k|i}}{\partial \beta_i} \\
=-r_{ij}^2 + \frac{1}{V_i} \sum_k v_{k|i} r_{ik}^2  \\
=-r_{ij}^2 + \sum_k p_{k|i} r_{ik}^2 
$$

Let's get ahead of something here and talk about expectation values, as we will
be having to notice them to get back to our old $\delta_{i,U}$ expression. For 
our finite discrete distribution, the expectation involves taking the 
probability-weighted sum over all possible values of $j$. You can already see
something like that in the score above where the second term is a sum over
$k$ of the squared distances, weighted by the respective probabilities.

Under some regularity conditions, the expectation of the score is zero (this is
true of all well-behaved probability distributions, not this specific example).
What are those conditions? They are that the log-likelihood must be 
differentiable with respect to $\beta_i$ (it is) and that the possible values of 
$j$ (the "support") do not depend on $\beta_i$ (they don't). So if we now
explicitly write out the expectation of the score it had better be zero:

$$
E\left[ s \right]
= \sum_j p_{j|i} \left( -r_{ij}^2 + \sum_k p_{k|i} r_{ik}^2 \right) \\
= \left(\sum_j p_{j|i} (-r_{ij}^2) \right) + \left(\sum_k p_{k|i} r_{ik}^2 \sum_j p_{j|i} \right) \\
= -\sum_j p_{j|i} r_{ij}^2 + \sum_k p_{k|i} r_{ik}^2 \\
= 0
$$

That's a relief.

Now we can define the Fisher Information. The Fisher Information, $I$, is the 
variance of the score, which we can write as the sum of two expectations:

$$
I = \mathrm{Var}(s) = E[s^2] - (E[s])^2
$$

But didn't we just say that $E[s] = 0$? We did! So the second term is zero and
we are left with:

$$
I = E[s^2]
$$

And I already prepared us for how to deal with an expectation: we'll square the
sum and multiply by the probabilities. So:

$$
I = \sum_j p_{j|i} \left( -r_{ij}^2 + \sum_k p_{k|i} r_{ik}^2 \right)^2 \\
= \sum_j p_{j|i} \left( (r_{ij}^2)^2 - 2 r_{ij}^2 \sum_k p_{k|i} r_{ik}^2 + \left( \sum_k p_{k|i}
 r_{ik}^2 \right)^2 \right) \\
= \sum_j p_{j|i} (r_{ij}^2)^2
 - 2 \left( \sum_k p_{k|i} r_{ik}^2 \right) \left( \sum_j p_{j|i} r_{ij}^2 \right)
 + \left( \sum_k p_{k|i} r_{ik}^2 \right)^2
$$

Now we don't have to worry about the $k$ and $j$ being different, so we can just
use $j$ everywhere and the second and third terms can therefore be combined,
leaving us with:

$$
I = \sum_j p_{j|i} (r_{ij}^2)^2
 - \left( \sum_j p_{j|i} r_{ij}^2 \right)^2
$$

This looks a lot like a variance. In fact, yes, now is the time to scroll back
up and see that this is exactly the variance expression we had for 
$\delta_{i,U}$, except for the factor of $2 \beta_i^2$ (this is also the point 
you say "oh, yeah").

So another way to express the correlation dimension is:

$$
\delta_{i,U} = 2 \beta_i^2 I
$$

I have searched the literature but was unable to find any Fisher
Information-based approach to intrinsic dimensionality estimation. The closest I
got was [Carter, Raich and Hero](https://doi.org/10.1109/TSP.2009.2031722) who
use a proxy to Fisher Information distance between probability distributions as
part of manifold learning, but when they estimate the intrinsic dimensionality,
they use the classic
[Levina and Bickel](https://papers.nips.cc/paper_files/paper/2004/hash/74934548253bcab8490ebd74afed7031-Abstract.html) 
maximum likelihood estimator. I would be interested to hear if anyone knows of
any other work that makes this connection.

*November 1 2025*: I should note that the equation I came up with here isn't
ideal for calculations of intrinsic dimensionality for large datasets. In t-SNE
calibration for large datasets, usually only $3U$ nearest neighbors are kept for
each point which has an effect on the calculated values of $\beta$ and the
affinities. The consequence of this is that the intrinsic dimensionality
estimate is biased downards. In practice, the Levina-Bickel method is more
robust, especially when combined with the modification suggested by
[Mackay and Ghahramani](https://www.inference.org.uk/mackay/dimension/). This
is less of an issue with the datasets used in smallvis, where there is no
truncation. This critique also applies to the finite difference method 
originated by Lee and co-workers.

For a practical application, see [part 2](https://jlmelville.github.io/smallvis/idp.html).

Up: [Documentation Home](https://jlmelville.github.io/smallvis/).
