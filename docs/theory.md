---
title: "Some Theory"
output:
  html_document:
    theme: cosmo
    toc: true
    toc_float:
      collapsed: false
---
Up: [Documentation Home](https://jlmelville.github.io/smallvis/).

Here are the t-SNE and LargeVis cost functions and gradients, making it a bit
easier to see how they are related. 

A bit of nomenclature first:

* $v_{ij}$ are the un-normalized input affinities. A large affinity means that
the data points are similar in the input space, so you'd probably like to see
them close together in the output configuration.
* $p_{ij}$ are the normalized input affinities, so that they sum to one. They
are symmetric, so that $p_{ij} = p_{ji}$ They are often interpreted as joint
probabilities. Probabilities of what exactly? It's not totally
obvious. Let's say it's the probability of observing an edge between the two
vertices $i$ and $j$ in the graph representing the neigborhood relationship 
of your input data. These are dependent only on the input data (and the choice 
of perplexity) and so are constant with respect to optimizing the output 
coordinates. If you decide to think of the problem in terms of graph laplacians
and $p_{ij}$ as "just" a normalized affinity, I won't blame you.
* $w_{ij}$ are the output configuration weights, aka un-normalized affinities 
or similarities.
* $q_{ij}$ are the output probabilities, where $q_{ij} = w_{ij} / Z$ and $Z$ is
the sum of all the weights: $Z = \sum_{kl} w_{kl}$.
* $w_{ij}$ in both t-SNE and LargeVis is defined as the Student's t-distribution 
with one degree of freedom (or the Cauchy distribution, if you prefer), 
$w_{ij} = 1 / \left(1 + d_{ij}^2 \right)$. $d_{ij}$ is the Euclidean distance
between point $i$ and $j$ in the output coordinates.
* $\mathbf{y_{i}}$ is the vector of coordinates of point $i$ (with $N$ points 
in total).

If you need more detail (or to be more confused) on how to derive the embedding 
gradients from their cost functions, see 
[here](http://jlmelville.github.io/sneer/gradients.html).

## Cost Functions

### t-SNE

The t-SNE cost function is:

$$
C_{tSNE} = 
\sum_{ij} p_{ij} \log \frac{p_{ij}}{q_{ij}}
=
\sum_{ij} \left( p_{ij} \log p_{ij} - p_{ij} \log q_{ij} \right)
$$
The first term in the sum has no dependence on the output 
coordinates, so is a constant we'll just mark as $C_{P}$.

Now let's write out $q_{ij}$ as $w_{ij} / Z$:

$$
C_{tSNE} = C_{P} - \sum_{ij} p_{ij} \log \left( \frac{w_{ij}}{Z} \right) = 
Cp - \sum_{ij} p_{ij} \log w_{ij} + \sum_{ij} p_{ij} \log Z
$$
Finally, we'll do some rearranging and re-write $Z$ back to a sum of weights:

$$
C_{tSNE} = 
Cp - \sum_{ij} p_{ij} \log w_{ij} + \log Z \sum_{ij} p_{ij} \\
= Cp - \sum_{ij} p_{ij} \log w_{ij} + \log \sum_{ij} w_{ij}
$$

Ignoring, the constant term, we can see that the SNE cost function consists
of an attractive term, where maximizing the $w_{ij}$ (which implies minimizing 
the distances) would minimize $-p_{ij} \log w_{ij}$; and a repulsive term, where 
minimizing the sum of $w_{ij}$ (and hence maximizing the distances), 
will minimize the log of the sum.

### LargeVis

The LargeVis paper describes partitioning the data into a set  of nearest
neighbors (we'll call that $E$), which only feel an attractive force; and
everything else ($\bar{E}$), which only feel repulsive forces. Perplexity
calibration and weight definition are carried out as in t-SNE, to give the
following (log) likelihood function:

$$
L_{LV} = \sum_{ \left(i, j\right) \in E} p_{ij} \log w_{ij} 
+\gamma \sum_{\left(i, j\right) \in \bar{E}} \log \left( 1 - w_{ij} \right)
$$
Likelihood functions are maximized, so the final cost function we're
interested in will be $-L_{LV}$, in order that we have a function to minimize.
Otherwise, it gets a bit confusing keeping track of the signs of the gradients 
compared to the other methods in `smallvis`.

Before we get to that though, I believe that the likelihood function as given
above is not a totally accurate representation of the current state of the
LargeVis implementation. The following is my understanding of how it works at
the moment.

First, although equation 2 of the LargeVis paper indicates that the input 
weights are normalized, I don't think the input weights actually are normalized, 
i.e. $v_{ij}$ is used rather than $p_{ij}$.

Second, I think that the repulsive terms are actually applied to both
neighborhood and non-neighborhood points, so that the repulsive (second) term in
the likelihood function is actually a sum over all points, not just $\bar{E}$. I
have experimented with restricting the repulsion to only the non-neighbors in
`smallvis`, but this leads to some horrible results involving lots of very small
well-separated clusters. On reflection it seems like this would be the expected
outcome of not allowing neighbors to repel each other: neighboring points that
happen to be initialized close to each other will feel a strong mutual force
that results in them reducing their distances to zero, which will overwhelm any
longer range attraction from more distant neighbors. I also can't find any 
evidence in the LargeVis source code that prevents neighborhood points from
turning up as part of the non-neighborhood sample: there *is* a check 
made in the code, but it only prevents one of the $k$ neighbors from 
contributing to the repulsive part of the gradient applied to any given point.

That leaves the attractive part. The partitioning works here and the
non-neighbors *are* excluded from contributing to the attractive part
of the cost, but that term is weighted by $p_{ij}$ anyway, so the perplexity
calibration guarantees that even if you did calculate probabilities for all
pairs of points, those that weren't part of the neighborhood would contribute
negligibly to that component of the cost function.

Taken together, you could therefore ignore the partitioning scheme, as long as
you were prepared to calculate the complete O(N^2) matrices, which is exactly
what `smallvis` does. 

Based on the above discussion, along with using the un-normalized input weights,
and wanting a cost function to minimize, the `smallvis` version of the LargeVis
cost function is:

$$
C_{LV} = 
-\sum_{ij} v_{ij} \log w_{ij} 
-\gamma \sum_{ij} \log \left( 1 - w_{ij} \right) \\
= -\frac{1}{N}\sum_{ij} p_{ij} \log w_{ij} 
-\gamma \sum_{ij} \log \left( 1 - w_{ij} \right)
$$
In this form, apart from the issue of whether the input weights are normalized
or not, the attractive terms of both t-SNE and LargeVis are very similar.

### UMAP

Note: the following is based on my examination of the UMAP source code, followed
up by some clarification from UMAP creator Leland McInnes, who kindly answered 
some of my questions (in quite a lot of detail) about the intent of the source 
code. If any of the following seems wrong or nonsensical, that is a reflection
on my understanding of UMAP, rather than UMAP itself.

[UMAP](https://github.com/lmcinnes/umap) (Uniform Manifold Approximation and 
Projection) attempts to model the underlying manifold of a dataset via a fuzzy
topological structure. 

UMAP treats the input and output data as two fuzzy sets, with strength of 
membership being equivalent to $v_{ij}$ and $w_{ij}$, respectively. The cost
function for UMAP is the cross entropy of the fuzzy sets:

$$
C_{UMAP} = 
\sum_{ij} \left[ v_{ij} \log \left( \frac{v_{ij}}{w_{ij}} \right) + 
(1 - v_{ij}) \log \left( \frac{1 - v_{ij}}{1 - w_{ij}} \right) \right]
$$

which can be expanded to:

$$
C_{UMAP} =
\sum_{ij} \left[ v_{ij} \log(v_{ij}) + 
(1 - v_{ij})\log(1 - v_{ij}) \right] - 
\sum_{ij} \left[ v_{ij} \log(w_{ij}) \right] - 
\sum_{ij} \left[ (1 - v_{ij}) \log(1 - w_{ij}) \right]
$$

Just like with t-SNE, the first term is a constant, which this time we'll
call $C_{V}$, leaving:

$$
C_{UMAP} =
C_{V} -
\sum_{ij} v_{ij} \log(w_{ij}) - 
\sum_{ij} (1 - v_{ij}) \log(1 - w_{ij})
$$

which looks a lot like LargeVis, except instead of a constant $\gamma$ term
in the repulsion, each pair is weighted according to (one minus) the input 
weight. Additionally, both the input and output weight terms are defined 
differently to LargeVis.

The UMAP input weights are given by:

$$v_{ij} = \exp \left[ -\left( r_{ij} - \rho_{i} \right) / \sigma_{i} \right]$$

where $r_{ij}$ are the input distances (not necessarily Euclidean), $\rho_{i}$ is
the distance to the nearest neighbor (ignoring zero distances where neighbors
are duplicates) and $\sigma_{i}$ is chosen by a binary search such that 
$\sum_{j} v_{ij} = \log_{2} k$ where $k$ is the size of the neighborhood. This
is similar in spirit to the perplexity calibration used by t-SNE and LargeVis.

As the use of $v_{ij}$ indicates, the input weight are not normalized in UMAP.
They are symmetrized, but in a different way to the arithmetic mean approach
in t-SNE:

$$V_{symm} = V + V^{T} - V \circ V^{T}$$
where $T$ indicates the transpose and $\circ$ is the Hadamard (i.e. entry-wise) 
product. This effectively carries out a fuzzy set union.

The output weights are given by:

$$w_{ij} = 1 / \left(1 + ad_{ij}^{2b}\right)$$

where $a$ and $b$ are determined by a non-linear least squares fit based on 
a couple of user-selected parameters that control the tightness of the squashing
function. By setting $a = 1$ and $b = 1$ you get the t-SNE style weighting back.
The current UMAP defaults result in $a = 1.929$ and $b = 0.7915$.

To get to the UMAP gradient, we need the derivative of the weight with respect
to the squared distance.

$$
\frac{\partial w_{ij}}{\partial d_{ij}^2} =  
\frac{-b a d_{ij}^{2\left(b - 1\right)}} 
{\left(a d_{ij}^2 + 1 \right)^2} = -b a d_{ij}^{2\left( b - 1\right)} w_{ij}^2 
=
-\frac{b \left( 1 - w_{ij} \right)}{d_{ij}^2}w_{ij}
$$
The derivative of the cost functon with respect to the weights is:

$$
\frac{\partial C}{\partial w_{ij}} =
\left[
-\frac{v_{ij}}{w_{ij}} +
\frac{ \left( 1 - v_{ij} \right)}{ \left( 1 - w_{ij} \right)}
\right]
= 
\left[
-\frac{v_{ij}}{w_{ij} \left(1 - w_{ij} \right)} +
\frac{w_{ij}}{ w_{ij} \left( 1 - w_{ij} \right)}
\right]
$$

## Gradients

Let's look at some gradients and see how t-SNE, LargeVis and UMAP compare.

### t-SNE

First, t-SNE written slightly differently to how its usually presented, in terms
of the un-normalized weights. This will hopefully illuminate the similarities
and differences with LargeVis and UMAP:

$$
\frac{\partial C_{tSNE}}{\partial \mathbf{y_i}} = 
4
\sum_j^N \left( v_{ij} -\frac{\sum_{kl} v_{kl}}{\sum_{kl} w_{kl}} w_{ij} \right) 
\frac{w_{ij}}{\sum_{kl} v_{kl}}
\left(\mathbf{y_i - y_j}\right)
$$
If you want a bit more detail in how this is derived, see [here](http://jlmelville.github.io/sneer/experimental-gradients.html). Note that
the sum $\sum_{kl}$ means to sum over all weights and the input weights 
$v_{ij}$ are not matrix-normalized but otherwise have had the other processing 
carried out on them as is usually done in perplexity-based calibration, i.e. 
row-normalization and symmetrization.

#### Digression

Brief aside which will hopefully find a more suitable place to live one day. The
input weight processing in the usual t-SNE treatment means that 
$\sum_{kl} v_{kl} = N$. An alternative form for the t-SNE gradient is therefore:

$$
\frac{\partial C_{tSNE}}{\partial \mathbf{y_i}} = 
\frac{4}{N}
\sum_j^N \left( v_{ij} -\frac{N}{Z} w_{ij} \right)
w_{ij}
\left(\mathbf{y_i - y_j}\right)
$$
In practice, because the usual t-SNE initialization starts with short distances,
the output weights are all 1 to begin with, so the $N/Z$ term is approximately
$1/N$ initially. Over the course of the optimization the value of $Z$ begins
to drop and hence $N/Z$ increases. For the usual parameter settings in t-SNE, 
it's usually less than 1, but for low perplexity values and quite converged
results (e.g. allowing `max_iter = 50000` at least), $N/Z$ can get
larger than 1 (e.g. around 4 for the `s1k` dataset at `perplexity = 5`).
End of digression.

### LargeVis

For LargeVis, the gradient can be written as:

$$
\frac{\partial C_{LV}}{\partial \mathbf{y_i}} = 
  4\sum_j^N \left(
    v_{ij}
    -\frac{\gamma}{1 - w_{ij}} w_{ij}
   \right)
   w_{ij}
   \left(\mathbf{y_i - y_j}\right)
$$

where I've tried to retain the structure of the expression that most resembles
that of t-SNE, but the following is probably more convenient for the purposes
of a stochastic gradient descent optimization:

$$
\frac{\partial C_{LV}}{\partial \mathbf{y_i}} = 
  4\sum_j^N \left(
    v_{ij} w_{ij}
    -\frac{\gamma w_{ij}}{d_{ij}^2 + \epsilon}
   \right)
   \left(\mathbf{y_i - y_j}\right)
$$

The $\epsilon$ term is needed computationally to avoid division by zero. Results
can be quite sensitive to changing this value. In `smallvis` you can see the
effect of it via the `lveps` parameter. It's set to `0.1` in LargeVis.

### UMAP

Likewise, here's one way to express the UMAP gradient:

$$
\frac{\partial C_{UMAP}}{\partial \mathbf{y_i}} = 
  4\sum_j^N 
\left(
   v_{ij}
-
  w_{ij}
\right)
   \frac{b}{d_{ij}^2}
   \left(\mathbf{y_i - y_j}\right)
= 
  4\sum_j^N 
\left(
   v_{ij}
-
  w_{ij}
\right)
   b\left(\frac{aw_{ij}}{1 - w_{ij}}\right)^{\frac{1}{b}}
   \left(\mathbf{y_i - y_j}\right)
$$
and more explicitly:

$$
\frac{\partial C_{UMAP}}{\partial \mathbf{y_i}} = 
  4\sum_j^N \left[
    abd_{ij}^{2\left(b - 1\right)} w_{ij} v_{ij} 
    -\frac{b \left(1 - v_{ij}\right) }{d_{ij}^2 + \epsilon} w_{ij}
   \right]
   \left(\mathbf{y_i - y_j}\right)
$$

Again we need a value for $\epsilon$, which is `0.001` in the UMAP source.

The t-SNE, LargeVis and UMAP gradients all have a similar form based around the
difference between the input and output weights, $v_{ij} - w_{ij}$, although
LargeVis and t-SNE re-weight the repulsion compared to UMAP. 

Up: [Documentation Home](https://jlmelville.github.io/smallvis/).
