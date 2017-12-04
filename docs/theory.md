---
title: "Some Theory"
output:
  html_document:
    theme: cosmo
    toc: true
    toc_float:
      collapsed: false
---

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
the sum of all the weights: $Z = \sum_{ij} w_{ij}$.
* $w_{ij}$ in both t-SNE and LargeVis is defined as the Student's t-distribution 
with one degree of freedom (or the Cauchy distribution, if you prefer), 
$w_{ij} = 1 / \left(1 + d_{ij}^2 \right)$. $d_{ij}$ is the Euclidean distance
between point $i$ and $j$ in the output coordinates.
* $\mathbf{y_{i}}$ is the vector of coordinates of point $i$ (with $N$ points 
in total).

## Cost Functions

### t-SNE

The t-SNE cost function is:

$$
C_{tSNE} = 
\sum_{ij} p_{ij} \ln \frac{p_{ij}}{q_{ij}}
=
\sum_{ij} \left( p_{ij} \ln p_{ij} - p_{ij} \ln q_{ij} \right)
$$
The first term in the sum has no dependence on the output 
coordinates, so is a constant we'll just mark as $C_{P}$.

Now let's write out $q_{ij}$ as $w_{ij} / Z$:

$$
C_{tSNE} = C_{P} - \sum_{ij} p_{ij} \ln \left( \frac{w_{ij}}{Z} \right) = 
Cp - \sum_{ij} p_{ij} \ln w_{ij} + \sum_{ij} p_{ij} \ln Z
$$
Finally, we'll do some rearranging and re-write $Z$ back to a sum of weights:

$$
C_{tSNE} = 
Cp - \sum_{ij} p_{ij} \ln w_{ij} + \ln Z \sum_{ij} p_{ij} =
Cp - \sum_{ij} p_{ij} \ln w_{ij} + \ln \sum_{ij} w_{ij}
$$

Ignoring, the constant term, we can see that the SNE cost function consists
of an attractive term, where maximizing the $w_{ij}$ (which implies minimizing 
the distances) would minimize $-p_{ij} \ln w_{ij}$; and a repulsive term, where 
minimizing the sum of $w_{ij}$ (and hence maximizing the distances), 
will minimize the log of the sum.

### LargeVis

LargeVis partitions the data into a set of nearest neighbors, which only feel 
an attractive force; and everything else, which only feel repulsive forces. 
When evaluated over an entire dataset, this seems to be difficult to optimize,
so for now, `smallvis` assigns both attractive and repulsive components to all
pairs:

$$
C_{LV} = 
-\sum_{ij} v_{ij} \ln w_{ij} 
-\gamma \sum_{ij} \ln \left( 1 - w_{ij} \right)
$$
Equation 2 of the LargeVis paper indicates that the input affinities
are normalized, i.e. they use $p_{ij}$, not $v_{ij}$ in the first term. But 
looking at the source code, the affinities seem to be un-normalized. So we're
doing the same here.

Nonetheless, you can see that the attractive terms of both t-SNE and LargeVis 
are very similar.

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

The output weights are given by:

$$w_{ij} = 1 / \left(1 + ad_{ij}^{2b}\right)$$

where $a$ and $b$ are determined by a non-linear least squares fit based on 
a couple of user-selected parameters that control the tightness of the squashing
function. By setting $a = 1$ and $b = 1$ you get the t-SNE style weighting back.
The current UMAP defaults result in $a = 1.929$ and $b = 0.7915$.

## Gradients

Here are the gradients for t-SNE, LargeVis and UMAP, with the t-SNE one written
slightly differently to how its usually presented to show the sum of attractive
and repulsive forces more clearly and for comparison with LargeVis:

$$
\frac{\partial C_{tSNE}}{\partial \mathbf{y_i}} = 
  4\sum_j^N \left(
    p_{ij} w_{ij}
    -
    q_{ij} w_{ij}
   \right)
   \left(\mathbf{y_i - y_j}\right)
$$

$$
\frac{\partial C_{LV}}{\partial \mathbf{y_i}} = 
  4\sum_j^N \left(
    v_{ij} w_{ij}
    -\frac{\gamma w_{ij}}{d_{ij}^2 + \epsilon}
   \right)
   \left(\mathbf{y_i - y_j}\right)
$$
$$
\frac{\partial C_{UMAP}}{\partial \mathbf{y_i}} = 
  4\sum_j^N \left(
    abd_{ij}^{2\left(b - 1\right)}v_{ij} w_{ij}
    -\frac{b \left(1 - v_{ij}\right) w_{ij}}{d_{ij}^2 + \epsilon}
   \right)
   \left(\mathbf{y_i - y_j}\right)
$$
The UMAP gradient is a good deal less simple than the t-SNE or LargeVis 
gradient, but it can be seen that it has a similar structure.

The $\epsilon$ term in the LargeVis and UMAP gradient is needed computationally
to avoid division by zero. Results can be quite sensitive to changing this value.
In `smallvis` you can see the effect of it via the `lveps` parameter. In the
stochastic gradient descent implementations of LargeVis and UMAP, the values
are set to `0.1` and `0.001`, respectively.

