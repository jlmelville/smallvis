---
title: "Notes on PaCMAP"
date: "November 14 2021"
output:
  html_document:
    theme: cosmo
    toc: true
    toc_float:
      collapsed: false
editor_options: 
  chunk_output_type: inline
---

PaCMAP (Pairwise Controlled Manifold Approximation:
[github](https://github.com/YingfanWang/PaCMAP),
[arXiv](https://arxiv.org/abs/2012.04456),
[JMLR](https://www.jmlr.org/papers/v22/20-1061.html)) is a UMAP/LargeVis-like
method with a focus on maintaining global as well as local structure. In my
opinion, its main innovation is advocating for the inclusion of "mid-near"
pairs in the optimization. These are formed from finding points which aren't 
the nearest neighbors but are closer on average than a random point from the
dataset. The cost function in PaCMAP is also fairly simple, while maintaining
some properties that the authors identify as important to avoid distorting local
manifold structure while maintaining the good clustering behavior of UMAP.

## The Usual Preliminaries

The output coordinates of points $i$ is the vector $\mathbf{y_i}$. The distance
between two output points $i$ and $j$ is the Euclidean distance given as 
$d_{ij}$. But we usually work with the squared distances:


$$
d_{ij}^2 = \| \mathbf{y_i} - \mathbf{y_j} \|_2^2
$$

Output similarities/weights are like those used in t-SNE, using a Cauchy-like
function:

$$
w_{ij} = \frac{1}{1 + d_{ij}^2} 
$$


Some derivatives:

$$
\begin{aligned}
&\frac{\partial{d_{ij}}}{\partial{\mathbf{y_i}}} = \frac{1}{d_{ij}} \left( \mathbf{y_i} - \mathbf{y_j} \right)
\\
&\frac{\partial{w_{ij}}}{\partial{d_{ij}}} = -2 d_{ij} w_{ij}^2
\\
&\frac{\partial{w_{ij}}}{\partial{\mathbf{y_i}}} = -2 w_{ij}^2 \left( \mathbf{y_i} - \mathbf{y_j} \right)
\end{aligned}
$$

and the gradient of the cost function is:

$$
\frac{\partial{C}}{\partial{\mathbf{y_i}}} = -2 \frac{\partial{C}}{\partial{w_{ij}}} w_{ij}^2 \left(\mathbf{y_i} - \mathbf{y_j}\right)
$$

As the cost function for this kind of dimensionality reduction method is usually
given in terms of $w_{ij}$, we can concentrate on the form of
$\partial{C}/\partial{w_{ij}}$ and plug it into the full gradient above if we
want to write an actual implementation.

## Cost Function

In the PaCMAP paper, the cost function is defined in terms of $\tilde{d_{ij}}$,
which is defined as:

$$
\tilde{d_{ij}} = 1 + d_{ij}^2
$$

I will be working with $w_{ij}$ rather than $\tilde{d_{ij}}$ where possible.

Near neighbor and mid-near attractive interactions are of the form:

$$
C^{+}_{PaCMAP} = \frac{a\tilde{d_{ij}}}{b + \tilde{d_{ij}}} = \frac{a}{1 + bw_{ij}}
$$

where $a$ is a positive scalar weight that changes over the course of the
optimization to put different emphasis on local vs global pairs, and $b$ is a
fixed positive value.

For non-neighbors (also referred to as "further" or "far" pairs), the repulsive
interactions are:

$$
C^{-}_{PaCMAP} = \frac{1}{1 + \tilde{d_{ij}}} = \frac{w_{ij}}{1 + w_{ij}}
$$

In the paper, the repulsive interaction does have a weight associated it (i.e.
the numerator can be > 1) but in practice the optimization procedure always
fixes the numerator as 1, so for simplicity's sake I have omitted it.

PaCMAP never gives a "full" cost function that we might expect: i.e. one that
is defined over all possible pairs of points.

## Gradient

The derivative of the attractive interaction is:

$$
\frac{\partial C_{PaCMAP}^+}{\partial{w_{ij}}} =
-\frac{ab}{\left(1 + bw_{ij}\right)^2}
$$

and therefore the gradient is:

$$
\frac{\partial C_{PaCMAP}^+}{\partial{\mathbf{y_i}}} =
2\frac{abw_{ij}^2}{\left(1 + bw_{ij}\right)^2} \left(\mathbf{y_i} - \mathbf{y_j}\right)
$$

The derivative of the repulsive interaction is:

$$
\frac{\partial C_{PaCMAP}^-}{\partial{w_{ij}}} = 
\frac{1}{\left(1 + w_{ij}\right)^2}
$$

hence the gradient is:

$$
\frac{\partial C_{PaCMAP}^-}{\partial{\mathbf{y_i}}} =
-2\frac{w_{ij}^2}{\left(1 + w_{ij}\right)^2} \left(\mathbf{y_i} - \mathbf{y_j}\right)
$$

(if you wanted to include an $a$ and $b$ in the repulsive contribution, the form
of the gradient would look identical to the attractive interaction, but with a
difference in sign).

### Comparison with UMAP

The equivalent gradients for UMAP are:

$$
\begin{aligned}
& \frac{\partial C_{UMAP}^+}{\partial{\mathbf{y_i}}} = 2 w_{ij} \left(\mathbf{y_i} - \mathbf{y_j}\right) \\
& \frac{\partial C_{UMAP}^-}{\partial{\mathbf{y_i}}} =
-2\frac{w_{ij}^2}{1 - w_{ij}} \left(\mathbf{y_i} - \mathbf{y_j}\right)
\end{aligned}
$$

I have made a simplifying assumption in the case of UMAP that the sampling of
positive edges is uniform and the output weight function parameters are set so
they match that of PaCMAP.

### Paper/Code form

Here is the gradient in terms of the expression in the paper and how it looks
in the PaCMAP code.

This requires the following derivatives:

$$
\begin{aligned}
&\frac{\partial{d_{ij}}}{\partial{\mathbf{y_i}}} = \frac{1}{d_{ij}} \left( \mathbf{y_i} - \mathbf{y_j} \right)
\\
&\frac{\partial{\tilde{d_{ij}}}}{\partial{d_{ij}}} = 2 d_{ij}
\\
&\frac{\partial{\tilde{d_{ij}}}}{\partial{\mathbf{y_i}}} = 2 \left( \mathbf{y_i} - \mathbf{y_j} \right)
\\
&\frac{\partial{C}}{\partial{\mathbf{y_i}}} = 2 \frac{\partial{C}}{\partial{\tilde{d_{ij}}}} \left(\mathbf{y_i} - \mathbf{y_j}\right)
\end{aligned}
$$

The attractive derivative is:

$$
\frac{\partial C_{PaCMAP}^+}{\partial{\tilde{d_{ij}}}} =
\frac{ab}{\left(b + \tilde{d_{ij}}\right)^2}
$$

leading to the gradient:

$$
\frac{\partial C_{PaCMAP}^-}{\partial{\mathbf{y_i}}} =
2\frac{ab}{\left(b + \tilde{d_{ij}}\right)^2} \left(\mathbf{y_i} - \mathbf{y_j}\right)
$$

The repulsive derivative is:

$$
\frac{\partial C_{PaCMAP}^-}{\partial{\tilde{d_{ij}}}} =
-\frac{1}{\left(1 + \tilde{d_{ij}}\right)^2}
$$

giving:

$$
\frac{\partial C_{PaCMAP}^+}{\partial{\mathbf{y_i}}} =
-2\frac{1}{\left(1 + \tilde{d_{ij}}\right)^2} \left(\mathbf{y_i} - \mathbf{y_j}\right)
$$

### Some R functions

To convince myself that have correctly expressed these gradients are equivalent,
here are two R functions that calculate the "force constant" part of the
gradients (i.e. everything that isn't the $\mathbf{y_i} - \mathbf{y_j}$ part):

```R
pacmap_kw <- function(dij, a = 1, b = 1) {
  dij2 <- dij * dij
  wij <- 1 / (1 + dij2)
  bw1 <- 1 + b * wij
  w1 <- 1 + wij
  list(
    attr = 2 * (a * b * wij * wij) / (bw1 * bw1),
    rep = -2 * wij * wij / (w1 * w1)
  )
}

pacmap_kd <- function(dij, a = 1, b = 1) {
  dij2 <- dij * dij
  dsquig <- 1 + dij2
  bd <- b + dsquig
  d1 <- 1 + dsquig
  list(
    attr = 2 * (a * b)  / (bd * bd),
    rep = -2 / (d1 * d1)
  )
}
```

## Weight Schedule

For neighbors, $b = 10$. For mid-near pairs, $b = 10000$. The $a$ weights
(referred to as $w_{NB}$, $w_{MN}$ and $w_{FP}$ for "near pairs", "mid-near
pairs" and "further pairs" respectively in the paper) have different values at
different points in the optimization schedule:

Iteration  | near | mid       | further
:----------|:-----|:----------|:-------
0 -> 100   | 2    | 1000 -> 3 | 1
101 -> 200 | 3    | 3         | 1
201 -> 450 | 1    | 0         | 1

I have kept the `further` weight column for consistency with the paper even
though as you can see, far pair weights are always 1. For the mid-near pairs
during the first 100 iterations, `1000 -> 3` means the weight is linearly
decayed from 1000 to 3.

Neither the $a$ or $b$ parameters of the cost function can be modified by the
user in the PaCMAP implementation. The number of iterations can be changed, but
the first 200 iterations always use the same weights and schedule, so if you set
`num_iters = 900`, this means the third stage of the optimization schedule will
last for 700 iterations, rather than doubling the duration of all three stages
proportionally.

## Force Constant Plots

To help compare with UMAP, I have plotted the attractive and repulsive "force
constants", $k_{attr}$ and $k_{rep}$, over the range of $w_{ij}$. The force
constant is the bit of the gradient that isn't the 
$\left(\mathbf{y_i} - \mathbf{y_j}\right)$ part, i.e. if you consider the
layout to be about finding an equilibrium of points attached to springs, this is
the stiffness of the spring between two points. Alternatively, consider this a
gradient plot, given a fixed displacement of 
$\left(\mathbf{y_i} - \mathbf{y_j}\right) = 1$.

First, the attractive interactions of near neighbors:

![Attractive Forces](../img/pacmap/pacmap-near.png)

The blue line is UMAP with default settings. The red line is UMAP with the
output function chosen to give the same $w_{ij}$ values as PaCMAP. This version
of UMAP I will refer to as t-UMAP as it shares the same output weights as t-SNE
as well as PaCMAP. The green line is PaCMAP with $a = w_{NB} = 1$, which occurs
during the final part of the optimization, and the yellow is with 
$a = w_{NB} = 3$, which occurs during the middle part of the optimization. The 
UMAP forces are not only much higher, at least for the default settings, the
attractive forces rapidly accelerate in magnitude for high $w$ (i.e. neighboring
points in the embedding which are supposed to be close). For t-UMAP, the
relationship between $w_{ij}$ and $k_{attr}$ is linear. For PaCMAP, there is a
much more gentle increase in force as $w_{ij}$ increases, and the rate of change
decreases rather than increases. For PaCMAP, a $w_{ij} >= 0.62$ feels at least
90% of the maximum attractive force whereas for UMAP, $w_{ij} >= 0.96$. To be
pedantic, if we take into account UMAP's gradient clipping, it's actually
$w_{ij} >= 0.92$, but the overall message stays the same.

We'll come back to the mid-pairs after we have looked at the repulsive forces:

![Repulsive Forces](../img/pacmap/pacmap-far.png)

The colors are the same as in the previous plot: blue for UMAP, red for t-UMAP
and green for PaCMAP (there is only ever one weighting for the repulsive force).
It's hard to put the full range of forces for these methods in one plot, so I
have zoomed in on the top portion of the plot. There is a very similar pattern
as with the attractive forces: UMAP repulsions steeply increase as non-similar
items get close together in the embedding. Nothing very interesting happens
in the rest of the UMAP curves that aren't shown here. Note that the repulsive
gradient in UMAP requires a small positive value to prevent division by zero
at very high $w_{ij}$, so repulsions can get very large and negative. The UMAP
implementation clips the gradient so it doesn't exceed `4`: this is the total
gradient, not just the force constant, but that also puts an upper bound on the
effective force constant as `4`. PaCMAP on the other hand has no singularities
in its repulsion, (the maximum repulsion is `0.5`) and so no need need for
gradient clipping. It also doesn't show an acceleration in repulsion at high
$w_{ij}$.

Finally, let's look at the mid pair attractions, a PaCMAP-only feature:

![Mid Pair Attractive Forces](../img/pacmap/pacmap-mid.png)

This plot doesn't have any contributions from UMAP. The yellow, green and grey
lines are the attractive interactions from the "near" pairs with 
$a = w_{NB} = 3$, $a = w_{NB} = 2$ and $a = w_{NB} = 1$ respectively. The green
and yellow lines are the same as those from the first plot. The mid near forces
are the cyan line with $w_{MN} = 1000$, which is the value at which the
optimization starts, and purple with $w_{MN} = 3$, the value at the end of the
first part of the optimization, the force constant being linear scaled in
between. The grey line is the force constant that applies to the near pairs over
the same part of the optimization. The force constant is effectively uniform
over the entire range of $w_{ij}$ and by the end of the first part of the
optimization, $k_{attr} \approx 0.0006$ and never increases thereafter, so it
doesn't seem like it plays a major organizing role after the first 100
iterations.

Overall, this seems to point to why UMAP can "tear" manifolds or promote
clustering vs PaCMAP: forces (attractive or repulsive) which are much higher
at high $w_{ij}$. Based on my experiments with the weighting parameters of
the UMAP output weight function, I don't think it's possible to generate curves
which look like that of PaCMAP (the t-UMAP curve is the best you can do).

## No Gradient Clipping

I discussed this in the section above, but it's worth calling out as another
point of difference with UMAP: the PaCMAP gradients don't require clipping.

## Uniform Near Pair Weights

Note that for a given type of pair (e.g. near pairs), they all get the same
weight. This is different to UMAP and LargeVis (and t-SNE), where there is an
input weight calibration stage that puts a higher weight on pairs which have a
higher similarity. In UMAP and LargeVis this manifests as certain edges being
sampled more regularly during optimization. This is not the case for PaCMAP,
which effectively assigns all pairs the same weight and samples every pair
during each epoch

Alternatively, you could lump the near and mid-near pairs together as the
"positive" edges, and consider PaCMAP as applying only 1 or 2 different possible
weight values to the positive edges during any given epoch.

## Input Preprocessing

By default, if the input dimensionality $d > 100$ then the data is centered (but
not scaled) and PCA is applied via the scikit-learn
[TruncatedSVD](https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.TruncatedSVD.html)
to reduce the input dimensionality to 100. If PCA is not asked for (or the data
is already low-dimensional), the data is range-scaled and then centered.

The Truncated SVD used in scikit-learn seems to use the same method by 
[Halko and co-workers](https://arxiv.org/abs/0909.4061) as used in
[irlba](https://cran.r-project.org/package=irlba)'s `svdr` function and cited by
the [rsvd package](https://cran.r-project.org/package=rsvd) (see also the 
[github](https://github.com/erichson/rSVD) and
[JSS publication](https://arxiv.org/abs/1608.02148)), but the scikit-learn
version is noticeably faster than either (even after installing OpenBLAS for 
faster linear algebra routines on Linux). Results for the later dimensions are
not the same between TruncatedSVD in Python and irlba, rsvd or even the unrelated
[bigstatsr](https://cran.r-project.org/package=bigstatsr) package so there may
be an accuracy/speed trade off to be investigated on the R side.

*16 November 2021* I have confirmed something along these lines for the 
`macosko2015` dataset (of dimensions `44808 x 3000`). For centered but unscaled
data, extracting 100 components explains around 37% of the variance, so
differences in results aren't due to numerical issues that arise when all the
variance in a dataset has been extracted and the last few columns are just
creative ways of expressing a vector that should be full of zeros. The
`TruncatedSVD` (in version 1.0 of scikit-learn) with default parameters produces
columns which are substantially less orthogonal than those produced by the R
packages mentioned above. After normalizing the columns, the dot products of
successive columns are fairly uniform across all 100 columns for all the R
methods, in the region of $10^{-14}$ to $10^{-16}$. For `TruncatedSVD`, the
columns start at a comparable level of orthogonality but quickly become less
orthogonal. The dot product of the tenth and eleventh column is around $10^{-7}$
and for the last ten columns the dot products are in the region of $10^{-4}$ to
$10^{-5}$. This suggests that for the purposes of approximate nearest neighbor
search, initial dimensionality reduction via SVD may be able to use
substantially looser tolerances than the defaults used in many R packages.

## Nearest Neighbors

The number of nearest neighbors is scaled according to dataset size. For 
$N < 10000$, the number of neighbors = 10. For $N > 10000$ the following formula 
is used:

$$
n_{nbrs} = \left \lfloor 10 + 15 \left(\log_{10}n - 4\right) \right \rceil
$$

where $\left \lfloor \cdot \right \rceil$ indicates my cumbersome attempt to
communicate "round to the nearest integer".

Perhaps some examples would help:

$N$       | $n_{nbrs}$
:---------|:--------
10 000    | 10
20 000    | 15
50 000    | 20
60 000    | 22
70 000    | 23
100 000   | 25
1 000 000 | 40

If you didn't want to think about it, it seems like setting `n_neighbors = 20`
would be fine.

### The Approximate Nearest Neighbor Calculation

[Annoy](https://github.com/spotify/annoy) is used for nearest neighbor
calculations. `n_trees = 20`. Instead of finding the `n_neighbors` nearest
neighbors however, actually `n_neighbors + 50` are searched for. This is to
allow finding the neighbors with the smallest scaled distances (see below).

PaCMAP supports all the distance metrics that Annoy does (i.e. Cosine, Manhattan
and Hamming in addition to Euclidean). There is one odd wrinkle to this (see
the mid-near section below).

### Scaled Distances

The closest `n_neighbors` of each point are not actually used. Instead the
`n_neighbors` closest points in terms of a "generalized" squared distance,
originally proposed by
[Zelnik-Manor and Perona](https://papers.nips.cc/paper/2004/hash/40173ea48d9567f1f393b20c855bb40b-Abstract.html),
is used:

$$
d_{ij}^{2,select} = \frac{r_{ij}^2}{\sigma_i \sigma_j}
$$

where $r_{ij}^2$ represents the distance in the input space. The local scaling
factor $\sigma_i$ is the mean average (unscaled) distance of the 4th-6th nearest
neighbors of $i$. This choice for $\sigma_i$ is the same as that used in the
Trimap method ([github](https://github.com/eamid/trimap), 
[arXiv](https://arxiv.org/abs/1910.00204)).

The closest scaled distance neighbors are searched from within the 
`n_neighbors + 50` with the smallest unscaled nearest neighbors. Scaled 
distances are not used in any other part of the method.

Roughly then we can expect to be looking for `60-75` nearest neighbors for
datasets up to $N = 100 000$.

## Mid Near Neighbors

The default number of mid-near neighbors is half that of the near neighbors
(rounded to the nearest integer). Each mid-near neighbor is chosen by picking
six items at random and keeping the one with the second-smallest distance. Note
that scaled distances are *not* used here.

Oddly, Euclidean distance is always used for measuring the distances at this
point no matter what distance metric was chosen for the near neighbors.

If we consider the near and mid-near neighbors together as making up the
"positive" edges that feel an attractive interaction during optimization, we
are using around `15-38` neighbors up to $N = 100 000$.

## Far Pairs

The default number of far pairs is twice that of the near neighbors. This is
slightly fewer negative examples than the UMAP default of `5` per point.
However the UMAP implementation does not apply the gradient descent update
to the negative points, which effectively halves its negative sampling rate,
so in practice these values are actually quite comparable.

As in UMAP, these are picked with at (uniform) random.

## Neighbors are Sampled Once and Re-used

In UMAP and LargeVis, the "negative" sampling (same as "far pairs" here) is
repeated during each epoch, so a given point $i$ always gets an entirely new
set of far pairs. PaCMAP samples both the mid-near and far pairs once before
the optimization begins and re-uses them at each epoch. This saves on repeated
random number generation, at the cost of increased memory usage for storing
the indices of the pairs. With default settings we can expect an extra storage
cost of $2.5 * N * n_{nbrs}$, which technically would make PaCMAP's memory costs
scale as $O(N \log N)$.

## Output Initialization

The default initialization uses PCA. The resulting coordinates are all shrunk
by a factor of `0.01`. Alternatively a random initialization is available using
the standard t-SNE initialization of random normal variables with a standard
deviation of `1e-4`. If a user-supplied matrix is provided, the coordinates
are centered and then shrunk by a factor of `1e-4`.

## Optimization

Unlike UMAP and LargeVis, PaCMAP updates the embedding using a full batch
gradient, i.e. coordinates are updated all at once at the end of each epoch.
The gradient is still stochastic however.

Optimization is carried out using [Adam](https://arxiv.org/abs/1412.6980),
popular in deep learning (or at least popular among people who write papers
about deep learning). The momentum parameters $\beta_1$ and $\beta_2$ are kept
at their default values ($0.9$ and $0.999$ respectively). The learning rate
$\alpha = 1$ and $\epsilon = 10^{-7}$. The latter was intended only to prevent
division by zero (although different values
[may also affect performance](https://papers.nips.cc/paper/2018/hash/90365351ccc7437a1309dc64e4db32a3-Abstract.html)).

## PaCMAP reports its Stochastic Cost

As can be seen in the consideration of the gradients, the per-point PaCMAP loss
can be easily calculated using values needed for the gradient. PaCMAP stores the
loss at the end of the gradient vector and can report it at different stages
during the optimization. Because the negative pairs are constant for the entire
optimization, we might expect there to be less variance than if a UMAP-style
approach to negative sampling was used. On the other hand, presumably this cost
function has a higher bias. And the variance issue may still bite you if you
try to compare different costs between different runs with different random
number seeds.

## See Also

Some notes on the related method
[NCVis](https://jlmelville.github.io/smallvis/ncvis.html).

## R Code for Plots

```R
# generate force constants over a range of w
pacmap_k <- function(a = 1, b = 10, n = 100, max_w = 0.99) {
  w <- (1:n) / n
  w <- w[w < max_w]
  ka <- (2 * a * b * w * w) / ((1 + b * w) ^ 2)
  kr <- (-2 * w * w) / ((1 + w) ^ 2)
  list(attr = ka, rep = kr, w = w)
}

umap_k <- function(a = 1.577, b = 0.895, n = 100, eps = 1e-3, max_w = 0.99) {
  w <- (1:n) / n
  w <- w[w < max_w]
  d2 <- ((1 - w) / (a * w)) ^ (1 / b)
  d2b1 <- d2 ^ (b - 1)

  ka <- 2 * a * b * d2b1 * w
  kr <- (-2 * b * w) / (eps + d2)
  list(attr = ka, rep = kr, w = w)
}

uk <- umap_k()
tk <- umap_k(a = 1, b = 1)
pk <- pacmap_k()
pk3 <- pacmap_k(a = 3)

# Colors from:
# inlmisc::GetColors(7, scheme = "bright")

lwd <- 2
plot(
  uk$w,
  uk$attr,
  type = "l",
  xlab = "w",
  ylab = "k-attr",
  main = "Near Pairs Attractive Force",
  sub = "blue: UMAP, red: t-UMAP, green: PaCMAP wNB=1, yellow: PaCMAP wNB=3",
  cex.sub = 0.75,
  lwd = lwd,
  col = "#4477AA"
)
lines(uk$w, tk$attr, col = "#EE6677", lwd = lwd)
lines(uk$w, pk$attr, col = "#228833", lwd = lwd)
lines(uk$w, pk3$attr, col = "#CCBB44", lwd = lwd)

plot(
  uk$w,
  uk$rep,
  type = "l",
  xlab = "w",
  ylab = "k-rep",
  main = "Far Pairs Repulsive Force",
  sub = "blue: UMAP, red: t-UMAP, green: PaCMAP",
  cex.sub = 0.75,
  lwd = lwd,
  col = "#4477AA",
  ylim = c(-2, 0)
)
lines(uk$w, tk$rep, col = "#EE6677", lwd = lwd)
lines(uk$w, pk$rep, col = "#228833", lwd = lwd)

pk2 <- pacmap_k(a = 2)
pkm0 <- pacmap_k(a = 1000, b = 10000)
pkm100 <- pacmap_k(a = 3, b = 10000)
plot(
  uk$w,
  pk3$attr,
  type = "l",
  xlab = "w",
  ylab = "k-attr",
  main = "PaCMAP Attractive Force",
  sub = "yellow: wNB = 3, grey: wNB = 2, green: wNB = 1, cyan: wMN = 1000, purple: wMN = 3",
  cex.sub = 0.75,
  lwd = lwd,
  col = "#CCBB44"
)
lines(uk$w, pk$attr, col = "#228833", lwd = lwd)
lines(uk$w, pk2$attr, col = "#BBBBBB", lwd = lwd)
lines(uk$w, pkm0$attr, col = "#66CCEE", lwd = lwd)
lines(uk$w, pkm100$attr, col = "#AA3377", lwd = lwd)
```

## Changelog

* November 14 2021
    * Created document.
* December 19 2021
    * Add equations for UMAP gradient.
    * Add plots of sample attractive and negative force constants.
