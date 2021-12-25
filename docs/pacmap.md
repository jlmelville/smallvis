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

## Pairwise Distance Comparison

To get a feel for the effect of the different neighbor calculations, here are
some histograms of pairwise distances for different sets of neighbors across
a variety of datasets. The different sets of neighbors are labelled as:

* `15`: The 15 nearest neighbors. Actually it's only the 14 nearest neighbors
because I omit the first nearest neighbor of each item which is always itself.
These are the exact nearest neighbors, not the result of an approximate nearest
neighbor calculation. This is the UMAP default value.
* `15s`: The 15 nearest neighbors using the scaled distances. These are selected
from the exact nearest 65 neighbors, as would be done in PaCMAP (i.e.
`n_neighbors + 50`).
* `150`: The 150 exact nearest neighbors. This is an order of magnitude higher
than the UMAP default, and would be the upper limit to the usual nearest
neighbor calculation in t-SNE, where often 3 * perplexity is used, and
perplexity is rarely set higher than 50.
* `mid near`: The mid near pairs, calculated as described above. Default in
PaCMAP is to calculate as many mid near pairs as `0.5 * n_neighbors`, but I have
calculated it using 15 mid near pairs for consistency (there is no effect on the
shape of the distribution of distances).
* `random`: 15 random neighbors.

For each of these set of neighbors, a histogram was generated of the distances.
I used the Euclidean distance in all cases. Overlaying all these histograms on
top of each other is way too confusing, so I have generated them using the same
break points on separate plots. Each column in the table below is one dataset,
so by looking down a column you can see how the pair-wise distance distribution
changes. I would expect the distances to increase as you look down the column.

The datasets are mostly discussed on the [uwot
examples](https://jlmelville.github.io/uwot/umap-examples.html) except for the
first three which are discussed in the PaCMAP paper. is less relevant here than
looking for overall trends.

The first three are low-dimensional (2 or 3 dimensions) quite highly structured
datasets. The `mammoth` dataset has 10,000 points and 3 features (X, Y, Z
coordinates: it's a 3D model of a mammoth). `sch10k` is the s-curve with a hole
also with 10,000 points and X, Y, Z coordinates. This is an S-shaped curve which
has a circular hole right in the middle (in the middle of the 'S'). The
`curve2d` is a 2D polynomial (1,450 points).

| mammoth | sch10k | curve2d |
|:----:|:----:|:-----:
![mammoth 15](../img/pacmap/mammoth15.png)|![sch10k 15](../img/pacmap/sch10k15.png)|![curve2d 15](../img/pacmap/curve2d15.png)
![mammoth 15s](../img/pacmap/mammoth15s.png)|![sch10k 15s](../img/pacmap/sch10k15s.png)|![curve2d 15s](../img/pacmap/curve2d15s.png)
![mammoth 150](../img/pacmap/mammoth150.png)|![sch10k 150](../img/pacmap/sch10k150.png)|![curve2d 150](../img/pacmap/curve2d150.png)
![mammoth mid](../img/pacmap/mammothmid.png)|![sch10k mid](../img/pacmap/sch10kmid.png)|![curve2d mid](../img/pacmap/curve2dmid.png)
![mammoth rand](../img/pacmap/mammothrand.png)|![sch10k rand](../img/pacmap/sch10krand.png)|![curve2d rand](../img/pacmap/curve2drand.png)

The next three datasets are all image datasets. They are fairly small in terms
of number of objects, but are quite high dimensional (features are pixels).
`coil20` is a set of 20 objects in 72 poses each (so 1440 items altogether) and
16384 features. `frey` is a set of 1965 video stills of Brendan Frey's face with
560 features. `oli` is the Olivetti faces (40 faces, 10 images each) with 4096
features.

| coil20 | frey | oli |
|:----:|:----:|:-----:
![coil20 15](../img/pacmap/coil2015.png)|![frey 15](../img/pacmap/frey15.png)|![oli 15](../img/pacmap/oli15.png)
![coil20 15s](../img/pacmap/coil2015s.png)|![frey 15s](../img/pacmap/frey15s.png)|![oli 15s](../img/pacmap/oli15s.png)
![coil20 150](../img/pacmap/coil20150.png)|![frey 150](../img/pacmap/frey150.png)|![oli 150](../img/pacmap/oli150.png)
![coil20 mid](../img/pacmap/coil20mid.png)|![frey mid](../img/pacmap/freymid.png)|![oli mid](../img/pacmap/olimid.png)
![coil20 rand](../img/pacmap/coil20rand.png)|![frey rand](../img/pacmap/freyrand.png)|![oli rand](../img/pacmap/olirand.png)

The next three datasets are of identical dimensions: 70000 images, with 784
features (pixels). `mnist` are handwritten digits, `fashion` are images of items
of clothing, `kuzushiji` is handwritten Japanese characters.

| mnist | fashion | kuzushiji |
|:----:|:----:|:-----:
![mnist 15](../img/pacmap/mnist15.png)|![fashion 15](../img/pacmap/fashion15.png)|![kuzushiji 15](../img/pacmap/kuzushiji15.png)
![mnist 15s](../img/pacmap/mnist15s.png)|![fashion 15s](../img/pacmap/fashion15s.png)|![kuzushiji 15s](../img/pacmap/kuzushiji15s.png)
![mnist 150](../img/pacmap/mnist150.png)|![fashion 150](../img/pacmap/fashion150.png)|![kuzushiji 150](../img/pacmap/kuzushiji150.png)
![mnist mid](../img/pacmap/mnistmid.png)|![fashion mid](../img/pacmap/fashionmid.png)|![kuzushiji mid](../img/pacmap/kuzushijimid.png)
![mnist rand](../img/pacmap/mnistrand.png)|![fashion rand](../img/pacmap/fashionrand.png)|![kuzushiji rand](../img/pacmap/kuzushijirand.png)

The last three datasets are: `cifar10`, which contains 60000 images
and 1024 features. `macosko2015` is an RNA-seq dataset with 44808 items and
3000 features. This dataset has a very high hubness (i.e. there is one item
which is the near neighbor of a large proportion of the dataset). `tasic2018`
is another RNA-seq dataset with 23822 items and 3000 features. `cifar10` 
and `tasic2018` have a lower hubness `macosko2015`, but still higher than 
`mnist`, `fashion`, and `kuzushiji`.

| cifar10 | macosko2015 | tasic2018 |
|:----:|:----:|:-----:
![cifar10 15](../img/pacmap/cifar1015.png)|![macosko2015 15](../img/pacmap/macosko201515.png)|![tasic2018 15](../img/pacmap/tasic201815.png)
![cifar10 15s](../img/pacmap/cifar1015s.png)|![macosko2015 15s](../img/pacmap/macosko201515s.png)|![tasic2018 15s](../img/pacmap/tasic201815s.png)
![cifar10 150](../img/pacmap/cifar10150.png)|![macosko2015 150](../img/pacmap/macosko2015150.png)|![tasic2018 150](../img/pacmap/tasic2018150.png)
![cifar10 mid](../img/pacmap/cifar10mid.png)|![macosko2015 mid](../img/pacmap/macosko2015mid.png)|![tasic2018 mid](../img/pacmap/tasic2018mid.png)
![cifar10 rand](../img/pacmap/cifar10rand.png)|![macosko2015 rand](../img/pacmap/macosko2015rand.png)|![tasic2018 rand](../img/pacmap/tasic2018rand.png)

My observations on this are:

* You have to look quite hard at the `15` and `15s` distributions to see any
effect of using scaled neighbors, but that's not hugely surprising: at the most
extreme, the procedure would select neighbors 51-65, and for `mammoth` data,
that distribution is over quite a small range of the possible distances so you
wouldn't be able to see much detail there anyway. But this is not the case for
most of the datasets though: distance concentration quickly kicks in (i.e. in
high dimensions, the difference between a close and a far distance shrinks). 
* Going from 15 nearest neighbors to 150 nearest neighbors makes little
difference in terms of the distribution.
* The mid-near distributions are noticeably different from the nearest neighbor
distributions. I think it would be reasonable to ask the question "instead of
this mid-near business, can't we just use more nearest neighbors?". I think the
answer here says no.
* The mid-near distributions are quite similar to the random distributions. This
is hardly surprising given how they are constructed. For most datasets, the
difference in distribution is pretty subtle though, due to distance
concentration. I would be interested to see whether you could get away with a
uniform random sample of distances and get the same effect as the mid-near
pairs. This might make it easier to introduce PaCMAP ideas into existing UMAP
code bases because you could use the random negative sampling phase as a source
of pairs to use with the mid-near attractive forces.

The `macosko2015` dataset, which I said was most affected by hubness shows the
least difference between the `15` distances and the `rand`. It's quite striking
how little the distribution changes. Unlike every other dataset I've worked with
this is the one that shows the biggest effect of a PCA pre-processing step. In
most cases, it's assumed that running PCA is a way to speed up nearest neighbor
distance calculations by reducing the input dimensionality, but with
`macosko2015`, the UMAP output is affected substantially. The PCA step produces
a much better separation of clusters. Personally, while it works out well for
that dataset, I also find it a bit worrying that a pre-processing step can do
that.

As applying PCA is the default action in PaCMAP, below I regenerated some of the
distributions where PCA was applied to reduce the data to 100 dimensions before
doing the nearest neighbor search or random distance calculations. The table
contains four rows: The non-PCA results for the `15s` results and then PCA
results below. Then I show the `mid near` distributions without PCA and then
after applying PCA. The non-PCA results are repeated from the tables above. The
trends are the same for `15`, `150` and `rand` so there's no point showing
those. Also, the first three datasets above aren't included, because they
already only had 2 or 3 features in them.

### Effect of PCA

| coil20 | frey | oli |
|:----:|:----:|:-----:
![coil20 15s](../img/pacmap/coil2015s.png)|![frey 15s](../img/pacmap/frey15s.png)|![oli 15s](../img/pacmap/oli15s.png)
![coil20 pca10015s](../img/pacmap/coil20pca10015s.png)|![frey pca10015s](../img/pacmap/freypca10015s.png)|![oli pca10015s](../img/pacmap/olipca10015s.png)
![coil20 mid](../img/pacmap/coil20mid.png)|![frey mid](../img/pacmap/freymid.png)|![oli mid](../img/pacmap/olimid.png)
![coil20 pca100mid](../img/pacmap/coil20pca100mid.png)|![frey pca100mid](../img/pacmap/freypca100mid.png)|![oli pca100mid](../img/pacmap/olipca100mid.png)

| mnist | fashion | kuzushiji |
|:----:|:----:|:-----:
![mnist 15s](../img/pacmap/mnist15s.png)|![fashion 15s](../img/pacmap/fashion15s.png)|![kuzushiji 15s](../img/pacmap/kuzushiji15s.png)
![mnistpca100 15s](../img/pacmap/mnistpca10015s.png)|![fashionpca100 15s](../img/pacmap/fashionpca10015s.png)|![kuzushijipca100 15s](../img/pacmap/kuzushijipca10015s.png)
![mnist mid](../img/pacmap/mnistmid.png)|![fashion mid](../img/pacmap/fashionmid.png)|![kuzushiji mid](../img/pacmap/kuzushijimid.png)
![mnistpca100 mid](../img/pacmap/mnistpca100mid.png)|![fashionpca100 mid](../img/pacmap/fashionpca100mid.png)|![kuzushijipca100 mid](../img/pacmap/kuzushijipca100mid.png)

| cifar10 | macosko2015 | tasic2018 |
|:----:|:----:|:-----:
![cifar10 15s](../img/pacmap/cifar1015s.png)|![macosko2015 15s](../img/pacmap/macosko201515s.png)|![tasic2018 15s](../img/pacmap/tasic201815s.png)
![cifar10pca100 15s](../img/pacmap/cifar10pca10015s.png)|![macosko2015pca100 15s](../img/pacmap/macosko2015pca10015s.png)|![tasic2018pca100 15s](../img/pacmap/tasic2018pca10015s.png)
![cifar10 mid](../img/pacmap/cifar10mid.png)|![macosko2015 mid](../img/pacmap/macosko2015mid.png)|![tasic2018 mid](../img/pacmap/tasic2018mid.png)
![cifar10pca100 mid](../img/pacmap/cifar10pca100mid.png)|![macosko2015pca100 mid](../img/pacmap/macosko2015pca100mid.png)|![tasic2018pca100 mid](../img/pacmap/tasic2018pca100mid.png)

`macosko2015` does indeed show a pretty dramatic change in its distribution,
particularly its nearest neighbors. `tasic2018` also seems to show this. Could
just be coincidence that those are both RNA-seq datasets, but probably a look at
how such datasets are normalized or otherwise pre-processed might be in order.
For the other datasets, PCA doesn't seem to distort the distributions as much.
However, as they are all image datasets, that doesn't necessarily mean that
there aren't problems lurking. Memo to self: find some non-image datasets, e.g.
the 20 Newsgroup dataset (available in the PaCMAP repo) for text categorization.
`macosko2015` does indeed show a pretty dramatic change in its distribution,
particularly its nearest neighbors. `tasic2018` also seems to show this. Could
just be coincidence that those are both RNA-seq datasets, but probably a look at
how such datasets are normalized or otherwise pre-processed might be in order.
For the other datasets, PCA doesn't seem to distort the distributions as much.
However, as they are all image datasets, that doesn't necessarily mean that
there aren't problems lurking. Memo to self: find some non-image datasets, e.g.
the 20 Newsgroup dataset (available in the PaCMAP repo) for text categorization.

## 20NG

*December 25 2021*: I followed up on my memo to myself and have now looked at
the [20 Newsgroups](http://qwone.com/~jason/20Newsgroups/) dataset aka `20NG`,
but I will refer to it as `ng20` because variables which begin with a digit are
difficult or impossible to work with in most programming languages.

Currently, the PaCMAP repo has [20NG data](https://github.com/YingfanWang/PaCMAP/blob/f309fb8faef55a80e819ab91c8f81ac7d3937d9a/data/20NG.npy)
available, but this is a dense numpy matrix with 100 features, so presumably it
has already been reduced down by the Truncated SVD. As I want to compare this
with a higher rank decomposition, I processed the 20NG data from scratch.
The 20NG dataset is sparse, and PaCMAP does not support sparse input data (Annoy
does not support sparse distance calculations), so some kind of pre-processing
of the data to produce a dense matrix is required.

I took the dataset preparation from this [UMAP
example](https://umap-learn.readthedocs.io/en/latest/sparse.html#a-text-analysis-example)
(note that UMAP can work on sparse inputs directly because it uses
[pynndescent](https://github.com/lmcinnes/pynndescent/) to find nearest neighbors).
The simple procedure I followed was:

```python
import sklearn.datasets
import sklearn.feature_extraction.text
from sklearn.decomposition import TruncatedSVD

ng20v =  sklearn.datasets.fetch_20newsgroups_vectorized(subset="all")
ng20tfidf = sklearn.feature_extraction.text.TfidfTransformer(norm="l2").fit_transform(ng20v.data)

svd = TruncatedSVD(3000)
ng20tfidf_svd3000 = svd.fit_transform(ng20tfidf)
```

For the data here, I used L2 normalization because L1 (as used in the UMAP
example) results in some documents having very few near neighbors, which leads
to most of the histograms being squashed to the left side of the plots due to a
small number of large distances. The L2 normalization also gives distance
distributions and UMAP and PaCMAP embeddings which are closer to the results I
get when using the PaCMAP data. So I am sticking with L2 normalization for
consistency and legibility (but I have no opinion on whether the L1 or L2
normalization is "better" for this dataset).

As noted above in the "Input Preprocessing" section, I have some doubts that the
default settings for Truncated SVD in sklearn always produces the truncated SVD.
However, for this data it did a fine job: the dot products of each column with
the next column was always below $10^{-14}$, which is good enough for me. I also
tried with [irlba](http://bwlewis.github.io/irlba/) in R, but on this occasion
it was irlba that had trouble converging with its default parameters. In terms
of the distributions for the 100 and 3000 dimension distances, results seemed
very comparable, so I am sticking with the Python TruncatedSVD results below.

Below is a comparison of the distance distributions using 3000 columns vs 100.

| ng20 15s | ng20 mid |
|:----:|:---:
![ng20 15s](../img/pacmap/ng2015s.png)|![ng20 mid](../img/pacmap/ng20mid.png)
![ng20pca100 15s](../img/pacmap/ng20pca10015s.png)|![ng20pca100 mid](../img/pacmap/ng20pca100mid.png)

Again, there is a noticeable shift in the distributions. The point here isn't
whether keeping 3000 columns vs 100 columns results in better embeddings (the
100D result actually seems better and you would probably use the cosine metric
instead of Euclidean for this dataset), but that the distributions of distances
is noticeably shifted for some datasets, which is something to be aware of. If
you only look at the MNIST digits for example, you *won't* see this shift which
might lull you into believing that reducing dimensionality to 100 PCs gives you
a speed increase for the nearest neighbor search with no other downstream
effects.

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

* December 25 2021
  * Add the 20NG dataset to the analysis of the effect of PCA.
* December 24 2021
  * Add some more datasets for the nearest neighbor histograms.
* December 23 2021
  * Add nearest neighbor histograms.
* December 19 2021
  * Add equations for UMAP gradient.
  * Add plots of sample attractive and negative force constants.
* November 16 2021
  * Added some comments on the Truncated SVD implementation.
* November 14 2021
    * Created document.
