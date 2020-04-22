---
title: "LargeVis"
date: "April 19, 2020"
output:
  html_document:
    theme: cosmo
    toc: true
    toc_float:
      collapsed: false
---
Up: [Documentation Home](https://jlmelville.github.io/smallvis/).

LargeVis ([paper](https://arxiv.org/abs/1602.00370) and 
[github repo](https://github.com/lferry007/LargeVis)) is a clever method that
allows a t-SNE like visualization using stochastic gradient descent, which 
scales much better than the Barnes-Hut approximation to t-SNE. LargeVis does
not minimize the same cost function as t-SNE, but it is similar. There is an 
[R package](https://github.com/elbamos/largeVis) that is sadly no longer on 
CRAN, but can be installed from github easily enough. The `lvish` function that
is part of [uwot](https://cran.r-project.org/package=uwot) is close enough to
the real LargeVis for my purposes. Both of these packages replicate the 
stochastic gradient descent approach to optimization, which is the point of the
method: to be scalable to large datasets. As a result, I'm not aware of other
code that has implemented the LargeVis method in an exact gradient form.

## Theory

The theory of LargeVis is covered in plenty of depth over on the
[theory page](https://jlmelville.github.io/smallvis/theory.html). Here's a
summary. An exact gradient cost function looks like:

$$
C_{LV} = 
-\sum_{ij} p_{ij} \log w_{ij} 
-\gamma \sum_{ij} \log \left( 1 - w_{ij} \right)
$$

where $p_{ij}$ are input probabilities generated in the same way as t-SNE
(perplexity, averaging, normalization) and $w_{ij}$ are output weights, also
defined in the same way as t-SNE. $\gamma$ is weight factor for the repulsive
contributions. This is something that needs to be experimented with. The
default value in the SGD implementation ($\gamma = 7$) is way too large once you
start considering all pairs of points.

The gradient with respect to the output coordinates, $\mathbf{y}$ is:

$$
\frac{\partial C_{LV}}{\partial \mathbf{y_i}} = 
  4\sum_j^N \left[
    \frac{p_{ij}}{ 1 + d_{ij}^2 }
    -\frac{\gamma }{ \left( \epsilon + d_{ij}^2 \right) \left( 1 + d_{ij}^2 \right) } 
   \right]
   \left(\mathbf{y_i - y_j}\right)
$$

where $d_{ij}$ is the Euclidean output distance between coordinates. $\epsilon$
is a small positive value chosen to prevent division by zero. In the LargeVis
source code, it's set to 0.1.

In an email discussion, Dmitry Kobak made the clever observation that if you set 
$\epsilon = 1$ (practically that's not a small value so you probably *shouldn't*
do this, but bear with us), and you can write the gradient as:

$$
\frac{\partial C_{LV}}{\partial \mathbf{y_i}} = 
  \frac{4}{N}\sum_j^N \left(
    v_{ij}
    -
    N\gamma w_{ij} 
   \right)
   w_{ij} \left(\mathbf{y_i - y_j}\right)
$$

I have also used the substitution $p_{ij} = v_{ij} / N$. This expression 
resembles other methods' gradients quite closely:

* [t-distributed Elastic Embedding](https://jlmelville.github.io/smallvis/tee.htm),
where the t-EE parameter $\lambda = N\gamma$.
* if you choose $\gamma = 1/Z$ at each iteration (where $Z = \sum w_{ij}$), you 
get the t-SNE gradient.

## Un-normalized LargeVis

Does it make sense to be matching $p_{ij}$ to $w_{ij}$ in the cost function,
rather than $v_{ij}$? As it happens, the LargeVis source code doesn't actually
generate the $p_{ij}$ from the $v_{ij}$. It doesn't need to, because the 
absolute values aren't needed in the stochastic gradient descent method it uses:
the relative values are used to generate the probability of sampling a point,
and that is unaffected by scaling by $N$ or not. In the gradient calculation,
$p_{ij}$ is set to 1, as the magnitude of the attractive force has already been
accounted for by how often it's sampled.

If we replace $p_{ij}$ with $v_{ij}$ to create an un-normalized LargeVis, the
gradient is now:

$$
\frac{\partial C_{ULV}}{\partial \mathbf{y_i}} = 
  4\sum_j^N \left(
    v_{ij}
    -\frac{\gamma }{ \epsilon + d_{ij}^2 } 
   \right)w_{ij}
   \left(\mathbf{y_i - y_j}\right)
$$

and if you set $\gamma = 1$ this now resembles yet another method: the t-UMAP
variant of UMAP, with the following minor differences:

* the $v_{ij}$ are generated in a different way in t-UMAP, although this doesn't
make a huge difference.
* the repulsive part of the gradient is scaled by $1 - v_{ij}$ in the exact
gradient version of t-UMAP, but outside the k-nearest neighbors of point $i$,
$v_{ij}$ is very close to zero, so that isn't a big deal either.

Note that the resemblance between un-normalized LargeVis and t-UMAP does not
require setting $\epsilon = 1$.

## How to set $\gamma$

There's no explanation in LargeVis for what $\gamma$ is supposed to represent
theoretically, beyond that it's a weight on the repulsion. One way to view it is
to see it as a way to scale up the contribution from the missing repulsions in
stochastic gradient descent, i.e. as the negative sampling rate increases,
$\gamma$ should decrease, in which case $\gamma = 1$ is the sensible choice for
the full gradient. But in the LargeVis paper the cost function is introduced
before any mention of stochastic gradient descent or sampling so it seems like
it is intended to be an intrinsic part of the cost function. 

Whatever the intention, we can look at how the SGD implementation balances the
weight of the repulsive and attractive parts of the gradient to work out a rough
value for what $\gamma$ should be in the exact gradient case to get a similar
balance.

In the SGD implementation the update for a given point $i$, is based around a
gradient using one other "positive" point $j$, considering the attractive
interaction with an effective $p_{ij} = 1$. There are $n$ negative samples
(default = 5) whose contributions are scaled by $\gamma$ (default = 7). We could
therefore say that for the default SGD LargeVis case that:

$$
\frac{\gamma n}{p} = 35 \implies \gamma = \frac{35p}{n}
$$

where $p$ is the total weight of the positive interaction (in the SGD case 1) 
and $n$ the total number of repulsive (negative) interactions (in the SGD case,
5). 

In the exact gradient we consider all pairs of points to be involved in both
contributions. Assuming we are using the normalized version of LargeVis, then
the sum of the attractive contributions is 1. Meanwhile, all points contribute
to the repulsive contribution with a weight independent of
$p_{ij}$. Therefore for the purposes of the exact gradient, $p = 1$ and $n = N$.

$$
\gamma = \frac{35}{N}
$$

this would suggest an upper range of $\gamma$ between ~0.1-1 for `iris` and
~0.001-0.01 for `mnist6k`. Setting $\gamma$ between 0.001 and 1 might be a
good starting range for the exact gradient.

The practical utility of that calculation is less important than what it says
about the relative weighting of attractions and repulsions in LargeVis. The
choice of $\gamma = 7$ in the SGD optimization *looks* like it's implying that
repulsions should be up-weighted, but if this rough calculation is anywhere
close to reality, this translates to a $\gamma < 1$ in the full gradient
calculation, and therefore also the LargeVis cost function.

Depending on the choice of $\gamma$, this means that repulsions could be being
down-weighted in LargeVis compared to t-SNE, where, for the datasets I use with
`smallvis` and a perplexity of 40, the values of $\gamma$ typically range from
$1/N$ at the start of the optimization to ~0.1 at the end (for more on this see
the [t-EE page](https://jlmelville.github.io/smallvis/tee.htm)). If that's the 
case, this would explain the more compressed clusters seen in LargeVis compared
to t-SNE.

## Datasets

See the [Datasets](https://jlmelville.github.io/smallvis/datasets.html) page.

## Settings

Here's an example invocation for LargeVis with `iris`.

```R
res <- smallvis(iris, method = list("largevis", normalize = FALSE, gamma = 0.01, gr_eps = 0.1), perplexity = 40, Y_init = "spca", eta = 0.5, g2tol = 1e-7, min_cost = -Inf)
```

The choice of `gamma` is dependent on whether `normalize = FALSE` or `TRUE`. If
`TRUE`, then gamma needs to be scaled down by a factor of $N$ (where $N$ is the 
number of points in the dataset). Because the suitable range of `gamma` values
already starts quite small with `normalize = FALSE` (I would look between $1/N$
and 1), I would suggest sticking with `normalize = FALSE`, even though that
isn't how the cost function is written in the paper.

On the other hand, if you don't normalize, then the learning rate becomes more
sensitive to $N$. I suggest starting with something like $100 / N$. For some
choices of `gamma` and `gr_eps`, the LargeVis cost can be negative (which has
no particular significance), so you will want to set `min_cost = -Inf` to stop
the optimization stopping. If you set `normalize = TRUE`, then based on my
experience with t-EE, then `eta = 100, gr_eps = 1` is an ok choice if you don't 
want the hassle of setting a data dependent learning rate.

Finally, for other choices of `gamma` and `gr_eps`, although they converge just 
fine, you may observe that the usual cost change tolerance doesn't trigger. This
is because the cost function can actually increase by a very very small amount
(although the formatting of the `error:` when `verbose = TRUE` won't show this).
This is entirely a numeric issue and not due to an actual divergence (the 
layout doesn't perceptibly change), and the gradient converges to a small 2-norm 
as shown by the `||G||2` part of the logging. Set the `g2tol` under these
circumstances to trigger early stopping due to the gradient 2-norm falling below
a certain value. This is dataset dependent but values between `1e-5` to `1e-7`
are probably good, but conservative choices.

## Raw results

Here are some results using the settings above, and trying out a few 
difference combinations of whether to normalize, values of `gamma` and `gr_eps`.
I used a dataset-dependent value of `eta` for all results, even when
`normalize = TRUE`. I increased `max_iter = 5000` for all the results to take
into account that the learning rate might not be optimal. 

The top left result is normalized LargeVis with `gamma = 1, gr_eps = 0.1`. This
is quite close to a naive transfer of the SGD settings. On the top right are
the same settings, but not normalizing. On the bottom left are the un-normalized
LargeVis results with `gr_eps` reduced to `1e-6`. On the bottom right is
the un-normalized LargeVis `gr_eps = 1e-6` and a substantially reduced $\gamma$,
`gamma = 1e-3`, a value I suspect is a reasonable choice for full-gradient
LargeVis.

### iris

|                             |                           |
:----------------------------:|:--------------------------:
![iris lvg1e0.1p](../img/lv/iris_lvg1e0.1p.png)|![iris lvg1e0.1v](../img/lv/iris_lvg1e0.1v.png)
![iris lvg1e1e_6v](../img/lv/iris_lvg1e1e_6v.png)|![iris lvg1e_31e_6v](../img/lv/iris_lvg1e_31e_6v.png)


### s1k

|                             |                           |
:----------------------------:|:--------------------------:
![s1k lvg1e0.1p](../img/lv/s1k_lvg1e0.1p.png)|![s1k lvg1e0.1v](../img/lv/s1k_lvg1e0.1v.png)
![s1k lvg1e1e_6v](../img/lv/s1k_lvg1e1e_6v.png)|![s1k lvg1e_31e_6v](../img/lv/s1k_lvg1e_31e_6v.png)


### oli

|                             |                           |
:----------------------------:|:--------------------------:
![oli lvg1e0.1p](../img/lv/oli_lvg1e0.1p.png)|![oli lvg1e0.1v](../img/lv/oli_lvg1e0.1v.png)
![oli lvg1e1e_6v](../img/lv/oli_lvg1e1e_6v.png)|![oli lvg1e_31e_6v](../img/lv/oli_lvg1e_31e_6v.png)


### frey

|                             |                           |
:----------------------------:|:--------------------------:
![frey lvg1e0.1p](../img/lv/frey_lvg1e0.1p.png)|![frey lvg1e0.1v](../img/lv/frey_lvg1e0.1v.png)
![frey lvg1e1e_6v](../img/lv/frey_lvg1e1e_6v.png)|![frey lvg1e_31e_6v](../img/lv/frey_lvg1e_31e_6v.png)


### coil20

|                             |                           |
:----------------------------:|:--------------------------:
![coil20 lvg1e0.1p](../img/lv/coil20_lvg1e0.1p.png)|![coil20 lvg1e0.1v](../img/lv/coil20_lvg1e0.1v.png)
![coil20 lvg1e1e_6v](../img/lv/coil20_lvg1e1e_6v.png)|![coil20 lvg1e_31e_6v](../img/lv/coil20_lvg1e_31e_6v.png)


### mnist6k

|                             |                           |
:----------------------------:|:--------------------------:
![mnist6k lvg1e0.1p](../img/lv/mnist6k_lvg1e0.1p.png)|![mnist6k lvg1e0.1v](../img/lv/mnist6k_lvg1e0.1v.png)
![mnist6k lvg1e1e_6v](../img/lv/mnist6k_lvg1e1e_6v.png)|![mnist6k lvg1e_31e_6v](../img/lv/mnist6k_lvg1e_31e_6v.png)


### fashion6k

|                             |                           |
:----------------------------:|:--------------------------:
![fashion6k lvg1e0.1p](../img/lv/fashion6k_lvg1e0.1p.png)|![fashion6k lvg1e0.1v](../img/lv/fashion6k_lvg1e0.1v.png)
![fashion6k lvg1e1e_6v](../img/lv/fashion6k_lvg1e1e_6v.png)|![fashion6k lvg1e_31e_6v](../img/lv/fashion6k_lvg1e_31e_6v.png)

Obviously `gamma = 1` with normalized LargeVis turns out to be a terrible
idea. As anticipated, the repulsions are way too high and a circle of fairly 
uniformly spaced points arises. If we use the un-normalized setting, then we
get results that are much more intelligible, although which hint at the same
problems that t-EE had for the larger datasets, which is that starting with
too high a repulsion leads to problems.

Reducing the value of $\epsilon$ by setting `gr_eps = 1e-6` also causes problems
with `gamma = 1`. All the results are the same as their initialized coordinates,
although expanded. This isn't unique to the scaled PCA initialization: I tried
this with random initialization as used by the LargeVis SGD implementation and
a spectral initialization and the same thing happens.

The cause of this issue is that the small distances in the initialization,
combined with the relatively large gamma, leads to exceptionally large
gradients. This leads to a very large change in the coordinates, which on
subsequent iterations means which very large distances which causes
exceptionally small gradients, so nothing happens after the initial explosion.

In the SGD implementation, clipping occurs to prevent too large a gradient.
`smallvis` does not clip gradients, so to get around this, we have a choice
between initially using a more dispersed initialization, a larger `gr_eps`, a
smaller `gamma`, a smaller `eta` or a more sophisticated adaptive step size
optimization method that is more conservative during initial iterations.

The results at the bottom right show the effect of using a reduced `gamma`, and
it gives the most pleasing output by far. So it seems to be the case that the
choice of `gamma` is quite critical to getting good results.

## Annealed-$\gamma$ results

Given the similarity of the LargeVis gradient to t-EE, and given that t-EE
results improve notably if you start with a small value of its $\lambda$
parameter which entirely analogous to $\gamma$ parameter of LargeVis, we can
probably expect to see an improvement by doing the same here.

These runs use the same settings as for the previous results, except that I
started with $\gamma$ set to a power of 10 closest to $1/N$ (using arguments I
made on the t-EE page to do with what the initial repulsion in t-SNE is like),
ran a full optimization with `max_iter = 1000` iterations, then used the final
coordinates as input to a new run with $\gamma$ set to be ten times larger, and
repeated until $\gamma = 1$. For the datasets studied here, this involves 3-5
rounds of optimization, so this doesn't involve more iterations than that used
in the previous section.

### iris

|                             |                           |
:----------------------------:|:--------------------------:
![iris lv1e-9](../img/lv/iris_lv1e-9.png)|![iris lv0.001](../img/lv/iris_lv0.001.png)
![iris lv0.1](../img/lv/iris_lv0.1.png)|![iris lv1](../img/lv/iris_lv1.png)


### s1k

|                             |                           |
:----------------------------:|:--------------------------:
![s1k lv1e-9](../img/lv/s1k_lv1e-9.png)|![s1k lv0.001](../img/lv/s1k_lv0.001.png)
![s1k lv0.1](../img/lv/s1k_lv0.1.png)|![s1k lv1](../img/lv/s1k_lv1.png)


### oli

|                             |                           |
:----------------------------:|:--------------------------:
![oli lv1e-9](../img/lv/oli_lv1e-9.png)|![oli lv0.001](../img/lv/oli_lv0.001.png)
![oli lv0.1](../img/lv/oli_lv0.1.png)|![oli lv1](../img/lv/oli_lv1.png)


### frey

|                             |                           |
:----------------------------:|:--------------------------:
![frey lv1e-9](../img/lv/frey_lv1e-9.png)|![frey lv0.001](../img/lv/frey_lv0.001.png)
![frey lv0.1](../img/lv/frey_lv0.1.png)|![frey lv1](../img/lv/frey_lv1.png)


### coil20

|                             |                           |
:----------------------------:|:--------------------------:
![coil20 lv1e-9](../img/lv/coil20_lv1e-9.png)|![coil20 lv0.001](../img/lv/coil20_lv0.001.png)
![coil20 lv0.1](../img/lv/coil20_lv0.1.png)|![coil20 lv1](../img/lv/coil20_lv1.png)


### mnist6k

|                             |                           |
:----------------------------:|:--------------------------:
![mnist6k lv1e-9](../img/lv/mnist6k_lv1e-9.png)|![mnist6k lv0.001](../img/lv/mnist6k_lv0.001.png)
![mnist6k lv0.1](../img/lv/mnist6k_lv0.1.png)|![mnist6k lv1](../img/lv/mnist6k_lv1.png)


### fashion6k

|                             |                           |
:----------------------------:|:--------------------------:
![fashion6k lv1e-9](../img/lv/fashion6k_lv1e-9.png)|![fashion6k lv0.001](../img/lv/fashion6k_lv0.001.png)
![fashion6k lv0.1](../img/lv/fashion6k_lv0.1.png)|![fashion6k lv1](../img/lv/fashion6k_lv1.png)


Visually, things are looking much better. And they are pretty consistent across
all values of `gr_eps` except for `gr_eps = 1`, where there seems to be more
effective repulsion or it's just that the gradient gets smaller earlier in the
optimization and progress halts sooner.

Note that there *is* a difference in the final costs, though, with larger
values of `gr_eps` causing the costs to be negative. Because you want `gr_eps`
to only be large enough to avoid division by zero and avoid optimization issues 
due to large changes in gradient value, I recommend setting `gr_eps = 1e-9` to 
avoid the potential for negative cost values. Of course, this is only an option
if you commit to using a small value of `gamma` initially, or you will run into
the exploding gradients problem we saw earlier.

## Early Exaggeration

The annealed-$\gamma$ results tell us that we can probably use whatever value
of $\gamma$ we want if we are prepared to get there by starting at a small value
and increasing every few iterations. Starting with a reduced repulsion is 
equivalent to a increased attraction, which is what early exaggeration does.
As this is already part of the `smallvis` interface, maybe this is worth trying.

Below are some experiments using `gamma = 1e-3, gr_eps = 1e-3`. The top left
image is what you get without any early exaggeration. Top right is with 
`exaggeration_factor = 10, stop_lying_iter = 100`, turning off exaggeration
fairly early, and then bottom left is with `stop_lying_iter = 250`, the more
typical exaggeration time used in BH t-SNE. The `exaggeration_factor = 10`, 
is the equivalent of starting with `gamma` set to an order of magnitude
lower than we end with. This is more aggressive than the typical t-SNE 
exaggeration factor of 4. The bottom right image sees if there is an advantage
from an even stronger `exaggeration_factor = 100`.

A typical setting is:

```R
iris_lvee <- smallvis(iris, method=list("largevis", gamma = 1e-3, gr_eps = 1e-3, normalize = FALSE), Y_init = "spca", eta = 0.15, perplexity = 40, g2tol = 1e-7, min_cost = -Inf, exaggeration_factor = 10, stop_lying_iter = 100)
```

### iris

|                             |                           |
:----------------------------:|:--------------------------:
![iris lvg1e_31e_3ee0](../img/lv/iris_lvg1e_31e_3ee0.png)|![iris lvg1e_31e_3ee100](../img/lv/iris_lvg1e_31e_3ee100.png)
![iris lvg1e_31e_3ee250](../img/lv/iris_lvg1e_31e_3ee250.png)|![iris lvg1e_31e_3ee100_100](../img/lv/iris_lvg1e_31e_3ee100_100.png)


### s1k

|                             |                           |
:----------------------------:|:--------------------------:
![s1k lvg1e_31e_3ee0](../img/lv/s1k_lvg1e_31e_3ee0.png)|![s1k lvg1e_31e_3ee100](../img/lv/s1k_lvg1e_31e_3ee100.png)
![s1k lvg1e_31e_3ee250](../img/lv/s1k_lvg1e_31e_3ee250.png)|![s1k lvg1e_31e_3ee100_100](../img/lv/s1k_lvg1e_31e_3ee100_100.png)


### oli

|                             |                           |
:----------------------------:|:--------------------------:
![oli lvg1e_31e_3ee0](../img/lv/oli_lvg1e_31e_3ee0.png)|![oli lvg1e_31e_3ee100](../img/lv/oli_lvg1e_31e_3ee100.png)
![oli lvg1e_31e_3ee250](../img/lv/oli_lvg1e_31e_3ee250.png)|![oli lvg1e_31e_3ee100_100](../img/lv/oli_lvg1e_31e_3ee100_100.png)


### frey

|                             |                           |
:----------------------------:|:--------------------------:
![frey lvg1e_31e_3ee0](../img/lv/frey_lvg1e_31e_3ee0.png)|![frey lvg1e_31e_3ee100](../img/lv/frey_lvg1e_31e_3ee100.png)
![frey lvg1e_31e_3ee250](../img/lv/frey_lvg1e_31e_3ee250.png)|![frey lvg1e_31e_3ee100_100](../img/lv/frey_lvg1e_31e_3ee100_100.png)


### coil20

|                             |                           |
:----------------------------:|:--------------------------:
![coil20 lvg1e_31e_3ee0](../img/lv/coil20_lvg1e_31e_3ee0.png)|![coil20 lvg1e_31e_3ee100](../img/lv/coil20_lvg1e_31e_3ee100.png)
![coil20 lvg1e_31e_3ee250](../img/lv/coil20_lvg1e_31e_3ee250.png)|![coil20 lvg1e_31e_3ee100_100](../img/lv/coil20_lvg1e_31e_3ee100_100.png)


### mnist6k

|                             |                           |
:----------------------------:|:--------------------------:
![mnist6k lvg1e_31e_3ee0](../img/lv/mnist6k_lvg1e_31e_3ee0.png)|![mnist6k lvg1e_31e_3ee100](../img/lv/mnist6k_lvg1e_31e_3ee100.png)
![mnist6k lvg1e_31e_3ee250](../img/lv/mnist6k_lvg1e_31e_3ee250.png)|![mnist6k lvg1e_31e_3ee100_100](../img/lv/mnist6k_lvg1e_31e_3ee100_100.png)


### fashion6k

|                             |                           |
:----------------------------:|:--------------------------:
![fashion6k lvg1e_31e_3ee0](../img/lv/fashion6k_lvg1e_31e_3ee0.png)|![fashion6k lvg1e_31e_3ee100](../img/lv/fashion6k_lvg1e_31e_3ee100.png)
![fashion6k lvg1e_31e_3ee250](../img/lv/fashion6k_lvg1e_31e_3ee250.png)|![fashion6k lvg1e_31e_3ee100_100](../img/lv/fashion6k_lvg1e_31e_3ee100_100.png)

For most datasets, early exaggeration doesn't have much effect, which is perhaps
not too surprising as small datasets have more space to arrange themselves. For
`mnist6k` and `fashion6k` there is a more obvious improvement to the final cost,
and in `mnist6k`, the visual output looks a bit better. Increasing the 
`exaggeration_factor = 100` doesn't seem to help. This means that the range
of effective `gamma` you can use even with early exaggeration is still quite 
small. For example, you won't get good results with `gamma = 1` this way 
(results are identical to not using early exaggeration). Increasing the 
`exaggeration_factor` further will lead to issues with the learning rate being
too large in the initial exaggeration phase, as discussed in
[Clustering with t-SNE, provably](https://doi.org/10.1137/18M1216134).

## Normalized Results vs t-SNE

I've used un-normalized LargeVis for most of the results above, but to
demonstrate that the normalized version of LargeVis, here are some results with
the following settings:

```R
iris_lv <- smallvis(iris, method = list("largevis", normalize = TRUE, gamma = 10 / nrow(iris) ^ 2, gr_eps = 0.1), eta = nrow(iris) / 100, perplexity = 40, g2tol = 1e-7, min_cost = -Inf, exaggeration_factor = 10, Y_init = "spca")
```

This uses a data-dependent $\gamma = 10 / N^2$ and learning rate of 
$\eta = N / 100$, although a fixed $\eta = 100$ works nearly as well for this
combination of $\epsilon$ and $\gamma$. With the `exaggeration_factor = 10`,
this means that the repulsion weight in the gradient will be $1 / N$ (remember
that the repulsion weight is $\gamma N$ in the normalized version of LargeVis)
during early exaggeration, which is approximately the initial t-SNE repulsion
weight. After early exaggeration the repulsion is reduced by a factor of 10, but
remains well below the repulsion weight of around 0.05-0.3 that t-SNE achieves
with these datasets and settings. So I anticipate that the LargeVis results will
be comparable to, but more compressed than, the t-SNE results.

### iris

|                             |                           |
:----------------------------:|:--------------------------:
![iris ee100P](../img/lv/iris_ee100P.png)|![iris tsne](../img/lv/iris_tsne.png)


### s1k

|                             |                           |
:----------------------------:|:--------------------------:
![s1k ee100P](../img/lv/s1k_ee100P.png)|![s1k tsne](../img/lv/s1k_tsne.png)


### oli

|                             |                           |
:----------------------------:|:--------------------------:
![oli ee100P](../img/lv/oli_ee100P.png)|![oli tsne](../img/lv/oli_tsne.png)


### frey

|                             |                           |
:----------------------------:|:--------------------------:
![frey ee100P](../img/lv/frey_ee100P.png)|![frey tsne](../img/lv/frey_tsne.png)


### coil20

|                             |                           |
:----------------------------:|:--------------------------:
![coil20 ee100P](../img/lv/coil20_ee100P.png)|![coil20 tsne](../img/lv/coil20_tsne.png)


### mnist6k

|                             |                           |
:----------------------------:|:--------------------------:
![mnist6k ee100P](../img/lv/mnist6k_ee100P.png)|![mnist6k tsne](../img/lv/mnist6k_tsne.png)


### fashion6k

|                             |                           |
:----------------------------:|:--------------------------:
![fashion6k ee100P](../img/lv/fashion6k_ee100P.png)|![fashion6k tsne](../img/lv/fashion6k_tsne.png)

This is more or less what I expected. Looking at the range of the axes, the
t-SNE results are more expanded but the relative location of the clusters are
quite comparable visually (with the exception of `iris` where the opposite seems
true). You can see that for `mnist6k` and `fashion6k`, the clusters are more
compressed, so it's probable that the degree of repulsions is size-dependent.

## Comparison with SGD results

Finally, let's compare the results in `smallvis` to the `lvish` function in
[uwot](https://cran.r-project.org/package=uwot), which should give results
fairly close to the "true" LargeVis results. I'm not aiming to reproduce the
`lvish` results in `smallvis`

In order to make things more similar, for both `smallvis` and `lvish` we will
use a knn kernel, where for a given perplexity $PP$, the nearest $PP$ neighbors
of point $i$ get an input affinity of $1/PP$ and 0 everywhere else. The
resulting matrix is symmetrized as usual to get the final $v_{ij}$ values.

For `smallvis`, the only difference between the previous set of results is 
adding `inp_kernel = "knn"` to the `method` list. For `lvish`, the settings
were:

```R
iris_lvish <- lvish(iris, n_epochs = 1000, kernel = "knn", nn_method = "FNN", perplexity = 40, init = "spca")
```

The choice of `n_epochs = 1000` was fairly arbitrary, except it's double what
the default UMAP choice is for small datasets, so hopefully the coordinates
are reasonably well converged. Note that `lvish` is a stochastic method, so
exact results are dependent on the random number seed.

The left-hand result is the `smallvis` result, the right-hand result is due to
`lvish`.

### iris

|                             |                           |
:----------------------------:|:--------------------------:
![iris lveek](../img/lv/iris_lveek.png)|![iris lvishk](../img/lv/iris_lvishk.png)


### s1k

|                             |                           |
:----------------------------:|:--------------------------:
![s1k lveek](../img/lv/s1k_lveek.png)|![s1k lvishk](../img/lv/s1k_lvishk.png)


### oli

|                             |                           |
:----------------------------:|:--------------------------:
![oli lveek](../img/lv/oli_lveek.png)|![oli lvishk](../img/lv/oli_lvishk.png)


### frey

|                             |                           |
:----------------------------:|:--------------------------:
![frey lveek](../img/lv/frey_lveek.png)|![frey lvishk](../img/lv/frey_lvishk.png)


### coil20

|                             |                           |
:----------------------------:|:--------------------------:
![coil20 lveek](../img/lv/coil20_lveek.png)|![coil20 lvishk](../img/lv/coil20_lvishk.png)


### mnist6k

|                             |                           |
:----------------------------:|:--------------------------:
![mnist6k lveek](../img/lv/mnist6k_lveek.png)|![mnist6k lvishk](../img/lv/mnist6k_lvishk.png)


### fashion6k

|                             |                           |
:----------------------------:|:--------------------------:
![fashion6k lveek](../img/lv/fashion6k_lveek.png)|![fashion6k lvishk](../img/lv/fashion6k_lvishk.png)

Results are most certainly not identical, but they're not completely different
from each other. The chosen value of `gamma` could perhaps be tweaked to make
some visualizations line up better, for example `oli` seems to be over-expanded
in `smallvis`.

Some differences may be irreconcilable: e.g. for `frey`, the separation between
the green cluster and the rest of the data is always much more pronounced with
`lvish` across multiple random seeds, different choices of `n_epochs` and
`learning_rate`. For that dataset, even when using `smallvis` to optimize the
output of `lvish`, I was unable to find a value of `gamma` which didn't result
in a large amount of re-arrangement of the coordinates. For `coil20` there seems
to be more expansion with the SGD results and the three clusters (black, cyan
and blue) on the right hand side of the plots are pushed further, similarly to
the situation with `frey`.

`mnist6k` and `fashion6k` alsoshow a slightly different relative arrangement of
clusters. This may be a limitation of the SGD implementation in terms of
representing the LargeVis cost function as implemented in `smallvis`, but this
is a mystery I will leave for another time. If it turns out to be due to a bug
in `uwot`, I will regenerate these results.

However, we see similar overall ranges of coordinates (based on looking at the
axis labels), so this is good enough evidence for me to confirm that the 
LargeVis repulsion is substantially down-weighted compared to t-SNE.

## Conclusions

LargeVis is not designed to be a practical visualization method outside of its
stochastic gradient descent implementation, so the fact that it's not as easy
to use in `smallvis` compared to t-SNE isn't very surprising.

Some general observations are:

* The balance of repulsions versus attractions in the t-SNE family of methods
is critical.
* Methods which don't start the optimization with attractions hugely domination
will have a hard time producing useful results.
* t-SNE automatically increases the repulsion as the optimization proceeds, 
which is why it can produce output with more expanded clusters successfully.
* Early exaggeration helps t-SNE with this. It's less useful for methods like
UMAP or LargeVis because the repulsions don't get dynamically re-weighted.
* As presented in the LargeVis paper, it's easy to be misled by the $\gamma$ 
value in the cost function. Its value > 1 in the stochastic gradient descent
optimization, suggests it is up-weighting repulsions relative to the 
attractive forces. 
* Once you recast the SGD settings in terms of the exact gradient, it's actually
**strongly down-weighting the repulsions**.
* This is still true even when the SGD $\gamma = 1$, which is effectively what
UMAP does.
* This implies that with UMAP's default settings, its cost function actually
contains an implicit $\gamma \ll 1$.
* In t-SNE, exaggeration is the same as decreasing the repulsive interactions,
which is equivalent to introducing $\gamma < 1$ into the repulsive part of the
t-SNE gradient (ignoring a constant scaling of the gradient which can be
absorbed into the learning rate). This is why t-SNE with late exaggeration looks
like UMAP and LargeVis.
* Therefore t-SNE with exaggeration permanently turned on should give you a
close equivalent to LargeVis or even UMAP visualizations, subject to which
implementation provides you with the scalability and features you need.

## Acknowledgement

Dmitry Kobak has been asking me to look at how the `smallvis` version LargeVis
behaves for ages, as well as asking "don't you think it's odd that t-SNE with
exaggerations looks so much like UMAP?", so I have him to thank for prodding
me into looking at this in a bit more detail.

Up: [Documentation Home](https://jlmelville.github.io/smallvis/).
