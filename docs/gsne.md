---
title: "g-SNE"
date: "February 28, 2019"
output:
  html_document:
    theme: cosmo
    toc: true
    toc_float:
      collapsed: false
---

Up: [Documentation Home](https://jlmelville.github.io/smallvis/).

The [g-SNE](https://doi.org/10.1101/331611) paper (with 
[matlab code on github](https://github.com/gyrheart/gsne)) suggests a
modification of t-SNE to add a second divergence that is more focused on
matching similarities due to longer distances. The idea is that a weighted
combination of the two will balance global and local structure.

The g-SNE cost function is:

$$C = \sum_{ij} p_{ij} \log \frac{p_{ij}}{q_{ij}} + \lambda \sum_{ij} \hat{p}_{ij} \log \frac{\hat{p_{ij}}}{\hat{q}_{ij}}$$

with $\lambda$ being a weighting factor as in NeRV, the un-hatted part of the
cost function being the same as t-SNE, and the hatted quantities referring to
using the reciprocal of the usual t-SNE kernel:

$$\hat{w}_{ij} = 1 + d_{ij}^2 = \frac{1}{w_{ij}}$$

The gradient is:

$$\frac{\partial C}{\partial \mathbf{y_i}} = 
  4
  \sum_j
  \left[
    \left(p_{ij} - q_{ij}\right) -\lambda \left( \hat{p}_{ij} - \hat{q}_{ij} \right)
  \right]w_{ij}
  \left(
   \mathbf{y_i - y_j}
  \right)
$$

As was also apparent from the cost function, you can see that this reduces back 
to t-SNE when $\lambda = 0$.

## Datasets

See the [Datasets](https://jlmelville.github.io/smallvis/datasets.html) page.

## Settings

To use g-SNE in smallvis, use `method = "gsne"`, which uses a default 
$\lambda = 1$, which is used in most of the datasets in the paper. Different
values of $\lambda$ can be specified by using e.g. 
`method = list("gsne", lambda = 5)`. In the g-SNE paper, values of `lambda` 
between 0 and 10 are considered, so we will look at `1`, `2.5`, `5` and `10`.

Other settings will be kept standard. Early exaggeration isn't mentioned in 
the paper, and in my experiments it made no difference to the results, but
it's on by default in the code on github, so I've also used it to here. The
exaggeration is *not* applied to $\hat{P}$ in either smallvis or the matlab
implementation.

```R
iris_gsne <- smallvis(iris, method = list("gsne", lambda = 5), perplexity = 40, eta = 100, Y_init = "spca", exaggeration_factor = 4)
```

An obvious point of comparison is t-SNE, which is the same as g-SNE with 
$\lambda = 0$. For preserving global distances, we can use metric MDS. We've
looked at [MMDS elsewhere](https://jlmelville.github.io/smallvis/mmds.html), 
but in this case we will use the t-SNE scaling (not that it makes much 
difference).

```R
iris_tsne <- smallvis(iris, perplexity = 40, eta = 100, Y_init = "spca", exaggeration_factor = 4)
iris_mmds <- smallvis(iris, eta = 1e-5, max_iter = 2000, Y_init = "spca", method = "mmds")

```

## Results

For each dataset, from left-to-right and top-to-bottom, `lambda` increased from
`0` (t-SNE) to `1`, `2.5`, `5` and `10`. The bottom right image is MMDS.

### iris

|                             |                           |
:----------------------------:|:--------------------------:
![iris tsne](../img/gsne/iris_tsne.png)|![iris gsnexl1](../img/gsne/iris_gsnexl1.png)
![iris gsnexl2.5](../img/gsne/iris_gsnexl2.5.png)|![iris gsnexl5](../img/gsne/iris_gsnexl5.png)
![iris gsnexl10](../img/gsne/iris_gsnexl10.png)|![iris mmds](../img/gsne/iris_mmds.png)

### s1k

|                             |                           |
:----------------------------:|:--------------------------:
![s1k tsne](../img/gsne/s1k_tsne.png)|![s1k gsnexl1](../img/gsne/s1k_gsnexl1.png)
![s1k gsnexl2.5](../img/gsne/s1k_gsnexl2.5.png)|![s1k gsnexl5](../img/gsne/s1k_gsnexl5.png)
![s1k gsnexl10](../img/gsne/s1k_gsnexl10.png)|![s1k mmds](../img/gsne/s1k_mmds.png)

### oli

|                             |                           |
:----------------------------:|:--------------------------:
![oli tsne](../img/gsne/oli_tsne.png)|![oli gsnexl1](../img/gsne/oli_gsnexl1.png)
![oli gsnexl2.5](../img/gsne/oli_gsnexl2.5.png)|![oli gsnexl5](../img/gsne/oli_gsnexl5.png)
![oli gsnexl10](../img/gsne/oli_gsnexl10.png)|![oli mmds](../img/gsne/oli_mmds.png)

### frey

|                             |                           |
:----------------------------:|:--------------------------:
![frey tsne](../img/gsne/frey_tsne.png)|![frey gsnexl1](../img/gsne/frey_gsnexl1.png)
![frey gsnexl2.5](../img/gsne/frey_gsnexl2.5.png)|![frey gsnexl5](../img/gsne/frey_gsnexl5.png)
![frey gsnexl10](../img/gsne/frey_gsnexl10.png)|![frey mmds](../img/gsne/frey_mmds.png)

### coil20

|                             |                           |
:----------------------------:|:--------------------------:
![coil20 tsne](../img/gsne/coil20_tsne.png)|![coil20 gsnexl1](../img/gsne/coil20_gsnexl1.png)
![coil20 gsnexl2.5](../img/gsne/coil20_gsnexl2.5.png)|![coil20 gsnexl5](../img/gsne/coil20_gsnexl5.png)
![coil20 gsnexl10](../img/gsne/coil20_gsnexl10.png)|![coil20 mmds](../img/gsne/coil20_mmds.png)

### mnist6k

|                             |                           |
:----------------------------:|:--------------------------:
![mnist6k tsne](../img/gsne/mnist6k_tsne.png)|![mnist6k gsnexl1](../img/gsne/mnist6k_gsnexl1.png)
![mnist6k gsnexl2.5](../img/gsne/mnist6k_gsnexl2.5.png)|![mnist6k gsnexl5](../img/gsne/mnist6k_gsnexl5.png)
![mnist6k gsnexl10](../img/gsne/mnist6k_gsnexl10.png)|![mnist6k mmds](../img/gsne/mnist6k_mmds.png)

### fashion6k

|                             |                           |
:----------------------------:|:--------------------------:
![fashion6k tsne](../img/gsne/fashion6k_tsne.png)|![fashion6k gsnexl1](../img/gsne/fashion6k_gsnexl1.png)
![fashion6k gsnexl2.5](../img/gsne/fashion6k_gsnexl2.5.png)|![fashion6k gsnexl5](../img/gsne/fashion6k_gsnexl5.png)
![fashion6k gsnexl10](../img/gsne/fashion6k_gsnexl10.png)|![fashion6k mmds](../img/gsne/fashion6k_mmds.png)

## Conclusions

Except for `iris`, it's clear that, in general, increasing $\lambda$ makes the
g-SNE progressively approach the sort of layout we'd expect with MMDS. I'm 
just not sure that I particularly like any of the intermediate results more
than just sticking with the t-SNE results. I think starting from the scaled PCA
initialization gives you a lot of the global layout you would want. `lambda`
values between 1 and 2.5 seem the best place to start exploring.

On the other hand, the authors use a 
[topological clique analysis involving Betti numbers](https://arxiv.org/abs/1502.06172) 
for some gene expression data and conclude that the g-SNE plot represents the
data better. Unfortunately, the plots are rather small, and for a second example
which required a value `lambda = 5` for optimal results, the plots aren't shown.

Up: [Documentation Home](https://jlmelville.github.io/smallvis/).
