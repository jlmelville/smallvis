---
title: "UMAP/t-SNE hybrids"
date: "February 17, 2018"
output:
  html_document:
    theme: cosmo
    toc: true
    toc_float:
      collapsed: false
---

This a companion to the experiments I did with 
[UMAP](https://jlmelville.github.io/smallvis/umap.html). Here I'm going to see
if borrowing some of UMAP's way of doing things affects t-SNE performance
measurably. See the previous UMAP page and the 
[theory](https://jlmelville.github.io/smallvis/theory.html) for details on UMAP.
This work was inspired by the release of the 
[UMAP paper](https://arxiv.org/abs/1802.03426), which shows a very nice result
for the COIL-20 data with UMAP. As of this writing, the settings used for that
result were not available, so this is an attempt to see if there's something
specific to one of the differences of UMAP from t-SNE that causes that 
difference, e.g. the method of generating the input our output weights. If so,
perhaps t-SNE could produce a result nearly as good, just by borrowing that 
feature.

By borrowing individual changes from UMAP we can define some new versions of 
t-SNE:

* SKD t-SNE: Instead of a perplexity-based calibration to generate the input 
weights, it uses smooth knn distances.
* u-SNE: Instead of the t-distribution to define the output kernel, the UMAP
function is $w_{ij} = 1 / \left(1 + ad_{ij}^{2b}\right)$. When using a function
like this, both LargeVis and UMAP add a value between `0.001` to `0.1` to the 
denominator of the resulting gradient to both avoid division by zero 
and to avoid rather odd-looking results (at least in the way `smallvis` 
calculates the non-stochastic gradient). But when used with SNE (perhaps due
to normalization?) I could get away with a much smaller value, which here
is `1e-10`.
* CE t-SNE: The fuzzy set cross-entropy is used instead of the KL divergence. 
In the way it's defined by UMAP, this adds an extra non-constant term to the 
cost function.
* UMAP also uses a normalized Laplacian for initialization, but as discussed
under [t-SNE initialization](https://jlmelville.github.io/smallvis/init.html),
this gives very similar results to using Laplacian Eigenmaps, which in turn
doesn't give very different results on most datasets to just using the scaled
PCA result, so we won't pursue that further here.

## Datasets

See the [Datasets](https://jlmelville.github.io/smallvis/datasets.html) page.

## Settings

For settings, we'll use two sets. First, the ones given in the 
[original t-SNE paper](http://www.jmlr.org/papers/v9/vandermaaten08a.html), 
except instead of random initialization, we'll use the scaled PCA approach. Also,
I've doubled the number of iterations allowed from 1000 to 2000 to avoid any
problems with convergence with these new methods. We'll also generate t-SNE 
results to compare with these variations.

For the second set of results, I'll set the perplexity to 15, which is closer
to the UMAP defaults.

```
iris_usne <- smallvis(iris, scale = FALSE, perplexity = 40, Y_init = "spca", method = list("usne", gr_eps = 1e-10), eta = 100, exaggeration_factor = 4, stop_lying_iter = 50, max_iter = 2000, epoch = 50, tol = 1e-5)

iris_skdtsne <- smallvis(iris, scale = FALSE, perplexity = 40, Y_init = "spca", method = "skdtsne", eta = 100, exaggeration_factor = 4, stop_lying_iter = 50, max_iter = 2000, epoch = 50, tol = 1e-5)

iris_cetsne <- smallvis(iris, scale = FALSE, perplexity = 40, Y_init = "spca", method = "cetsne", eta = 100, exaggeration_factor = 4, stop_lying_iter = 50, max_iter = 2000, epoch = 50, tol = 1e-5)

iris_tsne <- smallvis(iris, scale = FALSE, perplexity = 40, Y_init = "spca", method = "tsne", eta = 100, exaggeration_factor = 4, stop_lying_iter = 50, max_iter = 2000, epoch = 50, tol = 1e-5)
```

## Evaluation

Apart from visualizing the results, the mean neighbor preservation of the k
closest neighbors is used to provide a rough quantification of the quality of
the result. For the first set of results, with perplexity = 40, k = 40 and is 
labelled as `mnp@40`. The second set of results with perplexity = 15, k = 15
and is labelled as `mnp@15`.

## Results: Perplexity 40

### iris

|                             |                           |
:----------------------------:|:--------------------------:
![iris USNE](../img/umaptsne/iris_usne.png)|![iris SKD t-SNE](../img/umaptsne/iris_skdtsne.png)
![iris CE t-SNE DY](../img/umaptsne/iris_cetsne.png)|![iris t-SNE](../img/umaptsne/iris_tsne.png)

### s1k

|                             |                           |
:----------------------------:|:--------------------------:
![s1k USNE](../img/umaptsne/s1k_usne.png)|![s1k SKD t-SNE](../img/umaptsne/s1k_skdtsne.png)
![s1k CE t-SNE DY](../img/umaptsne/s1k_cetsne.png)|![s1k t-SNE](../img/umaptsne/s1k_tsne.png)

### oli

|                             |                           |
:----------------------------:|:--------------------------:
![oli USNE](../img/umaptsne/oli_usne.png)|![oli SKD t-SNE](../img/umaptsne/oli_skdtsne.png)
![oli CE t-SNE DY](../img/umaptsne/oli_cetsne.png)|![oli t-SNE](../img/umaptsne/oli_tsne.png)


### frey

|                             |                           |
:----------------------------:|:--------------------------:
![frey USNE](../img/umaptsne/frey_usne.png)|![frey SKD t-SNE](../img/umaptsne/frey_skdtsne.png)
![frey CE t-SNE DY](../img/umaptsne/frey_cetsne.png)|![frey t-SNE](../img/umaptsne/frey_tsne.png)

### coil20

|                             |                           |
:----------------------------:|:--------------------------:
![coil20 USNE](../img/umaptsne/coil20_usne.png)|![coil20 SKD t-SNE](../img/umaptsne/coil20_skdtsne.png)
![coil20 CE t-SNE DY](../img/umaptsne/coil20_cetsne.png)|![coil20 t-SNE](../img/umaptsne/coil20_tsne.png)

### mnist

|                             |                           |
:----------------------------:|:--------------------------:
![mnist USNE](../img/umaptsne/mnist_usne.png)|![mnist SKD t-SNE](../img/umaptsne/mnist_skdtsne.png)
![mnist CE t-SNE DY](../img/umaptsne/mnist_cetsne.png)|![mnist t-SNE](../img/umaptsne/mnist_tsne.png)

### fashion

|                             |                           |
:----------------------------:|:--------------------------:
![fashion USNE](../img/umaptsne/fashion_usne.png)|![fashion SKD t-SNE](../img/umaptsne/fashion_skdtsne.png)
![fashion CE t-SNE DY](../img/umaptsne/fashion_cetsne.png)|![fashion t-SNE](../img/umaptsne/fashion_tsne.png)

## Results: Perplexity 15

### iris

|                             |                           |
:----------------------------:|:--------------------------:
![iris USNE](../img/umaptsne/iris_usne15.png)|![iris SKD t-SNE](../img/umaptsne/iris_skdtsne15.png)
![iris CE t-SNE DY](../img/umaptsne/iris_cetsne15.png)|![iris t-SNE](../img/umaptsne/iris_tsne15.png)

### s1k

|                             |                           |
:----------------------------:|:--------------------------:
![s1k USNE](../img/umaptsne/s1k_usne15.png)|![s1k SKD t-SNE](../img/umaptsne/s1k_skdtsne15.png)
![s1k CE t-SNE DY](../img/umaptsne/s1k_cetsne15.png)|![s1k t-SNE](../img/umaptsne/s1k_tsne15.png)

### oli

|                             |                           |
:----------------------------:|:--------------------------:
![oli USNE](../img/umaptsne/oli_usne15.png)|![oli SKD t-SNE](../img/umaptsne/oli_skdtsne15.png)
![oli CE t-SNE DY](../img/umaptsne/oli_cetsne15.png)|![oli t-SNE](../img/umaptsne/oli_tsne15.png)


### frey

|                             |                           |
:----------------------------:|:--------------------------:
![frey USNE](../img/umaptsne/frey_usne15.png)|![frey SKD t-SNE](../img/umaptsne/frey_skdtsne15.png)
![frey CE t-SNE DY](../img/umaptsne/frey_cetsne15.png)|![frey t-SNE](../img/umaptsne/frey_tsne15.png)

### coil20

|                             |                           |
:----------------------------:|:--------------------------:
![coil20 USNE](../img/umaptsne/coil20_usne15.png)|![coil20 SKD t-SNE](../img/umaptsne/coil20_skdtsne15.png)
![coil20 CE t-SNE DY](../img/umaptsne/coil20_cetsne15.png)|![coil20 t-SNE](../img/umaptsne/coil20_tsne15.png)

### mnist

|                             |                           |
:----------------------------:|:--------------------------:
![mnist USNE](../img/umaptsne/mnist_usne15.png)|![mnist SKD t-SNE](../img/umaptsne/mnist_skdtsne15.png)
![mnist CE t-SNE DY](../img/umaptsne/mnist_cetsne15.png)|![mnist t-SNE](../img/umaptsne/mnist_tsne15.png)

### fashion

|                             |                           |
:----------------------------:|:--------------------------:
![fashion USNE](../img/umaptsne/fashion_usne15.png)|![fashion SKD t-SNE](../img/umaptsne/fashion_skdtsne15.png)
![fashion CE t-SNE DY](../img/umaptsne/fashion_cetsne15.png)|![fashion t-SNE](../img/umaptsne/fashion_tsne15.png)

## Conclusions

None of these changes have very much effect on t-SNE. The biggest effect is by
using the UMAP output kernel, but the overall shape of the embedding is 
unaffected with the default kernel settings. It also noticeably increases the
computational time due to the more complex gradient, so unless you really like
the more compressed clusters, t-SNE still wins out here.

The smoothed k-nearest neighbor distances is also quite interesting as an 
alternative to the gaussian kernel, producing slightly smaller, more 
well-separated clusters (compare the `fashion` and `mnist` dataset with "t-SNE"
and "SKD t-SNE"). This might be due to the fact that the input weights are set to
zero outside the k-nearest neighbors.

