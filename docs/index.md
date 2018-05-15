---
title: "smallvis"
output:
  html_document:
    theme: cosmo
---

Documentation related to the [smallvis](https://github.com/jlmelville/smallvis) 
R package:

## Installing

```R
install.packages("devtools")
devtools::install_github("jlmelville/smallvis/smallvis")
library(smallvis)
```

## Prerequisites

* The [vizier](https://github.com/jlmelville/vizier)
package to plot the coordinates during optimization. It's not on CRAN, and 
therefore requires a fairly new version of 
[devtools](https://cran.r-project.org/package=devtools) (1.9 or greater) to 
install this as a dependency from github. 
* Similarly, the [mize](https://github.com/jlmelville/mize) package (which is
on CRAN but not in a sufficiently advanced state currently).
* Optionally, [RSpectra](https://cran.r-project.org/package=RSpectra) package, 
which is used only if you want to initialize from a spectral method (`Y_init = "laplacian"` 
or `Y_init = "normlaplacian"`). If not present, then the standard R function 
`eigen` is used, but it's a *lot* slower.

## Explorations of `smallvis`

I have written some longer form documents that demonstrate the output of 
`smallvis` and its various options, with an emphasis on t-SNE.

* The [datasets](https://jlmelville.github.io/smallvis/datasets.html) used in my ruminations.
* Some [theory](https://jlmelville.github.io/smallvis/theory.html) showing a 
comparison of cost functions and gradients.
* Some material on various 
[spectral](https://jlmelville.github.io/smallvis/spectral.html)
methods, which justifies the use of Laplacian Eigenmaps for initialization (a bit).
* A look at the effect of the [scaling](https://jlmelville.github.io/smallvis/scale.html) and 
[PCA and whitening](https://jlmelville.github.io/smallvis/pcaw.html) preprocessing options on t-SNE.
* How [perplexity](https://jlmelville.github.io/smallvis/perplexity.html) affects t-SNE.
* Using [intrinsic dimensionality](https://jlmelville.github.io/smallvis/idp-theory.html)
to [automatically select perplexity](https://jlmelville.github.io/smallvis/idp.html). Also, 
[non-global perplexities](https://jlmelville.github.io/smallvis/idp-class.html).
* [Distance-preserving methods](https://jlmelville.github.io/smallvis/mmds.html) (geodesic as well as Euclidean).
* A comparison of various flavours of [Stochastic Neighbor Embedding](https://jlmelville.github.io/smallvis/sne.html) (t-distributed and otherwise).
* A comparison of [NeRV and JSE](https://jlmelville.github.io/smallvis/nervjse.html).
* [Heavy-Tailed SSNE (HSSNE)](https://jlmelville.github.io/smallvis/hssne.html).
* [Optimizing the Heavy-Tail parameter in HSSNE](https://jlmelville.github.io/smallvis/dhssne.html).
* Testing some different [initialization methods](https://jlmelville.github.io/smallvis/init.html)
for t-SNE. Related: does using an [average distance matrix](https://jlmelville.github.io/smallvis/averaging.html) 
based on the result of multiple random initializations do better than just picking the result with the lowest error? 
(It doesn't.)
* Some initial results using [UMAP](https://jlmelville.github.io/smallvis/umap.html)
and related methods. Also, looking at the effect of borrowing some features from UMAP
to see [the effect on t-SNE](https://jlmelville.github.io/smallvis/umaptsne.html) (sadly, not much).
* [t-SNE normalization](https://jlmelville.github.io/smallvis/norm.html).
* An evaluation of [late exaggeration](https://jlmelville.github.io/smallvis/lexagg.html) 
in t-SNE, as suggested in the [FIt-SNE](https://arxiv.org/abs/1712.09005) paper.
* [Elastic Embedding](https://jlmelville.github.io/smallvis/ee.html) and a 
[t-Distributed variant](https://jlmelville.github.io/smallvis/tee.html).
* [Optimizaton: L-BFGS](https://jlmelville.github.io/smallvis/opt.html).
* [Optimizaton: Spectral Direction](https://jlmelville.github.io/smallvis/specd.html).
* [Optimizaton: Conjugate Gradient](https://jlmelville.github.io/smallvis/cg.html).
* [Optimizaton: SGD methods](https://jlmelville.github.io/smallvis/sgd.html).
* [Perplexity and three clusters](https://jlmelville.github.io/smallvis/three-clusters.html)
  from from [How to Use t-SNE Effectively](https://distill.pub/2016/misread-tsne).
* How hard is it to unroll the [Swiss Roll with SNE](https://jlmelville.github.io/smallvis/swisssne.html)?
* Using input bandwidths in [SNE output functions](https://jlmelville.github.io/smallvis/bandwidths.html).

