---
title: "SNE"
date: "January 2, 2018"
output:
  html_document:
    theme: cosmo
    toc: true
    toc_float:
      collapsed: false
---

## Datasets

See the [Datasets](https://jlmelville.github.io/smallvis/datasets.html) page.

## Evaluation

Apart from visualizing the results, the mean neighbor preservation of the
40 closest neighbors is used to provide a rough quantification of the quality
of the result, labelled as `mnp@40` in the plots.

## Settings

Here's an example of generating the results using the `fashion` dataset. For
ASNE and SSNE, values of `eta` were chosen that showed decent behavior on the
`fashion` dataset: convergence within 1000 iterations without a hint of divergence,
e.g. the error increasing during an "epoch". These settings were applied to
the other datasets without any further fiddling, so to be on the safe side,
I allowed for double the number of iterations (`max_iter = 2000`) than for
t-SNE.

```
fashion_wtsne <- smallvis(fashion, scale = FALSE, perplexity = 40, Y_init = "spca", 
method = "wtsne", ret_extra = c("dx", "dy"))
fashion_asne <- smallvis(fashion, scale = FALSE, perplexity = 40, Y_init = "spca", method = "asne", ret_extra = c("dx", "dy"), eta = 0.1, max_iter = 2000)
fashion_ssne <- smallvis(fashion, scale = FALSE, perplexity = 40, Y_init = "spca", method = "ssne", ret_extra = c("dx", "dy"), eta = 10, max_iter = 2000)
```

## Results

The DBD results are on the left, the Adam results are on the right.

### iris

|                             |                           |
:----------------------------:|:--------------------------:
![iris t-SNE](../img/opt/iris_dbd.png)|![iris wt-SNE](../img/sne/iris_wtsne.png)
![iris ASNE](../img/sne/iris_asne.png)|![iris SSNE](../img/sne/iris_ssne.png)

### s1k

|                             |                           |
:----------------------------:|:--------------------------:
![s1k t-SNE](../img/opt/s1k_dbd.png)|![s1k wt-SNE](../img/sne/s1k_wtsne.png)
![s1k ASNE](../img/sne/s1k_asne.png)|![s1k SSNE](../img/sne/s1k_ssne.png)

### oli

|                             |                           |
:----------------------------:|:--------------------------:
![oli t-SNE](../img/opt/oli_dbd.png)|![oli wt-SNE](../img/sne/oli_wtsne.png)
![oli ASNE](../img/sne/oli_asne.png)|![oli SSNE](../img/sne/oli_ssne.png)

### frey

|                             |                           |
:----------------------------:|:--------------------------:
![frey t-SNE](../img/opt/frey_dbd.png)|![frey wt-SNE](../img/sne/frey_wtsne.png)
![frey ASNE](../img/sne/frey_asne.png)|![frey SSNE](../img/sne/frey_ssne.png)

### coil20

|                             |                           |
:----------------------------:|:--------------------------:
![coil20 t-SNE](../img/opt/coil20_dbd.png)|![coil20 wt-SNE](../img/sne/coil20_wtsne.png)
![coil20 ASNE](../img/sne/coil20_asne.png)|![coil20 SSNE](../img/sne/coil20_ssne.png)

### mnist

|                             |                           |
:----------------------------:|:--------------------------:
![mnist t-SNE](../img/opt/mnist_dbd.png)|![mnist wt-SNE](../img/sne/mnist_wtsne.png)
![mnist ASNE](../img/sne/mnist_asne.png)|![mnist SSNE](../img/sne/mnist_ssne.png)

### fashion

|                             |                           |
:----------------------------:|:--------------------------:
![fashion t-SNE](../img/opt/fashion_dbd.png)|![fashion wt-SNE](../img/sne/fashion_wtsne.png)
![fashion ASNE](../img/sne/fashion_asne.png)|![fashion SSNE](../img/sne/fashion_ssne.png)

## Conclusions

ASNE and SSNE results converged quicker than t-SNE.
