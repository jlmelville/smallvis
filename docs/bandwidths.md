---
title: "Bandwidths"
date: "January 27, 2018"
output:
  html_document:
    theme: cosmo
    toc: true
    toc_float:
      collapsed: false
---

Up: [Documentation Home](https://jlmelville.github.io/smallvis/).

The investigation of the 
[three clusters](https://jlmelville.github.io/smallvis/three-clusters.html) 
dataset found that using ASNE with the input bandwidths used in the output
kernel solved all the problems with strange distortions occuring at higher
perplexities. Unfortunately, there was reason to believe that this would
not generalize to other datasets. Let's see if that is the case.

Also, the original description of 
[NeRV](http://www.jmlr.org/papers/v11/venna10a.html) transferred bandwidths
directly from the input kernel to the output kernel. So there is a literature
precedent for this. However, later publications about NeRV don't mention this
step, so it may have turned out to be a bad idea. We shall see.

## Datasets

See the [Datasets](https://jlmelville.github.io/smallvis/datasets.html) page.

## Evaluation

Apart from visualizing the results, the mean neighbor preservation of the
40 closest neighbors is used to provide a rough quantification of the quality
of the result, labelled as `mnp@40` in the plots.

## Settings

Adding bandwidths massively changes the magnitude of the gradients compared
to non-bandwidth methods. This required extensive learning rate twiddling,
both on a per-method *and* per-dataset basis.

```
# eta = 1000 for basne
# eta = 20000 for btasne
mnist_bnerv <- smallvis(s1k, scale = FALSE, perplexity = 40, Y_init = "spca", method = "bnerv", eta = 100, max_iter = 2000, epoch = 100)
```

## Results

Given the lack of anything interesting happening with SSNE, t-SNE and their 
bandwidthed equivalents on the 
[Swiss Roll](https://jlmelville.github.io/smallvis/swisssne.html) data, to save
some effort, we'll only look at ASNE and t-ASNE, and their bandwidthed versions
below.

### iris

|                             |                           |
:----------------------------:|:--------------------------:
![iris ASNE](../img/sne/iris_asne.png)|![iris BASNE](../img/bandwidths/iris_basne.png)
![iris t-ASNE](../img/norm/iris_tasne.png)|![iris Bt-ASNE](../img/bandwidths/iris_btasne.png)
![iris NeRV (0.5)](../img/nerv-jse/iris_nerv0_5.png)|![iris BNeRV](../img/bandwidths/iris_bnerv.png)

### s1k

|                             |                           |
:----------------------------:|:--------------------------:
![s1k ASNE](../img/sne/s1k_asne.png)|![s1k BASNE](../img/bandwidths/s1k_basne.png)
![s1k t-ASNE](../img/norm/s1k_tasne.png)|![s1k Bt-ASNE](../img/bandwidths/s1k_btasne.png)
![s1k NeRV (0.5)](../img/nerv-jse/s1k_nerv0_5.png)|![s1k BNeRV](../img/bandwidths/s1k_bnerv.png)

### oli

|                             |                           |
:----------------------------:|:--------------------------:
![oli ASNE](../img/sne/oli_asne.png)|![oli BASNE](../img/bandwidths/oli_basne.png)
![oli t-ASNE](../img/norm/oli_tasne.png)|![oli Bt-ASNE](../img/bandwidths/oli_btasne.png)
![oli NeRV (0.5)](../img/nerv-jse/oli_nerv0_5.png)|![oli BNeRV](../img/bandwidths/oli_bnerv.png)

### frey

|                             |                           |
:----------------------------:|:--------------------------:
![frey ASNE](../img/sne/frey_asne.png)|![frey BASNE](../img/bandwidths/frey_basne.png)
![frey t-ASNE](../img/norm/frey_tasne.png)|![frey Bt-ASNE](../img/bandwidths/frey_btasne.png)
![frey NeRV (0.5)](../img/nerv-jse/frey_nerv0_5.png)|![frey BNeRV](../img/bandwidths/frey_bnerv.png)

### coil20

|                             |                           |
:----------------------------:|:--------------------------:
![coil20 ASNE](../img/sne/coil20_asne.png)|![coil20 BASNE](../img/bandwidths/coil20_basne.png)
![coil20 t-ASNE](../img/norm/coil20_tasne.png)|![coil20 Bt-ASNE](../img/bandwidths/coil20_btasne.png)
![coil20 NeRV (0.5)](../img/nerv-jse/coil20_nerv0_5.png)|![coil20 BNeRV](../img/bandwidths/coil20_bnerv.png)

### mnist

|                             |                           |
:----------------------------:|:--------------------------:
![mnist ASNE](../img/sne/mnist_asne.png)|![mnist BASNE](../img/bandwidths/mnist_basne.png)
![mnist t-ASNE](../img/norm/mnist_tasne.png)|![mnist Bt-ASNE](../img/bandwidths/mnist_btasne.png)
![mnist NeRV (0.5)](../img/nerv-jse/mnist_nerv0_5.png)|![mnist BNeRV](../img/bandwidths/mnist_bnerv.png)

### fashion

|                             |                           |
:----------------------------:|:--------------------------:
![fashion ASNE](../img/sne/fashion_asne.png)|![fashion BASNE](../img/bandwidths/fashion_basne.png)
![fashion t-ASNE](../img/norm/fashion_tasne.png)|![fashion Bt-ASNE](../img/bandwidths/fashion_btasne.png)
![fashion NeRV (0.5)](../img/nerv-jse/fashion_nerv0_5.png)|![fashion BNeRV](../img/bandwidths/fashion_bnerv.png)

## Conclusions

In terms of neighborhood preservation, there are a few datasets (`iris`, `oli`, 
`frey` and `coil20`) where adding bandwidths to ASNE improves the results by a 
tiny amount. But it's not a major improvement: t-ASNE always outperforms BASNE
except on the consistently anomalous `iris` dataset. So it's not hugely to
discover that t-ASNE and NeRV are improved by adding bandwidths to the output
kernel only for `iris`.

Adding bandwidths directly from the input kernel to the output kernel has
nothing to recommend it, so it's not surprising that the NeRV method stopped
doing it.

Up: [Documentation Home](https://jlmelville.github.io/smallvis/).
