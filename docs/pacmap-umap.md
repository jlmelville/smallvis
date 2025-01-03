---
title: "PaCMAP Compared to UMAP"
date: "January 21 2022"
output:
  html_document:
    theme: cosmo
    toc: true
    toc_float:
      collapsed: false
editor_options: 
  chunk_output_type: inline
---

Part 2 of many example images of PaCMAP. A continuation of [Notes on
PaCMAP](https://jlmelville.github.io/smallvis/pacmap.html) and [PaCMAP
Examples](https://jlmelville.github.io/smallvis/pacmap-examples.html).

Here I have tried to make as direct a comparison with UMAP as I can. I again
fiddled with the PaCMAP code to turn off some things: non-PCA data is centered
but not range-scaled and I changed the PCA initialization in PaCMAP to scale the
standard deviation to 1 for all columns. Back in [PaCMAP
Examples](https://jlmelville.github.io/smallvis/pacmap-examples.html) I
speculated that an initialization like this would avoid the poor results I saw
when I turned off mid-pair weights.

For UMAP, I used [uwot](https://cran.r-project.org/package=uwot) which I have
easy control over to mess about with for my needs. I used the exact nearest
neighbors for the UMAP case because I have pregenerated them for these datasets
which removes a large source of computational time. For initialization, I used
the coordinates from PaCMAP at iteration 0.

In the table below, the first row is PaCMAP. The first column uses mid-near
weights as usual, the second doesn't use them *and* sets the near-pairs to 1 for
all iterations. This is to test my theory that scaling the PCA initialization
should be good enough for PaCMAP. It also removes the influence of the mid-near
pairs vs UMAP. Additionally, if it works well enough without mid-near
interactions, that would suggest you could implement something PaCMAP-like
fairly easily in existing UMAP-style code bases (like `uwot`).

The second row is UMAP results. The left column is UMAP with default settings.
The right column is t-UMAP, which should be "gentler" in terms of its forces
(although not as gentle as PaCMAP).

The third row is UMAP results but I use the 100D PCA input like PaCMAP would (if
the input dimensionality > 100), and also the near pairs are the scaled nearest
neighbors. Previous results suggest that this won't have a noticeable effect
except on `macosko2015` and `ng20`. Again UMAP is on the left and t-UMAP is
on the right.

As we don't have enough images as it is, I have added two other datasets:
`coil100`, which is like `coil20` but has more images; and `isofaces`, a more
manifold-like dataset, beloved of many a spectral-based embedding paper.

Finally, sorry for not having useful titles on the images.

## UMAP, t-UMAP and PaCMAP

| iris |  | 
|:--|:--|
![iris-pacmap-15-spca1-it450](../img/pacmap/nomid/iris-pacmap-15-spca1-it450.png)|![iris-pacmap-15-spca1nomid1-it450](../img/pacmap/nomid/iris-pacmap-15-spca1nomid1-it450.png)|
![iris-umap-spca1](../img/pacmap/nomid/iris-umap-spca1.png)|![iris-tumap-spca1](../img/pacmap/nomid/iris-tumap-spca1.png)|
![iris-umap-spca1s](../img/pacmap/nomid/iris-umap-spca1s.png)|![iris-tumap-spca1s](../img/pacmap/nomid/iris-tumap-spca1s.png)|

| s1k |  | 
|:--|:--|
![s1k-pacmap-15-spca1-it450](../img/pacmap/nomid/s1k-pacmap-15-spca1-it450.png)|![s1k-pacmap-15-spca1nomid1-it450](../img/pacmap/nomid/s1k-pacmap-15-spca1nomid1-it450.png)|
![s1k-umap-spca1](../img/pacmap/nomid/s1k-umap-spca1.png)|![s1k-tumap-spca1](../img/pacmap/nomid/s1k-tumap-spca1.png)|
![s1k-umap-spca1s](../img/pacmap/nomid/s1k-umap-spca1s.png)|![s1k-tumap-spca1s](../img/pacmap/nomid/s1k-tumap-spca1s.png)|

| oli |  | 
|:--|:--|
![oli-pacmap-15-spca1-it450](../img/pacmap/nomid/oli-pacmap-15-spca1-it450.png)|![oli-pacmap-15-spca1nomid1-it450](../img/pacmap/nomid/oli-pacmap-15-spca1nomid1-it450.png)|
![oli-umap-spca1](../img/pacmap/nomid/oli-umap-spca1.png)|![oli-tumap-spca1](../img/pacmap/nomid/oli-tumap-spca1.png)|
![oli-umap-spca1s](../img/pacmap/nomid/oli-umap-spca1s.png)|![oli-tumap-spca1s](../img/pacmap/nomid/oli-tumap-spca1s.png)|

| frey | | 
|:--|:--|
![frey-pacmap-15-spca1-it450](../img/pacmap/nomid/frey-pacmap-15-spca1-it450.png)|![frey-pacmap-15-spca1nomid1-it450](../img/pacmap/nomid/frey-pacmap-15-spca1nomid1-it450.png)|
![frey-umap-spca1](../img/pacmap/nomid/frey-umap-spca1.png)|![frey-tumap-spca1](../img/pacmap/nomid/frey-tumap-spca1.png)|
![frey-umap-spca1s](../img/pacmap/nomid/frey-umap-spca1s.png)|![frey-tumap-spca1s](../img/pacmap/nomid/frey-tumap-spca1s.png)|

| coil20 |  | 
|:--|:--|
![coil20-pacmap-15-spca1-it450](../img/pacmap/nomid/coil20-pacmap-15-spca1-it450.png)|![coil20-pacmap-15-spca1nomid1-it450](../img/pacmap/nomid/coil20-pacmap-15-spca1nomid1-it450.png)|
![coil20-umap-spca1](../img/pacmap/nomid/coil20-umap-spca1.png)|![coil20-tumap-spca1](../img/pacmap/nomid/coil20-tumap-spca1.png)|
![coil20-umap-spca1s](../img/pacmap/nomid/coil20-umap-spca1s.png)|![coil20-tumap-spca1s](../img/pacmap/nomid/coil20-tumap-spca1s.png)|

This the one dataset where the `dcorr` is noticeably worse for the PaCMAP result
without the mid-pair forces vs with it.

| coil100 |  | 
|:--|:--|
![coil100-pacmap-15-spca1-it450](../img/pacmap/nomid/coil100-pacmap-15-spca1-it450.png)|![coil100-pacmap-15-spca1nomid1-it450](../img/pacmap/nomid/coil100-pacmap-15-spca1nomid1-it450.png)|
![coil100-umap-spca1](../img/pacmap/nomid/coil100-umap-spca1.png)|![coil100-tumap-spca1](../img/pacmap/nomid/coil100-tumap-spca1.png)|
![coil100-umap-spca1s](../img/pacmap/nomid/coil100-umap-spca1s.png)|![coil100-tumap-spca1s](../img/pacmap/nomid/coil100-tumap-spca1s.png)|

| mnist |  | 
|:--|:--|
![mnist-pacmap-15-spca1-it450](../img/pacmap/nomid/mnist-pacmap-15-spca1-it450.png)|![mnist-pacmap-15-spca1nomid1-it450](../img/pacmap/nomid/mnist-pacmap-15-spca1nomid1-it450.png)|
![mnist-umap-spca1](../img/pacmap/nomid/mnist-umap-spca1.png)|![mnist-tumap-spca1](../img/pacmap/nomid/mnist-tumap-spca1.png)|
![mnist-umap-spca1s](../img/pacmap/nomid/mnist-umap-spca1s.png)|![mnist-tumap-spca1s](../img/pacmap/nomid/mnist-tumap-spca1s.png)|

| fashion | | 
|:--|:--|
![fashion-pacmap-15-spca1-it450](../img/pacmap/nomid/fashion-pacmap-15-spca1-it450.png)|![fashion-pacmap-15-spca1nomid1-it450](../img/pacmap/nomid/fashion-pacmap-15-spca1nomid1-it450.png)|
![fashion-umap-spca1](../img/pacmap/nomid/fashion-umap-spca1.png)|![fashion-tumap-spca1](../img/pacmap/nomid/fashion-tumap-spca1.png)|
![fashion-umap-spca1s](../img/pacmap/nomid/fashion-umap-spca1s.png)|![fashion-tumap-spca1s](../img/pacmap/nomid/fashion-tumap-spca1s.png)|

| kuzushiji |  | 
|:--|:--|
![kuzushiji-pacmap-15-spca1-it450](../img/pacmap/nomid/kuzushiji-pacmap-15-spca1-it450.png)|![kuzushiji-pacmap-15-spca1nomid1-it450](../img/pacmap/nomid/kuzushiji-pacmap-15-spca1nomid1-it450.png)|
![kuzushiji-umap-spca1](../img/pacmap/nomid/kuzushiji-umap-spca1.png)|![kuzushiji-tumap-spca1](../img/pacmap/nomid/kuzushiji-tumap-spca1.png)|
![kuzushiji-umap-spca1s](../img/pacmap/nomid/kuzushiji-umap-spca1s.png)|![kuzushiji-tumap-spca1s](../img/pacmap/nomid/kuzushiji-tumap-spca1s.png)|

| cifar10 |  | 
|:--|:--|
![cifar10-pacmap-15-spca1-it450](../img/pacmap/nomid/cifar10-pacmap-15-spca1-it450.png)|![cifar10-pacmap-15-spca1nomid1-it450](../img/pacmap/nomid/cifar10-pacmap-15-spca1nomid1-it450.png)|
![cifar10-umap-spca1](../img/pacmap/nomid/cifar10-umap-spca1.png)|![cifar10-tumap-spca1](../img/pacmap/nomid/cifar10-tumap-spca1.png)|
![cifar10-umap-spca1s](../img/pacmap/nomid/cifar10-umap-spca1s.png)|![cifar10-tumap-spca1s](../img/pacmap/nomid/cifar10-tumap-spca1s.png)|

| norb |  | 
|:--|:--|
![norb-pacmap-15-spca1-it450](../img/pacmap/nomid/norb-pacmap-15-spca1-it450.png)|![norb-pacmap-15-spca1nomid1-it450](../img/pacmap/nomid/norb-pacmap-15-spca1nomid1-it450.png)|
![norb-umap-spca1](../img/pacmap/nomid/norb-umap-spca1.png)|![norb-tumap-spca1](../img/pacmap/nomid/norb-tumap-spca1.png)|
![norb-umap-spca1s](../img/pacmap/nomid/norb-umap-spca1s.png)|![norb-tumap-spca1s](../img/pacmap/nomid/norb-tumap-spca1s.png)|

| ng20 |  | 
|:--|:--|
![ng20-pacmap-15-spca1-it450](../img/pacmap/nomid/ng20-pacmap-15-spca1-it450.png)|![ng20-pacmap-15-spca1nomid1-it450](../img/pacmap/nomid/ng20-pacmap-15-spca1nomid1-it450.png)|
![ng20-umap-spca1](../img/pacmap/nomid/ng20-umap-spca1.png)|![ng20-tumap-spca1](../img/pacmap/nomid/ng20-tumap-spca1.png)|
![ng20-umap-spca1s](../img/pacmap/nomid/ng20-umap-spca1s.png)|![ng20-tumap-spca1s](../img/pacmap/nomid/ng20-tumap-spca1s.png)|

This dataset shows the biggest difference between PaCMAP and UMAP according to
both global and local metrics. This might really be saying that it's a bad idea
to use the Euclidean metric on this dataset no matter what you use. Hard to
see that any of the visualizations are incredibly insightful at a global level.

| mammoth |  | 
|:--|:--|
![mammoth-pacmap-15-spca1-it450](../img/pacmap/nomid/mammoth-pacmap-15-spca1-it450.png)|![mammoth-pacmap-15-spca1nomid1-it450](../img/pacmap/nomid/mammoth-pacmap-15-spca1nomid1-it450.png)|
![mammoth-umap-spca1](../img/pacmap/nomid/mammoth-umap-spca1.png)|![mammoth-tumap-spca1](../img/pacmap/nomid/mammoth-tumap-spca1.png)|
![mammoth-umap-spca1s](../img/pacmap/nomid/mammoth-umap-spca1s.png)|![mammoth-tumap-spca1s](../img/pacmap/nomid/mammoth-tumap-spca1s.png)|

| swissroll |  | 
|:--|:--|
![swissroll-pacmap-15-spca1-it450](../img/pacmap/nomid/swissroll-pacmap-15-spca1-it450.png)|![swissroll-pacmap-15-spca1nomid1-it450](../img/pacmap/nomid/swissroll-pacmap-15-spca1nomid1-it450.png)|
![swissroll-umap-spca1](../img/pacmap/nomid/swissroll-umap-spca1.png)|![swissroll-tumap-spca1](../img/pacmap/nomid/swissroll-tumap-spca1.png)|
![swissroll-umap-spca1s](../img/pacmap/nomid/swissroll-umap-spca1s.png)|![swissroll-tumap-spca1s](../img/pacmap/nomid/swissroll-tumap-spca1s.png)|

| s-curve with a hole | | 
|:--|:--|
![scurvehole-pacmap-15-spca1-it450](../img/pacmap/nomid/scurvehole-pacmap-15-spca1-it450.png)|![scurvehole-pacmap-15-spca1nomid1-it450](../img/pacmap/nomid/scurvehole-pacmap-15-spca1nomid1-it450.png)|
![scurvehole-umap-spca1](../img/pacmap/nomid/scurvehole-umap-spca1.png)|![scurvehole-tumap-spca1](../img/pacmap/nomid/scurvehole-tumap-spca1.png)|
![scurvehole-umap-spca1s](../img/pacmap/nomid/scurvehole-umap-spca1s.png)|![scurvehole-tumap-spca1s](../img/pacmap/nomid/scurvehole-tumap-spca1s.png)|

| isofaces |  | 
|:--|:--|
![isofaces-pacmap-15-spca1-it450](../img/pacmap/nomid/isofaces-pacmap-15-spca1-it450.png)|![isofaces-pacmap-15-spca1nomid1-it450](../img/pacmap/nomid/isofaces-pacmap-15-spca1nomid1-it450.png)|
![isofaces-umap-spca1](../img/pacmap/nomid/isofaces-umap-spca1.png)|![isofaces-tumap-spca1](../img/pacmap/nomid/isofaces-tumap-spca1.png)|
![isofaces-umap-spca1s](../img/pacmap/nomid/isofaces-umap-spca1s.png)|![isofaces-tumap-spca1s](../img/pacmap/nomid/isofaces-tumap-spca1s.png)|

PaCMAP reliably splits the isofaces dataset into two from a PCA initialization.
The `dcorr` and `rtp` metrics indicate that the UMAP results are to be preferred
but they are not a feast for the eyes that's for sure. The `emd` metric
marginally prefers the PaCMAP results though.

| macosko2015 | | 
|:--|:--|
![macosko2015-pacmap-15-spca1-it450](../img/pacmap/nomid/macosko2015-pacmap-15-spca1-it450.png)|![macosko2015-pacmap-15-spca1nomid1-it450](../img/pacmap/nomid/macosko2015-pacmap-15-spca1nomid1-it450.png)|
![macosko2015-umap-spca1](../img/pacmap/nomid/macosko2015-umap-spca1.png)|![macosko2015-tumap-spca1](../img/pacmap/nomid/macosko2015-tumap-spca1.png)|
![macosko2015-umap-spca1s](../img/pacmap/nomid/macosko2015-umap-spca1s.png)|![macosko2015-tumap-spca1s](../img/pacmap/nomid/macosko2015-tumap-spca1s.png)|

This the only dataset where the `rtp` metric suggests that UMAP and t-UMAP 
results on PCA-processed data do a noticeably better job of preserving the
global structure than not applying PCA.

| tasic2018 |  | 
|:--|:--|
![tasic2018-pacmap-15-spca1-it450](../img/pacmap/nomid/tasic2018-pacmap-15-spca1-it450.png)|![tasic2018-pacmap-15-spca1nomid1-it450](../img/pacmap/nomid/tasic2018-pacmap-15-spca1nomid1-it450.png)|
![tasic2018-umap-spca1](../img/pacmap/nomid/tasic2018-umap-spca1.png)|![tasic2018-tumap-spca1](../img/pacmap/nomid/tasic2018-tumap-spca1.png)|
![tasic2018-umap-spca1s](../img/pacmap/nomid/tasic2018-umap-spca1s.png)|![tasic2018-tumap-spca1s](../img/pacmap/nomid/tasic2018-tumap-spca1s.png)|

| lamanno2020 |  | 
|:--|:--|
![lamanno2020-pacmap-15-spca1-it450](../img/pacmap/nomid/lamanno2020-pacmap-15-spca1-it450.png)|![lamanno2020-pacmap-15-spca1nomid1-it450](../img/pacmap/nomid/lamanno2020-pacmap-15-spca1nomid1-it450.png)|
![lamanno2020-umap-spca1](../img/pacmap/nomid/lamanno2020-umap-spca1.png)|![lamanno2020-tumap-spca1](../img/pacmap/nomid/lamanno2020-tumap-spca1.png)
![lamanno2020-umap-spca1s](../img/pacmap/nomid/lamanno2020-umap-spca1s.png)|![lamanno2020-tumap-spca1s](../img/pacmap/nomid/lamanno2020-tumap-spca1s.png)

Overall, the global structure according the to `rtp` metric is pretty much the
same for all the methods. According to `dcorr`, UMAP might do a better job.
For PaCMAP, where there is a difference when not including the mid-pair forces,
once again the evidence suggests that even with a PCA initialization, having
the mid-pair forces helps.

For local preservation, UMAP seems better at retaining the 15 nearest neighbors,
but the difference disappears at 65 nearest neighbors. This would be in line
with the idea that UMAP's attractive forces are stronger than those of PaCMAP
and the consistently more torn and distorted manifolds for the
S-curve-with-a-hole, swiss roll and mammoth also show that. t-UMAP definitely
seems to tear manifolds less strongly than UMAP (but definitely more than
PaCMAP).

The effect of using the local scaling + PCA for finding nearest
neighbors seems slight for all methods, whether we look at global or local
measures of preservation. I couldn't see any particular visual trend either,
apart from the expected difference for `macosko2015` and `ng20` and it's the
effect of applying PCA that has the big effect there.

## Effect of Edge Weights

[Böhm and co-workers](https://arxiv.org/abs/2007.08902) note that replacing
UMAP's edge weights with uniform values didn't have a big effect on results
(see their figure 2 with the MNIST digits for example). But when I added
the PaCMAP near pair and far pair forces (no mid near forces) to `uwot`, I
got results that on a lot of occasions looked more like t-UMAP than the Python PaCMAP results, even
though I had gone through the various efforts outlined above to remove several
differing aspects between PaCMAP and UMAP.

Assuming it's not a bug (in which case this entire section will be stricken
from the record), I noticed a shift to a more PaCMAP-like behavior when I set
all the edge weights to 1. I also tested putting the t-UMAP forces into
PaCMAP and still got results that looked PaCMAP-like to me. So below is a more
in-depth investigation on each of the datasets with various ways of shoe-horning
t-UMAP into PaCMAP or vice versa.

For each set of four images, the methods are as follows.

Top left: UMAP with the attractive and repulsive interactions replaced with the
PaCMAP versions for near and far pairs respectively. The sampling is based on
the UMAP graph edges.

Top right: UMAP with the attractive and repulsive interactions replaced with the
PaCMAP versions for near and far pairs respectively. The sampling is based on
the PaCMAP graph edges (i.e. the same weight for all edges).

Center left: t-UMAP with sampling based on the PaCMAP graph edges (i.e. the same
weight for all edges).

Center right: PaCMAP Python code with the mid-pair interactions turned off and the
near and far pair interactions replaced with the t-UMAP attractive and repulsive
interactions respectively. Sampling uses PaCMAP graph edges (i.e. the same
weight for all edges).

Bottom left: PaCMAP with mid-pair forces turned off (repeated from the previous
section). This is what PaCMAP "should" look like.

Bottom right: t-UMAP with 100D PCA and local scaling (also repeated from above).
This is what t-UMAP "should" look like.

Lots of plots below with some commentary where I think there is an example that
illustrates the effect.

| iris |  |             
|:--|:--|
![iris-upacmap-spca1](../img/pacmap/tpacmap/iris-upacmap-spca1.png)|![iris-upacmapv1-spca1](../img/pacmap/tpacmap/iris-upacmapv1-spca1.png)|
![iris-tumapv1-spca1](../img/pacmap/tpacmap/iris-tumapv1-spca1.png)|![iris-pacmap-15-tpacmap-it450](../img/pacmap/tpacmap/iris-pacmap-15-tpacmap-it450.png)|
![iris-pacmap-15-spca1nomid-it450](../img/pacmap/nomid/iris-pacmap-15-spca1nomid-it450.png)|![iris-tumap-spca1s](../img/pacmap/nomid/iris-tumap-spca1s.png)|

| s1k |  | 
|:--|:--|
![s1k-upacmap-spca1](../img/pacmap/tpacmap/s1k-upacmap-spca1.png)|![s1k-upacmapv1-spca1](../img/pacmap/tpacmap/s1k-upacmapv1-spca1.png)|
![s1k-tumapv1-spca1](../img/pacmap/tpacmap/s1k-tumapv1-spca1.png)|![s1k-pacmap-15-tpacmap-it450](../img/pacmap/tpacmap/s1k-pacmap-15-tpacmap-it450.png)|
![s1k-pacmap-15-spca1nomid-it450](../img/pacmap/nomid/s1k-pacmap-15-spca1nomid-it450.png)|![s1k-tumap-spca1s](../img/pacmap/nomid/s1k-tumap-spca1s.png)|

Alright so in this case, the PaCMAP in `uwot` have the desired effect without
having to binarize the edge weights. Like adding the t-UMAP forces to `pacmap`
creates a very t-UMAP-like result.

| oli |  | 
|:--|:--|
![oli-upacmap-spca1](../img/pacmap/tpacmap/oli-upacmap-spca1.png)|![oli-upacmapv1-spca1](../img/pacmap/tpacmap/oli-upacmapv1-spca1.png)|
![oli-tumapv1-spca1](../img/pacmap/tpacmap/oli-tumapv1-spca1.png)|![oli-pacmap-15-tpacmap-it450](../img/pacmap/tpacmap/oli-pacmap-15-tpacmap-it450.png)|
![oli-pacmap-15-spca1nomid-it450](../img/pacmap/nomid/oli-pacmap-15-spca1nomid-it450.png)|![oli-tumap-spca1s](../img/pacmap/nomid/oli-tumap-spca1s.png)|

Certainly for `oli` just putting PaCMAP near and far forces into `uwot` (top 
left) does not produce a `PaCMAP` effect, but binarizing the edges (top right)
does then give something close to the original PaCMAP results (bottom left).

| frey | | 
|:--|:--|
![frey-upacmap-spca1](../img/pacmap/tpacmap/frey-upacmap-spca1.png)|![frey-upacmapv1-spca1](../img/pacmap/tpacmap/frey-upacmapv1-spca1.png)|
![frey-tumapv1-spca1](../img/pacmap/tpacmap/frey-tumapv1-spca1.png)|![frey-pacmap-15-tpacmap-it450](../img/pacmap/tpacmap/frey-pacmap-15-tpacmap-it450.png)|
![frey-pacmap-15-spca1nomid-it450](../img/pacmap/nomid/frey-pacmap-15-spca1nomid-it450.png)|![frey-tumap-spca1s](../img/pacmap/nomid/frey-tumap-spca1s.png)|

As usual, `frey` does its own thing. I wouldn't say this can tell us much
either way. 

| coil20 |  | 
|:--|:--|
![coil20-upacmap-spca1](../img/pacmap/tpacmap/coil20-upacmap-spca1.png)|![coil20-upacmapv1-spca1](../img/pacmap/tpacmap/coil20-upacmapv1-spca1.png)|
![coil20-tumapv1-spca1](../img/pacmap/tpacmap/coil20-tumapv1-spca1.png)|![coil20-pacmap-15-tpacmap-it450](../img/pacmap/tpacmap/coil20-pacmap-15-tpacmap-it450.png)|
![coil20-pacmap-15-spca1nomid-it450](../img/pacmap/nomid/coil20-pacmap-15-spca1nomid-it450.png)|![coil20-tumap-spca1s](../img/pacmap/nomid/coil20-tumap-spca1s.png)|

The PaCMAP-in-uwot result is rather odd: you can see one of the loops becomes
a straight line. This isn't a quirk of the random number generator, it happens
every time I run it. But looking at the bottom row you can see that t-UMAP
causes the loops to be less compact and more evenly spaced (bottom right)
versus PaCMAP (bottom left). The t-UMAP-in-PaCMAP result (center right) 
looks closer to the PaCMAP result than the t-UMAP result. 

| coil100 |  | 
|:--|:--|
![coil100-upacmap-spca1](../img/pacmap/tpacmap/coil100-upacmap-spca1.png)|![coil100-upacmapv1-spca1](../img/pacmap/tpacmap/coil100-upacmapv1-spca1.png)|
![coil100-tumapv1-spca1](../img/pacmap/tpacmap/coil100-tumapv1-spca1.png)|![coil100-pacmap-15-tpacmap-it450](../img/pacmap/tpacmap/coil100-pacmap-15-tpacmap-it450.png)|
![coil100-pacmap-15-spca1nomid-it450](../img/pacmap/nomid/coil100-pacmap-15-spca1nomid-it450.png)|![coil100-tumap-spca1s](../img/pacmap/nomid/coil100-tumap-spca1s.png)|

Ok, here it's pretty obvious, right? PaCMAP tends to create three clusters of
really bunched up loops with the others highly compressed and studded across the plot.
t-UMAP is much more even in its distribution. The PaCMAP-in-uwot version (top left)
looks nothing like the Python PaCMAP version (bottom left) until binarized edges
are used (top right). Meanwhile, the t-UMAP-in-pacmap results (center right) also
resembles the PaCMAP version.

| mnist |  | 
|:--|:--|
![mnist-upacmap-spca1](../img/pacmap/tpacmap/mnist-upacmap-spca1.png)|![mnist-upacmapv1-spca1](../img/pacmap/tpacmap/mnist-upacmapv1-spca1.png)|
![mnist-tumapv1-spca1](../img/pacmap/tpacmap/mnist-tumapv1-spca1.png)|![mnist-pacmap-15-tpacmap-it450](../img/pacmap/tpacmap/mnist-pacmap-15-tpacmap-it450.png)|
![mnist-pacmap-15-spca1nomid-it450](../img/pacmap/nomid/mnist-pacmap-15-spca1nomid-it450.png)|![mnist-tumap-spca1s](../img/pacmap/nomid/mnist-tumap-spca1s.png)|

| fashion | | 
|:--|:--|
![fashion-upacmap-spca1](../img/pacmap/tpacmap/fashion-upacmap-spca1.png)|![fashion-upacmapv1-spca1](../img/pacmap/tpacmap/fashion-upacmapv1-spca1.png)|
![fashion-tumapv1-spca1](../img/pacmap/tpacmap/fashion-tumapv1-spca1.png)|![fashion-pacmap-15-tpacmap-it450](../img/pacmap/tpacmap/fashion-pacmap-15-tpacmap-it450.png)|
![fashion-pacmap-15-spca1nomid-it450](../img/pacmap/nomid/fashion-pacmap-15-spca1nomid-it450.png)|![fashion-tumap-spca1s](../img/pacmap/nomid/fashion-tumap-spca1s.png)|

| kuzushiji |  | 
|:--|:--|
![kuzushiji-upacmap-spca1](../img/pacmap/tpacmap/kuzushiji-upacmap-spca1.png)|![kuzushiji-upacmapv1-spca1](../img/pacmap/tpacmap/kuzushiji-upacmapv1-spca1.png)|
![kuzushiji-tumapv1-spca1](../img/pacmap/tpacmap/kuzushiji-tumapv1-spca1.png)|![kuzushiji-pacmap-15-tpacmap-it450](../img/pacmap/tpacmap/kuzushiji-pacmap-15-tpacmap-it450.png)|
![kuzushiji-pacmap-15-spca1nomid-it450](../img/pacmap/nomid/kuzushiji-pacmap-15-spca1nomid-it450.png)|![kuzushiji-tumap-spca1s](../img/pacmap/nomid/kuzushiji-tumap-spca1s.png)|

| cifar10 |  | 
|:--|:--|
![cifar10-upacmap-spca1](../img/pacmap/tpacmap/cifar10-upacmap-spca1.png)|![cifar10-upacmapv1-spca1](../img/pacmap/tpacmap/cifar10-upacmapv1-spca1.png)|
![cifar10-tumapv1-spca1](../img/pacmap/tpacmap/cifar10-tumapv1-spca1.png)|![cifar10-pacmap-15-tpacmap-it450](../img/pacmap/tpacmap/cifar10-pacmap-15-tpacmap-it450.png)|
![cifar10-pacmap-15-spca1nomid-it450](../img/pacmap/nomid/cifar10-pacmap-15-spca1nomid-it450.png)|![cifar10-tumap-spca1s](../img/pacmap/nomid/cifar10-tumap-spca1s.png)|

| norb |  | 
|:--|:--|
![norb-upacmap-spca1](../img/pacmap/tpacmap/norb-upacmap-spca1.png)|![norb-upacmapv1-spca1](../img/pacmap/tpacmap/norb-upacmapv1-spca1.png)|
![norb-tumapv1-spca1](../img/pacmap/tpacmap/norb-tumapv1-spca1.png)|![norb-pacmap-15-tpacmap-it450](../img/pacmap/tpacmap/norb-pacmap-15-tpacmap-it450.png)|
![norb-pacmap-15-spca1nomid-it450](../img/pacmap/nomid/norb-pacmap-15-spca1nomid-it450.png)|![norb-tumap-spca1s](../img/pacmap/nomid/norb-tumap-spca1s.png)|

Another example where the PaCMAP result (bottom left) and t-UMAP results (bottom
right) are quite distinguishable, with the blue structures and especially the green
bits being a lot more compact for PaCMAP. The switch for PaCMAP-in-uwot (top
left) going to binarized edges (top right) is pretty clear.

| ng20 |  | 
|:--|:--|
![ng20-upacmap-spca1](../img/pacmap/tpacmap/ng20-upacmap-spca1.png)|![ng20-upacmapv1-spca1](../img/pacmap/tpacmap/ng20-upacmapv1-spca1.png)|
![ng20-tumapv1-spca1](../img/pacmap/tpacmap/ng20-tumapv1-spca1.png)|![ng20-pacmap-15-tpacmap-it450](../img/pacmap/tpacmap/ng20-pacmap-15-tpacmap-it450.png)|
![ng20-pacmap-15-spca1nomid-it450](../img/pacmap/nomid/ng20-pacmap-15-spca1nomid-it450.png)|![ng20-tumap-spca1s](../img/pacmap/nomid/ng20-tumap-spca1s.png)|

In the case of `ng20`, it seems the change in force types does most of the work,
not the effect of edge weights.

| mammoth |  | 
|:--|:--|
![mammoth-upacmap-spca1](../img/pacmap/tpacmap/mammoth-upacmap-spca1.png)|![mammoth-upacmapv1-spca1](../img/pacmap/tpacmap/mammoth-upacmapv1-spca1.png)|
![mammoth-tumapv1-spca1](../img/pacmap/tpacmap/mammoth-tumapv1-spca1.png)|![mammoth-pacmap-15-tpacmap-it450](../img/pacmap/tpacmap/mammoth-pacmap-15-tpacmap-it450.png)|
![mammoth-pacmap-15-spca1nomid-it450](../img/pacmap/nomid/mammoth-pacmap-15-spca1nomid-it450.png)|![mammoth-tumap-spca1s](../img/pacmap/nomid/mammoth-tumap-spca1s.png)|

`mammoth` proving its worth again as a way to interrogate these algorithms. The
effect is pretty stark here. t-UMAP tends to curve the legs and tusks while
PaCMAP keeps them straight. And yet it does seem a substantial effect comes
from binarizing edges: PaCMAP-in-uwot (top left) shows lots of curving until
edges are binarized (top right). The tusks still curve but the effect is
reduced. Meanwhile, t-UMAP-in-pacmap (center right) hsa straight lines. Using
binarized edges with t-UMAP (center left) also looks much "cleaner".

| swissroll |  | 
|:--|:--|
![swissroll-upacmap-spca1](../img/pacmap/tpacmap/swissroll-upacmap-spca1.png)|![swissroll-upacmapv1-spca1](../img/pacmap/tpacmap/swissroll-upacmapv1-spca1.png)|
![swissroll-tumapv1-spca1](../img/pacmap/tpacmap/swissroll-tumapv1-spca1.png)|![swissroll-pacmap-15-tpacmap-it450](../img/pacmap/tpacmap/swissroll-pacmap-15-tpacmap-it450.png)|
![swissroll-pacmap-15-spca1nomid-it450](../img/pacmap/nomid/swissroll-pacmap-15-spca1nomid-it450.png)|![swissroll-tumap-spca1s](../img/pacmap/nomid/swissroll-tumap-spca1s.png)|

None of the results here give as good results as plain PaCMAP (bottom right),
but adding binarized edges seems to always improve t-UMAP. Running out of differences
between UMAP and PaCMAP here, so the remaining difference in results is either due to the approximate nearest 
neighbors (unlikely as in 3 dimensions Annoy should do a very good job), or
because of the effect of UMAP always choosing different random neighbors.

| s-curve with a hole | | 
|:--|:--|
![scurvehole-upacmap-spca1](../img/pacmap/tpacmap/scurvehole-upacmap-spca1.png)|![scurvehole-upacmapv1-spca1](../img/pacmap/tpacmap/scurvehole-upacmapv1-spca1.png)|
![scurvehole-tumapv1-spca1](../img/pacmap/tpacmap/scurvehole-tumapv1-spca1.png)|![scurvehole-pacmap-15-tpacmap-it450](../img/pacmap/tpacmap/scurvehole-pacmap-15-tpacmap-it450.png)|
![scurvehole-pacmap-15-spca1nomid-it450](../img/pacmap/nomid/scurvehole-pacmap-15-spca1nomid-it450.png)|![scurvehole-tumap-spca1s](../img/pacmap/nomid/scurvehole-tumap-spca1s.png)|

The top left and bottom right results do stand out from the others in terms
of smoothness of the manifold.

| isofaces |  | 
|:--|:--|
![isofaces-upacmap-spca1](../img/pacmap/tpacmap/isofaces-upacmap-spca1.png)|![isofaces-upacmapv1-spca1](../img/pacmap/tpacmap/isofaces-upacmapv1-spca1.png)|
![isofaces-tumapv1-spca1](../img/pacmap/tpacmap/isofaces-tumapv1-spca1.png)|![isofaces-pacmap-15-tpacmap-it450](../img/pacmap/tpacmap/isofaces-pacmap-15-tpacmap-it450.png)|
![isofaces-pacmap-15-spca1nomid-it450](../img/pacmap/nomid/isofaces-pacmap-15-spca1nomid-it450.png)|![isofaces-tumap-spca1s](../img/pacmap/nomid/isofaces-tumap-spca1s.png)|

| macosko2015 | | 
|:--|:--|
![macosko2015-upacmap-spca1](../img/pacmap/tpacmap/macosko2015-upacmap-spca1.png)|![macosko2015-upacmapv1-spca1](../img/pacmap/tpacmap/macosko2015-upacmapv1-spca1.png)|
![macosko2015-tumapv1-spca1](../img/pacmap/tpacmap/macosko2015-tumapv1-spca1.png)|![macosko2015-pacmap-15-tpacmap-it450](../img/pacmap/tpacmap/macosko2015-pacmap-15-tpacmap-it450.png)|
![macosko2015-pacmap-15-spca1nomid-it450](../img/pacmap/nomid/macosko2015-pacmap-15-spca1nomid-it450.png)|![macosko2015-tumap-spca1s](../img/pacmap/nomid/macosko2015-tumap-spca1s.png)|

| tasic2018 |  | 
|:--|:--|
![tasic2018-upacmap-spca1](../img/pacmap/tpacmap/tasic2018-upacmap-spca1.png)|![tasic2018-upacmapv1-spca1](../img/pacmap/tpacmap/tasic2018-upacmapv1-spca1.png)|
![tasic2018-tumapv1-spca1](../img/pacmap/tpacmap/tasic2018-tumapv1-spca1.png)|![tasic2018-pacmap-15-tpacmap-it450](../img/pacmap/tpacmap/tasic2018-pacmap-15-tpacmap-it450.png)|
![tasic2018-pacmap-15-spca1nomid-it450](../img/pacmap/nomid/tasic2018-pacmap-15-spca1nomid-it450.png)|![tasic2018-tumap-spca1s](../img/pacmap/nomid/tasic2018-tumap-spca1s.png)|

| lamanno2020 |  | 
|:--|:--|
![lamanno2020-upacmap-spca1](../img/pacmap/tpacmap/lamanno2020-upacmap-spca1.png)|![lamanno2020-upacmapv1-spca1](../img/pacmap/tpacmap/lamanno2020-upacmapv1-spca1.png)
![lamanno2020-tumapv1-spca1](../img/pacmap/tpacmap/lamanno2020-tumapv1-spca1.png)|![lamanno2020-pacmap-15-tpacmap-it450](../img/pacmap/tpacmap/lamanno2020-pacmap-15-tpacmap-it450.png)
![lamanno2020-pacmap-15-spca1nomid-it450](../img/pacmap/nomid/lamanno2020-pacmap-15-spca1nomid-it450.png)|![lamanno2020-tumap-spca1s](../img/pacmap/nomid/lamanno2020-tumap-spca1s.png)|

The effect of binarizing edges is not always obvious or consistent across
datasets, but there are enough examples here to make me think that binarizing
the edges has quite a bit to do with PaCMAP's behavior on low-dimensional 
manifolds like the `mammoth` and even more structured data like `coil100`.
Although binarizing edge weights has been treated as an approximation that has no
effect on UMAP, I would caution against that, although it may be something that
only affects a small number of datasets in the same way preprocessing the data
with PCA may not have a big effect.

Without considering clusters or other more supervised methods of global
structure preservation, the global preservation metrics indicate that if you
have access to a sensibly-scaled PCA initialization for UMAP I am not convinced
that PaCMAP is a huge win over UMAP. However, first, the PaCMAP paper doesn't
say that, it is focused on being able to get good results from more random
initializations. Second, the Python UMAP implementation *doesn't* currently have
such a PCA implementation option. But `uwot` *does*, and that's what I use so this
is less of a concern for me. 

PaCMAP definitely seems to do a better job or ripping manifolds less (although
it certainly can't work miracles with `swiss roll` and will detach a leg on the
`mammoth`) so it would have a big advantage there. Unfortunately, apart from the
more synthetic simulation sets used here, it's not clear to me that for the
other datasets, if there is some hidden manifold structure within the
cluster-like blobs, that PaCMAP is any better at revealing that than UMAP.

Finally, I still have some qualms about whether the global preservation metrics
are helpful. In particular, I am suspicious of the `emd` metric which regularly
gives the opposite conclusion to `rtp` and `dcorr`. This may be down to the
stochastic sampling nature of my implementation, but I think there are deeper
issues about the normalization that is required. In particular, global metrics
could give entirely the opposite conclusion if there *is* the manifold structure
that PaCMAP is better at handling: successfully unrolling the `swiss roll` will
lead to very poor correlations with the original input distances. The spectral
initialization that UMAP does by default is usually much more successful with
these kinds of structures visually (`mammoth` is not a huge success, but
`isofaces` and `swiss roll` do well) but you wouldn't know it if you just looked
at the global preservation metrics (as usual the local preservation is
unchanged). Also, I think these metrics may reward large separations of roughly
clustered data. Arguably very best `norb` results are from when the mid-pair
weights are turned off and the initial PCA result caused large distances. The
different parts of the dataset merrily rearranged locally without shifting from
their locations due to PCA so both global and local results were pretty good.
But I don't know if that really represents the best possible visualization that
PaCMAP or UMAP could come up with for that dataset. However, the idea of
optimizing each neighborhood in isolation in a series of mini-UMAPs, then trying
to find a good alignment for them in the style of LTSA might be an interesting
way to leverage this observation.
