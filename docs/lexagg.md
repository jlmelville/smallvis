---
title: "Late Exaggeration"
date: "May 12, 2018"
output:
  html_document:
    theme: cosmo
    toc: true
    toc_float:
      collapsed: false
---

Up: [Documentation Home](https://jlmelville.github.io/smallvis/).

The [FIt-SNE](https://arxiv.org/abs/1712.09005) paper recommends the technique
of "late exaggeration". This is exactly the same as early exaggeration 
(multiply the input probabilities by a fixed constant), but doing it at the end
of the optimization, rather than the beginning. This seems to have the effect
of forcing any clusters to be tighter and further apart, which makes for better
visualizations. The effect seems reminiscent of the output of LargeVis and
UMAP, in fact.

In the paper, the authors set the exaggeration factor to 12 for the last 250
iterations of an embedding of the entire MNIST digits dataset, which mirrors the
recommended setting for early exaggeration with large datasets. 

## Datasets

See the [Datasets](https://jlmelville.github.io/smallvis/datasets.html) page.

## Settings

Settings will be the same as the 
[original t-SNE paper](http://www.jmlr.org/papers/v9/vandermaaten08a.html),
except with scaled PCA initialization rather than random initialization. Given
that the example in the FIt-SNE paper used late exaggeration settings that
matched the default early exaggeration settings, we shall do the same here,
except that the length and extent of the early exaggeration is lessened for the
sort of datasets that can be handled by `smallvis`: the `exaggeration_factor`
will be `4`, and will run for the last 100 iterations:

```
iris_late_ex <- smallvis(iris, perplexity = 40, Y_init = "spca", eta = 100,
                            exaggeration_factor = 4, stop_lying_iter = 100,
                            late_exaggeration_factor = 4, start_late_lying_iter = 900)
```

Initial results suggested that the exaggeration factor of 4 was a too high for
late exaggeration, so we will also look at reduced `late_exaggeration_factor` 
values of `2` and `1.5`. The early exaggeration factor is unchanged at 4 in
these cases. I looked at also reducing the early exaggeration factor when the
late exaggeration factor was reduced, but it doesn't have much of an effect, so
we needn't worry about that.

We'll also compare results with the standard t-SNE with early exaggeration only:

```
iris_tsne <- smallvis(iris, perplexity = 40, Y_init = "spca", eta = 100, 
                      exaggeration_factor = 4, stop_lying_iter = 100)
```

## Evaluation

Apart from visualizing the results, the mean neighbor preservation of the 40
closest neighbors (in line with the `perplexity`) is used to provide a rough 
quantification of the quality of the result. On the plots, this is labelled as 
`mnp@40`.

## Results

For each dataset, the top left image is the result of using late exaggeration
with an exaggeration factor of 4, the top right is with an exaggeration factor
of 2, and the bottom left is with an exaggeration factor of 1.5. The bottom
right is standard t-SNE with no late exaggeration. All results used early
exaggeration with an exaggeration factor of 4.

### iris

|                             |                           |
:----------------------------:|:--------------------------:
![iris lexagg 4](../img/lexagg/iris_lexagg4.png)|![iris lexagg 2](../img/lexagg/iris_lexagg2.png)
![iris lexagg 1.5](../img/lexagg/iris_lexagg1_5.png)|![iris exagg 4](../img/lexagg/iris_exagg4.png)

### s1k

|                             |                           |
:----------------------------:|:--------------------------:
![s1k lexagg 4](../img/lexagg/s1k_lexagg4.png)|![s1k lexagg 2](../img/lexagg/s1k_lexagg2.png)
![s1k lexagg 1.5](../img/lexagg/s1k_lexagg1_5.png)|![s1k exagg 4](../img/lexagg/s1k_exagg4.png)

### oli

|                             |                           |
:----------------------------:|:--------------------------:
![oli lexagg 4](../img/lexagg/oli_lexagg4.png)|![oli lexagg 2](../img/lexagg/oli_lexagg2.png)
![oli lexagg 1.5](../img/lexagg/oli_lexagg1_5.png)|![oli exagg 4](../img/lexagg/oli_exagg4.png)

### frey

|                             |                           |
:----------------------------:|:--------------------------:
![frey lexagg 4](../img/lexagg/frey_lexagg4.png)|![frey lexagg 2](../img/lexagg/frey_lexagg2.png)
![frey lexagg 1.5](../img/lexagg/frey_lexagg1_5.png)|![frey exagg 4](../img/lexagg/frey_exagg4.png)


### coil20

|                             |                           |
:----------------------------:|:--------------------------:
![coil20 lexagg 4](../img/lexagg/coil20_lexagg4.png)|![coil20 lexagg 2](../img/lexagg/coil20_lexagg2.png)
![coil20 lexagg 1.5](../img/lexagg/coil20_lexagg1_5.png)|![coil20 exagg 4](../img/lexagg/coil20_exagg4.png)

### mnist

|                             |                           |
:----------------------------:|:--------------------------:
![mnist lexagg 4](../img/lexagg/mnist_lexagg4.png)|![mnist lexagg 2](../img/lexagg/mnist_lexagg2.png)
![mnist lexagg 1.5](../img/lexagg/mnist_lexagg1_5.png)|![mnist exagg 4](../img/lexagg/mnist_exagg4.png)

### fashion

|                             |                           |
:----------------------------:|:--------------------------:
![fashion lexagg 4](../img/lexagg/fashion_lexagg4.png)|![fashion lexagg 2](../img/lexagg/fashion_lexagg2.png)
![fashion lexagg 1.5](../img/lexagg/fashion_lexagg1_5.png)|![fashion exagg 4](../img/lexagg/fashion_exagg4.png)


## Conclusions

A late exaggeration of 4 is clearly way too large for the datasets looked at
here. Either 1.5 or 2 seems much more appropriate. For these settings, we see
the expected result of much smaller, better-separated clusters with little
effect on the neighbor preservation values. The most suitable value seems
dataset-dependent (as usual): to my eyes, `iris` and `oli` don't benefit much
from any late exaggeration, although maybe an exaggeration even more mild than
1.5 would help.

For the other datasets, these results are very reminiscent of the
[UMAP](https://jlmelville.github.io/smallvis/umap.html) and 
[u-SNE](https://jlmelville.github.io/smallvis/umaptsne.html) results, but have
the advantage of being faster and easier to optimize.

Finally, if late exaggeration and early exaggeration are useful, why not always
be exaggerating?  Below are results from applying a constant exaggeration
with `exaggeration_factor = 1.5` by setting `stop_lying_iter = 1000`:

|                             |                           |
:----------------------------:|:--------------------------:
![iris cexagg 1.5](../img/lexagg/iris_cexagg1_5.png)|![s1k cexagg 1.5](../img/lexagg/s1k_cexagg1_5.png)
![oli cexagg 1.5](../img/lexagg/oli_cexagg1_5.png)|![frey cexagg 1.5](../img/lexagg/frey_cexagg1_5.png)
![coil20 cexagg 1.5](../img/lexagg/coil20_cexagg1_5.png)|![mnist cexagg 1.5](../img/lexagg/mnist_cexagg1_5.png)
![fashion cexagg 1.5](../img/lexagg/fashion_cexagg1_5.png)|

It seems that allowing a non-exaggerated phase of the optimization allows the 
clusters to re-arrange and spread out a bit, because the late-exaggeration
results are more spread out (particularly noticeable if you look at the range
on the axes), but these constant-exaggeration results look pretty good and 
still resemble a more-compressed version of the t-SNE results, although probably
starting from the usually good scaled PCA initialization may have something to 
do with that.

If you find that t-SNE results are not separating clusters in your data in a
good way, then it does seem that late exaggeration helps, but for the datasets
that `smallvis` can deal with, I recommend an `exaggeration_factor` of around
1.5, much lower than the value of 12 mentioned in the FIt-SNE paper, and even 4
is much too large based on the results shown here. 

Finally, even with the reduced exaggeration, for all examples I looked at,
the cost does increase over the last 100 iterations, suggesting that the 
learning rate or momentum values are too high for the late exaggeration phase.
Controlling this would introduce coupling between the optimizer and the 
exaggeration code, which I am not inclined to introduce, so it's just something
to be aware of.

## Acknowledgement

Thank you to Dmitry Kobak, who drew my attention to the similarity between
UMAP/LargeVis results and the output of late exaggeration, prompting me to
implement it in `smallvis`.

Up: [Documentation Home](https://jlmelville.github.io/smallvis/).

