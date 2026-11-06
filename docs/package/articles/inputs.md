# Inputs and initialization

Scaling, neighborhood construction, and initialization can change an
embedding as much as the choice of cost function. `smallvis` exposes
these choices so that you can hold them fixed while comparing methods.
This guide covers the current interface; the linked experiments explore
their effects.

## Choose the input representation

Supply a numeric matrix, a data frame, or an R `dist` object. Rows are
observations and columns are features; nonnumeric columns of a data
frame are filtered out. Use a `dist` object for precomputed distances:
an ordinary square matrix is still treated as observations and features.

``` r

set.seed(42)
distance_fit <- smallvis(
  dist(iris[, 1:4]),
  scale = FALSE,
  Y_init = "rand",
  perplexity = 30,
  epoch_callback = FALSE,
  n_threads = 1,
  verbose = FALSE
)
dim(distance_fit)
```

    ## [1] 150   2

Feature scaling and PCA need the original features, so choose those
before calculating distances. PCA initialization also needs features;
random or suitable spectral initialization can be used with distance
inputs.

## Scale and reduce the features

The default `scale = "absmax"` applies one global rescaling. Use
`scale = TRUE` to center and scale each feature, or `scale = FALSE` to
retain its original units. The [input
scaling](https://jlmelville.github.io/smallvis/package/articles/scale.md)
experiment illustrates why these choices matter.

Input PCA is separate from PCA initialization. With `pca = TRUE`,
`initial_dims` controls how many input components are retained before
affinities are calculated. Here the four iris features are reduced to
three:

``` r

pca_fit <- smallvis(
  iris,
  scale = TRUE,
  pca = TRUE,
  initial_dims = 3,
  Y_init = "spca",
  perplexity = 30,
  ret_extra = "X",
  epoch_callback = FALSE,
  n_threads = 1,
  verbose = FALSE
)
dim(pca_fit$X)
```

    ## [1] 150   3

Use `pca = "whiten"` to additionally give the retained components unit
variance. Whitening changes the relative importance of components, so it
changes the distances used by the embedding. See the [PCA and whitening
experiments](https://jlmelville.github.io/smallvis/package/articles/pcaw.md)
before treating it as a routine speed optimization.

## Set the neighborhood size

For Gaussian affinities, `perplexity` controls the effective
neighborhood size and must lie between 1 and the number of observations
minus 1. Try a few values when exploring a dataset; changing perplexity
changes the affinities being fitted, so final costs at different
perplexities are not a common score for selecting the best picture.

For UMAP’s smoothed-neighbor affinities, the same argument controls
neighborhood size rather than Gaussian perplexity. The [function
reference](https://jlmelville.github.io/smallvis/package/reference/smallvis.html#arguments)
gives its precise meaning for each input choice.

See
[perplexity](https://jlmelville.github.io/smallvis/package/articles/perplexity.md),
[intrinsic
dimensionality](https://jlmelville.github.io/smallvis/package/articles/idp.md),
and [multiscale
perplexities](https://jlmelville.github.io/smallvis/package/articles/msperp.md)
for experiments with alternatives to a single global setting. The
supported forms are documented in the [intrinsic
dimensionality](https://jlmelville.github.io/smallvis/package/reference/smallvis.html#intrinsic-dimensionality-perplexity)
and
[multiscale](https://jlmelville.github.io/smallvis/package/reference/smallvis.html#multiscale-perplexities)
reference sections.

## Choose the starting coordinates

The default `Y_init = "rand"` draws a small random Gaussian layout. Set
the random seed immediately before the call when you want to repeat a
run in the same environment. `Y_init = "spca"` starts from scaled PCA
scores, avoiding the variation from random starting coordinates. These
are different uses of PCA: `pca` changes the input distances; `Y_init`
changes where optimization begins.

`Y_init = "laplacian"` and `Y_init = "normlaplacian"` use spectral
coordinates from the input graph. They need a suitable affinity matrix.
The optional `RSpectra` package accelerates these initializations;
without it, smallvis uses base R’s full eigendecomposition. You can also
pass an explicit matrix with one row per observation and one column per
output dimension to `Y_init`.

The [initialization
experiments](https://jlmelville.github.io/smallvis/package/articles/init.md)
compare starting layouts, and [Spectral
Methods](https://jlmelville.github.io/smallvis/package/articles/spectral.md)
explains the different graph constructions and eigenvectors. For
repeated random starts, continue to [optimization and repeated
runs](https://jlmelville.github.io/smallvis/package/articles/optimization.md).
