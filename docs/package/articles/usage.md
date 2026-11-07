# Using smallvis

`smallvis` lets you experiment with neighborhood embedding methods,
including t-SNE and the LargeVis and UMAP cost functions. You can change
the input scaling, affinities, initialization, and optimizer
independently. Most methods use quadratic time and memory and suit small
datasets; [Barnes-Hut
t-SNE](https://jlmelville.github.io/smallvis/package/articles/barnes-hut.md)
provides an approximate route for larger inputs.

## Install

The package lives in the `smallvis` subdirectory of the repository.
Install it and its dependencies, including the GitHub plotting package
[vizier](https://github.com/jlmelville/vizier), with:

``` r

install.packages("pak")
pak::pak("jlmelville/smallvis/smallvis")
```

Installation from source requires an R-compatible C++ toolchain, such as
Rtools on Windows.

## Make an embedding

[`smallvis()`](https://jlmelville.github.io/smallvis/package/reference/smallvis.md)
selects the numeric columns of a data frame, so the species column in
`iris` stays out of the calculation. Here t-SNE starts from scaled PCA
coordinates and returns a two-column matrix. Plotting during
optimization is disabled so that the example also works in scripts.

``` r

library(smallvis)

embedding <- smallvis(
  iris,
  Y_init = "spca",
  perplexity = 30,
  epoch_callback = FALSE,
  n_threads = 1,
  verbose = FALSE
)
dim(embedding)
```

    ## [1] 150   2

``` r

species <- as.integer(iris$Species)
plot(
  embedding,
  col = c("#0072B2", "#D55E00", "#009E73")[species],
  pch = species,
  xlab = "t-SNE 1",
  ylab = "t-SNE 2"
)
```

![Two-dimensional t-SNE embedding of iris, with species distinguished by
color and plotting symbol.](usage_files/figure-html/plot-1.png)

The coordinates describe a neighborhood embedding: nearby points are the
relationships the method tries to preserve. The [perplexity
experiments](https://jlmelville.github.io/smallvis/package/articles/perplexity.md)
show how changing the neighborhood size changes the picture.

## Inspect the result

Set `ret_extra = TRUE` to return a list containing the coordinates in
`Y`, final objective contributions in `costs`, and an optimization
history in `itercosts`.

``` r

fit <- smallvis(
  iris,
  Y_init = "spca",
  perplexity = 30,
  ret_extra = TRUE,
  epoch_callback = FALSE,
  n_threads = 1,
  verbose = FALSE
)
c(iterations = fit$iter, final_cost = sum(fit$costs))
```

    ##   iterations   final_cost 
    ## 1000.0000000    0.1231099

The final cost is evaluated at the returned coordinates after removing
exaggeration. Values in the history were recorded during optimization
and can include exaggeration. The [function
reference](https://jlmelville.github.io/smallvis/package/reference/smallvis.html#value)
describes the available diagnostics, including optional matrices that
can be expensive to calculate or retain.

## Plot during optimization

The default `epoch_callback = TRUE` uses `vizier` to plot progress. You
can also supply your own callback. This example colors points by species
and puts the iteration number in the title:

``` r

iris_plot <- function(Y, iter, cost = NULL) {
  plot(Y, col = as.integer(iris$Species), main = paste("Iteration", iter))
}
smallvis(iris, Y_init = "spca", epoch_callback = iris_plot)
```

## Where next?

- [Inputs and
  initialization](https://jlmelville.github.io/smallvis/package/articles/inputs.md):
  scaling, PCA, distance inputs, and starting coordinates.
- [Optimization and repeated
  runs](https://jlmelville.github.io/smallvis/package/articles/optimization.md):
  methods, optimizers, exaggeration, and selecting among random starts.
- [Barnes-Hut and native
  acceleration](https://jlmelville.github.io/smallvis/package/articles/barnes-hut.md):
  approximate t-SNE and thread controls.
- [Theory and
  experiments](https://jlmelville.github.io/smallvis/package/articles/index.md):
  explanations and comparisons behind the package.
- [About
  smallvis](https://jlmelville.github.io/smallvis/package/articles/background.md):
  motivation and earlier development notes.
