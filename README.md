# smallvis

An R package for experimenting with dimensionality reduction using neighborhood embedding methods,
including t-SNE, LargeVis, and UMAP cost functions. It gives you control over input scaling, nearest
neighbors, initialization, and optimization so you can compare methods and try out ideas.

Most methods use quadratic time and memory and suit small datasets. Barnes-Hut supports larger
datasets for t-SNE, HSSNE, t-distributed elastic embedding, LargeVis, UMAP, t-UMAP, and Jensen-Shannon
SNE by combining sparse neighbor affinities with multithreaded C++ code.

## Install

The R package lives in the `smallvis` subdirectory. Install it and its dependencies, including the
GitHub package `vizier`, with

```r
install.packages("pak")
pak::pak("jlmelville/smallvis/smallvis")
```

Installation from source requires an R-compatible C++ toolchain, such as Rtools on Windows.

## Quick start

`smallvis()` uses the numeric columns of a data frame. This example initializes t-SNE with scaled
PCA coordinates and plots the final embedding.

```r
library(smallvis)

embedding <- smallvis(
  iris,
  Y_init = "spca",
  perplexity = 30,
  epoch_callback = FALSE,
  n_threads = 1,
  verbose = FALSE
)
plot(embedding, col = as.integer(iris$Species), pch = 19)
```

Set `ret_extra = TRUE` to return a list containing the embedding in `Y`, costs, and other diagnostics.
Use `smallvis_rep()` to select the lowest-cost embedding from repeated random starts, or
`smallvis_perpstep()` to optimize through a series of decreasing perplexities.

## Documentation

- [Package reference](https://jlmelville.github.io/smallvis/package/reference/index.html)
- [Getting started](https://jlmelville.github.io/smallvis/package/articles/usage.html)
- [Guides, theory, and experiments](https://jlmelville.github.io/smallvis/package/articles/index.html),
  including inputs, initialization, perplexity, methods, and optimization

## License

GPL version 2 or later. LargeVis-specific code is also available under
[Apache 2.0](https://www.apache.org/licenses/LICENSE-2.0), and Barnes-Hut and UMAP-related code under
the BSD 3-clause license. See the [code origins and copyright notices](smallvis/inst/COPYRIGHTS).
