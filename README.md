# smallvis

An R package for small-scale dimensionality reduction using 
neighborhood-preservation
dimensionality reduction methods, including [t-Distributed Stochastic Neighbor Embedding](https://lvdmaaten.github.io/tsne/), 
[LargeVis](https://arxiv.org/abs/1602.00370) and 
[UMAP](https://github.com/lmcinnes/umap). 

LargeVis and UMAP are of particular interest because they seem to give 
visualizations which are very competitive with t-SNE, while being more amenable
to scaling up to large datasets. Why do these methods work so well? This package
is designed to make it easier to experiment with and compare these methods.

One way it does this is by abandoning the more advanced nearest-neighbor 
methods, distance approximations, sampling, and multi-threaded stochastic 
gradient descent techniques. The price paid for this simplification is that the 
algorithms are back to being O(N^2) in storage and computation costs (and being 
in pure R). Unlike UMAP, 
[the official implementation of LargeVis](https://github.com/lferry007/LargeVis),
and the 
[Barnes-Hut implementation of t-SNE](https://github.com/lvdmaaten/bhtsne),
this package is therefore *not* suitable for large scale visualization.
Hence the name smallvis.

## Prerequisites

By default, `smallvis` uses the [vizier](https://github.com/jlmelville/vizier)
package to plot the coordinates during optimization. It's not on CRAN, and 
therefore requires a fairly new version of 
[devtools](https://cran.r-project.org/package=devtools) (1.9 or greater) to 
install this as a dependency from github.

There is also an optional dependency on the 
[RSpectra](https://cran.r-project.org/package=RSpectra) package, which is used
only if you want to initialize from a spectral method 
(set `Y_init = "laplacian"` to do this and see the paper by 
[Linderman and Steinerberger](https://arxiv.org/abs/1706.02582) for details on
why you might want to). If not present, then the standard R function `eigen`
is used, but this is much slower (because we only need the first few eigenvectors,
and `eigen` calculates all of them). On my Sandy Bridge-era laptop running
R 3.4.2 on Windows 10, using `Rspectra::eigs` to fetch the top 3 eigenvectors 
from a 6,000 x 6,000 affinity matrix takes about 6 seconds; using `eigen` takes 
around 25 minutes.

## Installing

```R
install.packages("devtools")
devtools::install_github("jlmelville/smallvis")
library(smallvis)
```

## Using

```R
# By default, we use all numeric columns found in a data frame, so you don't need to filter out factor or strings
# set verbose = TRUE to log progress to the console
# Automatically plots the results during optimization
tsne_iris <- smallvis(iris, perplexity = 25, verbose = TRUE)

# Using a custom epoch_callback
uniq_spec <- unique(iris$Species)
colors <- rainbow(length(uniq_spec))
names(colors) <- uniq_spec
iris_plot <- function(x) {
  plot(x, col = colors[iris$Species])
}

tsne_iris <- smallvis(iris, perplexity = 25, epoch_callback = iris_plot, verbose = TRUE)

# Default method is t-SNE, use largevis cost function instead
# Also needs a gamma value specified, and not as easy to optimize as t-SNE:
# reduce learning rate (eta) and increase maximum iterations
largevis_iris <- smallvis(iris, method = "largevis", gamma = 1, perplexity = 25, epoch_callback = iris_plot, 
                          eta = 0.1, max_iter = 5000, verbose = TRUE)

# UMAP: see https://github.com/lmcinnes/umap
umap_iris <- smallvis(iris, method = "umap", perplexity = 25, eta = 0.01)

# use (scaled) PCA initialization so embedding is repeatable
tsne_iris_spca <- smallvis(iris, perplexity = 25, epoch_callback = iris_plot, Y_init = "spca")

# or initialize from Laplacian Eigenmap (similar to UMAP initialization)
tsne_iris_lap <- smallvis(iris, perplexity = 25, epoch_callback = iris_plot, Y_init = "lap")

# return extra information in a list, like with Rtsne
tsne_iris_extra <- smallvis(iris, perplexity = 25, epoch_callback = iris_plot, ret_extra = TRUE)

# more (potentially large and expensive to calculate) return values, but you have to ask for them specifically
tsne_iris_extra_extra <- smallvis(iris, perplexity = 25, epoch_callback = iris_plot,
                              ret_extra = c("P", "Q", "DX", "DY", "X"))

# Repeat embedding 10 times and keep the one with the best cost
tsne_iris_best <- smallvis_rep(nrep = 10, iris, perplexity = 25, ret_extra = TRUE)
plot(tsne_iris_best$Y)

# Classical momentum optimization instead of delta-bar-delta
umap_iris_mom <- smallvis(iris, scale = FALSE, opt = list("mom", eta = 1e-2, mu = 0.8),
                          method = "umap", Y_init = "spca")

# L-BFGS optimization via the mize package
umap_iris_lbfgs <- smallvis(iris, scale = FALSE, opt = list("l-bfgs", c1 = 1e-4, c2 = 0.9),
                            method = "umap", Y_init = "spca", max_iter = 300)
```

## Things To Be Aware Of

* There isn't a UMAP publication yet, so the implementation in `smallvis` relies
on my examination of the UMAP source code, with some much-appreciated clarification from UMAP creator
[Leland McInnes](https://github.com/lmcinnes). Expect some bugs, and any horrifically
bogus results should be double-checked with the output of the 
[official UMAP implementation](https://github.com/lmcinnes/umap) before casting calumnies
on the quality of UMAP itself.
* LargeVis requires the use of a `gamma` parameter, which weights the contribution
of attractive and repulsive contributions to the cost function. In the real LargeVis,
it is recommended to set this value to `7`, but this relies on the specifics of
the stochastic gradient descent method. In `smallvis`, this value is very 
dataset dependent: the more data you have, the smaller gamma should be to avoid
over-emphasising repulsive interactions.
* LargeVis partitions each pairwise interaction into either an attractive or repulsive
contribution. In `smallvis`, each interaction is a combination of both.
* Both UMAP and LargeVis use a classic stochastic gradient descent approach 
with a decaying learning rate. The implementation in `smallvis`
uses the same delta-bar-delta optimization method used in t-SNE. It works well
in my experience, but may require some tuning and more iterations compared to 
optimizing t-SNE.
* In this setting, LargeVis and UMAP gradient requires quite a large value of 
epsilon to avoid division by zero and get decent results. It's hard-coded to 
0.1 in the LargeVis source code, so I have used the same value by default in 
`smallvis`. It can be controlled by the `lveps` parameter.
* Mainly for my own benefit, there is also a 
[theory](https://jlmelville.github.io/smallvis/theory.html) page 
showing a comparison of cost functions and gradients. Also, some material on
the various [spectral](https://jlmelville.github.io/smallvis/spectral.html)
methods, which justifies the use of Laplacian Eigenmaps (a bit).
* As an example of `smallvis` results, you can read about the 
[effect of initialization methods](https://jlmelville.github.io/smallvis/init.html)
on t-SNE.

## See Also

For large scale visualization in R see:

* The Barnes-Hut t-SNE package [Rtsne](https://cran.r-project.org/package=Rtsne)
* The [largeVis](https://cran.r-project.org/package=largeVis) package.
* RSpectra [https://cran.r-project.org/package=RSpectra] can help spectral methods
scale.

Also of relevance are:

* [UMAP](https://github.com/lmcinnes/umap) (in Python)
* [LargeVis](https://github.com/lferry007/LargeVis) (in C++)
* [Spectra](http://spectralib.org/), the C++ library that RSpectra wraps.

Much of the code here is based on my [fork](https://github.com/jlmelville/rtsne) 
of Justin Donaldson's R package for [t-SNE](https://cran.r-project.org/package=tsne).

## License

[GPLv2 or later](https://www.gnu.org/licenses/gpl-2.0.txt). Any 
LargeVis-specific code (e.g. cost and gradient calculation) can also be 
considered [Apache 2.0](https://www.apache.org/licenses/LICENSE-2.0).
Similary, UMAP-related code is also licensed as 
[BSD 3-clause](https://opensource.org/licenses/BSD-3-Clause).
