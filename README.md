# smallvis

An R package for small-scale dimensionality reduction using neighbor embedding 
methods.

A package to compare the cost functions used in [t-Distributed Stochastic Neighbor Embedding](https://lvdmaaten.github.io/tsne/) and 
[LargeVis](https://arxiv.org/abs/1602.00370). The latter method seems to work
very similarly to t-SNE despite not normalizing the output weights. The
[official implementation](https://github.com/lferry007/LargeVis) is focused 
heavily on providing O(N) scaling, and as a result, it's not that convenient for
studying why it works. 

This package is designed to make it easier to compare
the output of embedding using the t-SNE and LargeVis cost functions, but at the
cost of them being O(N^2) in storage and computation costs (and being in pure R).
Unlike LargeVis (and the 
[Barnes-Hut implementation of t-SNE](https://github.com/lvdmaaten/bhtsne)) it is
therefore *not* suitable for large scale visualization. Hence the name smallvis.

## Prerequisites

By default, `smallvis` uses the [vizier](https://github.com/jlmelville/vizier)
package to plot the coordinates during optimization. It's not on CRAN, and 
therefore requires a fairly new version of 
[devtools](https://cran.r-project.org/package=devtools) (1.9 or greater) to 
install this as a dependency from github.

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

# By default, we use all numeric columns found in a data frame, so you don't need to filter out factor or strings
# set verbose = TRUE to log progress to the console
tsne_iris <- smallvis(iris, perplexity = 25, epoch_callback = iris_plot, verbose = TRUE)

# Default method is t-SNE, use largevis cost function instead
# Also needs a gamma value specified, and not as easy to optimize as t-SNE:
# reduce learning rate (eta) and increase maximum iterations
largevis_iris <- smallvis(iris, method = "largevis", gamma = 1, perplexity = 25, epoch_callback = iris_plot, 
                          eta = 10, max_iter = 5000, verbose = TRUE)

# use (scaled) PCA initialization so embedding is repeatable
tsne_iris_spca <- smallvis(iris, perplexity = 25, epoch_callback = iris_plot, Y_init = "spca")

# return extra information in a list, like with Rtsne
tsne_iris_extra <- smallvis(iris, perplexity = 25, epoch_callback = iris_plot, ret_extra = TRUE)

# more (potentially large and expensive to calculate) return values, but you have to ask for them specifically
tsne_iris_extra_extra <- smallvis(iris, perplexity = 25, epoch_callback = iris_plot,
                              ret_extra = c("P", "Q", "DX", "DY", "X"))

# Repeat embedding 10 times and keep the one with the best cost
tsne_iris_best <- smallvis_rep(nrep = 10, iris, perplexity = 25, ret_extra = TRUE)
plot(tsne_iris_best$Y)
```

## Things To Be Aware Of

* LargeVis requires the use of a `gamma` parameter, which weights the contribution
of attractive and repulsive contributions to the cost function. In the real LargeVis,
it is recommended to set this value to `7`, but this relies on the specifics of
the stochastic gradient descent method. In `smallvis`, this value is very 
dataset dependent: the more data you have, the smaller gamma should be to avoid
over-emphasising repulsive interactions.
* As mentioned above, the real LargeVis uses a classic stochastic gradient 
descent approach with a decaying learning rate. The implementation in `smallvis`
uses the same delta-bar-delta optimization method used in t-SNE. It works well
in my experience, but may require some tuning and seems to require more 
iterations compared to optimizing t-SNE.
* The LargeVis gradient requires quite a large value of epsilon to
avoid division by zero and get results. It's hard-coded to 0.1 in the LargeVis
source code, so I have used the same value by default in `smallvis`.
* Mainly for my own benefit, there is also a 
[theory](https://jlmelville.github.io/smallvis/theory.html) page 
showing a comparison of the cost function and gradient for t-SNE and LargeVis.

## See Also

For large scale visualization in R see:

* The Barnes Hut t-SNE package [Rtsne](https://cran.r-project.org/package=Rtsne)
* The [largeVis](https://cran.r-project.org/package=largeVis) package.

Much of the code here is based on my [fork](https://github.com/jlmelville/rtsne) 
of Justin Donaldson's R package for [t-SNE](https://cran.r-project.org/package=tsne).

## License

[GPLv2 or later](https://www.gnu.org/licenses/gpl-2.0.txt), although any 
LargeVis-specific code (e.g. cost and gradient calculation) can also be 
considered [Apache 2.0](https://www.apache.org/licenses/LICENSE-2.0).
