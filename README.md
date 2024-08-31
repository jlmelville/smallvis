# smallvis

An R package for small-scale dimensionality reduction using 
neighborhood-preservation
dimensionality reduction methods, including 
[t-Distributed Stochastic Neighbor Embedding](https://lvdmaaten.github.io/tsne/), 
and a non-stochastic version of the [LargeVis](https://arxiv.org/abs/1602.00370)
and [UMAP](https://arxiv.org/abs/1802.03426) cost functions.

The purpose of this package is to make it easier to experiment with different
dimensionality reduction methods while having more control over things like
input scaling, nearest neighbor calculations, initialization and optimization, 
which can make comparisons between different packages difficult. Be warned, 
most implementations are not optimized for speed and scale like O(N^2), but see
below for my dream to upgrade `smallvis` to more of a mediumvis.

*August 11 2024* **The Turbo Championship Edition Update**. I have briefly
brought `smallvis` back from the dead to speed it up a bit. I have added:

* Barnes-Hut t-SNE. This will scale up to larger datasets and it is feasible
to run it on the full MNIST digits dataset (i.e. 70,000 items). This uses a
(2D only) C++ translation of the cython implementation in the Python
[openTSNE](https://github.com/pavlin-policar/openTSNE) package originally
authored by [Pavlin Poličar](https://github.com/pavlin-policar). It's BSD
3-clause licensed (and can be found in `src/bh.h`). Use it with 
`method = "bhtsne"`. The degree of approximation can be controlled with 
`theta`. Be aware that it's not as fast as e.g.
[Rtsne](https://github.com/jkrijthe/Rtsne), at least during the optimization
step. That package uses Laurens van der Maaten's original C++ code which I am
very unsure can be redistributed with R code due to its BSD 4-clause license. I
would love to be wrong about that though! The Quad Tree implementation could
be used with other embedding methods, but I haven't got round to implementing
that yet.
* A C++ multi-threaded perplexity search using only the nearest neighbors of
each point (3 times the perplexity) is used with BH t-SNE.
* My own [rnndescent](https://cran.r-project.org/package=rnndescent) package
replaces FNN for nearest neighbor search. Apart from being a monument to my ego,
it can be faster than FNN for brute force search because it can be
multi-threaded (use `n_threads` to control this). Also, approximate nearest
neighbor search becomes quite important with larger datasets.
* For exact search, I have started adding multi-threaded C++ code to calculate
the gradient. Set `use_cpp = TRUE` to use this. It's not as big a win in speed
up as you might hope because the R code is using some very optimized linear
algebra for some steps which will blow my puny C++ code out of the water.
However in many cases the linear algebra libraries won't be mult-threaded so
sheer brute force threads can overcome this. Just don't expect setting `n` 
threads to give you `n` times the speed. Like Barnes-Hut, this requires me to
implement the gradients in C++ for each method and I haven't done that yet.
* [irlba](https://cran.r-project.org/package=irlba) is now a dependency for 
doing PCA on larger datasets.

I will probably at least attempt to apply Barnes-Hut and multi-threading to some
other methods. On the other hand, it's taken me five years to get back to this,
so don't hold your breath.

## Prerequisites

`smallvis` uses the [vizier](https://github.com/jlmelville/vizier)
package to plot the coordinates during optimization. It's not on CRAN, and 
therefore requires a fairly new version of 
[devtools](https://cran.r-project.org/package=devtools) (1.9 or greater) to 
install this as a dependency from github.

There is also an optional dependency on the
[RSpectra](https://cran.r-project.org/package=RSpectra) package, which is used
only if you want to initialize from a spectral method (set `Y_init =
"laplacian"` or `Y_init = "normlaplacian"` to do this and see the paper by
[Linderman and Steinerberger](https://arxiv.org/abs/1706.02582) for details on
why you might want to). If not present, then the standard R function `eigen` is
used, but this is much slower (because we only need the first few eigenvectors,
and `eigen` calculates all of them). On my Sandy Bridge-era laptop running R
3.4.2 on Windows 10, using `Rspectra::eigs` to fetch the top 3 eigenvectors from
a 6,000 x 6,000 affinity matrix takes about 6 seconds; using `eigen` takes
around 25 minutes.

## Installing

```R
install.packages("devtools")
devtools::install_github("jlmelville/smallvis", subdir = "smallvis")
library(smallvis)
```

## Using

```R
# By default, we use all numeric columns found in a data frame, so you don't need to filter out factor or strings
# set verbose = TRUE to log progress to the console
# Automatically plots the results during optimization
tsne_iris <- smallvis(iris, perplexity = 25, verbose = TRUE)

# Barnes-Hut recommended settings:
bhtsne_iris <- smallvis(iris, bh = TRUE, n_threads = 4, perplexity = 30, 
                        nn = "approximate", inp_kernel = "knn", theta = 1.0,
                        exaggeration_factor = 12, stop_lying_iter = 250,
                        Y_init = "spca")

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
largevis_iris <- smallvis(iris, method = "largevis", perplexity = 25, epoch_callback = iris_plot, 
                          eta = 0.1, max_iter = 5000, verbose = TRUE)
                          
# For extra control over method-specific parameters pass a list as the "method" parameter:
# In largevis, gamma controls the balance of repulsion vs attraction
# The smallvis man page lists the method-specific parameters which can be controlled in this way
largevis_iris <- smallvis(iris, method = list("largevis", gamma = 1), perplexity = 25, epoch_callback = iris_plot, 
                          eta = 0.1, max_iter = 5000, verbose = TRUE)

# UMAP: see https://github.com/lmcinnes/umap
# UMAP also has extra parameters, but we use the defaults here
umap_iris <- smallvis(iris, method = "umap", perplexity = 25, eta = 0.01)

# use (scaled) PCA initialization so embedding is repeatable
tsne_iris_spca <- smallvis(iris, perplexity = 25, epoch_callback = iris_plot, Y_init = "spca")

# or initialize from Laplacian Eigenmap (similar to UMAP initialization)
tsne_iris_lap <- smallvis(iris, perplexity = 25, epoch_callback = iris_plot, Y_init = "lap")

# or initialize from normalized Laplacian eigenvectors (even closer to UMAP initialization)
tsne_iris_nlap <- smallvis(iris, perplexity = 25, epoch_callback = iris_plot, Y_init = "normlap")

# return extra information in a list, like with Rtsne
tsne_iris_extra <- smallvis(iris, perplexity = 25, epoch_callback = iris_plot, ret_extra = TRUE)

# more (potentially large and expensive to calculate) return values, but you have to ask for them specifically
tsne_iris_extra_extra <- smallvis(iris, perplexity = 25, epoch_callback = iris_plot,
                              ret_extra = c("P", "Q", "DX", "DY", "X"))

# Repeat embedding 10 times and keep the one with the best cost
tsne_iris_best <- smallvis_rep(nrep = 10, X = iris, perplexity = 25, ret_extra = TRUE)
iris_plot(tsne_iris_best$Y)

# Let smallvis pick a perplexity for you, using the Intrinsic Dimensionality Perplexity
tsne_iris_idp <- smallvis(iris, epoch_callback = ecb, perplexity = "idp", Y_init = "spca",
                          exaggeration_factor = 4)
                          
# Classical momentum optimization instead of delta-bar-delta
umap_iris_mom <- smallvis(iris, scale = FALSE, opt = list("mom", eta = 1e-2, mu = 0.8),
                          method = "umap", Y_init = "spca")

# L-BFGS optimization via the mize package
umap_iris_lbfgs <- smallvis(iris, scale = FALSE, opt = list("l-bfgs", c1 = 1e-4, c2 = 0.9),
                            method = "umap", Y_init = "spca", max_iter = 300)
                            
# Early Exaggeration
tsne_iris_ex <- smallvis(iris, eta = 100, exaggeration_factor = 4, stop_lying_iter = 100)

# and Late Exaggeration as suggested by Linderman and co-workers
tsne_iris_lex <- smallvis(iris, eta = 100, exaggeration_factor = 4, stop_lying_iter = 100,
                          late_exaggeration_factor = 1.5, start_late_lying_iter = 900) 
```

## Available Embedding Methods 

* [t-Distributed Stochastic Neighbor Embedding](https://lvdmaaten.github.io/tsne/).
* Sammon mapping and metric Multidimensional Scaling.
* Metric MDS with geodesic distances, like an iterative version of [Isomap](https://dx.doi.org/10.1126/science.290.5500.2319).
* [SNE](https://papers.nips.cc/paper/2276-stochastic-neighbor-embedding) and [Symmetric SNE (PDF)](https://www.cs.toronto.edu/~amnih/papers/sne_am.pdf)
* [Heavy-Tailed Symmetric SNE](http://papers.nips.cc/paper/3770-heavy-tailed-symmetric-stochastic-neighbor-embedding) (HSSNE).
* [Neighbor Retrieval Visualizer](http://www.jmlr.org/papers/v11/venna10a.html) (NeRV).
* [Jensen-Shannon Embedding](http://www.sciencedirect.com/science/article/pii/S0925231213001471) (JSE).
* [Weighted SNE using degree centrality](http://www.jmlr.org/proceedings/papers/v32/yange14.html) (wt-SSNE).
* [Elastic Embedding (PDF)](http://faculty.ucmerced.edu/mcarreira-perpinan/papers/icml10.pdf).
* [LargeVis](https://arxiv.org/abs/1602.00370) (the cost function, not the stochastic gradient descent part).
* [UMAP](https://arxiv.org/abs/1802.03426) (the cost function and calibration method): see [uwot](https://github.com/jlmelville/uwot) for a more functional version.
* [Alpha-Beta-SNE (ABSNE)](http://proceedings.mlr.press/v37/narayan15.html) and [ft-SNE](https://github.com/jiwoongim/ft-SNE) which
generalize t-SNE to a family of divergences.
* [Global t-SNE (GSNE)](https://github.com/gyrheart/gsne).

## Things To Be Aware Of

* March 23 2019: Methods that use the exponential function (e.g. NeRV, JSE, SSNE, 
ASNE) are now more robust, but sadly a lot slower, due to me implementing the
[log-sum-exp "trick"](https://statmodeling.stat.columbia.edu/2016/06/11/log-sum-of-exponentials/)
to avoid numeric underflow. This mainly helps JSE, which showed a tendency to
have its gradients suddenly explode. It's still difficult to optimize, though.
Perhaps [symmetric JSE](https://jlmelville.github.io/smallvis/nervjse.html#addendum_2:_jse_revisited)
can help under those circumstances.
* Feb 13 2018: the [UMAP paper](https://arxiv.org/abs/1802.03426) is out, but I 
have yet to read and understand it fully, so the implementation in `smallvis` 
currently relies
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

## My Idle Thoughts

`smallvis` exists mainly to satisfy my urge to answer the various, minor, 
[stamp-collecting](https://quoteinvestigator.com/2015/05/08/stamp/) questions 
that have occurred to me as I have read the dimensionality reduction literature. 
Those that I have cobbled together into something that demonstrates the use of
`smallvis` can be found at the 
[documentation](https://jlmelville.github.io/smallvis/) page.

## See Also

For large scale visualization in R see:

* The Barnes-Hut t-SNE package [Rtsne](https://cran.r-project.org/package=Rtsne)
* The [largeVis](https://cran.r-project.org/package=largeVis) package.
* [RSpectra](https://cran.r-project.org/package=RSpectra) can help spectral methods
scale.

Also of relevance are:

* [UMAP](https://github.com/lmcinnes/umap) (in Python)
* [uwot](https://github.com/jlmelville/uwot) a package implementing LargeVis 
and UMAP.
* [LargeVis](https://github.com/lferry007/LargeVis) (in C++)
* [Spectra](http://spectralib.org/), the C++ library that RSpectra wraps.
* [FIt-SNE](https://github.com/KlugerLab/FIt-SNE), an FFT-based t-SNE library. 
I have implemented the "late exaggeration" method that it uses. See 
[their paper](https://arxiv.org/abs/1712.09005) for more details.

Much of the code here is based on my [fork](https://github.com/jlmelville/rtsne) 
of Justin Donaldson's R package for [t-SNE](https://cran.r-project.org/package=tsne).

## License

[GPLv2 or later](https://www.gnu.org/licenses/gpl-2.0.txt). Any 
LargeVis-specific code (e.g. cost and gradient calculation) can also be 
considered [Apache 2.0](https://www.apache.org/licenses/LICENSE-2.0).
Similarly, Barnes-Hut and UMAP-related code is also licensed as 
[BSD 3-clause](https://opensource.org/licenses/BSD-3-Clause).
