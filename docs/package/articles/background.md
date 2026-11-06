# About smallvis

`smallvis` exists to make it easier to experiment with neighborhood
embedding methods and to understand how their component parts affect the
results. The [usage
guide](https://jlmelville.github.io/smallvis/package/articles/usage.md)
cover the current interface; this page collects the background and dated
development notes.

## Earlier implementation notes

- March 23 2019: Methods that use the exponential function (e.g. NeRV,
  JSE, SSNE, ASNE) are now more robust, but sadly a lot slower, due to
  me implementing the [log-sum-exp
  “trick”](https://statmodeling.stat.columbia.edu/2016/06/11/log-sum-of-exponentials/)
  to avoid numeric underflow. This mainly helps JSE, which showed a
  tendency to have its gradients suddenly explode. It’s still difficult
  to optimize, though. Perhaps [symmetric
  JSE](https://jlmelville.github.io/smallvis/package/articles/nervjse.html#addendum-2-jse-revisited)
  can help under those circumstances.
- Feb 13 2018: the [UMAP paper](https://arxiv.org/abs/1802.03426) is
  out, but I have yet to read and understand it fully, so the
  implementation in `smallvis` currently relies on my examination of the
  UMAP source code, with some much-appreciated clarification from UMAP
  creator [Leland McInnes](https://github.com/lmcinnes). Expect some
  bugs, and any horrifically bogus results should be double-checked with
  the output of the [official UMAP
  implementation](https://github.com/lmcinnes/umap) before casting
  calumnies on the quality of UMAP itself.
- LargeVis requires the use of a `gamma` parameter, which weights the
  contribution of attractive and repulsive contributions to the cost
  function. In the real LargeVis, it is recommended to set this value to
  `7`, but this relies on the specifics of the stochastic gradient
  descent method. In `smallvis`, this value is very dataset dependent:
  the more data you have, the smaller gamma should be to avoid
  over-emphasising repulsive interactions.
- LargeVis partitions each pairwise interaction into either an
  attractive or repulsive contribution. In `smallvis`, each interaction
  is a combination of both.
- Both UMAP and LargeVis use a classic stochastic gradient descent
  approach with a decaying learning rate. The implementation in
  `smallvis` uses the same delta-bar-delta optimization method used in
  t-SNE. It works well in my experience, but may require some tuning and
  more iterations compared to optimizing t-SNE.
- In this setting, LargeVis and UMAP gradient requires quite a large
  value of epsilon to avoid division by zero and get decent results.
  It’s hard-coded to 0.1 in the LargeVis source code, so I have used the
  same value by default in `smallvis`. It can be controlled with the
  method option `gr_eps`.
- Mainly for my own benefit, there is also a
  [theory](https://jlmelville.github.io/smallvis/package/articles/theory.md)
  page showing a comparison of cost functions and gradients. Also, some
  material on the various
  [spectral](https://jlmelville.github.io/smallvis/package/articles/spectral.md)
  methods, which justifies the use of Laplacian Eigenmaps (a bit).

## Background

*August 11 2024* **The Turbo Championship Edition Update**. I have
briefly brought `smallvis` back from the dead to speed it up a bit. I
have added:

- Barnes-Hut t-SNE. This will scale up to larger datasets and it is
  feasible to run it on the full MNIST digits dataset (i.e. 70,000
  items). This uses a (2D only) C++ translation of the cython
  implementation in the Python
  [openTSNE](https://github.com/pavlin-policar/openTSNE) package
  originally authored by [Pavlin
  Poličar](https://github.com/pavlin-policar). It’s BSD 3-clause
  licensed (and can be found in `src/bh.h`). Use it with
  `method = "tsne", bh = TRUE`. The degree of approximation can be
  controlled with `theta`. Be aware that it’s not as fast as e.g.
  [Rtsne](https://github.com/jkrijthe/Rtsne), at least during the
  optimization step. That package uses Laurens van der Maaten’s original
  C++ code which I am very unsure can be redistributed with R code due
  to its BSD 4-clause license. I would love to be wrong about that
  though! The Quad Tree implementation could be used with other
  embedding methods, but I haven’t got round to implementing that yet.
- A C++ multi-threaded perplexity search using only the nearest
  neighbors of each point (3 times the perplexity) is used with BH
  t-SNE.
- My own [rnndescent](https://cran.r-project.org/package=rnndescent)
  package replaces FNN for nearest neighbor search. Apart from being a
  monument to my ego, it can be faster than FNN for brute force search
  because it can be multi-threaded (use `n_threads` to control this).
  Also, approximate nearest neighbor search becomes quite important with
  larger datasets.
- For exact search, I have started adding multi-threaded C++ code to
  calculate the gradient. Set `use_cpp = TRUE` to use this. It’s not as
  big a win in speed up as you might hope because the R code is using
  some very optimized linear algebra for some steps which will blow my
  puny C++ code out of the water. However in many cases the linear
  algebra libraries won’t be mult-threaded so sheer brute force threads
  can overcome this. Just don’t expect setting `n` threads to give you
  `n` times the speed. Like Barnes-Hut, this requires me to implement
  the gradients in C++ for each method and I haven’t done that yet.
- [irlba](https://cran.r-project.org/package=irlba) is now a dependency
  for doing PCA on larger datasets.

I will probably at least attempt to apply Barnes-Hut and multi-threading
to some other methods. On the other hand, it’s taken me five years to
get back to this, so don’t hold your breath.

## My Idle Thoughts

`smallvis` exists mainly to satisfy my urge to answer the various,
minor,
[stamp-collecting](https://quoteinvestigator.com/2015/05/08/stamp/)
questions that have occurred to me as I have read the dimensionality
reduction literature. Those that I have cobbled together into something
that demonstrates the use of `smallvis` can be found at the
[articles](https://jlmelville.github.io/smallvis/package/articles/index.md)
page.

## See Also

For large scale visualization in R see:

- The Barnes-Hut t-SNE package
  [Rtsne](https://cran.r-project.org/package=Rtsne)
- The [largeVis](https://cran.r-project.org/package=largeVis) package.
- [RSpectra](https://cran.r-project.org/package=RSpectra) can help
  spectral methods scale.

Also of relevance are:

- [UMAP](https://github.com/lmcinnes/umap) (in Python)
- [uwot](https://github.com/jlmelville/uwot) a package implementing
  LargeVis and UMAP.
- [LargeVis](https://github.com/lferry007/LargeVis) (in C++)
- [Spectra](http://spectralib.org/), the C++ library that RSpectra
  wraps.
- [FIt-SNE](https://github.com/KlugerLab/FIt-SNE), an FFT-based t-SNE
  library. I have implemented the “late exaggeration” method that it
  uses. See [their paper](https://arxiv.org/abs/1712.09005) for more
  details.

Much of the code here is based on my
[fork](https://github.com/jlmelville/rtsne) of Justin Donaldson’s R
package for [t-SNE](https://cran.r-project.org/package=tsne).

## License

[GPLv2 or later](https://www.gnu.org/licenses/gpl-2.0.txt). Any
LargeVis-specific code (e.g. cost and gradient calculation) can also be
considered [Apache 2.0](https://www.apache.org/licenses/LICENSE-2.0).
Similarly, Barnes-Hut and UMAP-related code is also licensed as [BSD
3-clause](https://opensource.org/licenses/BSD-3-Clause).
