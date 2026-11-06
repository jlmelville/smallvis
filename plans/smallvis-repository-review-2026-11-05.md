# smallvis repository review

**Fix native memory safety and numerical correctness before adding methods.** The three exported
functions, the readable R implementations of the objectives, and the long-form experiments make a
coherent project. Its biggest liability is that changing an input representation, backend, or
diagnostic request can silently change the mathematical problem. The existing R implementations
also give the project a relatively cheap way to test the accelerated paths.

I reviewed the working checkout at `3382ecf`. It included a local deletion of a shadowed t-ASNE
definition, new t-ASNE tests, and untracked experiment plans. The deletion removes an earlier
definition that R was already replacing. The new tests of normalization are worth keeping. I
treated the plans as development context and the research articles as part of the public project.

## Fix the native boundaries first

1. **Default perplexity can crash R on a small dataset.** With 12 rows, `use_cpp = TRUE`, and the
   default perplexity of 30, Gaussian calibration fails. The fallback passes `begin() + 30` to
   `partial_sort` on a vector of length 12, then indexes it up to 30. I reproduced invalid reads,
   an invalid write, and a segfault under Valgrind through `smallvis()`, even at `max_iter = 0`.
   Validate each Gaussian target against the number of available non-self neighbors before
   calibration, including per-point and multiscale targets where supported. Bound the native
   selection independently. Kernel bandwidths and UMAP neighbor counts have different contracts.
   [Dense fallback][fallback] · [R caller][calibration]

2. **Public output dimensions reach kernels that assume exactly two coordinates.** On the same
   input, `method = "mmds", k = 1, use_cpp = TRUE` caused invalid reads and writes. The kernel
   indexes pairs even though its coordinate and gradient storage have one element per point.
   Exact t-SNE weights and gradients also stride by two. Barnes–Hut infers the number of points
   from the flattened embedding length divided by two and uses it to index the sparse graph.
   Route dense embeddings of other dimensions to the working general R path, or generalize the
   dense loops. Reject non-2D Barnes–Hut at entry and check the actual matrix shape at the native
   boundary. An octree is a separate feature. [Dense kernels][kernels] · [BH entry][bhentry]

3. **Threaded calibration races on an ordinary integer.** Helgrind reported conflicting accesses
   to `bad_perp++` from two workers calibrating 24 equidistant observations. The dense and kNN
   versions both share the counter by reference. R then discards the count. Failure handling has
   a second defect even with one thread. On that fixture, all 24 dense rows took the fallback and
   all 24 acquired positive self-weight. Use per-worker counts reduced after joining, exclude self
   before selecting neighbors, and report failed rows once. Define finite fallback weights for ties
   and all-zero distances. [Dense counter and selection][fallback] · [kNN counter][knncounter]

4. **Barnes–Hut loses probability mass even at `theta = 0`.** On one 24-point fixture with an
   all-coincident initialization, the reported KL was **−41.6190**. Direct evaluation using the
   very same sparse P and exact Student weights gave **1.5359**. With distinct generic coordinates
   on the same fixture, the two calculations agreed to rounding. The tree skips an entire
   coincident leaf including the other observations, reinserts only one point when splitting a
   leaf with multiplicity, and can absorb a point at an internal centroid without passing it to a
   child. A four-point example with distinct coordinates also exposed the internal-centroid case.
   Preserve subtree counts and exclude exactly the query observation from the normalizer. Compare
   both the force and the cost to dense sums, using the same P and several insertion orders.
   [Tree insertion and self-interactions][tree] · [Force and normalizer][bhforce]

   The same `eps` serves as a coordinate-merging tolerance and a numerical safeguard, and the cost
   and gradient calls pass it differently. With the documented choice `eps = 1e-3` and initial
   coordinates at the default `1e-4` scale, the computed repulsion on a 24-point probe disappeared
   to rounding. Keep tree membership separate from stabilization, use consistent settings in cost
   and gradient, and check finite coordinates before tree construction. [BH cost and gradient][bhr]

## Make comparable experiments compute comparable quantities

**Coordinates and precomputed distances.** With scaling disabled, `X` and `dist(X)` produced
different Gaussian P matrices. The distance branch squares values for `"gauss"`, while the default
caller sends `"gaussian"`. The kNN branch is worse. It computes a thresholded column and then
writes the constant 1 to the entire column. On tie-free data with 24 observations and four
neighbors, the precomputed route connected all 24 entries of every row, including self. Fix both
branches, define the direction of neighbor edges and handling of ties, and compare P before
optimization. The dense native Gaussian route also needs either a path for `dist` input or a
documented fallback. [Gaussian distances][affdist] · [Distance kNN graph][knngraph]

**Method dispatch.** For exact t-SNE, the top-level `inp_kernel = "knn"` gave the default Gaussian
P exactly. Putting the same setting inside `method = list("tsne", inp_kernel = "knn")` changed P.
In the reverse direction, `bh = TRUE` with the bare string produced a sparse P, while
`method = list("tsne")` silently ran dense t-SNE under the same controls. Explicitly selecting a
Gaussian input with Barnes–Hut reached `@x` on a base matrix and errored. These are particularly
costly mistakes in a package used to compare algorithms. Normalize string and list specifications
once, keep one allowlist and constructor route, explicitly set precedence, and reject or implement
unsupported combinations before expensive work. Preserve the custom function route and method
defaults. Also correct the usage article's historical instruction to use `method = "bhtsne"`.
The public selector accepts `method = "tsne", bh = TRUE`. [Character dispatch][dispatch] ·
[List dispatch][listdispatch] · [Usage article][usage]

**Repeated starts.** `smallvis_rep()` re-runs input preparation for every start. For two starts on
a fixed 600 × 30 dataset with Barnes–Hut and approximate neighbors, 7,098 entries of P differed
by more than `1e-10`. With exact neighbors, the two P matrices agreed. Evaluating two embeddings
against different P matrices defeats the promise to select the best initialization. Prepare one
input graph and one set of affinities for the built-in method, then initialize each optimization
from fresh coordinates and fresh optimizer state. Check that all repeats share P and that reported
costs can be reproduced from it. This also saves repeated neighbor search and calibration.
[Repeat loop][repeat] · [Neighbor selection][neighbors]

**Multiscale normalization.** The multiscale helper accepts the caller's normalization and
symmetrization settings, then hard-codes symmetric globally normalized t-SNE affinities. On ASNE,
ordinary P rows summed to one. Wrapping the same perplexity as a one-scale multiscale request
changed the row sums to 0.021–0.077. Forward the method's settings or explicitly limit multiscale
to supported methods. A one-scale request should reproduce the corresponding single-scale P.
[Multiscale averaging][msp]

**Spectral initialization.** The Laplacian eigenmap requires the largest algebraic eigenvalues of
the random-walk matrix after discarding its stationary mode. The implementation selects by
magnitude and assumes the first returned vector is the one to discard. On a ten-node cycle, both
installed eigensolver paths returned modes with eigenvalues 1 and about −0.809. The two desired
modes both have eigenvalue about +0.809. The existing fixture checks agreement with coordinates
but does not establish this ordering. Select and sort the intended modes explicitly. A small
regular graph gives an exact subspace test unaffected by eigenvector signs or rotations.
[Laplacian eigenmap][spectral] · [Spectral tests][spectraltests]

There is a separate public side effect here. Holding X and the seed fixed, adding `"V"` to
`ret_extra` changed pairwise distances in the initial Laplacian embedding by up to 1.117. The
initializer prefers V when the diagnostic request happens to retain it and otherwise uses P.
Choose the intended graph independently of which objects the user requests for inspection.
[Graph selection][spectralinput]

**Preprocessing.** These small fixes affect common inputs and deserve their own regression cases.

- `scale = TRUE` filters a matrix element by element. The 24 × 4 matrix in the probes became
  0 × 1 before scaling. Select columns by variance with dimension-preserving indexing. Default
  `"absmax"` also divides by `abs(max(X))`. A centered vector with entries −9 and nine copies of
  +1 kept maximum absolute magnitude 9. Use the documented maximum absolute value and handle a
  zero denominator deliberately. The global scale factor matters for distance-based methods even
  when Gaussian perplexity calibration can absorb it. [Input scaling][scale]
- Skipping PCA leaves the outer `pca` flag set. A 24 × 4 input with `pca = TRUE, Y_init = "pca"`
  and the default 50 initial dimensions returned the first two input columns and reported 50 PCs.
  Even when PCA runs, the initializer always slices two columns, including for `k = 3`. Carry the
  applied preprocessing state through to initialization and metadata, and honor `k`.
  [Preprocessing decision][pcaprep] · [PCA initializer][pcainit]
- The truncated-SVD branch returns `lambda` without assigning it. A 24 × 14 matrix with
  `pca = "whiten", initial_dims = 3` failed with `object 'lambda' not found`. Derive the retained
  eigenvalues from the returned standard deviations and check whitened variances against the
  full-SVD path. [PCA scores and whitening][whiten]

## Optimizers and diagnostics

- **AdaMax changes its recurrence across dimensions.** Its comment correctly describes an
  elementwise maximum, but the code takes a row maximum across all embedding dimensions. With
  first gradient `(1, 4)` and unit learning rate, the update was `(−0.25, −1)` instead of
  approximately `(−1, −1)`. The stored epsilon is also unused. Use an elementwise maximum and the
  intended denominator safeguard. A two-step test with unequal coordinate gradients and a zero
  gradient is enough. [AdaMax update][adamax]
- **Some stabilized forces are offered to line-search optimizers as exact gradients.** The finite
  difference tests for UMAP and LargeVis deliberately set `gr_eps = 0`. At their ordinary nonzero
  settings, the directional derivative of the reported scalar differs from the force handed to
  `mize`. For LargeVis, I measured 132.64 versus 76.26 on the same coordinates and direction.
  Stabilized forces are a legitimate experiment with a gradient-only optimizer. Give a
  function-based optimizer a matched scalar and derivative, or reject that combination with a
  precise explanation, as the package already does for exaggeration. [LargeVis force][largevis] ·
  [UMAP force][umap] · [Optimizer guard][optguard]
- **Final costs can be missing or stale when an optimizer stops.** Termination before the first
  epoch with `opt = list("l-bfgs", max_gr = 1)` returned `iter = 1` and an empty cost history.
  `smallvis_rep()` then errored with `argument is of length zero`. After a prior epoch the same
  ranking expression can consume the earlier cost. Finalize the objective at the coordinates
  being returned and rank repeats using that value. Preserve the epoch history for its own
  purpose. [Termination branch][terminate] · [Ranking][repeat]
- **Perplexity stepping fails when there is nothing to step through.** A valid target of 10 on
  24 observations failed with `object 'perps' not found`. Also, temporary `tol = 0` and `epoch`
  settings survive into the final run if the caller omitted those arguments. Use an empty step
  sequence and build preliminary calls from a copy of the final arguments. Check the final run
  against a direct call from the same coordinates. [Perplexity stepping][steps]
- **Several published diagnostics have the wrong meaning.** The dense C++ route discards the
  calibration result, then reads the undefined `x2ares` for `ret_extra = "beta"`. The degree
  diagnostics reverse the documented row and column sums and halve only one term of their
  advertised mean. In squared metric MDS, the DX/DY export computes square roots but returns the
  squared matrices. Each has a cheap independent reference. [Affinity diagnostics][calibration] ·
  [Degree definitions][degdocs] · [Squared MDS export][mdsexport]
- **ASNE per-observation KL contributions mix axes.** The cached P log P term uses column sums,
  while the P log Q term uses row sums. The total still agrees, which lets aggregate gradient
  checks pass. In the 24-point probe the returned vector contained −0.466, while every
  independently computed row KL was positive. Align the reductions and document what `costs`
  returns for each family, including the single scalar currently returned by Barnes–Hut.
  [Input entropy cache][plogp] · [ASNE cost][asnecost]
- **The Gaussian output stabilization still permits all-zero rows.** The self-distance sets the
  maximum log weight to zero, and the diagonal is discarded only after exponentiation. For ASNE
  initialized at equally spaced coordinates 100 units apart, the result contained NaN costs.
  Exclude self when computing normalization shifts and use the axis appropriate to the objective.
  [Gaussian output probabilities][expq]
- **t-EE exaggeration changes a field its force never reads.** Its hook multiplies V, while the
  gradient reads P. At a fixed initial layout, learning rate, momentum schedule, and five
  iterations, factors 1 and 12 gave identical coordinates. Decide whether to support exaggeration
  for t-EE, then either scale the force term it uses or clearly reject the option. [t-EE][tee]

## C++ efficiency and maintenance

1. **Start with avoidable quadratic copies.** Every dense t-SNE gradient call copies P and W into
   new vectors. At 5,000 points those two copies alone require 400 MB of payload per call, before
   the live R matrices and other allocations. Dense calibration also creates an N × N vector, an
   R matrix, and a transposed R result. Allocate outputs and acquire pointers to materialized R
   inputs on the main thread, then let workers read those buffers and write disjoint ranges.
   Measure peak memory and iteration time before pursuing more elaborate changes. The 400 MB is
   calculated from the allocations, with no timing claim. [Gradient copies][copies] ·
   [Calibration output][denseoutput]
2. **Use one thread launcher.** There are two near-identical launchers specialized by lambda
   arity and two more bespoke launch loops in perplexity calibration. With more threads than
   points the shared helper can assign all work to the last worker. If thread creation fails
   partway, destruction of already joinable threads can terminate R. A bounded worker count,
   one ordinary worker signature with a chunk ID, and exception-safe joining address all three
   issues. Use owning smart pointers for quadtree children when touching insertion. A persistent
   thread pool or allocator should wait for a measured bottleneck. [Parallel loops][threads] ·
   [Calibration launcher][perpthreads]
3. **Remove a broken unused conversion.** The private squared-distance-to-Student-weight
   conversion launches N scalar transforms for an N² matrix, leaving almost all output zero.
   The only R reference is its private wrapper. Searching package code, tests, docs sources and
   experiment scripts found no consumer. Remove the R helper and native implementation together
   and regenerate the bridge. By contrast, the 13 apparently unreachable built-in optimizers are
   selected dynamically, including DBD. [Conversion][conversion] · [Dispatch][optdispatch]

## Scope and order of work

Keep the three exports, the small method-specific cost lists, the long-form articles, and the R
versions of the equations. The R source is valuable both to a reader and as a reference
implementation. The most productive structural changes are a shared method-construction path,
input preparation that can be reused across repeats, and one final-cost calculation. Each fixes
an observed failure. The approximately 600-line main driver is worth extracting along those
boundaries as the fixes land. A new object system would add another interface to learn.

1. Land the two bounds fixes, race fix and tree correction. Run the small failing calls in
   isolated processes under a memory or race checker and compare exact forces with dense R.
2. Restore input, graph, initializer and optimizer equivalence. Use small fixtures that assert
   which mathematical problem is being solved. Keep the finite-difference tests and a few
   end-to-end trajectories. Add tests for P, zero self-weight, actual output dimension, and
   the final cost of the returned coordinates before changing historical coordinate fixtures.
3. For users, add a short capability table to the existing usage article covering input forms,
   dimensions, affinity kernels, exact/BH behavior, and diagnostic availability. An effective
   method/backend field in `ret_extra` would make experiments easier to interpret. If demand for
   comparing objectives on the same graph extends beyond repeated starts, consider a small
   precomputed-affinity or neighbor-input interface with an example. First prove the internal
   reuse with `smallvis_rep()`.
4. Keep the current Linux/Windows checks, formatting, lint and package site. A focused native
   diagnostic on changes to these kernels would earn its cost. For now I would defer a plugin
   system, persistent experiment store, distributed execution, a broad benchmark service, and
   extending Barnes–Hut to every method. The local proposal for staged mize optimization is a
   separate numerical feature and can follow these fixes.

Two small maintenance details affect the experience of a stranger finding the repo. The README
points readers to per-file licensing notices, but the LargeVis and UMAP R files have no such
notices. Give readers a direct pointer to the applicable terms and attribution. Also, the sole
runtime use of the GitHub dependency `vizier` is the default plotting callback. If installation
friction matters, making that dependency optional with a base-graphics fallback would preserve the
interactive experience. Both changes can be handled independently of the numerical work.
[README][readme] · [Plot callback][plot]

## Verification

I built and installed the working package from source in a temporary directory on Linux with R
4.5.2 and GCC 15.2. `R CMD check --no-manual` finished with status OK using installed dependencies.
Its attempted CRAN index lookup was unavailable. The separate installed test run passed 579
expectations, with seven partial-argument-matching warnings and no failures, errors, or skips.
The additional results above come from public API probes, direct comparisons of internal kernels
where indicated, Valgrind Memcheck, and Helgrind. The review did not run Windows or a complete
sanitizer matrix. The generated Rcpp bridges were used for tracing calls and excluded from
handwritten C++ style judgments.

For a compact reproduction of the Gaussian input mismatch, compare P before taking any steps.

```r
set.seed(42)
x <- matrix(rnorm(24 * 4), 24)
affinity <- function(input) {
  smallvis::smallvis(
    input, scale = FALSE, perplexity = 4, max_iter = 0,
    ret_extra = "P", epoch_callback = FALSE, verbose = FALSE
  )$P
}
max(abs(affinity(x) - affinity(dist(x))))
# 0.002929218
```

The architecture map found three public roots and 162 R functions in this worktree. Manual review
of dynamic method dispatch and the untracked experiment scripts mattered for interpreting
reachability. No checkout files were modified by this review.

[fallback]: ../smallvis/src/perplexity.cpp#L177
[calibration]: ../smallvis/R/perplexity.R#L121
[kernels]: ../smallvis/src/gradients.cpp#L11
[bhentry]: ../smallvis/src/bh.cpp#L193
[knncounter]: ../smallvis/src/perplexity.cpp#L82
[tree]: ../smallvis/src/bh.h#L75
[bhforce]: ../smallvis/src/bh.cpp#L6
[bhr]: ../smallvis/R/sne.R#L112
[affdist]: ../smallvis/R/perplexity.R#L234
[knngraph]: ../smallvis/R/perplexity.R#L469
[dispatch]: ../smallvis/R/smallvis.R#L1559
[listdispatch]: ../smallvis/R/smallvis.R#L1862
[repeat]: ../smallvis/R/iterated.R#L75
[neighbors]: ../smallvis/R/perplexity.R#L511
[msp]: ../smallvis/R/perplexity.R#L577
[spectral]: ../smallvis/R/init.R#L46
[spectraltests]: ../smallvis/tests/testthat/test_spectral.R#L26
[spectralinput]: ../smallvis/R/smallvis.R#L1101
[scale]: ../smallvis/R/input.R#L25
[pcaprep]: ../smallvis/R/input.R#L44
[pcainit]: ../smallvis/R/init.R#L18
[whiten]: ../smallvis/R/init.R#L137
[adamax]: ../smallvis/R/opt.R#L515
[largevis]: ../smallvis/R/largevis.R#L40
[umap]: ../smallvis/R/umap.R#L54
[optguard]: ../smallvis/R/smallvis.R#L1051
[terminate]: ../smallvis/R/smallvis.R#L1270
[steps]: ../smallvis/R/iterated.R#L224
[degdocs]: ../smallvis/R/smallvis.R#L628
[mdsexport]: ../smallvis/R/mmds.R#L115
[plogp]: ../smallvis/R/sne.R#L31
[asnecost]: ../smallvis/R/cost.R#L99
[expq]: ../smallvis/R/cost.R#L43
[tee]: ../smallvis/R/ee.R#L123
[copies]: ../smallvis/src/gradients.cpp#L184
[denseoutput]: ../smallvis/src/perplexity.cpp#L268
[threads]: ../smallvis/src/threads.h#L40
[perpthreads]: ../smallvis/src/perplexity.cpp#L216
[conversion]: ../smallvis/src/gradients.cpp#L166
[optdispatch]: ../smallvis/R/opt.R#L6
[usage]: ../smallvis/vignettes/articles/usage.Rmd#L278
[readme]: ../README.md#L52
[plot]: ../smallvis/R/epoch.R#L32
