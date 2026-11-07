# smallvis repository review, 2026-11-06

The previous cleanup substantially improved the repo, but **there are still correctness issues
worth fixing before adding methods**. I found one native termination failure and several
reproducible numerical bugs, mostly in less-used methods and combinations of options.

The project's architecture is broadly appropriate: three public entry points, readable method
implementations, a small native acceleration layer, and substantial explanatory articles. The next
improvements can preserve that character.

I reviewed `e65bcbaf6b71669c92efe75cfdd595739e894a31`. No repository files were changed during the
review. This document preserves the review delivered in chat, with repository-relative source
links. The [follow-up implementation plan](smallvis-correctness-followup-plan.md) records the
subsequent request to plan correctness fixes, cleanup, and allocation improvements. Its optional
features are excluded from implementation scope.

1. **High priority: Barnes-Hut subdivision can fail to terminate.**

   In [bh.h][tree], distinct points always trigger subdivision, without checking whether
   floating-point arithmetic can still separate them.

   Two x-coordinates, `1 + .Machine$double.eps` and `1 + 2*.Machine$double.eps`, with equal
   y-coordinates reproduce this. The standalone native probe produced an AddressSanitizer stack
   overflow; a public `smallvis()` call with that initialization timed out.

   This is an uncommon input, but finite-resolution termination is a basic requirement for this
   tree. Add a terminal representation for points that further subdivision cannot separate,
   preserving counts and self-exclusion. Test adjacent representable coordinates alongside exact
   duplicates.

2. **High priority: two AB-SNE branches calculate incorrect objective values, and this changes
   stopping behavior.**

   These errors survive gradient checks because they affect terms independent of the embedding.

   In [absneb0][absneb0], the cached term uses `P * log(P^alpha)` where the documented formula
   requires `P^alpha * log(P^alpha)`. With `alpha = 2, lambda = 2`, one probe reported **-2.350247**
   against an independently calculated **0.01215411**. A public optimization stopped after
   **one iteration** because the incorrect negative cost crossed the default `min_cost`.

   In [absneamb][absneamb], each column subtracts `N`, although there are only `N - 1` non-self
   terms. An eight-point public example with matching P and Q returned **-8**, instead of zero.

   Correct the constants and add scalar-reference tests, including zero divergence at matching
   distributions. This is more valuable than adding further trajectory snapshots.

3. **Coordinate input and precomputed distances still disagree for EE and g-SNE.**

   I swept all 48 built-in selectors and confirmed two separate representation bugs:

   - [Elastic Embedding][ee] retains a compact `dist` object where subsequent arithmetic expects
     a full matrix. R recycles its entries, and its normalization also differs. At identical
     output coordinates, one comparison gave costs **26.24927 versus 51.12120**, with a recycling
     warning.
   - [g-SNE][gsne] uses squared distances for coordinate input but ordinary distances for `dist`
     input when constructing `Phat`.

   Convert EE's distances to a matrix before normalization, and square g-SNE's precomputed
   distances. Extend the existing input-equivalence tests to these method-specific quantities,
   gradients, and costs.

4. **Gaussian calibration still mishandles zero distances and underflow.**

   With `scale = FALSE`, all-coincident observations cause the [R calibrator][r-calibration] to
   fail with "missing value where TRUE/FALSE needed"; the native route instead produces its
   fallback.

   More seriously, `X = diag(24), perplexity = 1` produces non-finite affinities and a NaN cost in
   **both** backends. When all weights underflow, the [native code][cpp-calibration] sets entropy
   to zero. That falsely satisfies the target `log(1)`, allowing an all-zero row into
   normalization.

   Treat zero or non-finite normalization as calibration failure independently of the entropy
   residual. Give the R and C++ paths consistent behavior for unattainable targets and coincident
   inputs.

5. **The line-search compatibility checks remain incomplete.**

   [USNE][usne] and [pseudo-separable t-SNE][pstsne] accept L-BFGS while supplying forces that do
   not differentiate their reported scalar.

   Across four finite-difference step sizes, one directional comparison gave:

   | Method | Supplied force | Scalar derivative |
   | --- | ---: | ---: |
   | USNE | 0.03768988 | 0.03326312 |
   | Pseudo-separable t-SNE | 0.1411800 | 0.1915116 |

   Extend the existing optimizer restrictions to these methods, or implement matched objectives.
   USNE also needs explicit diagonal handling: setting `gr_eps = 0` currently produces non-finite
   gradients, so that setting alone is not a working escape route.

6. **Dynamic HSSNE has two state-transition bugs.**

   Its [epoch hook][dhssne] checks `alpha + 0.1` before admitting the candidate `alpha - 0.1`.
   A public run starting at `alpha = 0.1` selected zero, despite the constructor's positive lower
   bound.

   The hook also changes `alpha` while retaining weights calculated using the previous value.
   In a probe changing 0.5 to 0.6, the cached and freshly recalculated gradients differed by
   **0.01045**.

   Check the actual candidate against the bound and invalidate dependent weights when accepting
   it. Test the first gradient after the transition against a fresh calculation.

7. **Unsupported multiscale combinations should fail deliberately.**

   The [multiscale branch][multiscale] returns only P, bypassing quantities that some methods
   require.

   Even a single-scale request such as `perplexity = list("multiscale", 4)` breaks weighted t-SNE,
   weighted SSNE, unnormalized t-SNE, and pseudo-separable t-SNE. Bandwidth-transfer methods can
   initialize without their required `beta` values and then fail during gradient calculation.

   Support combinations whose meaning is clear and inexpensive to implement. Reject the others
   early with a useful message. Defining a transferred bandwidth for a mixture of scales would
   be a separate methodological addition.

8. **Neighbor extraction still assumes self occupies the first position in two places.**

   [UMAP calibration][umap-neighbors] and [geodesic neighborhood construction][geodesic-neighbors]
   retain positional self-removal assumptions that the corrected Gaussian and kNN paths have
   already addressed.

   With duplicated observations and UMAP `k = 4`, the distance-input route produced four non-self
   neighbors for some rows, while other rows had three. Separately,
   `method = list("gmmds", k = 1)` fails because subsetting drops matrix dimensions.

   Filter self by observation index and preserve dimensions. Reuse the established approach in
   the corrected neighbor paths.

There are also several smaller cleanup opportunities that fit the existing scope:

- **Finish forwarding backend controls.** [HSSNE][hssne-controls] omits `use_cpp` when
  constructing its inherited initializer. [AB-SNE's special branches][absne-controls] omit both
  `use_cpp` and `n_threads`, while its general branch hard-codes native calibration. These controls
  should behave consistently across parameter values.
- **Remove confirmed leftovers.** The first [t-ASNE definition][shadowed-tasne] is still shadowed
  by the later definition. `user_idp_perps()` and `graph_to_sparse()` have no consumers in the
  inspected package, tests, or documentation sources. Three calls in [fdiv.R][argument-typo] also
  contain `use_cpp == use_cpp` where a named argument was evidently intended; replacement
  initializers currently conceal the mistake.
- **Fix one misleading test loop.** The [native warning test][warning-test] iterates over `X` but
  always calls `smallvis(diag(24), ...)`. Its second iteration does not exercise the duplicate-input
  case it appears to cover.

**The clearest remaining efficiency improvement is eliminating more quadratic copies.** The
previous t-SNE gradient change establishes a suitable pattern. [MDS gradients][mds-copies] still
copy both input distance matrices into vectors on every call: 400 MB of additional payload at
5,000 observations. [Student weights][weight-copy] still stage an entire matrix before copying
it into R storage. Borrowing protected input buffers and writing directly into allocated R output
would be focused improvements. Those are allocation calculations, not measured speedup claims;
use the existing benchmark machinery to establish the benefit.

**The test suite needs stronger mathematical coverage across method families.** It already has
useful finite-difference tests, native ownership checks, and public regressions. The AB-SNE
findings demonstrate the missing complement: independently evaluated scalar objectives,
normalization identities, and boundary cases. Extend those checks to related methods when fixing
a family-level bug.

Given the new tree failure, a small native CI job now earns its maintenance cost. The existing
standalone tree and threading harnesses are useful, but the configured workflows do not execute
them. Compile and run those harnesses with sanitizers, including the new boundary regression.

For architecture, I would retain the method-specific functions and ordinary cost lists. They make
the equations accessible. The main complexity concentrations are the optimization lifecycle and
affinity preparation; extract small helpers where a fix establishes a clear responsibility. There
is little evidence that an object-system rewrite, plugin registry, or generalized execution
framework would improve this project.

Two optional additions would fit its purpose:

- Return the **effective backend and stopping reason** with diagnostics, making experiments easier
  to interpret.
- Eventually expose **precomputed affinities or neighbors** for controlled comparisons. This is
  useful, but it adds a public contract and should follow the correctness work.

A persistent experiment database, distributed execution, broad benchmark service, or expanded
Barnes-Hut dimensionality would be materially larger commitments. None is needed to meet the
quality bar described in the review request.

Validation used a fresh temporary source build on Linux with R 4.5.2. `R CMD check --no-manual`
finished **OK**, with **968 passing expectations and no failures, warnings, or skips**. All
handwritten C++ translation units passed strict GCC 15.2 and Clang 21.1 checks. The existing native
harnesses passed focused ASan/UBSan runs; the additional adjacent-coordinate probe exposed the
tree failure. Leak detection was disabled for those sanitizer runs. Windows checks and a full
article-site rebuild were not run.

[tree]: ../smallvis/src/bh.h#L81
[absneb0]: ../smallvis/R/absne.R#L196
[absneamb]: ../smallvis/R/absne.R#L304
[ee]: ../smallvis/R/ee.R#L24
[gsne]: ../smallvis/R/gsne.R#L34
[r-calibration]: ../smallvis/R/perplexity.R#L339
[cpp-calibration]: ../smallvis/src/perplexity.cpp#L176
[usne]: ../smallvis/R/sne.R#L768
[pstsne]: ../smallvis/R/sne.R#L669
[dhssne]: ../smallvis/R/sne.R#L400
[multiscale]: ../smallvis/R/perplexity.R#L77
[umap-neighbors]: ../smallvis/R/umap.R#L276
[geodesic-neighbors]: ../smallvis/R/perplexity.R#L478
[hssne-controls]: ../smallvis/R/sne.R#L280
[absne-controls]: ../smallvis/R/absne.R#L17
[shadowed-tasne]: ../smallvis/R/sne.R#L1065
[argument-typo]: ../smallvis/R/fdiv.R#L83
[warning-test]: ../smallvis/tests/testthat/test_native_perplexity.R#L75
[mds-copies]: ../smallvis/src/gradients.cpp#L198
[weight-copy]: ../smallvis/src/gradients.cpp#L151
