# smallvis correctness and maintenance implementation plan

**Status** Phases 1 through 3 completed on 2026-11-05. Implementation stops before Phase 4, as
requested.

**Next work** Phase 4 remains pending. Resume its optimizer and diagnostic work only when requested.

The [full repository review](smallvis-repository-review-2026-11-05.md) records the observations and
source locations at `3382ecf`. Its clean Linux package check and 579 passing test expectations are
the baseline. Existing edits to `smallvis/R/sne.R`, the new t-ASNE tests, and the other files under
`plans/` belong to ongoing work. Recheck the worktree before starting a phase.

The objective is to make the existing embedding experiments safe and comparable. Keep the three
public entry points, the readable R equations, and the method-specific cost lists. Ship each fix
with its named regression. Separate changes that alter scientific results from performance-only
changes so future readers can identify why a result moved.

## Phase 1. Make native calibration, dimensions, and worker lifetime safe

The affected sources are the [calibration kernels](../smallvis/src/perplexity.cpp),
[dense gradients](../smallvis/src/gradients.cpp), [thread helper](../smallvis/src/threads.h),
and their R callers. Use the existing testthat suite for public behavior. Run the dangerous
pre-fix cases in fresh processes against a freshly installed copy of the candidate package.

1. **Bound Gaussian perplexity before selecting neighbors.** Apply the available non-self neighbor
   limit to every requested Gaussian target, including vector and multiscale forms where supported.
   Check the native fallback's iterator bounds as well. Preserve the separate meanings of sigma,
   kNN degree, and UMAP neighborhood size.
   **Regression:** With `X = as.matrix(iris[1:12, 1:4])`, `use_cpp = TRUE`, and `max_iter = 0`,
   the default 30 must produce an informative R error and leave the process usable. On 12 rows,
   test 11 as valid and 12 as invalid. Include one invalid element in a supported vector or
   multiscale request. Exercise the native bound directly in an isolated process under Memcheck.
2. **Enforce the actual output dimension.** For dense MDS and t-SNE, prefer the existing general R
   calculation when a 2D specialization cannot handle the requested shape. For Barnes–Hut, reject
   dimensions other than two before building affinities and guard the actual matrix at native
   entry. Cover both cost and gradient routes.
   **Regression:** Fix the starting coordinates and compare the first dense step for MDS and t-SNE
   at `k = 1` and `k = 3` with the R reference, including the returned shape. Keep a 2D accelerated
   control. Barnes–Hut must give a clear dimension error at 1 and 3. Repeat the original 1D MDS
   call under Memcheck and require zero invalid reads or writes attributable to the package.
3. **Preserve zero self-weight on calibration failure.** Exclude the query by index before selecting
   fallback neighbors. Define a consistent tie policy. Handle an all-zero distance row before
   estimating beta, so its fallback weights remain finite.
   **Regression:** For `diag(24)`, perplexity 4, and no input scaling, inspect the dense fallback
   before symmetrization. Each row must put weight on four other observations, with zero diagonal
   and positive row mass. Exercise duplicate observations with all-zero neighbor distances in the
   dense and kNN paths. Require finite normalized probabilities, no self-edge, and unit row sums
   where row normalization is requested.
4. **Make failed-row counts race-free and observable.** Accumulate worker-local counts and reduce
   them after joining. Keep the native calibration result long enough to return available
   diagnostics and issue one warning per calibration call.
   **Regression:** Force failures on 24 equidistant rows in each native calibrator, first with one
   worker and then with two. Both must report exactly 24 failures and produce equivalent weights.
   Check the single public warning and its count. Run the two-worker cases under Helgrind. Exact
   counts alone can pass in the presence of a race. On a separate tie-free fixture, requesting
   dense native `ret_extra = "beta"` must yield one finite precision per row whose weights attain
   the target perplexity within calibration tolerance.
5. **Make worker startup and cleanup safe.** Consolidate the duplicate launch loops while preserving
   serial mode and per-worker reductions. Cap workers at available work, split rows fairly, join
   successfully started threads if a later launch fails, and propagate worker exceptions after
   joining. Keep R API calls on the main thread.
   **Regression:** Cover zero-work, serial, and more-workers-than-rows cases in a focused launcher
   test. Record row visits and require exactly one visit per row. In a small standalone native
   test, force a failure after the first worker starts and a separate exception inside a worker.
   Both must reach the caller after joining. Keep fault injection out of the user-facing API.

## Phase 1 implementation and validation

All five Phase 1 items are implemented. Gaussian targets are checked in R and at native selection
boundaries. Dense MDS and t-SNE use general R expressions outside the 2D specializations.
Barnes–Hut checks its required dimension before calibration and in both native entries.

Native fallbacks exclude self, break distance ties by observation index, and handle zero-distance
rows with finite weights. The nearest-neighbor caller also removes self by index. A public
all-duplicate input showed that self occupied the first slot in only one of 24 returned neighbor
rows, which made positional removal insufficient. Each native calibration now reduces worker-local
failure counts and retains its available diagnostics.

One shared launcher caps workers at available rows, balances chunks, joins started workers on
startup failure, and propagates worker errors after joining. Its replaceable launcher is an internal
test seam for the requested partial-startup regression.

The regressions are in [perplexity tests](../smallvis/tests/testthat/test_native_perplexity.R),
[dimension tests](../smallvis/tests/testthat/test_native_dimensions.R), and the
[standalone worker test](../smallvis/tests/native/threads.cpp). The original 12-row crash was
reproduced under Memcheck against the pre-fix installed package before adding the first fix.

Validation used R 4.5.2, GCC 15.2, Clang 21.1.8, and Valgrind 3.26.0 on Linux. The base
revision is `3382ecf` with the uncommitted Phase 1 changes described above.

- The installed package suite passed 660 expectations. Its seven partial-argument warnings are
  the same seven recorded in the original review. The 81 new expectations had no warnings,
  failures, or skips.
- Memcheck over the installed native regressions and serial, two-worker, and excess-worker
  calibration probes reported zero errors and zero definitely, indirectly, or possibly lost bytes.
- Helgrind over the same installed probes reported zero errors.
- The standalone worker regression passed GCC and Clang builds with
  `-std=c++17 -pthread -Wall -Wextra -Werror`. Memcheck reported zero errors and all allocations
  freed. Helgrind also reported zero errors, including startup and worker exception cases.
- Both compilers accepted the three handwritten C++ translation units with
  `-Wall -Wextra -Werror -fsyntax-only`. Air 0.12.0, clang-format, lintr 3.4.0, and
  `git diff --check` passed.
- A clean source build and `R CMD check --no-manual` passed. The checks used installed dependencies.
  Windows checks were not run locally.

The temporary evidence directory is `/tmp/smallvis-phase1.8ia8acQ9`. Relevant logs are
`suite.log`, `final.memcheck`, `final.helgrind`, `threads.memcheck`, `threads.helgrind`,
and `check/smallvis.Rcheck/00check.log`. The commands below allow the checks to be repeated
after those temporary files have expired.

The candidate was built and installed outside the checkout, and the applied files were compared
byte for byte with that source. The original t-ASNE changes and tests were retained. Exported C++
signatures did not change, so the Rcpp bridges required no regeneration. Package help was
regenerated from roxygen. The full review remains the original audit.

### Commands for repeating the native checks

Run these from the repository root. The installed test lane checks the installed candidate.

```sh
mkdir -p /tmp/smallvis-phase1-lib
R CMD INSTALL --preclean --install-tests --library=/tmp/smallvis-phase1-lib smallvis
R_LIBS=/tmp/smallvis-phase1-lib OPENBLAS_NUM_THREADS=1 \
  R -d "valgrind --vgdb=no --tool=memcheck --error-exitcode=99 --leak-check=full" \
  --vanilla --slave -e \
  'testthat::test_package("smallvis", filter="native_", stop_on_failure=TRUE)'
R_LIBS=/tmp/smallvis-phase1-lib OPENBLAS_NUM_THREADS=1 \
  R -d "valgrind --vgdb=no --tool=helgrind --error-exitcode=99" \
  --vanilla --slave -e \
  'testthat::test_package("smallvis", filter="native_", stop_on_failure=TRUE)'
c++ -std=c++17 -pthread -Wall -Wextra -Werror smallvis/tests/native/threads.cpp \
  -o /tmp/smallvis-threads
valgrind --vgdb=no --tool=memcheck --error-exitcode=99 --leak-check=full \
  /tmp/smallvis-threads
valgrind --vgdb=no --tool=helgrind --error-exitcode=99 /tmp/smallvis-threads
```

## Phase 2. Restore Barnes–Hut cost and force correctness

Use [the tree](../smallvis/src/bh.h), [its traversal](../smallvis/src/bh.cpp), and the
[R objective](../smallvis/R/sne.R). Use one independent dense reference across these tests.
Fix a symmetric, zero-diagonal P of total mass one, set `theta = 0` and disable exaggeration.
From Y compute off-diagonal Student weights, normalize them to Q, and evaluate both KL and the
exact four-times t-SNE gradient. Reuse the same P for both calculations. This isolates tree errors
from neighbor-search differences.

1. **Exclude only the query from a coincident leaf.** For an all-coincident 2D layout on N points,
   calculate expected KL using uniform Q over the N(N−1) ordered non-self pairs. The reported KL
   must match that sum and the force must be zero. A force-only test misses the broken normalizer.
   The review's 24-point case returned −41.6190 versus the direct 1.5359.
2. **Preserve multiplicity when a leaf splits.** Insert several identical points first, then a
   distinct point that forces a split. Assert dense-reference cost and gradient, repeat with the
   rows permuted, and undo the permutation before comparing gradients. Use nonuniform P so the
   force comparison can distinguish missing mass.
3. **Descend past an internal centroid.** Insert the four distinct rows `(-1,0)`, `(1,0)`,
   `(0,0)`, and `(3,1)` in that order. The third row equals the current internal centroid.
   Require exact-reference KL and gradient for this order and another insertion order. Keep a
   generic distinct-point control.
4. **Separate coordinate identity from numerical stabilization.** Use the same documented
   stabilization convention in cost and gradient. On distinct initial coordinates at scale
   `1e-4`, run with machine epsilon and with `eps = 1e-3`. Compare to the dense reference under
   that convention. In particular, subtract the attractive component and compare the repulsive
   component directly. A loose absolute tolerance on the total force can hide its disappearance.
5. **Reject non-finite coordinates before tree construction.** Supply a correctly shaped
   user initialization containing one NaN, then one Inf. Check that the public call returns an
   informative error in a subprocess with a timeout. Exercise both native tree entry points so an
   invalid layout produced during optimization also fails safely.

Close the native phases with the committed public regressions, the focused Memcheck and Helgrind
runs, and the full package suite. Record package path, compiler, and checker results when the fixes
land. Keep exact-force assertions at `theta = 0`. Existing approximate-mode examples should still
produce finite costs and embeddings. Change a historical numerical fixture only after the direct
calculation explains the difference.

## Phase 2 implementation and validation

All five Phase 2 items are implemented against base revision `3ae8e4f`. The tree tracks which
subtree contains the query. A coincident leaf contributes its full multiplicity minus that one
query, and an internal node containing the query is expanded even at a large theta. Splitting
transfers the complete count into the appropriate child. Points at an internal centroid also
descend into a child. Child nodes now use owning smart pointers.

Only identical coordinates share a leaf. The Student weight is `1 / (1 + squared_distance)`,
and `eps` is added to the sum of non-self weights in both the cost and gradient. The R cost
callback now forwards the configured epsilon. Coordinate validation runs after initialization in
R and during native tree construction, covering layouts produced later by optimization.

The [Barnes–Hut regressions](../smallvis/tests/testthat/test_native_bh.R) use a fixed nonuniform
symmetric P and independent dense cost, force, and repulsion calculations. They cover coincident
points, duplicate-leaf splits, internal-centroid insertion, point permutations, one and two
workers, and epsilon at machine precision and `1e-3`. A two-point case at large theta checks
self-exclusion when summarizing would otherwise include the query. No historical coordinate
fixture needed updating.

The NA, NaN, and positive/negative Inf probes run through the public API and both native entries
in a subprocess with a 30-second timeout. The child verifies that it loaded the candidate package.
The test dependencies `callr` and `pkgload` were added to Suggests to support this in source and
installed test runs. The earlier dimension regression also now passes its resolved method list,
so its accelerated control exercises the selected backend.

The [standalone tree regression](../smallvis/tests/native/bh.cpp) checks subtree counts after
duplicate and centroid insertion, plus invalid-coordinate exceptions. The pre-fix numerical
probes disagreed with the dense reference and the invalid native probe reached a C stack overflow.
The finite-coordinate regression is now bounded to keep that failure contained if it recurs.

Validation used R 4.5.2, GCC 15.2, Clang 21.1.8, and Valgrind 3.26.0 on Linux.

- The installed package suite passed 710 expectations, including 53 new Barnes–Hut expectations.
  The seven partial-argument warnings are the existing baseline warnings.
- Memcheck over all 134 native regression expectations and direct native error probes reported
  zero errors and zero definitely, indirectly, or possibly lost bytes.
- Helgrind over the same native checks reported zero errors.
- Strict GCC and Clang builds passed for all three handwritten C++ translation units and for
  the standalone tree regression.
- The standalone tree regression passed AddressSanitizer, UndefinedBehaviorSanitizer, and
  LeakSanitizer. This run used the host environment because LeakSanitizer rejects the sandbox's
  ptrace environment. Memcheck separately reported all 31 allocations freed and zero errors.
- A clean source build and `R CMD check --no-manual` passed. Air, clang-format, lintr, and
  `git diff --check` passed. Windows validation remains outstanding.

The temporary evidence directory is `/tmp/smallvis-phase2.qZssJcbi`. It contains `suite.log`,
`native.memcheck`, `native.helgrind`, `tree.memcheck`, `tree.helgrind`, and
`check/smallvis.Rcheck/00check.log`. The Phase 1 native commands above also run the new R tests.
To repeat the standalone tree checks from the repository root, use the following commands.

```sh
c++ -std=c++17 -pthread -Wall -Wextra -Werror smallvis/tests/native/bh.cpp \
  -o /tmp/smallvis-bh
valgrind --vgdb=no --tool=memcheck --error-exitcode=99 --leak-check=full /tmp/smallvis-bh
valgrind --vgdb=no --tool=helgrind --error-exitcode=99 /tmp/smallvis-bh
g++ -std=c++17 -pthread -fsanitize=address,undefined -fno-omit-frame-pointer -g \
  smallvis/tests/native/bh.cpp -o /tmp/smallvis-bh-sanitize
ASAN_OPTIONS=detect_leaks=1 /tmp/smallvis-bh-sanitize
```

The applied files were compared byte for byte with the tested candidate. Exported Rcpp signatures
are unchanged. Package help was regenerated from roxygen. The pre-cleanup files remain in stash
`1bc0e8f`, including the t-ASNE changes and earlier experiment plans. Phase 3 completion is
recorded below.

## Phase 3. Make equivalent calls use the same inputs and initialization

- **Dispatch and affinity representation.** Resolve the method name, method list, and top-level
  controls once. Preserve custom method functions and existing method-specific defaults. Shared
  public defaults apply to both built-in method forms, explicit list options override shared
  controls, and other method options retain constructor defaults.
  Check that string/list t-SNE with the same controls produce identical P, kNN differs from the
  Gaussian control, and requesting BH selects a sparse route in either form. Handle or reject
  explicit dense Gaussian input for BH before it reaches sparse-slot access.
- **Coordinates, distances, and multiscale inputs.** On tie-free X, compare P from X and
  `dist(X)` with scaling disabled. Cover Gaussian distances, the distance-input fallback when
  `use_cpp = TRUE`, and a directed four-neighbor graph before symmetrization. Independently check
  each row's four indices and zero diagonal. A one-scale multiscale request must reproduce a
  single scale for both globally normalized t-SNE and row-normalized ASNE.
- **Scaling and PCA.** Compare `scale = TRUE` on a matrix and its data-frame equivalent to
  columnwise standardization after dropping a constant column. Check `"absmax"` on −9 and nine
  values of +1 and define the all-constant case. Test skipped PCA on a 24 × 4 matrix requesting
  50 PCs, actual PCA initialization at `k = 3`, and whitening from 14 features to three. Check
  output dimensions, PCA metadata, the retained subspace, and whitened variances.
- **Spectral initialization.** On a ten-node cycle, both eigensolver paths must span the
  two-dimensional eigenspace at +cos(2π/10) and exclude the stationary mode. Compare subspaces
  to allow signs and rotations. With X and seed held fixed, requesting V as extra output must
  leave pairwise distances in the starting layout unchanged.

## Phase 3 implementation and validation

All four Phase 3 groups are implemented against base revision `8134a12`. No Phase 4 optimizer,
repeated-start, stepping, or diagnostic-formula changes are included.

Built-in string and list dispatch now use one constructor path with the same shared public defaults.
Explicit non-NULL list options override shared controls. NULL kernel and neighbor controls preserve
constructor defaults, as do method-specific options absent from the public interface. Custom
constructors retain their defaults and receive declared, explicitly supplied shared controls, plus
the existing perplexity/neighbor-count argument. Invalid list options fail before input processing.
Resolved perplexity and native controls also reach result metadata and optional distance output.
Both t-SNE forms select Barnes–Hut when requested, and its initializer converts dense Gaussian
weights to a general sparse matrix before native sparse-slot access.

Coordinate and precomputed-distance Gaussian affinities now use the same squared distances.
Precomputed distances, vector targets, and exponential kernels use the supported R calibrator
when the dense native calibrator cannot handle them. Exact distance neighbor lookup supports dense
and sparse output, removes self by index, and selects exactly k directed neighbors, including k = 1
and tied distances. Nearest-neighbor Gaussian calibration accepts distance inputs and fractional
perplexities. Multiscale calibration preserves each method's normalization and symmetrization.

Matrix and data-frame scaling now operate on columns and retain a one-column matrix after numeric
column selection. Abs-max scaling uses the largest absolute centered entry. Constant input raises
a clear error for range, abs-max, and variance scaling because the requested normalization is
undefined. PCA initialization reuses preprocessing scores only when preprocessing actually ran,
respects the requested output dimension, and preserves original feature counts in metadata.
Truncated whitening retains its eigenvalues, scalar full-SVD scores retain matrix shape, and
unavailable PCA initialization dimensions fail clearly instead of returning invalid coordinates.

Laplacian eigenmaps solve the symmetric normalized adjacency problem, request the largest algebraic
eigenvalues, remove the stationary mode, and map back to the generalized eigenvectors. Tests cover
the ten-cycle Fourier subspace and a graph with unequal degrees. Both spectral initializers handle
small inputs with a dense solve. Graph selection and symmetrization are independent of requested
optional output, so asking for V no longer changes starting coordinates.

The regressions are in [dispatch tests](../smallvis/tests/testthat/test_dispatch.R),
[input equivalence tests](../smallvis/tests/testthat/test_input_equivalence.R),
[preprocessing tests](../smallvis/tests/testthat/test_preprocessing.R), and
[spectral contracts](../smallvis/tests/testthat/test_spectral_contracts.R). The existing API helper
now leaves the input kernel unset unless a test selects it, preserving each method's own default.
The existing Laplacian coordinate fixture was replaced by eigen-equation and weighted-orthogonality
assertions. On its graph, the old fixture's eigenvalues were -0.3935862 and -0.1759269, whereas the
leading nonstationary algebraic eigenvalues are -0.02654955 and -0.03492391. This independently
explains the result change. The unchanged normalized-Laplacian fixture still passes.

Validation used R 4.5.2 and GCC 15.2 on Linux.

- The installed full suite passed 802 expectations with no failures, errors, or skips. The seven
  partial-argument warnings are the previously recorded baseline warnings.
- Running the new regression groups against archived HEAD reproduced failures in dispatch, input
  equivalence, scaling/PCA, and spectral initialization.
- Coordinate and distance neighbor indices agree on the tie-free fixture. Their maximum distance
  difference is 2.45e-7 because the coordinate neighbor backend uses native precision. The
  nearest-neighbor Gaussian comparisons therefore use tolerance 1e-6, with separate index and
  distance checks. Dense Gaussian and exponential equivalence use tolerance 1e-10.
- A clean source build and `R CMD check --no-manual` passed with status OK. Air, lintr, and
  `git diff --check` passed. The check used installed dependencies. Windows checks were not run.
- Roxygen help was regenerated. The usage article and method help now show the supported
  `method = "tsne", bh = TRUE` selector. No C++ source or generated Rcpp bridge changed.

The temporary evidence directory is `/tmp/smallvis-phase3.PL3ldwET`. It contains `suite.log`,
`before-tests.log`, `before-init-tests.log`, `before-spectral-tests.log`, `proof.log`,
`format.log`, `lint.log`, and `check/smallvis.Rcheck/00check.log`. To repeat the package checks,
run the full test suite from `smallvis/`, then build that directory and check the source tarball
with `R CMD check --no-manual` as specified below. The applied files were compared byte for byte
with the tested source. The pre-cleanup stash remains untouched.

For Phase 4, perplexity stepping must account for the now-consistent precedence of an explicit
`method = list(..., perplexity = ...)` over the top-level perplexity. Preliminary and final runs
must resolve the intended target together rather than silently retaining a nested override.

## Phase 4. Make optimizer results, repeated starts, and diagnostics trustworthy

- **Final scoring and repeated starts.** Compute the final objective at the returned coordinates,
  including a stop before the first epoch and a stop after an earlier epoch. Preserve the
  separate epoch history. For built-in repeated starts, prepare the input affinities once and
  reset coordinates and optimizer state for each start. Verify that neighbor search runs once,
  returned P is identical across approximate-neighbor repeats, and independently rescoring each
  returned Y against that P selects the reported best run. Cover `keep_all`, matrix returns,
  and a first-step L-BFGS termination.
- **Perplexity stepping.** Use separate preliminary and final arguments. On 24 observations with
  target 10, the empty step schedule must match a direct run with the same seed. On a nonempty
  schedule, check the final call against a direct call from the last preliminary layout, including
  default and explicitly supplied epoch and tolerance settings.
- **Optimizer arithmetic.** Test two AdaMax updates against an elementwise recurrence with
  unequal gradients in two dimensions and one zero gradient. For UMAP and LargeVis, compare the
  scalar's directional derivative with the supplied gradient when using a function-based
  optimizer. Retain stabilized forces for gradient-only optimization. Until a matching objective
  is supplied, reject combinations with inconsistent derivatives. Keep the existing guard on
  staged mize exaggeration.
- **Exaggeration and Gaussian output weights.** For t-EE, compare the attractive force against
  its defining formula at factors 1 and 12, and check that removing exaggeration restores the
  original force. If the intended method has no exaggerated form, reject the option explicitly.
  For ASNE at coordinates spaced 100 units apart, compare the row probabilities and finite KL
  with independently shifted off-diagonal log weights.
- **Diagnostic meanings.** Use an asymmetric weight matrix to distinguish incoming and outgoing
  degrees and to check their average. Compare squared-MDS DX and DY with ordinary Euclidean
  distances. Check each ASNE observation's KL against its own row of P and Q as well as the total.
  Document the current scalar BH cost alongside the per-observation results of other methods.

## Phase 5. Reduce allocations and document supported behavior

- Measure dense t-SNE peak memory and iteration time in fresh processes on the same fixed input,
  compiler, and thread count. Remove redundant P/W copies and dense calibration transposes while
  preserving output values and input ownership. At 5,000 points those P/W copies account for
  400 MB of payload per gradient evaluation. Record measured peak memory separately. Compare
  one- and two-thread results with the dense reference before claiming an improvement.
- Recheck consumers of the private squared-distance-to-Student-weight conversion. If still
  unused, remove it with its R wrapper and regenerate Rcpp bindings. Verify clean installation
  and absence of stale registrations or references. If a supported caller emerges, fix the
  conversion and test every entry of a 3 × 3 matrix against `1 / (1 + D2)`. Keep the dynamically
  selected optimizers.
- Add a capability table to the existing usage article. Its BH selector example was corrected
  during Phase 3.
  Cover input forms, output dimensions, kernel/backend combinations, and diagnostics exercised
  above. Point the README directly to the applicable LargeVis and UMAP terms and attribution
  after verifying their provenance.

## Completion and deferred decisions

Run focused testthat cases during each fix. From the package directory run the full suite with
`Rscript -e 'testthat::test_local()'`, build the package, and run `R CMD check --no-manual` on the
tarball. Use the existing Linux and Windows checks, formatting, and lint. Regenerate bindings only
when their declarations change. At each phase boundary record the revision, completed items, test
results, and next open item here. Close the plan by reconciling it with the saved full review.

The following are separate choices after the fixes:

- Consider an effective-backend field in `ret_extra` if users need a machine-readable record
  beyond the capability table. Consider a public precomputed-affinity interface after the
  internal reuse in repeated starts has demonstrated its value.
- Move `vizier` to an optional plotting dependency if installation friction warrants a fallback.
  Add a recurring native diagnostic CI job only if the focused commands justify its upkeep.
- Keep the staged mize proposal, `smallvis-mize-staged-optimization-plan.md` in the pre-cleanup
  stash, independent. Defer a thread pool, octree, plugin system, persistent experiment store,
  distributed execution, and broad benchmark infrastructure until a concrete workload calls for
  them.

**Progress on 2026-11-05** Phases 1 through 3, their regressions, native validation for the
first two phases, and Phase 3 package checks are complete. Phases 4 and 5 are pending.
