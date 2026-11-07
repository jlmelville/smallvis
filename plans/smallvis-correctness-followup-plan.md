# smallvis correctness, cleanup, and allocation follow-up

**Status: phases 1–5 complete; phases 6–7 pending.** The user requested that the
[second repository review](smallvis-repository-review-2026-11-06.md) be saved and converted into
a phased plan on 2026-11-06, then explicitly authorized phases 1, 2, and 3 in sequence. Source and
test changes for phases 1–3 are committed as `a81ef3e`. The user requested implementation of phase 4
with commits after completed changes. USNE/pseudo-separable optimizer contracts are committed as
`817a608`; the dynamic HSSNE fix is committed as `1a5bed0`. The user has authorized phase 5.
The backend-forwarding correction is committed as `e9fc703`; the unused-code cleanup is validated
and recorded below. Next: await direction on phase 6, native allocations and measurements.

**Baseline:** `e65bcbaf6b71669c92efe75cfdd595739e894a31`, with a clean worktree before saving these
documents. The earlier [correctness plan](smallvis-correctness-implementation-plan.md) is complete
and remains the historical record of the first cleanup. This plan addresses the newly saved review.

## Purpose and boundaries

Make the existing embedding methods safe, numerically trustworthy, and consistent across their
supported inputs and controls. Reduce the demonstrated allocation overhead after establishing
correctness. Preserve the three public entry points, readable R equations, method-specific cost
lists, current licensing, and the project's experimental and educational character.

The plan includes the eight correctness findings, backend-control cleanup, confirmed unused code,
the misleading test loop, focused mathematical regressions, native test automation, and the MDS
gradient and Student-weight allocation improvements. Update existing help or articles only where
the corrected behavior requires it. Keep scientific corrections separate from performance-only
changes so changed numerical results have an identifiable cause.

New result fields, stopping-reason reporting, public affinity or neighbor inputs, new methods,
expanded Barnes-Hut dimensions, a thread pool, an object-system rewrite, and experiment or
benchmark services are outside this plan. The previous capability-reference proposal and staged
mize optimization proposal remain separate. There is no CRAN-release work, broad documentation
rewrite, or rerun of historical experiments in scope.

These two documents live in the existing repository-level `plans/` directory, outside the
`smallvis/` package root. The user subsequently authorized committing completed changes, including
phases 1–3 and each subsequent authorized change. Preserve unrelated work. Pushing and publishing
remain outside scope.

## Review disposition

The correctness observations were reproduced during the review against the baseline above.
Cleanup observations came from source and consumer inspection. Allocation sizes are calculated
from live buffers; performance gains remain unmeasured. Phase completion requires fresh evidence
against the implementation revision.

| Review item | Disposition | Completion evidence |
| --- | --- | --- |
| 1. Barnes-Hut nontermination | Phase 1 | Bounded public probe, exact reference, native sanitizers |
| 2. AB-SNE scalar objectives | Phase 2 | Independent scalar formulas, zero identity, stopping regression |
| 3. EE and g-SNE input equivalence | Phase 3 | Equal constructed quantities, objectives, and gradients |
| 4. Gaussian calibration boundaries | Phase 1 | Finite non-self fallback and warnings in R and native paths |
| 5. USNE and pseudo-separable line searches | Phase 4 | Matched derivative or early rejection; finite USNE diagonal |
| 6. Dynamic HSSNE transitions | Phase 4 | Valid candidates and fresh gradient after an accepted change |
| 7. Multiscale combinations | Phase 3 | Supported results or deliberate early rejection |
| 8. Self-neighbors and geodesic k = 1 | Phase 3 | Identity-based self exclusion and preserved dimensions |
| Backend-control forwarding | Phase 5 | Public routing checks at ordinary and special parameter values |
| Shadowed definition, unused helpers, argument typos | Phase 5 | Confirmed consumers, focused diff, preserved public behavior |
| Native-warning test loop | Phase 1 | Each fixture reaches the public call and has its own assertions |
| MDS and Student-weight copies | Phase 6 | Numerical equivalence, input ownership, comparative memory evidence |
| Mathematical coverage across affected families | Phases 1–5 | Independent invariants attached to each correction |
| Native harness CI | Phase 7 | Existing harnesses and new boundary case execute in a bounded job |
| Retain the current architecture | All phases | Helpers serve demonstrated fixes; public surface stays stable |
| Effective backend and stopping-reason fields | Excluded by user | Revisit only on a separate feature request |
| Public affinities or neighbors | Excluded by user | Revisit only on a separate feature request |
| Larger infrastructure and additional dimensions | Outside scope | Separate objective and evidence would be required |

## Context and evidence to retain

The package lives in `smallvis/`. Public execution and optimizer compatibility checks are in
`R/smallvis.R`; calibration is in `R/perplexity.R` and `src/perplexity.cpp`; the tree is in
`src/bh.h`; method formulas are in their named R files. Existing tests are in `tests/testthat/`
and standalone native harnesses in `tests/native/`. Allocation tooling is under the repository's
`benchmarks/` directory. Paths in phase descriptions are relative to the package unless stated
otherwise.

The review used R 4.5.2, testthat 3.3.2, GCC 15.2, and Clang 21.1 on Linux. A clean source
`R CMD check --no-manual` passed with 968 expectations and no test failures, warnings, or skips.
Strict compiler checks and the two existing native harnesses also passed. An additional adjacent-
coordinate tree probe failed under ASan. Those results establish the baseline and its coverage
gap; they do not establish that any planned fix works. Windows, a full-package sanitizer run,
and leak detection were not exercised by that review.

For the scalar, input-equivalence, and directional-derivative reproductions, use this fixed setup:

```r
set.seed(42)
X <- matrix(rnorm(24 * 4), 24)
Y <- matrix(rnorm(48), 24)
```

Use `scale = FALSE`, `perplexity = 4`, `Y_init = Y`, `epoch_callback = FALSE`, `verbose = FALSE`,
and `n_threads = 1` unless a phase specifies otherwise. Compare at fixed coordinates before
optimization. For scalar output request `ret_extra = TRUE`; add matrices only where an independent
reference needs them. Keep new executable probes in temporary storage until implementation is
authorized, then put the useful regressions in the existing suite.

Temporary review logs and scripts were under `/tmp/smallvis-review-20261106.bmx9dw` and
`/tmp/smallvis-review-*.R`. They are expendable supporting evidence, not required inputs to this
plan. Preserve meaningful results in this document as phases finish; store raw outputs outside Git.

## Phase 1. Bound native tree construction and Gaussian calibration

Work in `src/bh.h`, its traversal callers in `src/bh.cpp`, `src/perplexity.cpp`, and
`R/perplexity.R`. Extend `test_native_bh.R`, `test_native_perplexity.R`, and `tests/native/bh.cpp`.
Keep tree and calibration corrections separately reviewable within the phase.

1. Reproduce the tree failure in a timed subprocess before editing. Use the native two-point
   input `(1 + eps, 0), (1 + 2 * eps, 0)`, where `eps` is double machine epsilon. For the public
   route, use `X = cbind(1:24, sin(1:24))` and an initializer with 23 copies of the first point
   and one copy of the second; set `bh = TRUE`, `nn = "exact"`, `max_iter = 0`, and request costs.
   The baseline public process exceeded a 15-second timeout; the native probe overflowed its stack.
2. Detect when subdivision cannot separate distinct representable coordinates. Choose the smallest
   terminal-leaf treatment that preserves point multiplicity, excludes exactly one query, and
   permits exact pair contributions at `theta = 0`. Preserve ordinary approximate traversal and
   exact duplicate behavior. A broad coordinate-merging epsilon would alter the corrected contract.
3. Make Gaussian calibration distinguish a valid entropy result from zero or non-finite total
   weight. Handle zero mean distance in the R initialization. Preserve the existing finite
   fallback intent and warning ownership for each dense and neighbor route; their graph supports
   need not become identical. Stable weight evaluation may be used where needed.
4. Exercise `diag(24)` with perplexity 1 and 4, and all-zero `24 x 4` input with scaling disabled.
   The baseline gives NaN affinities at perplexity 1 in both dense backends and a condition error
   for the all-zero R route. Include successful calibration on ordinary data as a control.
5. Fix the native warning test to pass the loop's `X` into `smallvis()`. Assert positive row mass,
   zero self-weight, finite normalized P, and one correctly counted warning for fallback cases.

**Acceptance:** timed public cases finish; tree costs and gradients agree with dense sums at
`theta = 0`, including insertion-order variants and one/two threads. Existing duplicate and
centroid tests remain valid. Calibration yields valid affinities or an intentional, documented
domain error, never silent NaNs or an R condition error. Run the native tree harness under
ASan/UBSan and obtain active leak evidence for the changed ownership path, using a supported
leak-sanitizer run or Valgrind if needed. Recheck threaded calibration with the existing focused
memory/race checks if its worker state changes. Record any unavailable lane explicitly.

### Phase 1 outcome and evidence, 2026-11-06

Implemented against the baseline above, without changes to native signatures or generated bridges:

- Tree subdivision now requires representable progress along a differing coordinate. A terminal
  leaf retains exact coordinate/count groups when subdivision stalls or its extent overflows.
  Cost and gradient traversal evaluate those groups directly, excluding one query observation.
  The bucket is allocated only for these exceptional leaves; ordinary nodes carry one extra pointer.
- The R calibrator initializes zero mean distance with beta zero. R and native calibrators reject
  zero or non-finite total weight as an invalid entropy evaluation and use their existing fallback.
  Dense fallbacks and neighbor fallbacks retain their distinct support policies.
- The public warning regression now uses each actual input fixture, including coincident input,
  and covers perplexities 1 and 4 through dense R, dense native, and native neighbor calibration.
  It checks warnings, positive row mass, finite normalized P, raw support, and zero self-weight.
  Successful dense R and native calibration still reconstructs Gaussian weights and target entropy.
- A timed subprocess covers the original public tree failure plus native costs and gradients under
  three insertion orders and one/two threads. Forces are compared relative to their tiny magnitude,
  which distinguishes exact pair contributions from accidentally merging the coordinates.
  The standalone harness checks bucket multiplicity, overflowing extent, and exception unwinding.

Validation of the final phase 1 source on Linux, R 4.5.2, Rcpp 1.1.1.1.1, testthat 3.3.2:

| Check | Result |
| --- | --- |
| `Rscript --vanilla -e 'testthat::test_local()'` from `smallvis/` | 1073 passing expectations; no failures, warnings, or skips |
| Fresh source `R CMD build`, then `R CMD check --no-manual` | Status OK; no counted errors, warnings, or notes |
| `air format smallvis --check`, Air 0.12.0 | Passed |
| Loaded-namespace `lintr::lint_package("smallvis")` from the repository root | No lints |
| Strict GCC 15.2 and Clang 21.1 C++17 syntax checks, all three handwritten `.cpp` files | Passed with `-Wall -Wextra -Wpedantic -Werror` |
| Tree harness under both compilers with ASan/UBSan and a 30-second timeout | Passed; leak detection disabled in these runs |
| Tree harness under Valgrind 3.26 with full leak checks and error exit status | Zero errors; all 49 allocations freed; zero bytes retained |

The initial full-suite run inherited `warnPartialMatchArgs = TRUE` from the local R profile and
reported seven existing partial-argument warnings in unchanged tests. `Rscript --vanilla` confirms
the default is FALSE and the full suite is clean without that profile. These test cleanups remain
outside phase 1. LeakSanitizer itself failed because ptrace is active in this sandbox; Valgrind
provides the active leak evidence for the same final tree harness, including exception unwinding.
Calibration changes leave worker buffers and shared bookkeeping unchanged, so no new race-check
lane was required. Windows and a full-package sanitizer run were not exercised.
The package check could not refresh the CRAN dependency index in this restricted environment;
it checked against installed dependencies and still completed with Status OK.

Temporary evidence is under `/tmp/smallvis-phase1-check-31x4u_0c` and
`/tmp/smallvis-phase1-{tests-vanilla,native,calibration,lint}.log`. These are supporting artifacts,
not required inputs. The durable regressions are in the existing package tests.

## Phase 2. Correct AB-SNE scalar objectives and stopping behavior

Work in `R/absne.R`, with focused scalar regressions alongside `test_gradients.R` and
`test_final_results.R`. Use the formulas already documented in `vignettes/articles/absne.Rmd`.
The existing AB-SNE coverage checks finite-difference gradients across ordinary and special branches
and one rounded public trajectory for ordinary parameters. Missing invariants are independently
computed absolute scalar values, zero divergence at equal distributions, and stopping governed by
the corrected scalar. Phase 2 adds those checks without changing gradient formulas.

1. Correct `absneb0()`'s cached input term to use `P^alpha * log(P^alpha)`. Correct
   `absneamb()`'s constant for the number of non-self pairs. Preserve the existing gradients.
2. Evaluate off-diagonal scalar formulas independently from production helpers for representative
   ordinary and special branches, including `alpha = lambda = 2` and `alpha = 1, lambda = 0`.
   Choose finite, well-conditioned probabilities so epsilon clipping does not obscure the identity.
3. Reproduce the public zero-divergence witness with `X = diag(8)`, `k = 8`, `Y_init = diag(8)`,
   perplexity 7, scaling disabled, `method = list("absne", alpha = 1, lambda = 0)`, and zero
   iterations. P and Q match to roundoff; the baseline cost is -8.
4. Check the ordinary fixed fixture with `method = list("absne", alpha = 2, lambda = 2)`,
   `opt = list("adam", eta = 0.001)`, `max_iter = 20`, `epoch = 1`, and `tol = 0`. The baseline
   stops after one iteration with a negative cost. Verify the corrected stopping decision using
   the independent scalar; do not replace this with a new unexplained coordinate snapshot.

**Acceptance:** scalar references agree, equal distributions have zero divergence within numerical
tolerance, and nonnegative branches do not trigger premature stopping from a spurious constant.
Existing finite-difference gradients still pass. Treat altered costs and stopping trajectories as
intended corrections and explain them where current user documentation needs an update.

### Phase 2 outcome and evidence, 2026-11-06

Corrected two scalar expressions in `R/absne.R`: `absneb0()` now weights its cached input log term
by `P^alpha`, and `absneamb()` subtracts `N - 1` non-self terms per column. Gradient expressions
and optimizer stopping logic are unchanged. Corrected costs can change relative-tolerance and
minimum-cost stopping decisions, as intended.

`test_absne_objectives.R` checks ten parameter combinations spanning the ordinary formula and all
four special forms against an independent off-diagonal reference in `helper_objectives.R`.
The scalar fixture explicitly stays above clipping thresholds. Four configurations of the two
corrected branches, including negative alpha, also give zero column contributions when P equals Q
on the eight-point identity-matrix witness. These tests failed on the old scalar expressions.

The new regression in `test_final_results.R` reproduces the fixed 24-point public run. Before the
fix it stopped at iteration 1 with a negative cost. It now reaches all 20 requested Adam iterations
with the default `min_cost = 0`, `epoch = 1`, and `tol = 0`. Every recorded epoch agrees with the
independent scalar, all epoch costs are finite and positive, and the final score decreases and
matches the returned coordinates. No coordinate snapshot or disabled stopping check is needed.

The AB-SNE article now explains that both sums exclude self pairs, removing the incorrect
unresolved claim about different diagonal counts. It also explains why input-only terms matter
for reported costs and cost-based stopping even though they leave gradients unchanged. Its
historical experiments and static example remain unchanged.

Validation on the same Linux/R toolchain as phase 1:

| Check | Result |
| --- | --- |
| Focused `test_local(filter = "absne_objectives\|gradients\|final_results")` in a clean R session | 459 passing expectations; no failures, warnings, or skips |
| Full `Rscript --vanilla -e 'testthat::test_local()'` from `smallvis/` | 1132 passing expectations; no failures, warnings, or skips |
| Fresh source `R CMD build`, then `R CMD check --no-manual` | Status OK; checked archive contains the exact changed source and new tests |
| Air 0.12.0 and loaded-namespace lintr | Passed; no lints |
| R Markdown formatting helper on `absne.Rmd` | Passed; no code formatting changes |
| Temporary article render and static HTML inspection | Passed; 82 well-formed MathML expressions, 17 display equations, four valid local links |

The package check again used installed dependencies because the CRAN index was unavailable.
The standalone preview includes copies of the four unchanged neighboring article targets for link
inspection. Browser layout inspection was not performed; rendering and static checks do not
establish visual readability. No native source changed in phase 2, so phase 1's native safety
evidence remains applicable without repeating those lanes.

Supporting artifacts are `/tmp/smallvis-phase2-check-tl00heb4`,
`/tmp/smallvis-phase2-{focused,tests,lint,render}.log`, and
`/tmp/smallvis-phase2-article/absne.html`. Package tests retain the executable regressions; these
temporary artifacts are not prerequisites for future phases.

## Phase 3. Restore input, neighbor, and multiscale contracts

Work in `R/ee.R`, `R/gsne.R`, `R/perplexity.R`, `R/umap.R`, and the relevant method initializers.
Extend `test_input_equivalence.R` and the existing smooth-kNN and API tests.

1. Convert EE's compact distances to a full matrix before normalization. Square g-SNE's
   precomputed distances when constructing `Phat`. Compare X and `dist(X)` under the fixed setup
   for the constructed quantities, scalar objective, gradient, and a short deterministic run.
2. Filter self-neighbors by identity in UMAP and geodesic construction, using the established
   neighboring code as the pattern. Preserve dimensions when only one neighbor remains. Keep
   UMAP's current meaning of k, which reserves a slot for self, distinct from geodesic k's
   non-self count; this correction should not silently increase every UMAP neighborhood.
3. Use six identical zero rows followed by `cbind(1:10, (1:10)^2)` for the duplicated-input UMAP
   witness at k = 4. Check non-self counts and diagonal zero for coordinate and distance inputs.
   Tied neighbor identities and approximate searches need not produce identical graphs. Test
   `method = list("gmmds", k = 1)` through the public API with both input representations.
4. Inventory multiscale consumers of P, degree weights, raw affinities, and bandwidths. Derive
   quantities from the averaged P where that preserves the method's defined mathematics, such
   as degree weights for weighted methods. Explicitly reject combinations requiring unavailable
   or undefined quantities before optimization. Do not invent a bandwidth-transfer rule for a
   multiscale mixture. Retain working t-SNE and ASNE multiscale behavior.
5. Cover the observed failures in weighted t-SNE, weighted SSNE, unnormalized t-SNE,
   pseudo-separable t-SNE, and bandwidth-transfer methods. For each, record whether it is
   supported by a verified construction or rejected with a method-specific explanation.

**Acceptance:** EE and g-SNE representation checks agree without recycling warnings; neighbor
counts honor each contract; geodesic k = 1 works; multiscale requests produce meaningful results
or clear errors before a missing-field failure. Preserve single-scale equivalence for supported
combinations, including `list("multiscale", 4)`. Update existing help for intentional restrictions.

### Phase 3 outcome and evidence, 2026-11-06

EE now expands precomputed distances before normalizing its negative weights; g-SNE squares them
when constructing `Phat`. Public regressions compare independent target matrices, initial costs
and gradients, and four deterministic Adam iterations for coordinate and distance input. Callback
inspection exposes the constructed targets without adding public result fields.

UMAP and geodesic graph construction now use the existing neighbor helper and remove self by
identity before truncation. UMAP reserves one self slot; geodesic k counts only non-self neighbors.
Geodesic distance input now selects exactly k neighbors before undirected symmetrization instead
of including all ties at the cutoff. Tied identities across representations remain unspecified.
Avoiding matrix slices also fixes the geodesic k = 1 dimension drop.

UMAP retains a continuous cardinality target while rounding only the search width down. The old
distance route accepted fractional k, whereas the coordinate backend rejected it; both now handle
k = 4.5 consistently. Its supported range is 2 through n. Regressions cover duplicated input,
a controlled neighbor pool with self first, later, or absent, the public range endpoints, and a
connected one-neighbor chain through geodesic MDS. Ordinary UMAP reference values remain unchanged.

The multiscale consumer inventory led to these dispositions:

| Methods or data | Disposition | Evidence |
| --- | --- | --- |
| t-SNE, ASNE | Retain averaged Gaussian P with each method's normalization | Existing and new singleton/mixture comparisons |
| `wtsne`, `wssne` | Derive `pdeg = colSums(P)` and the existing degree-product matrix | Independent degree/product assertions and finite runs |
| `tsneu`, `pstsne` | Derive symmetric conditional weights as `n * P` | Singleton costs/trajectories and mixture weights agree |
| `btsne`, `bssne`, `basne`, `btasne` | Reject automatic bandwidth transfer; accept explicit output beta | Rejection before calibration, plus scalar-beta singleton and mixture runs |
| `bnerv` | Reject multiscale because it requires an automatically transferred bandwidth | Method-specific public error before calibration |
| `bhssne` | Retain its existing fixed output beta construction | Consumer inspection; it overwrites the requested input beta |
| Other P-only methods, including EE and g-SNE | Retain their existing averaged-P construction | No missing auxiliary quantity in consumer inspection |
| Per-scale beta, raw-affinity degrees, intrinsic dimensionality | Remain unavailable for mixtures | No invented bandwidth or diagnostic aggregation |

Automatic-transfer rejection is tested for top-level and nested method perplexities. Supported
affected methods complete short gradient-based runs, agree with ordinary calibration for
`list("multiscale", 4)`, and use the average of the independently calibrated P matrices for
`list("multiscale", c(4, 6))`. Optimizer compatibility changes remain in phase 4.

Validation on Linux with the same R toolchain as earlier phases:

| Check | Result |
| --- | --- |
| Full `Rscript --vanilla -e 'testthat::test_local()'` from `smallvis/` | 1368 passing expectations; no failures, warnings, or skips |
| Fresh source `R CMD build`, then `R CMD check --no-manual` | Status OK; installed suite also passes all 1368 expectations |
| Checked archive compared against the changed R files, help, and new tests | Exact source included, including both untracked test files |
| Air 0.12.0 and loaded-namespace lintr | Passed; no lints |
| Roxygen 8.1.0 generation and a second pass | Only `man/smallvis.Rd` changes; second pass leaves help, DESCRIPTION, and NAMESPACE unchanged |

The function help now describes UMAP's count/range, geodesic neighbor semantics, degree diagnostics,
and the bandwidth-transfer restrictions. Existing guides already link to these reference sections;
no article or site rebuild was needed. The package check again completed using installed
dependencies because the CRAN index was unavailable. Windows was not exercised. No native source
or bridge signature changed in phase 3, so earlier native safety evidence still maps to that code.

Temporary logs are `/tmp/smallvis-phase3-{focused,tests,lint,roxygen,roxygen-idempotence}.log`;
the fresh archive and check are under `/tmp/smallvis-phase3-check-xd5gpp9h`. The old defects were
observed before implementation in the input/neighbor and multiscale regressions; their logs are
supporting evidence, while the maintained tests retain the acceptance witnesses.

## Phase 4. Match optimizer contracts and dynamic HSSNE state

Work in `R/sne.R`, the existing compatibility checks in `R/smallvis.R`, and, only if needed,
`R/cost.R` or `R/epoch.R`. Extend `test_optimizer_contracts.R` and focused epoch regressions.

1. Mark stabilized USNE forces and the current pseudo-separable t-SNE method as unsuitable for
   function-based optimizers using the existing mechanism. Keep their gradient-only experiments
   available. A new surrogate objective for pseudo-separable t-SNE is outside this correction.
2. Correct USNE's diagonal handling at `gr_eps = 0`. Permit function-based optimization for an
   unstabilized configuration only after its finite gradient matches its scalar derivative.
   Account for coincident coordinates only as required by the kernel's actual differentiability;
   do not claim an exact finite derivative where the mathematical limit is undefined.
3. Repeat the directional check using `U = Y / sqrt(sum(Y * Y))` and steps `1e-3`, `1e-4`,
   `1e-5`, and `1e-6`. The saved review gives the baseline discrepancies. Verify deliberate early
   rejection for incompatible public calls and continued operation of gradient-only optimizers.
4. Check dynamic HSSNE's proposed `alpha - 0.1` against the lower bound. When accepting a new
   alpha, invalidate weights and all dependent cached quantities before the next evaluation.
   Inspect any optimizer-held values reused across this objective change; reset or reject an
   incompatible combination through the existing contract if necessary.
5. Test the first gradient after an accepted 0.5-to-0.6 transition against a fresh evaluation.
   For the public lower-bound witness, initialize with `Y * 0.1`, start alpha at 0.1, use Adam
   with eta `1e-5`, six iterations, epoch 1, `tol = 0`, and `min_cost = -Inf`. Observe alpha
   through the existing callback; the baseline admits zero.

**Acceptance:** every newly admitted function-based configuration has a matched derivative;
incompatible combinations fail before expensive work. Dynamic HSSNE never admits the invalid
candidate and never uses weights from a prior alpha. Preserve ordinary optimizer and exaggeration
regressions. Keep any helper extraction limited to the lifecycle responsibility it resolves.

Implementation is split into two commits: force/objective compatibility, then dynamic HSSNE
transitions. Existing USNE coverage uses a tiny positive stabilizer rather than zero, and existing
compatibility coverage omits USNE and pseudo-separable t-SNE. New regressions must distinguish
both missing flags and diagonal NaNs. For the UMAP kernel, coincident-pair forces tend to zero
only when `b > 0.5`; negative `min_dist = -0.5`, currently accepted by curve fitting, produces
`b` near 0.426 and provides a real non-differentiable public configuration to reject at coincidence.

### Phase 4 force/objective compatibility checkpoint, 2026-11-06

USNE now marks stabilized forces as gradient-only, and pseudo-separable t-SNE always uses that
restriction. The common rejection diagnostic recommends an applicable gradient-only optimizer
without suggesting the nonexistent `gr_eps` option for pseudo-separable t-SNE. Rejection occurs
before input calibration. Both methods retain finite Adam runs.

Unstabilized USNE zeros self-force coefficients explicitly. Coincident non-self pairs receive their
zero-force limit only for `b > 0.5`; otherwise the derivative is undefined and the gradient reports
an explanatory error. Independent normalized-kernel scalar checks agree with the finite gradient
at steps `1e-3`, `1e-4`, `1e-5`, and `1e-6`, including a fitted `b < 0.5` away from coincidence.
Differentiable coincident configurations also pass derivative checks and public optimization.

The unchanged default mismatches were repeated at all four steps: USNE force 0.03768988 versus
scalar derivative 0.03326312; pseudo-separable t-SNE force 0.1411800 versus 0.1915116. These remain
intentional gradient-only experiments; no surrogate objective was added.

Validation: 442 focused expectations and 1400 full-suite expectations passed without failures,
warnings, or skips. A fresh `R CMD check --no-manual` completed with Status OK. Air and lintr pass.
Function help was regenerated with roxygen 8.1.0; the second pass was idempotent. The package check
used installed dependencies because the CRAN index was unavailable, as in earlier phases.
Artifacts are `/tmp/smallvis-phase4a-{focused,tests,lint,directional,roxygen}.log` and
`/tmp/smallvis-phase4a-check-1lujw9_q`. This checkpoint is committed as `817a608`; dynamic HSSNE
was still pending at that checkpoint. Inspection of installed mize 0.2.5.9001 confirms function
and gradient caches
tagged for the next iteration persist across steps. Restrict dynamic HSSNE to gradient-only
optimizers instead of adding a mize reinitialization protocol; fixed-alpha HSSNE remains available
for function-based optimization. Adaptive gradient histories remain intentional optimizer state.

### Phase 4 dynamic HSSNE outcome, 2026-11-06

The decrement candidate is now checked against the positive lower bound using `alpha - 0.1`.
When a different alpha is selected, the epoch hook invalidates W through `cost_clear()` and clears
Z, the gradient, and column costs. Input probabilities and their alpha-independent cache remain
intact. An unchanged alpha retains its coherent state. Candidate ordering and tie selection are
preserved; no new adaptation algorithm or optimizer reset mechanism was introduced.

Dynamic HSSNE now uses the existing gradient-only compatibility restriction, before input
calibration. This prevents mize from reusing objective data across alpha changes. Fixed-alpha
HSSNE still works with L-BFGS; internal adaptive-gradient and momentum histories remain available
for the intended dynamic experiments.

Regressions reproduce the invalid-zero-alpha case and compare the very first cost/gradient after
an accepted 0.5-to-0.6 transition against freshly initialized fixed-alpha HSSNE, without manually
clearing the transition's cache. A constructed current-alpha optimum checks the no-change case.
The earlier rounded DHSSNE public trajectory still passes without rewriting its expected values.

Final phase 4 validation on Linux:

| Check | Result |
| --- | --- |
| Focused dynamic-HSSNE and optimizer-contract tests | 85 passing expectations; no failures, warnings, or skips |
| Full `Rscript --vanilla -e 'testthat::test_local()'` from `smallvis/` | 1421 passing expectations; no failures, warnings, or skips |
| Fresh source `R CMD build`, then `R CMD check --no-manual` | Status OK; installed suite also passes all 1421 expectations |
| Checked archive compared with changed source, help, and new tests | Exact candidate included, including untracked `test_dynamic_hssne.R` |
| Air, lintr, R Markdown formatting | Passed; no lints, code-block formatting changes, or residual wide code |
| Roxygen 8.1.0 and second-pass checksums | Help regeneration is idempotent; DESCRIPTION and NAMESPACE unchanged |
| Current optimization guide and DHSSNE article renders | Passed against the freshly checked installed package |
| Static rendered HTML checks | Passed: ten optimization-guide targets, four DHSSNE targets, thirty valid MathML expressions |

The help and optimization guide describe the method restrictions; the DHSSNE article now states
the positive candidate bound. The article's historical experiments remain unchanged. Rendered
copies are under `/tmp/smallvis-phase4b-articles`; unchanged neighboring pages were copied only
to resolve local links. Browser layout inspection was not performed. No native source or bridge
changed in phase 4, and Windows was not exercised.

Logs use `/tmp/smallvis-phase4b-*.log`; the fresh archive/check is under
`/tmp/smallvis-phase4b-check-6jor1jdq`. The CRAN index remained unavailable, so the successful check
used installed dependencies. Temporary check staging now inventories both tracked and unignored
untracked package files, avoiding manual copy entries for newly added tests. The final archive's
relevant files were compared byte-for-byte with the checkout.

## Phase 5. Finish control forwarding and remove confirmed leftovers

Work in `R/sne.R`, `R/absne.R`, `R/fdiv.R`, `R/nerv.R`, `R/perplexity.R`, and `R/util.R`.

1. Forward `use_cpp` through HSSNE's inherited initializer. Forward `use_cpp` and `n_threads`
   through AB-SNE's special branches and replace its hard-coded native selection with the caller's
   setting. Inspect directly related constructors for the same omission while touching this code.
2. Fix the three `use_cpp == use_cpp` argument typos in `fdiv.R`. Check the effective initializer
   and update routes rather than assuming a corrected argument changes the overridden closures.
3. Remove the earlier shadowed `tasne()` definition. Reconfirm that `user_idp_perps()` and
   `graph_to_sparse()` have no supported consumers, including dynamic lookup, before removing
   them. Retain dynamically selected methods and optimizers; static reachability alone is
   insufficient evidence for deletion.

**Acceptance:** targeted public routing tests demonstrate native and R selection at both ordinary
and special parameter values, including thread forwarding. Existing numerical equivalence tests
pass within their documented tolerances. No exports are removed or added, and the diff contains
only the identified cleanup and required tests. Regenerate bridge files only if declarations change.

The checkout started clean at `1a5bed0`. Dynamic HSSNE had the same omitted inherited-initializer
`use_cpp` as fixed HSSNE and belongs in this correction. The f-divergence typos occupied
the positional symmetrize slot of the base constructor, but each method replaces that initializer
and supplies its own update closure; inspect effective routing rather than claiming a new numerical
effect from the typo alone. Commit the control fixes separately from removal of unused code.
An AST scan of related constructors also found the same missing `use_cpp` in SNERV's inherited
SSNE initializer. Include that one-line forwarding fix and its public routing case.

### Phase 5 backend-control checkpoint, 2026-11-06

HSSNE, dynamic HSSNE, and SNERV now forward `use_cpp` to their inherited initializer. All four
AB-SNE special branches forward both backend and thread controls, and the general branch honors
the selected backend rather than forcing native calibration. The three f-divergence comparisons
are corrected to named arguments; their overridden initializers and updates already honored the
controls, so no separate numerical effect is claimed for those typos.

A public regression covers eleven method configurations across top-level and nested controls,
R/native selection, and one/two threads. Wrappers record real native-calibrator calls and output
distance controls while calling the actual implementations. Equal input distances isolate routing
from legitimate calibration-precision differences; P, coordinates, and scalar costs agree across
the selected routes. The original test exposed 38 failures before the fixes; the added SNERV
case covers the further omission found by the constructor scan.

Validation: 751 focused expectations and 1696 full-suite expectations passed without failures,
warnings, or skips. A fresh `R CMD check --no-manual` completed with Status OK. Air and lintr pass;
existing rounded trajectories needed no updates. A repeat AST scan found no omitted named backend
or thread controls in the inspected SNE constructor families. The checked archive matches the
changed source and regression file byte-for-byte. No signatures, help, exports, or native sources
changed. The CRAN index was unavailable, so checking used installed dependencies.

Logs are `/tmp/smallvis-phase5a-{focused,tests,lint,red}.log` and
`/tmp/smallvis-phase5-constructor-audit{,-final}.log`; the fresh archive/check is under
`/tmp/smallvis-phase5a-check-n7uoxs53`. This correction is committed as `e9fc703`; unused-code
removal follows in its own commit.

### Phase 5 cleanup outcome, 2026-11-06

Removed 98 lines: the earlier shadowed `tasne()` definition, unused `user_idp_perps()`, and unused
`graph_to_sparse()`. Repository-wide searches found only the helper definitions and historical
review/plan mentions. The runtime method lookup uses the explicit method selector list, and the
optimizer lookup selects optimizer constructors; neither is a consumer of these helpers.
The method registry and all supported methods remain intact.

A parsed-expression snapshot from before cleanup matches the remaining active `tasne()` definition
exactly. Existing gradient and public API tests continue to exercise that method. No tests were
added solely to assert the absence of unused internal functions.

Final phase 5 validation on the same Linux/R toolchain:

| Check | Result |
| --- | --- |
| Full `Rscript --vanilla -e 'testthat::test_local()'` from `smallvis/` | 1696 passing expectations; no failures, warnings, or skips |
| Fresh source `R CMD build`, then `R CMD check --no-manual` | Status OK; installed suite also passes all 1696 expectations |
| Final archive compared with modified source and routing regression | Exact current source included |
| Installed public exports | Exactly `smallvis`, `smallvis_rep`, and `smallvis_perpstep` |
| Air and loaded-namespace lintr | Passed; no lints |
| Roxygen synchronization | Recorded and installed versions are 8.1.0; `needs_roxygenize()` is FALSE |

The control fixes restore the documented options, so help changes were unnecessary in this phase.
No signatures, generated bridges, help, native sources, or exports changed. Windows was not
exercised. The package check again used installed dependencies because the CRAN index was
unavailable. Logs are `/tmp/smallvis-phase5b-{tests,lint,retained-tasne}.log`;
the final archive/check is under `/tmp/smallvis-phase5b-check-9jvq2y2m`. Raw temporary artifacts
are not prerequisites for continuing the plan.

## Phase 6. Remove the identified native allocations and measure the result

Work in `src/gradients.cpp`, extend `test_native_allocations.R`, and adapt the repository's
`benchmarks/native-allocations.R` and `benchmarks/README.md` only as needed for these measurements.
Use the completed correctness work as the comparison baseline so numerical fixes do not confound
the allocation comparison.

1. Borrow the protected R input buffers for `mmds_grad_cpp()` on the main thread, following the
   existing `tsne_grad_cpp()` pattern. Workers must not call the R API or modify input storage.
2. Allocate Student-weight output in R once and fill its buffer directly. Preserve the zero
   diagonal policy of the owning caller and the existing row/column conventions. The distance
   kernels have similar staging, but expanding to them requires separate measured justification;
   they are not additional commitments in this phase.
3. Reuse ownership tests with aliased inputs and independent gradient/weight references. Compare
   one and two threads and complete short MDS and t-SNE runs against the pre-allocation revision.
4. First update the benchmark's variant check: its existing `baseline` mode requires the old
   `_smallvis_d2_to_tweight_cpp` routine, which is absent from both versions in this comparison.
   Record the two actual source revisions, loaded libraries, DLLs, and common build settings.
   Add only the MDS-gradient and Student-weight cases needed here; retain useful historical cases.
5. Smoke the harness at a small size, then compare separately installed builds in fresh processes.
   Hold synthetic inputs, BLAS threads, compiler flags, and timing boundaries constant. Use at
   least three alternating repetitions at a practical size and one/two worker threads. Include
   full embedding measurements so copy removal is assessed in its actual consumer.

**Acceptance:** objectives and gradients remain equivalent, aliased inputs remain unchanged, and
native checks find no new invalid access. Static inspection confirms removal of the two MDS input
copies and the Student-weight staging matrix. Comparative process peak RSS and timings are
recorded with limitations; at 5,000 points the removed payloads are 400 MB and 200 MB respectively,
not promises about RSS or elapsed time. Investigate an unexplained material regression before
retaining a change. Add no timing assertions to routine unit tests and claim speed gains only when
the measurements support them.

## Phase 7. Automate native regressions and close the review

1. Run the existing tree and threading harnesses, plus the new tree boundary case, in a small
   bounded Linux CI job using ASan/UBSan. Reuse the repository's pinned actions, read-only default
   permissions, and existing hosted-runner conventions. Ensure native source, harness, and
   workflow changes exercise the job. Do not introduce a service, periodic benchmark job, or
   broad sanitizer matrix.
2. Complete the full package, formatting, lint, compiler, and native verification below. Use the
   existing Linux/Windows package workflow for platform evidence when hosted execution is
   authorized; local-only checks must not be reported as a Windows pass.
3. Re-read the saved review and reconcile every row of the disposition table with the final
   changes and evidence. Correct relevant help and generated documentation through their source
   generators. Preserve optional-feature exclusions and the earlier plan's completed status.
4. Record the final revision, results, remaining limitations, and any changed dispositions here.
   Mark implementation complete only when the admitted work and required evidence are complete.

The phase 5 documentation cross-check found that the Barnes-Hut epsilon help still says points
share a leaf only at identical coordinates. Clarify this during the final documentation pass:
epsilon does not merge distinct coordinates; finite-resolution terminal leaves preserve separate
coordinate groups. This reconciles the existing help with the phase 1 tree correction.

## Validation and recovery

After implementation is authorized, start each phase by checking Git status and the current
revision. Reproduce its distinguishing failures before changing behavior. Keep regressions with
their fixes, and use the package's existing tests rather than adding a separate test framework.

Run focused tests from `smallvis/`, selecting the affected files, for example:

```sh
Rscript -e 'testthat::test_local(filter = "native_bh|native_perplexity")'
```

At package-wide integration, run the full suite and loaded-namespace lint check from `smallvis/`:

```sh
Rscript -e 'testthat::test_local()'
Rscript -e 'pkgload::load_all(); lints <- lintr::lint_package(); print(lints); quit(status = if (length(lints)) 1L else 0L)'
```

Use the configured tool versions and `air format smallvis --check` from the repository root.
Build the source archive and run `R CMD check --no-manual` in a temporary directory. Run strict
GCC/Clang checks on handwritten C++ and the focused native harnesses after native changes.
Exclude generated Rcpp exports from handwritten style judgments. Regenerate roxygen or Rcpp
outputs only when their source declarations change, and verify a second generation is unchanged.

The existing native harnesses can be compiled from `smallvis/` with the commands in their source
headers. Add sanitizer flags for diagnostic runs and a process timeout for the dangerous tree
case. Keep raw logs, installed comparison packages, build archives, and benchmark CSVs outside
the checkout. If a required diagnostic is unavailable, name the limitation and the remaining
check; do not infer success from an older run.

Keep phases separately reviewable and preserve the baseline for the allocation comparison.
On resumption, reconcile completed source changes and evidence with the progress list. Do not
reset unrelated files or apply old temporary scripts blindly. If a proposed fix requires a new
methodological contract or feature, leave that expansion out and record the smallest supported
resolution of the existing failure.

## Decisions and open implementation choices

- Decision: save the source review separately from the execution plan. Rationale: findings and
  baseline evidence remain available while implementation status evolves. Date/author:
  2026-11-06 / Codex, following the user's save-and-plan request.
- Decision: plan corrections, cleanup, and allocation work without starting implementation.
  Rationale: the user explicitly requested this boundary. Date/author: 2026-11-06 / user.
- Decision: exclude optional features. Rationale: the user explicitly deferred them; they are
  preserved only in the source review and exclusion table. Date/author: 2026-11-06 / user.
- Decision: implement phase 1 only. Rationale: the user requested starting the first part of the
  plan. This boundary was subsequently extended to phase 2 by explicit instruction.
  Date/author: 2026-11-06 / user.
- Decision: implement phase 2, including scalar-reference and stopping regressions. Rationale:
  the user explicitly requested the AB-SNE scalar fixes and stopping behavior.
  Date/author: 2026-11-06 / user.
- Decision: implement phase 3 only after completing phases 1 and 2. Rationale: the user explicitly
  requested input equivalence, neighbor handling, and multiscale contracts.
  Date/author: 2026-11-06 / user.
- Decision: commit completed phases 1–3, then implement phase 4 and commit completed changes.
  Rationale: the user explicitly requested this sequence and ongoing commit practice. No push or
  publication was requested. Date/author: 2026-11-06 / user.
- Decision: implement phase 5 and commit backend controls separately from unused-code cleanup.
  Rationale: the user authorized phase 5; the two changes have different behavior and validation
  responsibilities. Phases 6–7 remain pending. Date/author: 2026-11-06 / user and Codex.
- Decision: derive degrees and symmetric conditional weights from averaged P, and require an
  explicit output beta where automatic bandwidth transfer would be undefined. Rationale: these
  retain the methods' mathematics without defining a bandwidth for a mixture.
  Date/author: 2026-11-06 / Codex.
- Decision: retain exact coordinate/count groups only when subdivision cannot make progress.
  Rationale: this bounds construction without merging distinct points or changing normal approximate
  traversal. Date/author: 2026-11-06 / Codex.
- Proposed method: use the existing compatibility mechanism to restrict unmatched line searches
  and reject undefined multiscale combinations. Rationale: these restore meaningful behavior
  without adding objectives or APIs. Date/author: 2026-11-06 / Codex.
- Open implementation choices: the native CI job's placement. Resolve this within phase 7
  using the stated witnesses; choosing a routine implementation does not expand the scope.

## Progress and outcomes

- [x] Save the second review and create this phased plan.
- [x] Phase 1: native tree and Gaussian calibration boundaries.
- [x] Phase 2: AB-SNE scalar objectives.
- [x] Phase 3: input, neighbor, and multiscale contracts.
- [x] Phase 4: optimizer compatibility and dynamic HSSNE state.
- [x] Phase 5: backend controls and cleanup.
- [ ] Phase 6: native allocations and comparative measurements.
- [ ] Phase 7: native CI, integration checks, and review reconciliation.

Phases 1–5 are implemented and validated as recorded above. Their regressions distinguish the
original native, calibration, scalar, stopping, input, neighbor, multiscale, optimizer, and dynamic
state and backend-control failures. Confirmed unused internals have been removed. Completed
changes follow the user's commit-after-each-change instruction. Phases 6–7, benchmarks, workflows,
generated bridges, and site output are untouched.
