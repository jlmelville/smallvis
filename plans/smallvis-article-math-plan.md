# Article mathematics presentation

## Purpose and current state

The owner requested a complete pass over article mathematics, using the clearer layouts already
adopted by mize and flotsam. Help readers follow definitions, derivations, and comparisons without
changing the methods or the author's explanatory voice. This follows the completed pkgdown
migration; it is a separate presentation pass.

Status: implementation and validation complete. The owner reviewed the preview and approved
committing, pushing, and following CI. The checkout started clean at `bb2e9df`, matching upstream.
Preserve the experiment stashes, especially `1bc0e8f888be8064874928e86dde501b8ae41958`.

## Progress

- [x] Inventory all articles using parsed Markdown, excluding code from math detection.
- [x] Check the R documentation skill and inspect mize/flotsam's adopted equation layouts.
- [x] Improve layouts and mathematical typography where they help the reader.
- [x] Check formula preservation and separately validate any genuine formula correction.
- [x] Render the collection; inspect every math article at desktop and mobile widths.
- [x] Run formatting and package checks.
- [x] Complete final scope review.

## Scope and decisions

The 45 articles contain 30 pages with math: 1,249 inline expressions and 256 display equations.
The inventory is in scratch state, including each parsed expression and its article. Simple
parameter mentions do not need embellishment. Prefer standard KaTeX syntax: aligned equalities
and updates, grouped definitions, readable fractions, upright descriptive labels and operators,
and consistent vector indices. Use boxes or annotations only if they carry explanatory value.

Preserve derivation steps, notation choices, historical results, code, figures, links, heading IDs,
and all package behavior. Do not invent a universal notation system, add custom macros or plugins,
or change the renderer. A mathematical typo discovered during this work must be distinguished
from typography and supported by a local derivation or numerical witness before correction.
Generated `docs/package/` remains owned by publishing CI; use temporary preview output.

## Work and acceptance

1. Review all 30 mathematical articles, including short inline-only pages. Align multi-step
   derivations and related updates; break long displays at mathematical boundaries. Keep short
   formulas short and retain prose between reasoning steps.
2. Compare parsed article structures against the starting revision. Verify unchanged code,
   figures, headings, and links; account for each mathematical edit. Separate notation and layout
   from substantive corrections, with a numerical witness for the latter.
3. Build pkgdown in a temporary directory. In Chromium, inspect all 30 math pages at ordinary
   article width and mobile width for KaTeX errors, missing expressions, page overflow, clipped
   lines, misleading alignment, and readability. Visually inspect the densest derivations and
   representative inline typography. Long formulas may scroll within their own mobile container.
4. Run the R Markdown formatter on changed sources, configured Air and lintr, full package tests,
   and `R CMD check --no-manual`. Confirm website-only edits stay outside the package archive.
   Do not rerun the historical embedding experiments or claim new evidence for their results.

## Context, evidence, and recovery

Sources live in `smallvis/vignettes/articles/`; `_pkgdown.yml` selects KaTeX and
`smallvis/pkgdown/extra.css` already contains display-equation overflow. The relevant skill is
`r-docs-pkgdown`, especially `references/rmarkdown-editing.md` and its mathematical-edit section.
Mize's momentum articles and flotsam's LTSA articles demonstrate useful aligned derivations and
grouped constraints. Supported syntax is documented at <https://katex.org/docs/supported.html>.

Scratch evidence is under `/tmp/smallvis-math-review/`. The baseline parsed inventory is
`before-math.json`. Reconcile source changes with this plan on resumption; do not regenerate
baseline evidence from an edited tree. Keep validation scripts and generated previews in scratch
state. Git diffs and the clean starting revision are the recovery boundary; never reset unrelated
work or apply/drop a stash.

## Outcomes

Reviewed all 30 mathematical articles and edited 21. The nine unchanged articles use simple
parameters or short expressions that already read clearly. Retained every original display;
the Swiss Roll neighbor rule is now one cases display rather than several inline fragments.
The final collection contains 257 displays and 1,246 inline expressions (1,503 total).

Presentation improvements include aligned equalities and paired coordinate updates; line breaks
at mathematical boundaries; matrix and cases environments; upright method names, descriptive
subscripts, and operators; italic
indices outside bold vectors; a proper bold Greek parameter; conditional-probability bars;
upright differentials; and correctly scoped hats and tildes. The final t-SNE cost labels its
attractive and repulsive terms with underbraces. One boxed identity makes the central variance
result easy to find again in the long intrinsic-dimension derivation. The changes use standard
KaTeX syntax without macros, plugins, renderer changes, or CSS changes.

### Corrections distinct from typography

- `theory.Rmd` and `ncvis.Rmd`: the denominator in the derivative of the output kernel must contain
  `a d_{ij}^{2b}`, matching the kernel definition, rather than `a d_{ij}^2`. This matters when
  `b != 1`. The NCVis prose also now calls `d_{ij}^2` a squared Euclidean distance.
- `ftsne.Rmd`: the CH-SNE cost divides squared probability differences by `q_{ij}`, matching the
  implementation and the displayed gradient; its former denominator `q_{ij}^2` was incorrect.
- `absne.Rmd` and `ftsne.Rmd`: added brackets to put all AB-divergence summands and both
  Jensen-Shannon terms under their respective sums and prefactors.
- `theory.Rmd`: replaced the inconsistent `Cp` with the previously defined `C_P`.
- `idp-class.Rmd`: the local dimension in the introduction is `delta`, as in the equations,
  rather than `D`. `idp-theory.Rmd` now keeps the probability subscript outside the variance
  operator. `specd.Rmd` consistently uses subscripted line-search constants.

### Validation evidence

- Independent formula checks in scratch `check-formulas.R` use 24 deterministic observations,
  nontrivial coordinates, and a normalized direction. CH-SNE, Jensen-Shannon, and AB-divergence
  (`alpha = 0.7`, `beta = 1.2`) costs match the package: maximum absolute discrepancy
  `3.47e-18`. Their gradients agree with directional differences at steps `1e-4`, `1e-5`, and
  `1e-6`: maximum discrepancy `7.72e-10`. The CH-SNE check also directly transcribes the displayed
  gradient and confirms that the former cost denominator gives a different result.
- The kernel derivative was checked at 18 combinations of `a = {0.3, 1.7}`, `b = {0.6, 1, 2.1}`,
  and squared distance `{0.2, 1.3, 3.5}`, using the same three difference steps. Maximum error
  `4.36e-8`; the former expression differs by as much as `0.501`. Both simplified derivative
  identities agree within `1e-12`.
- Parsed Markdown comparison across all 45 articles preserves code, figures, headings, links, and
  prose structure. The only exceptions are the single NCVis word correction and the explicitly
  reviewed Swiss Roll cases paragraph. All expression braces balance; all multiline displays
  use an equation environment. The mathematical diff was separately reconciled with the
  typography changes and the corrections above.
- R Markdown formatting: 45 files, 107 R blocks, zero changes needed and zero residual wide code
  lines. Air and loaded-namespace lintr pass.
- Full development tests: 968 passes, zero failures/skips, seven existing partial-argument
  warnings. Source build and `R CMD check --no-manual`: `Status: OK`, with installed tests passing.
  Website articles remain excluded from the package archive. No implementation files changed.
- Complete pkgdown build passed. Initial Chromium review covered all 30 mathematical articles at
  1280 and 390 pixels, including all 1,505 expressions then present: no missing math, KaTeX errors,
  broken images, or whole-page overflow. It identified one desktop UMAP display with 17 pixels of
  internal scrolling, subsequently split over two lines.
- Final browser checks cover all 30 mathematical articles at both widths, with a focused rerun
  after the last spectral matrix adjustment. All 1,503 expressions and 257 displays render at
  each width. No KaTeX errors, missing expressions, broken images, whole-page overflow, or desktop
  equation scrolling. Forty-seven longer displays scroll inside their mobile containers.
  Screenshots were visually reviewed for the central theory, spectral, intrinsic-dimension,
  divergence, NCVis, PaCMAP, and Swiss Roll layouts. Static checks pass for all 45 article pages,
  including local asset references, links, and section targets.
- Final scope is 21 article sources and this plan. No package code, generated site output,
  dependencies, workflows, or stashes changed. The preview is served locally at
  <http://localhost:8768/package/articles/theory.html>.

### Reviewed article inventory

Names correspond to sources under `smallvis/vignettes/articles/`. Counts are for the final source.

| Article | Inline | Display | Edited |
| --- | ---: | ---: | --- |
| `absne` | 66 | 17 | Yes |
| `cg` | 1 | 0 | No |
| `dhssne` | 30 | 0 | No |
| `ee` | 15 | 6 | Yes |
| `ftsne` | 4 | 9 | Yes |
| `gsne` | 7 | 3 | Yes |
| `hssne` | 13 | 2 | Yes |
| `idp-class` | 17 | 5 | Yes |
| `idp-theory` | 133 | 56 | Yes |
| `idp` | 17 | 4 | Yes |
| `largevis` | 97 | 6 | Yes |
| `mmds` | 7 | 0 | No |
| `msperp` | 2 | 0 | No |
| `ncvis` | 130 | 43 | Yes |
| `nervjse` | 20 | 6 | Yes |
| `norm` | 16 | 0 | No |
| `opt` | 5 | 0 | No |
| `optsne` | 6 | 0 | Yes |
| `pacmap-examples` | 3 | 0 | Yes |
| `pacmap` | 94 | 19 | Yes |
| `perplexity` | 6 | 1 | Yes |
| `sne` | 19 | 0 | No |
| `specd` | 8 | 0 | Yes |
| `spectral` | 331 | 43 | Yes |
| `swisssne` | 9 | 2 | Yes |
| `tee` | 98 | 7 | Yes |
| `theory` | 84 | 25 | Yes |
| `three-clusters` | 1 | 0 | No |
| `umap` | 0 | 3 | Yes |
| `umaptsne` | 7 | 0 | No |
