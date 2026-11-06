# Consolidate the smallvis documentation

## Purpose and current state

On 2026-11-06, the owner approved moving the complete article collection into pkgdown rather
than migrating a few examples. Make the existing practical guidance, theory, and experiments
discoverable through one site. Preserve the author's explanations, experimental results, and
existing public article URLs. The owner reviewed the local preview and approved committing,
pushing, and linking the repository's About website to the pkgdown site.

Status: migration, the approved KaTeX follow-up, local validation, and owner review are complete.
The clean checkout was fast-forwarded to `73a0408`, the generated documentation update following
the last accepted usage correction. The experiment stash remains
`1bc0e8f888be8064874928e86dde501b8ae41958` and must remain untouched.

## Scope and decisions

- Move all 40 research Markdown sources from `docs/` to website-only R Markdown articles under
  `smallvis/vignettes/articles/`. Keep the original dates, equations, experiments, and figures.
- Move the shared image tree into the website article sources, retaining every image byte.
  The existing `^vignettes$` package-build exclusion keeps these assets out of installations.
- Organize the article index by topic. Provide focused practical guides for getting started,
  inputs and initialization, optimization, and Barnes-Hut/native acceleration. Retain the old
  usage guide's background material on a separate page.
- Keep pkgdown at `/smallvis/package/`, with static redirects from the old root article URLs.
  Preserve fragment navigation and update internal article links to their canonical destinations.
- Keep historical examples and numerical results precomputed. Run a bounded set of current
  practical examples during site builds; do not claim to reproduce the historical experiments.
- Preserve the function reference as the detailed API contract. Do not add the rejected
  capability table, rewrite mathematical derivations, change package behavior, or expand licensing.
- The existing workflow owns generated `docs/package/` output. Build a separate temporary preview
  for validation and leave deployment output to CI unless publication requires a source change.

## Progress

- [x] Inventory the collection and inspect the publishing boundary.
- [x] Refresh upstream without modifying the experiment stash.
- [x] Move sources/assets and connect navigation and legacy redirects.
- [x] Restructure practical guidance and exercise its examples.
- [x] Build and inspect all articles, links, equations, and representative layouts.
- [x] Run package, formatting, and lint checks; reconcile the final diff and stash identity.
- [x] Switch to KaTeX, correct the inline exponent grouping, and validate the complete site.

## Discoveries

The 40 articles contain no executable R chunks, iframes, or embedded scripts. They reference
2,471 distinct existing images (about 66 MiB); the full shared tree contains 2,826 images (about
82 MiB). Migration does not require rerunning these experiments. Their many equations and figure
grids require rendered inspection, while changing rendering does not establish numerical validity.

The root home page's old installation, prerequisites, and article-list fragments also redirect to
matching current destinations. The existing CI publishes only to `docs/package/`, leaving
root-level compatibility redirects intact. Its `smallvis/**` path filter already covers the
relocated sources and assets.

## Validation and acceptance

1. Inventory every original article and image before moving it. Confirm all 40 destination
   articles exist and image hashes remain identical after the move.
2. Build the complete pkgdown site from the current package in a temporary destination. Require
   every article in the index and search, resolving internal links and local assets, and redirects
   that reach the corresponding page and preserve fragments.
3. Inspect rendered pages in Chromium at desktop and mobile widths. Check navigation, equations,
   ordinary examples, long code, and dense figure grids. Inspect all articles for broken images,
   math-rendering errors, and whole-page overflow; inspect representative screenshots visually.
4. Run practical examples in fresh article processes. Keep installation snippets display-only.
   Check the guide's observable result shapes and finite outputs without freezing random layouts.
5. Run the full installed test suite, build and check the source package with `--no-manual`, and
   inspect its archive to confirm website articles and images remain excluded. Run configured Air
   and lintr checks, R Markdown formatting, and whitespace checks. Audit Actions if workflows
   change.
6. Confirm no package implementation or generated bindings changed, no temporary build artifacts
   entered the checkout, and the preserved experiment stash has the same identity.

## Validation results so far

- All 40 original articles were migrated, together with four current practical guides and one
  background page (45 articles total). All 2,826 original image hashes match after relocation.
- All 256 display equations are unchanged apart from whitespace. Mathematical headings were given
  readable sidebar labels with their original explicit IDs; no formulas or results were revised.
- The article index and search include all 45 articles. Every original section ID survives.
- The complete source build and `R CMD check --no-manual` passed: 968 expectations, zero errors,
  warnings, or notes. The archive contains no website articles, images, pkgdown config, or output.
- The full development suite passed with its seven existing partial-argument warnings. Air and
  loaded-namespace lintr passed; roxygen reports no synchronization needed.
- All 107 labeled R fences pass the formatter with no residual wide lines. Unlabeled historical
  examples were labeled only after parsing them as R; historical code remains display-only.
- Initial automated inspection found no broken images or explicit math errors. Visual review and
  comparison with the source equation inventory caught 19 display equations in six articles that
  Pandoc's MathML conversion left as raw TeX. MathJax, also used by the original articles, rendered
  all 25 display equations in the theory pilot, including its three skipped MathML equations.
  Final checks now compare rendered display counts with each article's independent source count
  and search visible prose for unrendered delimiters.
- MathJax equations scroll within their containers, preserving page width on mobile. Figure grids
  also scroll on small screens to retain readable labels. One long inline expression in the
  spectral article was split into equivalent expressions to fit a mobile line.
- The initial MathJax browser evidence covers all 45 articles and the index at both widths: 92
  current page/width cases, with no broken images, raw TeX, math errors, sidebar math artifacts,
  whole-page overflow, or desktop code overflow. All 256 original display equations render.
- All 40 legacy article URLs preserve query parameters and section fragments in Chromium. The old
  home page and its three section destinations also redirect correctly. Interactive search finds
  migrated content. Static checks report no broken article links or missing local assets.
- Each practical guide also ran in a fresh R session outside pkgdown. All promised coordinate
  shapes and finite values passed, as did the repeated-start cost selection and PCA input shape.

## Renderer decision and 404 follow-up

During MathJax validation, pkgdown's generated `404.html` added a nonexistent `/NA` script URL
when converting its inline configuration to absolute links. The same `tweak_link_absolute()` code
occurs in the installed 2.2.0 and the inspected CRAN 2.2.1 source. The page still displayed
its message and navigation normally. A custom build wrapper or copied template was omitted.

The owner subsequently chose KaTeX, matching the local mize, funconstrain, flotsam, and rnndescent
configurations. The renderer is now `katex`, and the CSS uses their `.katex-display` overflow rule.
The spectral article's inline diffusion-map exponent now has the missing grouping braces, with no
change to its intended meaning. The built page renders the corrected expression at both widths.

The complete KaTeX site build passed with pkgdown 2.2.0. Browser checks covered all 45 articles,
the article index, five reference pages, the homepage, and the generated 404 page at 1280 and
390 pixels: 106 cases, with no failed pages, broken images, raw TeX, math errors, sidebar artifacts,
or whole-page overflow. All 256 original display equations render. Representative equations and
figure grids were visually inspected at both widths. One unchanged reference example has 15 pixels
of internal horizontal scrolling; the article examples fit at desktop width.

Static checks cover all 99 generated HTML files, all 45 indexed articles, all 876 search entries,
all preserved section IDs, and all 2,826 image hashes. There are no missing local assets, including
on the generated 404 page, and no `/NA` script URL. All legacy redirects and interactive search
passed against the KaTeX build. R Markdown formatting passed for all 107 R blocks with no residual
wide lines. No 404 workaround or workflow change was needed.

The previously checked source archive still matches all 61 checked package files. DESCRIPTION has
no source changes, and the changed renderer configuration, CSS, and article are excluded from the
archive, so the earlier successful package checks remain applicable. Historical experiments were
preserved, not numerically rerun. Evidence is in `/tmp/smallvis-docs-migration.wtZVAp`, particularly
`katex-site-build.log`, `katex-site-browser-results.json`, `katex-static-results.json`, and
`katex-legacy.log`. The existing local preview now serves the KaTeX build.

## Outcome and publication

The complete article collection now has one source and one pkgdown index. Practical guidance is
separate from dated explorations, and existing public article links remain usable. No package
implementation, generated bindings, function help, dependencies, licensing, or workflows changed.
The existing CI builds and publishes `docs/package/` from master; generated output was not copied
into the checkout. On 2026-11-06 the owner approved publishing the reviewed migration and setting
the repository's About website to `https://jlmelville.github.io/smallvis/package/`.

## Recovery and evidence

Keep scratch scripts, manifests, preview HTML, screenshots, and command logs in a dedicated `/tmp`
directory (`/tmp/smallvis-docs-migration.wtZVAp`). Source moves can be compared
against `73a0408`; regenerate the preview from the final working tree after substantive corrections.
Record completion evidence and any unverified boundaries here before handing the result back.
