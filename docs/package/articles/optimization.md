# Optimization and repeated runs

`smallvis` makes it possible to change the objective and optimizer
independently. This is useful for experiments, but a setting that works
well for t-SNE need not work well for another cost function. Start by
changing one choice at a time, keeping preprocessing and initialization
fixed.

## Configure a method

Use a string for a method’s defaults, or a list whose first element
names the method and whose other elements override its options. For
example, this configures the balance of repulsion and attraction in the
LargeVis cost function:

``` r

largevis_fit <- smallvis(
  iris,
  method = list("largevis", gamma = 1),
  Y_init = "spca",
  perplexity = 30,
  eta = 0.1,
  max_iter = 1000,
  epoch_callback = FALSE,
  n_threads = 1,
  verbose = FALSE
)
dim(largevis_fit)
```

    ## [1] 150   2

The package implements the LargeVis and UMAP cost functions with its own
optimization machinery. These runs are different from the original
algorithms’ sampled stochastic optimization procedures. The
[LargeVis](https://jlmelville.github.io/smallvis/package/articles/largevis.md)
and
[UMAP](https://jlmelville.github.io/smallvis/package/articles/umap.md)
articles examine those differences; the [method-specific
options](https://jlmelville.github.io/smallvis/package/reference/smallvis.html#method-specific-options)
reference lists the available controls.

## Change the optimizer

The default is the delta-bar-delta scheme familiar from t-SNE. The `opt`
argument also accepts optimizers provided through
[mize](https://github.com/jlmelville/mize), with options supplied in a
list. For UMAP, t-UMAP, and LargeVis, function-based line searches
require `gr_eps = 0` so the gradient matches the objective used by the
search:

``` r

umap_fit <- smallvis(
  iris,
  scale = FALSE,
  method = list("umap", gr_eps = 0),
  opt = list("l-bfgs", c1 = 1e-4, c2 = 0.9),
  Y_init = "spca",
  max_iter = 300,
  epoch_callback = FALSE,
  n_threads = 1,
  verbose = FALSE
)
dim(umap_fit)
```

    ## [1] 150   2

The [alternative
optimizers](https://jlmelville.github.io/smallvis/package/reference/smallvis.html#alternative-optimizers)
reference gives the supported settings. Longer experiments compare
[L-BFGS](https://jlmelville.github.io/smallvis/package/articles/opt.md),
[conjugate
gradient](https://jlmelville.github.io/smallvis/package/articles/cg.md),
[spectral
directions](https://jlmelville.github.io/smallvis/package/articles/specd.md),
and [gradient-based adaptive
updates](https://jlmelville.github.io/smallvis/package/articles/sgd.md).

## Use exaggeration

Early exaggeration strengthens attractive interactions for an initial
phase. In this t-SNE run, input probabilities are multiplied by 4 for
the first 100 iterations:

``` r

exaggerated_fit <- smallvis(
  iris,
  Y_init = "spca",
  perplexity = 30,
  eta = 100,
  exaggeration_factor = 4,
  stop_lying_iter = 100,
  epoch_callback = FALSE,
  n_threads = 1,
  verbose = FALSE
)
```

`late_exaggeration_factor` and `start_late_lying_iter` enable a later
phase; the [late exaggeration
experiments](https://jlmelville.github.io/smallvis/package/articles/lexagg.md)
explore its effect. Function-based line searches are restricted to
configurations whose objective and gradient agree; smallvis reports
unsupported combinations rather than silently treating them as ordinary
optimization.

## Compare random starts

[`smallvis_rep()`](https://jlmelville.github.io/smallvis/package/reference/smallvis_rep.md)
runs from random initial coordinates and returns the result with the
lowest final cost. It overrides `Y_init`; built-in methods share one
prepared input graph across the starts. With `ret_extra = TRUE`,
`all_costs` contains the final cost from each run:

``` r

set.seed(42)
best <- smallvis_rep(
  nrep = 3,
  X = iris,
  perplexity = 30,
  ret_extra = TRUE,
  epoch_callback = FALSE,
  n_threads = 1,
  verbose = FALSE
)
best$all_costs
```

    ## [1] 0.1257062 0.1293343 0.1258536

``` r

sum(best$costs)
```

    ## [1] 0.1257062

This comparison holds the objective and affinities fixed. Costs from
different methods or input graphs generally have different meanings and
should not be used as a universal embedding score. Use `keep_all = TRUE`
to inspect every layout; see the [repeated-run
reference](https://jlmelville.github.io/smallvis/package/reference/smallvis_rep.md).

## Step down through perplexities

[`smallvis_perpstep()`](https://jlmelville.github.io/smallvis/package/reference/smallvis_perpstep.md)
starts with larger perplexities, then optimizes at the requested target.
The `step_iter` budget is taken out of `max_iter`, rather than added to
it:

``` r

stepped_fit <- smallvis_perpstep(
  X = iris,
  Y_init = "spca",
  perplexity = 30,
  step_iter = 50,
  max_iter = 300,
  epoch_callback = FALSE,
  n_threads = 1,
  verbose = FALSE
)
dim(stepped_fit)
```

    ## [1] 150   2

The [perplexity-stepping
reference](https://jlmelville.github.io/smallvis/package/reference/smallvis_perpstep.md)
describes the schedule and which controls apply to the final phase.
