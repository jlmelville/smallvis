# Best t-SNE Result From Multiple Initializations

Run t-SNE multiple times from a random initialization, and return the
embedding with the lowest cost.

## Usage

``` r
smallvis_rep(nrep = 10, keep_all = FALSE, ...)
```

## Arguments

- nrep:

  Number of repeats.

- keep_all:

  If `TRUE`, then the return value is a list of lists, indexed from 1 ..
  `nrep`, with each entry the result from each
  [`smallvis`](https://jlmelville.github.io/smallvis/package/reference/smallvis.md)
  run. Otherwise just the result with the lowest error is returned.

- ...:

  Arguments to apply to each
  [`smallvis`](https://jlmelville.github.io/smallvis/package/reference/smallvis.md)
  run.

## Value

The
[`smallvis`](https://jlmelville.github.io/smallvis/package/reference/smallvis.md)
result with the lowest final cost, or if `keep_all` is `TRUE` all
results as a list, indexed as 1 .. `nrep`. If `ret_extra` is not
`FALSE`, then the final costs for all `nrep` runs are also included in
the return value list as a vector called `all_costs`. In this case, if
`keep_all` is `TRUE`, then `all_costs` appears as an extra item on all
results. Additionally, each result will have an extra entry `best_rep`,
giving the index of the result with the lowest cost.

## Details

This function ignores any value of `Y_init` you set, and uses
`Y_init = "rand"`.

## Examples

``` r
set.seed(42)
best <- smallvis_rep(
  nrep = 2, X = iris, perplexity = 30, max_iter = 100,
  ret_extra = TRUE, epoch_callback = FALSE, n_threads = 1,
  verbose = FALSE
)
range(best$all_costs)
#> [1] 0.1421794 0.2018069

if (FALSE) { # \dontrun{
# Return best result out of five random initializations
tsne_iris_best <- smallvis_rep(
  nrep = 5, X = iris, perplexity = 50, method = "tsne",
  ret_extra = TRUE
)
# How much do the costs vary between runs?
range(tsne_iris_best$all_costs)
# Display best embedding found
plot(tsne_iris_best$Y)

# Keep all results
# First result is in tsne_iris_rep[[1]], second in tsne_iris_rep[[2]] etc.
tsne_iris_rep <- smallvis_rep(
  nrep = 5, X = iris, perplexity = 50, method = "tsne",
  ret_extra = TRUE, keep_all = TRUE
)
# Index of result with smallest error is in special list item 'best_rep'
best_iris <- tsne_iris_rep[[tsne_iris_rep[[1]]$best_rep]]
} # }
```
