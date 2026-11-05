# Dimensionality Reduction With Perplexity Stepping

Carry out dimensionality reduction of a (small) dataset using one of a
variety of neighbor embedding methods, using a decreasing value of
perplexity to avoid bad local minima.

## Usage

``` r
smallvis_perpstep(step_iter = 250, ...)
```

## Arguments

- step_iter:

  Number of iterations to carry out the perplexity stepping. Must be \<
  the value of `max_iter`.

- ...:

  Arguments to be passed to
  [`smallvis`](https://jlmelville.github.io/smallvis/package/reference/smallvis.md).
  See 'Details' for information on which arguments may be modified or
  ignored during certain parts of the embedding.

## Value

The result of the final run of
[`smallvis`](https://jlmelville.github.io/smallvis/package/reference/smallvis.md)
at the target perplexity.

## Details

This function uses ideas similar to those in the NeRV (Venna et al.,
2010) and JSE (Lee et al., 2013), where to avoid local minima, the
initial optimization steps use affinities with larger bandwidths (NeRV)
or larger perplexity values (JSE). This implementation uses a series of
decreasing perplexity values, as in JSE.

For details on the arguments that can be passed to the dimensionality
reduction routine, see the help text for
[`smallvis`](https://jlmelville.github.io/smallvis/package/reference/smallvis.md).

To avoid spending too much extra time in perplexity calibrations, the
extra perplexities start at the power of 2 closest to, but not greater
than, half the dataset size (in terms of number of objects). Further
calibrations are then carried out halving the perplexity each time,
until the perplexity specified by the user is reached.

The number of iterations spent in the larger perplexity values is
specified by the `step_iter` parameter. This determines the total number
of iterations, e.g. if `step_iter = 250` and extra optimizations at a
perplexity of 1024, 512, 256, 128 and 64 will be carried out, these will
run for 50 iterations each. To keep the number of iterations equivalent
to that used by a single run of
[`smallvis`](https://jlmelville.github.io/smallvis/package/reference/smallvis.md),
the value of `step_iter` is subtracted from the value of `max_iter`
before the optimization at the target perplexity is carried out, e.g. if
`max_iter = 1000` and `step_iter = 250`, the final optimization will run
for 750 iterations only.

Any value of `tol`, `exaggeration_factor` and `stop_lying_iter` provided
is used only with the final optimization.

## References

Venna, J., Peltonen, J., Nybo, K., Aidos, H., & Kaski, S. (2010).
Information retrieval perspective to nonlinear dimensionality reduction
for data visualization. *Journal of Machine Learning Research*, *11*,
451-490.

Lee, J. A., Renard, E., Bernard, G., Dupont, P., & Verleysen, M. (2013).
Type 1 and 2 mixtures of Kullback-Leibler divergences as cost functions
in dimensionality reduction based on similarity preservation.
*Neurocomputing*, *112*, 92-108.

## Examples

``` r
embedding <- smallvis_perpstep(
  X = iris, perplexity = 30, Y_init = "spca",
  step_iter = 25, max_iter = 100, epoch_callback = FALSE,
  n_threads = 1, verbose = FALSE
)
dim(embedding)
#> [1] 150   2

if (FALSE) { # \dontrun{
# t-SNE on the iris with L-BFGS optimization
# The 1000 max_iter is split between 250 iterations at perplexity = 64
# and then 750 iterations at perplexity = 40.
iris_lbfgs_pstep <- smallvis_perpstep(
  step_iter = 250, X = iris, scale = FALSE, verbose = TRUE, Y_init = "spca",
  ret_extra = c("DX", "DY"), perplexity = 40, max_iter = 1000, opt = list("l-bfgs")
)
} # }
```
