# smallvis

`smallvis` is an R package for experimenting with neighborhood embedding methods, including t-SNE,
LargeVis, and UMAP cost functions. You can control input scaling, nearest neighbors, initialization,
and optimization to compare methods and try out ideas.

Install from GitHub with `pak::pak("jlmelville/smallvis/smallvis")`. The R package lives in the
`smallvis` subdirectory of the repository.

Start with [Using smallvis](articles/usage.html) to make and plot an embedding. Continue with the
guides to [inputs and initialization](articles/inputs.html),
[optimization and repeated runs](articles/optimization.html), or
[Barnes-Hut and native acceleration](articles/barnes-hut.html).

The [article collection](articles/index.html) brings together the theory and experiments behind
the package, including spectral methods, perplexity, initialization, optimization, and comparisons
with other embedding methods. Dated explorations retain their original results and examples.
The [function reference](reference/index.html) documents the current arguments and return values.
