# smallvis

`smallvis` is an R package for experimenting with neighborhood embedding
methods, including t-SNE, LargeVis, and UMAP cost functions. You can
control input scaling, nearest neighbors, initialization, and
optimization to compare methods and try out ideas.

Install from GitHub with `pak::pak("jlmelville/smallvis/smallvis")`. The
R package lives in the `smallvis` subdirectory of the repository.

Start with [Using
smallvis](https://jlmelville.github.io/smallvis/package/articles/usage.md)
to make and plot an embedding. Continue with the guides to [inputs and
initialization](https://jlmelville.github.io/smallvis/package/articles/inputs.md),
[optimization and repeated
runs](https://jlmelville.github.io/smallvis/package/articles/optimization.md),
or [Barnes-Hut and native
acceleration](https://jlmelville.github.io/smallvis/package/articles/barnes-hut.md).

The [article
collection](https://jlmelville.github.io/smallvis/package/articles/index.md)
brings together the theory and experiments behind the package, including
spectral methods, perplexity, initialization, optimization, and
comparisons with other embedding methods. Dated explorations retain
their original results and examples. The [function
reference](https://jlmelville.github.io/smallvis/package/reference/index.md)
documents the current arguments and return values.
