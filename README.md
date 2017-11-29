# smallvis

An R package for small-scale dimensionality reduction using neighbor embedding 
methods.

A package to compare the cost functions used in [t-Distributed Stochastic Neighbor Embedding](https://lvdmaaten.github.io/tsne/) and 
[LargeVis](https://arxiv.org/abs/1602.00370). The latter method seems to work
very similarly to t-SNE despite not normalizing the output weights. The
[official implementation](https://github.com/lferry007/LargeVis) is focused 
heavily on providing O(N) scaling, and as a result, it's not that convenient for
studying why it works. 

This package is designed to make it easier to compare
the output of embedding using the t-SNE and LargeVis cost functions, but at the
cost of them being O(N^2) in storage and computation costs (and being in pure R).
Unlike LargeVis (and the 
[Barnes-Hut implementation of t-SNE](https://github.com/lvdmaaten/bhtsne)) it is
therefore *not* suitable for large scale visualization. Hence the name smallvis.

## Installing

```R
install.packages("devtools")
devtools::install_github("jlmelville/smallvis")
library(smallvis)
```

## Using

```R
uniq_spec <- unique(iris$Species)
colors <- rainbow(length(uniq_spec))
names(colors) <- uniq_spec
iris_plot <- function(x) {
  plot(x, col = colors[iris$Species])
}

# By default, we use all numeric columns found in a data frame, so you don't need to filter out factor or strings
# set verbose = TRUE to log progress to the console
tsne_iris <- smallvis(iris, perplexity = 25, epoch_callback = iris_plot, verbose = TRUE)

# Default method is t-SNE, use largevis cost function instead
# Also needs a gamma value specified, and not as easy to optimize as t-SNE:
# reduce learning rate (eta) and increase maximum iterations
largevis_iris <- smallvis(iris, method = "largevis", gamma = 1, perplexity = 25, epoch_callback = iris_plot, 
                          eta = 10, max_iter = 5000, verbose = TRUE)

# use (scaled) PCA initialization so embedding is repeatable
tsne_iris_spca <- smallvis(iris, perplexity = 25, epoch_callback = iris_plot, Y_init = "spca")

# return extra information in a list, like with Rtsne
tsne_iris_extra <- smallvis(iris, perplexity = 25, epoch_callback = iris_plot, ret_extra = TRUE)

# more (potentially large and expensive to calculate) return values, but you have to ask for them specifically
tsne_iris_extra_extra <- smallvis(iris, perplexity = 25, epoch_callback = iris_plot,
                              ret_extra = c("P", "Q", "DX", "DY", "X"))

# Repeat embedding 10 times and keep the one with the best cost
tsne_iris_best <- smallvis_rep(nrep = 10, iris, perplexity = 25, ret_extra = TRUE)
plot(tsne_iris_best$Y)
```

## Things To Be Aware Of

* LargeVis requires the use of a `gamma` parameter, which weights the contribution
of attractive and repulsive contributions to the cost function. In the real LargeVis,
it is recommended to set this value to `7`, but this relies on the specifics of
the stochastic gradient descent method. In `smallvis`, this value is very 
dataset dependent: the more data you have, the smaller gamma should be to avoid
over-emphasising repulsive interactions.
* As mentioned above, the real LargeVis uses a classic stochastic gradient 
descent approach with a decaying learning rate. The implementation in `smallvis`
uses the same delta-bar-delta optimization method used in t-SNE. It works well
in my experience, but may require some tuning and seems to require more 
iterations compared to optimizing t-SNE.
* The LargeVis gradient requires quite a large value of epsilon to
avoid division by zero and get results. It's hard-coded to 0.1 in the LargeVis
source code, so I have used the same value by default in smallvis.

## See Also

For large scale visualization in R see:

* The Barnes Hut t-SNE package [Rtsne](https://cran.r-project.org/package=Rtsne)
* The [largeVis](https://cran.r-project.org/package=largeVis) package.

Much of the code here is based on my [fork](https://github.com/jlmelville/rtsne) 
of Justin Donaldson's R package for [t-SNE](https://cran.r-project.org/package=tsne).

## License

[GPLv2 or later](https://www.gnu.org/licenses/gpl-2.0.txt), although any 
LargeVis-specific code (e.g. cost and gradient calculation) can also be 
considered [Apache 2.0](https://www.apache.org/licenses/LICENSE-2.0).

## Some Theory

Here are the t-SNE and LargeVis cost functions and gradients, making it a bit
easier to see how they are related. Some nomenclature:

* $p_{ij}$ is the input (joint) probability of... um, well, it's not totally
obvious. Let's say it's the probability of observing an edge between the two
vertices $i$ and $j$ in the graph representing the neigborhood relationship 
of your input data. These are dependent only on the input data (and the choice 
of perplexity) and so are constant with respect to optimizing the output 
coordinates. If you decide to think of the problem in terms of graph laplacians
and $p_{ij}$ as a normalized affinity, I won't blame you.
* $q_{ij}$ is the equivalent "probability" based on the output coordinates.
* $w_{ij}$ are the weights, aka un-normalized affinities or similarities. 
$q_{ij} = w_{ij} / Z$ where $Z$ is the sum of all the weights: 
$Z = \sum_{ij} w_{ij}$.
* $w_{ij}$ in both t-SNE is the Student's t-distribution with one degree of
freedom (or the Cauchy distribution, if you prefer), 
$w_{ij} = 1 / \left(1 + d_{ij}^2 \right)$. $d_{ij}$ is the Euclidean distance
between point $i$ and $j$ in the output coordinates.
* $\mathbf{y_{i}}$ is the vector of coordinates of point $i$ (with $N$ points 
in total).

### t-SNE

The t-SNE cost function is:

$$
C_{tSNE} = 
\sum_{ij} p_{ij} \ln \frac{p_{ij}}{q_{ij}}
=
\sum_{ij} \left( p_{ij} \ln p_{ij} - p_{ij} \ln q_{ij} \right)
$$
The first term in the sum has no dependence on the output 
coordinates, so is a constant we'll just mark as $C_{P}$. 

Now let's write out $q_{ij}$ as $w_{ij} / Z$:

$$
C_{tSNE} = C_{P} - \sum_{ij} p_{ij} \ln \left( \frac{w_{ij}}{Z} \right) = 
Cp - \sum_{ij} p_{ij} \ln w_{ij} + \sum_{ij} p_{ij} \ln Z
$$
Finally, we'll do some rearranging and re-write $Z$ back to a sum of weights:

$$
C_{tSNE} = 
Cp - \sum_{ij} p_{ij} \ln w_{ij} + \ln Z \sum_{ij} p_{ij} =
Cp - \sum_{ij} p_{ij} \ln w_{ij} + \ln \sum_{ij} w_{ij}
$$

Ignoring, the constant term, we can see that the SNE cost function consists
of an attractive term, where maximizing the $w_{ij}$ (which implies minimizing 
the distances) would minimize $-p_{ij} \ln w_{ij}$; and a repulsive term, where 
minimizing the sum of $w_{ij}$ (and hence maximizing the distances), 
will minimize the log of the sum.

### LargeVis

For computational efficiency, the LargeVis cost function partitions the data 
into a set of nearest neighbors, which only feel an attractive force; and 
everything else, which only feel repulsive forces. But as `smallvis` has 
abandoned such efficiencies, we can write it similarly to the t-SNE cost 
function:

$$
C_{LV} = 
-\sum_{ij} p_{ij} \ln w_{ij} 
-\gamma \sum_{ij} \ln \left( 1 - w_{ij} \right)
$$

If we want, we can expand the right hand sum further to:

$$
C_{LV} = 
-\sum_{ij} p_{ij} \ln w_{ij} 
-\gamma \sum_{ij} \ln w_{ij}
-2 \gamma \sum_{ij} \ln d_{ij}
$$

Don't know if that makes anything much clearer, though.

You can see that in both t-SNE and LargeVis, the attractive terms are identical. 
Any difference in behavior we see must therefore be due to how repulsion 
between points is handled.

## Gradients

Here are the gradients for t-SNE and LargeVis, with the t-SNE one written
slightly differently to how its usually presented to show the sum of attractive
and repulsive forces more clearly and for comparison with LargeVis:

$$
\frac{\partial C_{tSNE}}{\partial \mathbf{y_i}} = 
  4\sum_j^N \left(
    p_{ij} w_{ij}
    -
    q_{ij} w_{ij}
   \right)
   \left(\mathbf{y_i - y_j}\right)
$$

$$
\frac{\partial C_{LV}}{\partial \mathbf{y_i}} = 
  4\sum_j^N \left(
    p_{ij} w_{ij}
    -\frac{\gamma w_{ij}}{d_{ij}^2 + \epsilon}
   \right)
   \left(\mathbf{y_i - y_j}\right)
$$

The $\epsilon$ term in the LargeVis gradient is needed computationally to avoid
division by zero. Results are quite sensitive to changing this value.
In `smallvis` you can see the effect of it via the `lveps` parameter.
