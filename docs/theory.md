# Some Theory

Here are the t-SNE and LargeVis cost functions and gradients, making it a bit
easier to see how they are related. 

A bit of nomenclature first:

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
* $w_{ij}$ in both t-SNE and LargeVis is defined as the Student's t-distribution 
with one degree of freedom (or the Cauchy distribution, if you prefer), 
$w_{ij} = 1 / \left(1 + d_{ij}^2 \right)$. $d_{ij}$ is the Euclidean distance
between point $i$ and $j$ in the output coordinates.
* $\mathbf{y_{i}}$ is the vector of coordinates of point $i$ (with $N$ points 
in total).

## Cost Functions

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
