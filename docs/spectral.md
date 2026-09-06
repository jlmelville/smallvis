---
title: "Spectral Methods"
output:
  html_document:
    theme: cosmo
    toc: true
    toc_float:
      collapsed: false
editor_options: 
  markdown: 
    wrap: 100
---

Up: [Documentation Home](https://jlmelville.github.io/smallvis/).

Originally written in December 2017, with some minor corrections and typos over the years. Some
substantial new additions were added in December 2021 (these are called out in the text).

Spectral methods are only tangentially related to `smallvis`, in that a spectral method is available
for initialization (`Y_init = "laplacian"` or `Y_init = "normlaplacian"`) and that the attractive
part of most of the cost functions approximates a spectral method. Various methods are related to
each other, but the nomenclature is not always consistent between papers, and a lot of the more
accessible material seems to be a bit unclear, so I need something to remind myself without having
to work it all out over and over again every time. So it's going here.

There is also a substantial overlap between graph Laplacians and diffusion maps, so this page now
discusses those connections also.

This page is just to remind me of the definitions of matrices and the algorithmic procedures, not
anything to do with theoretical properties. For that, see the further reading section.

## The Usual Preliminaries

We have a matrix of input data, $X$, with $N$ observations of $K$ features.

Now let's represent it as a fully-connected undirected graph, where each object in the dataset is a
vertex, and the strength of the connection between each vertex is given by a weighted edge.

The graph is represented as a matrix, $W$: an $N$ x $N$ matrix of edge weights (or "similarities" or
"affinities": you may also see the matrix written as $A$), where if $w_{ij}$ is large, that means
object $i$ and $j$ are considered similar. $W$ is sometimes also called a kernel matrix and the
function that generated the similarities as a kernel.

Let's start with a Gaussian function of the Euclidean distances (which is a common choice):

$$w_{ij} = \exp\left(-r_{ij}^2 / \sigma\right)$$

where $r_{ij}$ is the Euclidean distance between point $i$ and $j$ and $\sigma$ is a bandwidth
parameter of some kind (you can set it to one and forget it exists if you prefer).

This function is sometimes referred to as a radial basis function (RBF) or the heat kernel.

If we use the same positive bandwidth for every pair and keep all the entries, including the
diagonal ones, this gives a positive-semidefinite (PSD) matrix $W$: all its eigenvalues are
non-negative. This makes life easier for the diffusion maps discussion later. For the graph
Laplacians, we only need $W$ to be:

- Symmetric.
- Contain all non-negative values.

If $N$ is large then it is usual to sparsify $W$ by only keeping the $k$-largest off-diagonal values
in each column. This creates the k-nearest neighbor similarity graph, which must then be
resymmetrized (e.g. by adding its transpose). There are several ways to do this, which give slightly
different graphs.

It would be natural to try a different kernel, let the bandwidth vary from point to point, or use a
sparse graph like those used in t-SNE and UMAP. These are all reasonable choices, but they may lose
the PSD property of our full Gaussian matrix. Density alone doesn't settle this: a dense affinity
matrix can have negative eigenvalues, while kernels can be designed to remain PSD with variable
bandwidths, exact zeros, or both. For examples and ways to combine kernels, see chapter 4 of
[Rasmussen and Williams (PDF)](https://gaussianprocess.org/gpml/chapters/RW4.pdf). We can leave the
consequences until we get to diffusion maps.

If your data is naturally a graph, then you skip all of the above, you already have the data you
need to create $W$, which is now an adjacency matrix. It's likely that in that case $W$ is a sparse
matrix, but the graph Laplacian construction is the same.

Some of the theory did assume a fixed bandwidth across your data, but there is recent work on
variable-bandwidth kernels. See the work by [Berry and Harlim](https://arxiv.org/abs/1406.5064) if
you are of a more theoretical inclination than me. For newer convergence results covering kernelized
k-nearest neighbor graphs and self-tuned affinities, see [Cheng and
co-workers](https://arxiv.org/abs/2410.23212).

### Self-loops: Dealing with the Diagonal

If your $W$ matrix is derived from a graph and you don't have loops (i.e. no edges from a node to
itself), then the diagonal of $W$ will be all zeros. Kernel matrices usually have ones on the
diagonal (an obvious consequence of a Gaussian kernel). What ends up on the diagonal doesn't affect
the relationship between the different graph Laplacians but can affect the numerical values
themselves.

It also affects the PSD guarantee. If we set the entire diagonal of a nonzero symmetric $W$ to zero,
its eigenvalues sum to zero, so they cannot all be non-negative. The random walk matrix $P$
introduced below must then have negative eigenvalues too. This has significance for diffusion maps,
although it doesn't stop us using the graph for Laplacian Eigenmaps.

*December 10 2021*: I recently stumbled upon an example of where this difference lead to some
confusion for me. The [sklearn's SpectralEmbedding
class](https://scikit-learn.org/stable/modules/generated/sklearn.manifold.SpectralEmbedding.html)
generates an affinity matrix (i.e. the diagonal contains all 1s), but to create the graph Laplacian
derived from that affinity matrix, uses [scipy's csgraph.laplacian
function](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csgraph.laplacian.html).
This returns a symmetrized normalized graph Laplacian (see below) calculated under the assumption
that there are all zeros on the diagonal. This doesn't have an enormous numerical effect on the
output, but caused me some bewilderment. *March 7 2026*: At some point the documentation for
`scipy.sparse.csgraph.laplacian` was updated to mention that the diagonal is replaced by zeroes so
no more bewilderment for future readers.

*December 11 2021*: I just noticed that the spectral clustering method of Ng, Jordan and Weiss
explicitly calls for setting the diagonal of the affinity matrix to zero. The Shi and Malik spectral
clustering paper does not mention this (and the weight heat map figures in their paper suggest the
diagonals are all 1s). I'm pretty sure that the diffusion map literature doesn't mention doing
anything special to the diagonal of the affinity matrix either. For the full Gaussian construction
discussed here, keeping the diagonal ones preserves the PSD property.

## The Degree Matrix

The degree matrix is a diagonal matrix where each value in the diagonal is the sum of each edge
associated with a vertex. In other words, sum the rows of $W$ and put that in the diagonal of an $N$
x $N$ matrix. Or sum the columns; $W$ is, after all, symmetric.

$$d_{ii} = \sum_{j} w_{ij}$$

We will be inverting $D$ below, so we need $d_{ii}>0$: if you have isolated vertices, you will need
to remove them or decide how to handle them separately. I will also assume that the graph is
connected. Otherwise, with positive degrees, you get one zero Laplacian eigenvalue per connected
component, so there is more than one trivial eigenvector to discard.

## Some Graph Laplacians

Now that we have $W$ and $D$, we can create some Laplacians, using the naming scheme given by von
Luxburg:

### Unnormalized Graph Laplacian

$$L = D - W$$

Also referred to as the *combinatorial* graph Laplacian by Tremblay and Loukas.

This is the one graph Laplacian in this document which is unaffected by the values on the diagonal.

### Symmetrized Normalized Laplacian

$$L_{sym} = D^{-1/2} L D^{-1/2} = I - D^{-1/2} W D^{-1/2}$$

also sometimes called just the *normalized* Laplacian and referred to by the nomenclature $L_{n}$,
but see below for another normalized Laplacian, which makes this naming ambiguous.

$D^{-1/2} W D^{-1/2}$ is sometimes called the (symmetrized) normalized adjacency matrix.

For $i \ne j$, this normalization can be written entrywise as:

$$\left(L_{sym}\right)_{ij} = -\frac{w_{ij}}{\sqrt{d_{ii} d_{jj}}}$$

On the diagonal,

$$\left(L_{sym}\right)_{ii} = 1 - \frac{w_{ii}}{d_{ii}},$$

which equals 1 only when $w_{ii}=0$.

### Random Walk Normalized Laplacian

$$L_{rw} = D^{-1} L = I - D^{-1} W$$

von Luxburg notes that this is *also* sometimes referred to as the normalized Laplacian, so it's
best to use the longer names von Luxburg uses to avoid confusion. And yes, $D^{-1} W$ is *also*
sometimes referred to as a normalized adjacency matrix.

### Random Walk Transition Matrix

People really like giving $D^{-1} W$ lots of different names:

$$P = D^{-1} W$$

Also called the *Diffusion Operator* in Socher's report. This also means you can write the Random
Walk Normalized Laplacian as:

$$L_{rw} = I - P$$

$P$ is row-normalized (i.e. all rows add up to 1).

It's also pretty easy to see that:

$$L_{sym} = I - D^{1/2} P D^{-1/2}$$

The matrix

$$P_{sym} = D^{1/2} P D^{-1/2} = D^{-1/2} W D^{-1/2}$$

shows up in the discussion of diffusion maps. There doesn't seem to be a fixed name or symbol for
it, so I will go with $P_{sym}$ here in analogy with $L_{sym}$ to indicate it's a symmetrized
version of $P$.

## Eigenvectors

As the "spectral" bit indicates, a lot of eigenvectors are going to be mentioned below. Usually the
eigenvectors are sorted according to the associated eigenvalues.

### Smallest and Largest

Most of the time, when I talk about "largest" and "smallest" eigenvectors, I am referring to the
eigenvectors associated with the largest and smallest eigenvalues, respectively. For graph
Laplacians, the eigenvalues of interest are all non-negative and real, so there shouldn't be any
ambiguity about whether I mean a big negative or positive eigenvalue (none of them are negative
except close to zero due to numerical issues). However, I do discuss the eigenvalues of $P$ a bit as
well, and those *can* be negative.

### "Top" Eigenvectors

The "top" eigenvectors also means the eigenvectors associated with the largest eigenvalues. This is
commonly found in discussions of singular value decomposition, especially when used as part of
principal component analysis, but can also find its way into discussions of spectral clustering and
graph laplacians. The "top" eigenvector is sometimes referred to as the "dominant" eigenvector.

### "First" Eigenvectors

I will also borrow the nomenclature of von Luxburg which refers to the "first" $k$ eigenvectors as
being the eigenvectors associated with the $k$ *smallest* eigenvalues. There is a bit of a clash of
naming here, as you might think that "top" and "first" eigenvectors refer to the same thing, but
they don't. In fact they have opposite meanings. Great.

### Eigenvalues

Both versions of the normalized adjacency matrix $D^{-1} W$ and $D^{-1/2} W D^{-1/2}$, have
eigenvalues that vary between -1 and 1.

For graph Laplacians, the smallest eigenvalue is 0. For the normalized Laplacians $L_{sym}$ and
$L_{rw}$, the maximum value an eigenvalue can attain is 2. *December 9 2021*: For proof, see
[Spectral Graph Theory](http://www.math.ucsd.edu/~fan/research/revised.html), Part 5 of Lemma 1.7 in
Section 1.3, 'Basic facts about the spectrum of a graph'. At the time I write this, chapter 1 is
freely available at the link above.

To keep the notation straight, I'll number the normalized Laplacian eigenvalues from smallest to
largest:

$$0=\lambda_1 \leq \lambda_2 \leq \dots \leq \lambda_N,$$

and use the same index for the corresponding eigenvalue of $P$ or $P_{sym}$:

$$\mu_i = 1-\lambda_i$$

This means $1=\mu_1 \geq \mu_2 \geq \dots \geq \mu_N$. When discussing $P$, "largest" means most
positive unless I explicitly say largest in magnitude. In particular, $P$ and $P_{sym}$ can have
negative eigenvalues: any $\lambda_i>1$ gives $\mu_i<0$, with $\lambda_i=2$ corresponding to
$\mu_i=-1$. This will matter when we get to diffusion maps.

For our full Gaussian matrix with its diagonal retained, the intervals are narrower:
$0\leq\mu_i\leq1$ and $0\leq\lambda_i\leq1$. More generally, this holds whenever $W$ is PSD: the
normalization $D^{-1/2}WD^{-1/2}$ preserves that property, and $P$ has the same eigenvalues as this
symmetric matrix.

## Connections Between Laplacian Eigenvectors

*December 8 2021* While I have made a few typo corrections and clarifications since originally
writing this, I am calling out this brief section because it adds some relationships that may be of
practical interest: it's better to choose the Laplacian format and spectral decomposition method
that meets your robustness and speed needs, and then convert the result into what you want. You
probably aren't writing the software to do the linear algebra yourself so having flexibility over
the choice of package you use for that can be very helpful.

The eigenvalues of $L_{sym}$ and $L_{rw}$ are the same. So if you just need the eigenvalues it
doesn't matter which of those Laplacians you use. We can therefore refer to $\lambda_i$ as the ith
eigenvalue in most cases of interest without having to worry about which of the two Laplacians we
are referring to.

The eigenvectors are also related. If $v_{rw}$ is an eigenvector of $L_{rw}$ and $v_{sym}$ is an
eigenvector of $L_{sym}$, then:

$$v_{rw} = D^{-1/2} v_{sym}$$

As a spoiler for later, the eigenvectors of $P$ are the same as $L_{rw}$ but the eigenvalue ordering
is reversed, because the eigenvalues of $P$ are equal to $1 - \lambda$. So the smallest eigenvalues
of $L_{rw}$ correspond to the largest eigenvalues of $P$.

Further, $P$ and $P_{sym}$ have the same eigenvalues and the relationship between eigenvectors is
the same as that between $L_{rw}$ and $L_{sym}$:

$$v_{P} = D^{-1/2} v_{Psym}$$

I mainly discuss the eigenvectors of $L_{rw}$ below as these are most closely related to Laplacian
Eigenmaps and Diffusion Maps. Assume if you see reference to a vector $v$ without any subscript that
this refers to the eigenvector of $L_{rw}$.

### Scaling the Eigenvectors

Eigenvectors aren't defined to have any particular length. Different software will return the
eigenvectors with different lengths, so you will need to make a decision about their length. The
Laplacian Eigenmaps literature often refers to the smallest eigenvector $v_{rw,1}$ as a vector of
1s, so you could do that. The actual length will then depend on the dimensions of your matrices. The
other obvious choice is to scale all the vectors to unit length.

Because of the relationship between the eigenvectors of $L_{rw}$ and $L_{sym}$ described above, if
you go with all 1s for $L_{rw}$ then the smallest eigenvector of $L_{sym}$ is $D^{1/2}\mathbf{1}$.
Knowing this may have some practical value: most eigenvector software routines let you supply a
guess for at least one of the eigenvectors so this is an easy choice involving numbers you already
have to hand.

Forgetting to scale eigenvector output has been a common cause of temporary panic and wild goose
chasing on my part when writing this document.

## Laplacian Eigenmaps

Solve the generalized eigenvalue problem:

$$Lv = \lambda D v$$

The Laplacian Eigenmap uses the smallest eigenvectors. But not the very smallest eigenvector, $v_1$,
which is constant (we can scale it to be a vector of 1s), and corresponds to an eigenvalue of zero.
So if you want to reduce to two dimensions, use the second-smallest and third-smallest eigenvectors.

Let's do a brief bit of rearranging:

$$Lv = \lambda D v \\
D^{-1} L v = \lambda v \\
\left(D^{-1} D - D^{-1} W \right) v = \lambda v \\
\left(I - P \right) v = \lambda v \\
L_{rw} v = \lambda v$$

So it turns out that the standard eigenvalue problem with $L_{rw}$ will produce the same results as
the generalized eigenvalue problem with $L$ and $D$. A non-generalized eigenvalue problem is
preferable to the generalized problem, at least in R, because generalized problems require
installing the CRAN package [`geigen`](https://cran.r-project.org/package=geigen). You could even
use the eigenvectors of $P$, although you have to bear in mind that the eigenvalues of $P$ differ
from $L_{rw}$. Compared to $L_{rw}$, when using $P$ the order of the eigenvectors are reversed, i.e.
you want those associated with the *largest* eigenvalues, ignoring the uninformative top
eigenvector. $P$ might even be preferable because it's ever-so-slightly less work to calculate than
$L_{rw}$. We'll revisit the relationship between $L_{rw}$ and $P$ when we talk about diffusion maps.

### Output

Now that you have $k$ nontrivial eigenvectors, stack them columnwise to form an $N$ x $k$ matrix.
I'll label it $Y$.

$$Y = \left[v_2 | v_3 | \dots | v_{k+1} \right]$$

where, as noted above, we have discarded $v_1$. The rows of that matrix are the coordinates of the
graph vertices in the reduced dimension, i.e. the ith row of the 2D Laplacian Eigenmap representing
vertex i would be:

$$y_i = \left(v_{i,2}, v_{i,3} \right)$$

To fix the axis scales, the usual Laplacian Eigenmaps normalization is $Y^{\mathsf{T}}DY=I$. If you
start with orthogonal unit-length eigenvectors of $L_{sym}$, multiplying them by $D^{-1/2}$ already
takes care of this.

### The Connection with Locally Linear Embedding

The Laplacian Eigenmap paper demonstrates a connection between LE and LLE, in that LLE is
approximately computing the eigenvectors of $L^2$, which has the same eigenvectors as $L$ (and the
square of the eigenvalues).

## Spectral Clustering and Normalization

von Luxburg describes three different spectral clustering algorithms, which all involve forming a
Laplacian matrix, calculating some eigenvectors, and then forming the reduced-dimension matrix $Y$
from column-stacking the eigenvectors.

1.  Un-normalized: compute the first $k$ eigenvectors of $L$.
2.  Normalized ([Shi and Malik](https://ieeexplore.ieee.org/document/868688)): compute the first $k$
    *generalized* eigenvectors of $L$. This is just what Laplacian Eigenmaps do, so from the above
    discussion we know that it is equivalent to computing the first $k$ eigenvectors of $L_{rw}$
    (hence justifying the term "normalized").
3.  Normalized ([Ng, Jordan and
    Weiss](https://papers.nips.cc/paper/2001/hash/801272ee79cfde7fa5960571fee36b9b-Abstract.html)):
    compute the first $k$ eigenvectors of $L_{sym}$. This version requires an additional row
    normalization step of the output matrix, $Y$, before you can do clustering: normalize the rows
    so each row has length 1 (normalization to unit $l_2$ norm).

After some additional theoretical discussions, von Luxburg concludes that clustering on the
un-normalized graph Laplacian has some undesirable properties, so you definitely want to use one of
the normalized Laplacians for clustering. Of the two normalized Laplacians, the Shi-Malik approach
(once cast in terms of using $L_{rw}$) is the least effort. A slight downside to $L_{rw}$ is that
unlike $L$ and $L_{sym}$, it is not symmetric, and symmetric matrices usually have access to
slightly more methods (or more efficient methods) for solving the eigenproblem.

Conclusion: use $L_{rw}$ (Laplacian Eigenmaps). *December 8 2021* Probably I should have said: use
the *eigenvectors* of $L_{rw}$, but given the straightforward conversion between the eigenvectors of
$L_{sym}$ and those of $L_{rw}$, you don't need to form $L_{rw}$ directly for computational
purposes. In fact, see the 'Using Truncated SVD' section below for why you might want to operate on
a matrix related to $L_{sym}$ instead of $L_{rw}$.

## Diffusion Maps

From the perspective of what sort of matrices you build and what you do to them, diffusion maps
share a lot of that machinery with Laplacian Eigenmaps. But as we get into it, we'll see that their
application has a very different emphasis. Not much of it is germane to spectral embedding or
spectral clustering, so feel free to skip to the "Using Truncated SVD" or "Repeated Eigendirections
and Other Problems" section.

Rather than visualize two or three eigenvectors, the motivation is to use distances on the map to
describe a random walk. Two vertices should be close if walks started there have similar
probabilities of ending up elsewhere on the graph.

The connection to our affinity matrix is straightforward: $p_{ij}=w_{ij}/d_{ii}$ is the probability
of moving from vertex $i$ to vertex $j$ in one step. So choosing the affinities also chooses the
walk. A Gaussian favors short moves; a local bandwidth changes what counts as short at each point; a
neighborhood mask rules out some moves entirely.

Which eigenvectors are selected and how the eigenvalues scale them do become more important for
diffusion maps than for spectral embedding and spectral clustering. For now, we'll assume we are
working with the full fixed-bandwidth Gaussian with its diagonal ones, so we know the eigenvalues of
$P$ are non-negative.

The calculation starts by solving the eigenproblem for $P$ rather than $L_{rw}$:

$$P v = \mu v$$

The eigenvectors are the same whether you use $L_{rw}$ or $P$, but the eigenvalues are different, so
I am using $\mu$ instead of $\lambda$ to differentiate them from the eigenvalues associated with
Laplacian Eigenmaps and the usual spectral clustering algorithms. As it happens, the eigenvalues are
related by:

$$\mu = 1 - \lambda$$

With our Gaussian construction, we keep the eigenvectors with the largest $\mu$, leaving out the
constant eigenvector. These are the same ones we used for Laplacian Eigenmaps. Once we've seen how
the distances work, we can come back to what changes when a different affinity construction gives
negative eigenvalues.

### Diffusion Time and Coordinates

Once you have selected the eigenvectors, the eigenvalues are used to scale them when forming the $Y$
matrix. We can use $v_2$ and $v_3$ just as in Laplacian Eigenmaps. Here's what one row of $Y$ would
look like in this 2D diffusion map case:

$$y_i = \left(\mu_2 v_{i,2}, \mu_3 v_{i,3}\right)$$

where $v_1$ is the uninteresting top eigenvector of all 1s (and $\mu_1 = 1$).

Where does the diffusion come in? Because you can evaluate the transition probabilities after a
positive integer number of steps $t$ with the iterated matrix $P^{t}$, you can get a sense of the
geometry of the data at different scales by seeing how the probability changes over time. And
there's not even that much extra work to do: the eigenvectors of the iterated matrix are the same as
the original matrix $P$, and the eigenvalues are $\mu_i^t$.

For a given value of $t$, the 2D diffusion map at time $t$ is therefore:

$$y_i = \left(\mu_2^t v_{i,2}, \mu_3^t v_{i,3}\right)$$

If you want an even more obvious connection to Laplacian Eigenmaps:

$$y_i = \left[\left(1 - \lambda_{2}\right)^{t} v_{i,2}, \left(1 - \lambda_{3}\right)^{t} v_{i,3}\right]$$

It's now easier to see the effect of $t$: it relatively upweights the contribution of earlier
eigenvectors in determining the distance between points on the map. So larger diffusion times
increasingly emphasize the slowest, most global modes.

### Diffusion Distances

I was a bit dismissive about scaling the eigenvectors when we were talking about Laplacian
Eigenmaps. With diffusion maps we want to use the coordinates to calculate diffusion distances, so
we need to be a bit more principled. We will now briefly get a little bit into the weeds about the
stationary distribution of the random walk: a set of probabilities for being at each vertex which
stays the same after another step. Writing these probabilities as a column vector $\pi$, this means
$\pi^{\mathsf{T}}P=\pi^{\mathsf{T}}$.

For our connected undirected graph, it can be shown that

$$\pi_i=\frac{d_{ii}}{\sum_j d_{jj}}.$$

Remember that $d_{ii}=\sum_j w_{ij}$ is the total edge weight at vertex $i$. So we are normalizing
the vertex degrees to sum to one. $p_{ii}$ is the probability of staying at vertex $i$ for one step.
You can check the stationary distribution using the symmetry of $W$:

$$\sum_i\pi_i p_{ij}
=\frac{\sum_i w_{ij}}{\sum_k d_{kk}}
=\frac{d_{jj}}{\sum_k d_{kk}}=\pi_j.$$

Put the $\pi_i$ on the diagonal of a matrix $\Pi=D/\sum_j d_{jj}$. The normalization we need is
$V^{\mathsf{T}}\Pi V=I$, where the columns of $V$ are the eigenvectors of $P$. This is much like the
Laplacian Eigenmaps normalization. If $U$ contains orthogonal unit-length eigenvectors of $P_{sym}$,
the conversion is

$$V=\Pi^{-1/2}U=\sqrt{\sum_j d_{jj}}\,D^{-1/2}U.$$

This also makes the trivial eigenvector a vector of all 1s, as promised.

The diffusion distance compares the transition probabilities from two starting vertices after $t$
steps, weighted by the stationary probabilities:

$$D_t^2(i,j)=\sum_{k=1}^{N}
\frac{\left((P^t)_{ik}-(P^t)_{jk}\right)^2}{\pi_k}.$$

Having done the normalization above, [Nadler and co-workers](https://arxiv.org/abs/math/0506090)
showed that this can also be written in terms of the eigenvectors:
$$D_t^2(i,j)=\sum_{\ell=2}^{N}\mu_\ell^{2t}
\left(v_{i,\ell}-v_{j,\ell}\right)^2.$$

So with all the nontrivial eigenvectors retained, ordinary Euclidean distances on the map are
exactly the diffusion distances. This is why the normalization matters: rescaling the axes changes
those distances.

### Why Do All This?

With our Gaussian construction, Diffusion Maps are just like Laplacian Eigenmaps but with some
stretching of the axes. Big deal! And if all we wanted was a visualization, that would be a fair
summary. But as we have seen above, scaling is what gives the distances their random walk
interpretation. We need both the eigenvectors and their weights, not just a couple of directions to
plot.

It's really more about using the distances in the embedded space, which you can think of as having
been "smoothed" or "denoised". The diffusion distance effectively averages over multiple paths
between points, which makes it more robust to noise. This may be better than trying to measure the
distance between points via geodesics using shortest paths in the k-nearest neighbor graph in the
original space. You can also tune $t$ to increase emphasis from local to more global structure. Once
you have these distances, you can re-use them to e.g. create a new k-nearest neighbor graph and so
on.

The trick is to make sure you include enough dimensions for these distances to be meaningful.
Probably two or three dimensions isn't enough. I'll come back to this point in the "A Truncated SVD
recipe" section below.

### What If the Affinities Give Negative Eigenvalues?

There are perfectly reasonable choices of kernel to build an affinity matrix that lose the PSD
guarantee. For example, maybe you want to use the self-tuning kernel for spectral clustering
described by [Zelnik-Manor and
Perona](https://proceedings.neurips.cc/paper/2004/hash/40173ea48d9567f1f393b20c855bb40b-Abstract.html).
Maybe you want a sparse kernel. Maybe you have decided to remove all self-loops and hence the
affinity matrix has a diagonal of all zeros. All of this can lead to PSD being lost. As an aside,
there is research in the field of Gaussian Processes (GP) about the task of designing kernels that
are PSD, sparse and with variable bandwidths, but it seems non-trivial.

These choices still define a valid random walk, provided $W$ is symmetric and non-negative and the
degrees are positive. If the resulting $W$ is no longer PSD, then $P$ has negative eigenvalues too.
At this point we can no longer assume that diffusion maps select the same eigenvectors as Laplacian
Eigenmaps.

The distance formula tells us what to do: each contribution is weighted by $\mu_\ell^{2t}$, so
select the largest eigenvalue *magnitudes*. For example, a Laplacian eigenvalue $\lambda=1.99$ gives
$\mu=-0.99$. Laplacian Eigenmaps would put this eigenvector near the back of the queue, but its
contribution to diffusion distances decays just as slowly as that of a mode with $\mu=0.99$.

For an example where a negative eigenvalue matters, consider a bipartite graph: the vertices fall
into two sets, with affinities only between the sets. A walk switches sets at every step. Here $P$
has an eigenvalue of $-1$ associated with that alternation, and the diffusion distances retain the
distinction between the two sets.

### Lazy Walks

Sometimes when you have negative eigenvalues, you don't want diffusion distances to be so sensitive
to their contribution. In that case you could just use $P$ as normal and order the eigenvalues
algebraically, back to the Laplacian-Eigenmap-like ordering. But although that pushes those negative
eigenvalues down the list, if you retain enough eigenvectors they will eventually start showing up
with the same contribution to the distances. Another option is a *lazy walk*: at each step, toss a
fair coin to decide whether to stay put or take a step using $P$. Its matrix is

$$P_{lazy}=\frac{I+P}{2}.$$

This keeps the eigenvectors and their normalization, but changes each eigenvalue to $(1+\mu_i)/2$,
which is between 0 and 1. The coordinate weights are now $((1+\mu_i)/2)^t$, so the largest weights
select the same directions as Laplacian Eigenmaps. In particular, eigenvalues near $-1$, which had
large weights in the original walk, now have weights near zero. There are more details on lazy
walks in diffusion maps in section 28.4 of [Matthew Hirn's lecture notes
(PDF)](https://matthewhirn.com/wp-content/uploads/2021/04/cmse890_spring2021_lecture21.pdf).

### Book-keeper Beware!

When checking a diffusion maps implementation, look at both the eigenvalue ordering and the
coordinate weights. For the original walk, use the largest $|\mu_i|$ and weights $\mu_i^t$; for the
lazy walk, use the largest $\mu_i$ and weights $((1+\mu_i)/2)^t$. I'll use the original walk with an
integer number of steps $t$ below unless stated otherwise.

### Anisotropic Diffusion

And actually, there's a bit more to them than even that. You can also define different diffusion
operators. At this point, sadly, I've never actually seen the procedure written down unambiguously.
Both the Socher report and Wikipedia page seem to leave out some steps. Here's my attempt. I'll try
and stick with the notation I've seen used elsewhere, which is truly unfortunate.

*December 8 2021*: The [pyDiffMap](https://github.com/DiffusionMapsAcademics/pyDiffMap) package
seems to carry out a procedure very close to what follows, so this should be reasonably reliable.

*October 12 2025*: Having spent more time with the Coifman and Lafon diffusion maps paper, I have
corrected the procedure below. The only change is that originally I used the degree matrix $D$ to
form the symmetrized diffusion operator matrix and I should have used the degree matrix formed from
the alpha-normalized kernel matrix. To be generous to my past self, I later mentioned this making
conceptually more sense to me, but it's actually just the correct thing to do.

The full procedure is something like:

- Form $W$ and $D$ as usual.
- Specify $\alpha$, the anisotropic diffusion parameter, a value between 0 and 1.
- Normalize $W$ according to $\alpha$:

$$W^{\left( \alpha \right)} = D^{-\alpha} W D^{-\alpha}$$

- Form a new diagonal degree matrix, $D^{\left( \alpha \right)}$, based on the new kernel matrix,
  $W^{\left( \alpha \right)}$:

$$d^{\left( \alpha \right)}_{ii} = \sum_{j} w^{\left(\alpha\right)}_{ij}$$

- Form a new diffusion operator using the new kernel matrix and the inverse of the new degree matrix
  (notation is not great here):

$$P^{\left( \alpha \right)} = \left(D^{\left( \alpha \right)}\right)^{-1} W^{\left( \alpha \right)}$$

- Generate the diffusion map using $P^\left( \alpha \right)$, and scale the eigenvectors with the
  eigenvalues as described above.

- If you want to use a symmetrized version of $P^{\left( \alpha \right)}$, then form:

$$P^{\left( \alpha \right)}_{sym} = \left(D^{\left(\alpha\right)}\right)^{1/2} P^{\left( \alpha \right)} \left(D^{\left(\alpha\right)}\right)^{-1/2}$$

and remember to convert the eigenvectors back via:

$$v_{P\left(\alpha\right)} = \left(D^{\left(\alpha\right)}\right)^{-1/2} v_{P\left(\alpha\right)sym}$$

It's particularly lamentable that you have to deal with reading both
$\left(D^{\left( \alpha \right)}\right)^{-1}$ (the inverse of the diagonal matrix
$D^{\left( \alpha \right)}$) *and* the entirely different $D^{-\alpha}$ (invert the diagonal matrix
$D$, then raise the resulting diagonal values to the power of $\alpha$). For the stationary
probabilities used in the normalization above, use the new degrees $d^{\left(\alpha\right)}_{ii}$
too.

This normalization preserves PSD when we start with a PSD $W$, because it scales matching rows and
columns by positive factors. So with our full Gaussian kernel, changing $\alpha$ doesn't introduce
negative eigenvalues.

### Choice of $\alpha$

*October 12 2025*: A new section that goes into a bit more detail about why being able to choose
$\alpha$ is useful.

There are three canonical choices for $\alpha$ which correspond to interpretations of partial
differential equations. These are limiting results, so they come with assumptions about smoothness,
sampling and how the bandwidth changes as the sample size grows:

- $\alpha = 0$, you get back the diffusion map based on the random walk-style diffusion operator
  (and Laplacian Eigenmaps). This includes the *density* of sampling in the embedding. Probably
  useful for clustering.
- $\alpha = 1$ approximates the Laplace-Beltrami operator. This "corrects" for the density of
  sampling, so in the limit the embedding reflects the *geometry* of the underlying manifold without
  the effect of sampling density.
- $\alpha = 0.5$ gives a weighted, Fokker-Planck-type operator in the limit, under an appropriate
  sampling model.

In terms of *why* you want to do this, you may or may not want the effect of density to affect the
embedding. Imagine your manifold is a nice simple 2D sheet embedded in a higher-dimensional space,
but for some reason your data is sampled from one part of the sheet more densely than another, then
$\alpha = 0$ would reflect that density difference in the embedding, causing a distortion.
$\alpha = 1$ would show the geometry of the sheet without the effect of the sampling density. Values
of $\alpha$ between 0 and 1 give a compromise between these two extremes.

What about the Fokker-Planck business? This is a rather specific case where you have data that is
generated by some dynamic process. Imagine a molecular dynamics simulation of something like a
protein where each item in the data is a "snapshot" of the protein structure at a particular time.
The protein will spend most of its time in a few low-energy states, but occasionally jump between
those states via transition states which are much higher in energy and hence not sampled as often.
So a density-based approach is going to capture the low-energy states well but not the transition
states. A geometry-based approach will better capture the low energy states joined by these
"valley-like" transition states but you would lose the knowledge that the protein spends most of its
time in the low-energy states. Intermediate values of $\alpha$ give a compromise. What would be
special about using $\alpha = 0.5$ is that you could have a more physical interpretation of the
distances between points in the embedding space, e.g. the kinetics of the transitions between states
can be revealed by the random walk. But that requires the random walk to represent the actual
dynamics, or to be calibrated against them. Geometric proximity alone doesn't tell you how quickly
the protein moves between states.

## Using Truncated SVD

*December 8 2021*. For symmetric graph Laplacians, it should be possible to use a truncated SVD to
get the eigenvectors and eigenvalues. In terms of packages and routines available, it's always
seemed to me that there are a lot more options for SVD than for eigenvalue problems.

For this approach to be feasible for large matrices, you still need access to both a truncated SVD
routine (i.e. one that doesn't need to calculate all the singular vectors at once) and one that
works with sparse matrices. It's likely you won't want to code all that up from scratch. And it's
probably not going to be as fast in all cases as a dedicated eigenvalue library.

But if you find yourself in a situation where you do have access to fast SVD routines, but not to
eigenvalue problem solvers, this might be worth considering.

### Connection Between SVD and Spectral Decomposition

For symmetric real positive semidefinite matrices like $L$ and $L_{sym}$, the singular vectors can
be chosen to be the same as the eigenvectors. Also, the singular values and the eigenvalues
coincide.

For a symmetric matrix with negative eigenvalues, the singular values are the *absolute* values of
the eigenvalues. That changes the ordering we need for Laplacian Eigenmaps, but gives exactly the
magnitude ordering we need for diffusion distances. So truncated SVD can be useful for both; we just
need to choose the right matrix. Let's deal with Laplacian Eigenmaps first.

Can we use this to our advantage? To a certain extent, yes. An SVD of $L_{sym}$ can replace the
spectral decomposition. However, as we want the eigenvectors with the smallest eigenvalues of
$L_{sym}$, the equivalent singular vectors are the bottom singular vectors, whereas the truncated
SVD routines we want to use target the top ones.

It would be preferable to use the top singular vectors rather than the bottom singular vectors. Then
we can make use of truncated SVD approaches like [irlba](https://cran.r-project.org/package=irlba)
which run in far less time and memory and hence scale to larger matrices.

Something like the $P$ matrix would be ideal because we want the top eigenvalues from that matrix. A
couple of minor problems to overcome though:

1.  $P$ isn't symmetric, an obvious deal-breaker. However, $P_{sym}$ *is* symmetric, has the same
    eigenvalues and it's easy to convert the eigenvectors.

2.  The eigenvalues of $P$ and $P_{sym}$ range from -1 to 1. This is a slightly subtler problem.
    Singular values are always non-negative, i.e. we lose track of the sign of the eigenvalues and
    only get their magnitudes through SVD. Unfortunately, this ruins the Laplacian ordering of the
    singular vectors. For $P$, we want the most positive eigenvalues here. We can get around that by
    taking advantage of the fact that shifting a matrix by $cI$ ($c$ being a scalar value) shifts
    the eigenvalues by $c$. The eigenvalues of $cI + P_{sym}$ are therefore $c + \mu$. Picking $c=1$
    means the eigenvalues of $I + P_{sym}$ are in the range of 0-2, which is safe for us to run SVD
    on: the matrix is positive semidefinite, so the singular values now have the order we want.
    Another way to write this is

    $$I+P_{sym}=2I-L_{sym}.$$

3.  Ordering the vectors is now very confusing: for graph Laplacians, we were talking about the
    $k$th eigenvector as that associated with the $k$th *smallest* eigenvalue, but when using $P$,
    it's actually the $k$th *largest* eigenvalue. And now we have brought SVD into the mix where the
    singular vectors are ordered by decreasing singular value: the $k$th singular vector is that
    associated with the $k$th *largest* singular value.

### A Truncated SVD recipe

*December 30 2021*: Thanks to [a reddit comment on spectral
embeddings](https://www.reddit.com/r/MachineLearning/comments/rrrjrz/comment/hqirwix) I have
discovered that several network embedding methods (e.g. [NetSMF](https://arxiv.org/abs/1906.11156))
use truncated SVD to get eigenvectors of graph Laplacians, similar to what is described below.

Here's a procedure for calculating the Laplacian eigenvectors via truncated SVD. As noted above,
indexing can get confusing, but it *is* consistent in the sense that if you stick with the
convention of ordering eigenvectors from smallest to largest, and ordering singular vectors from
largest to smallest, the $k$th eigenvector and $k$th singular vector are the same (and the
eigenvalues and singular values are converted easily).

1.  Form $W$ and $D$.
    - If you want an anisotropic diffusion map, form $W^{\left( \alpha \right)}$ and its degree
      matrix $D^{\left( \alpha \right)}$.
2.  Form the symmetric matrix $I + D^{-1/2} W D^{-1/2} = I + P_{sym}=2I-L_{sym}$, or
    $I + P^{\left(\alpha\right)}_{sym}$ for diffusion maps.
3.  Via truncated SVD find the top $k + 1$ singular vectors.
    - The top singular vectors correspond to the smallest eigenvectors of $L_{sym}$.
    - To put it another way: the kth *largest* singular vector is the same as the kth *smallest*
      eigenvector.
    - The kth *smallest* eigenvalue $\lambda_k$ and the kth *largest* singular value, $d_k$ are
      related by $\lambda_{k} = 2 - d_{k}$.
    - Just as the case with the *smallest* eigenvector, the singular vector associated with the
      *largest* singular value is the trivial (degree-weighted) vector.
4.  If you want the eigenvectors of $L_{rw}$, convert in the normal way, i.e.
    $v_{rw} = D^{-1/2} v_{sym}$.
5.  The diffusion map eigenvalues are $\mu_{k} = 1 - \lambda_{k} = d_{k} - 1$.

This recipe also selects the diffusion map modes when $P$ has non-negative eigenvalues. For the lazy
walk, use $d_k/2$ as the eigenvalues instead, giving coordinate weights $(d_k/2)^t$. In either case,
convert the vectors using $\Pi^{-1/2}$ for the diffusion-distance normalization described above.

#### Truncated SVD for Diffusion Distances

For the original walk, there is an even simpler option: run truncated SVD on $P_{sym}$ without the
identity shift. This works whether or not $W$ is PSD. Keep the leading left singular vectors
$u_\ell$, with singular values $s_\ell$, and form coordinates

$$y_{i,\ell}=s_\ell^t\frac{u_{i,\ell}}{\sqrt{\pi_i}}.$$

The signs of the eigenvalues don't matter if we just want the diffusion distances.

There is one wrinkle: eigenvalues with opposite signs but the same magnitude can give mixed singular
vectors. They have the same distance weight, so retaining the whole tied subspace still gives the
same distances. If we cut through a tie, the truncated map is not unique.

What about the constant eigenvector we normally leave out? It contributes nothing to distances, so
we can leave it in for this calculation. We may carry an extra coordinate, but we save ourselves the
trouble of identifying and removing it from the SVD.

So we can use $I+P_{sym}$ for the Laplacian ordering, and $P_{sym}$ for the original walk's
magnitude ordering.

### How Many Coordinates?

If you use this sort of embedding as initialization for something like UMAP or t-SNE, then you
always have a fixed number of dimensions to work with (very likely 2 or 3). But if you decide to
work with diffusion maps more directly and want the diffusion distances to be as useful as possible
you will have to decide how many dimensions to keep.

Probably something in the style of the "fraction of variance explained" or thresholding as used in
PCA is appropriate, but using the eigenvalues rather than the variance. But I would argue that the
correct quantity to look at is the square of each eigenvalue because the distance calculation in the
embedded space at $t=1$ is:

$$d(i,j)= \left( \sum_{\ell=2}^{N} \mu_\ell^2 (v_{i,\ell} - v_{j,\ell})^2 \right)^{1/2}$$

so the contribution of each dimension to the squared distance is weighted by the square of the
eigenvalue. If $S$ is the set of nontrivial eigenvectors we keep, I would suggest choosing it so
that:

$$\frac{\sum_{\ell\in S} \mu_\ell^2}{\sum_{\ell=2}^{N} \mu_\ell^2} \geq \text{threshold}$$

Here we want the largest magnitudes, which might include negative eigenvalues if the walk isn't
lazy. This measures how much of the average squared diffusion distance we retain, weighting each
pair of starting vertices by $\pi_i\pi_j$, in the spirit of explained variance.

But we need to know how to calculate the total sum of the squares of the eigenvalues. We know that
the sum of the eigenvalues of a square matrix is the trace, but alas we need the sum of the squares
of the eigenvalues. We don't want to be calculating the square explicitly, but fortunately it's also
true that the squared Frobenius norm of a symmetric matrix equals the trace of the matrix's square,
i.e. we can just sum the squares of the elements of $P_{sym}$.

*But* we also have to take into account that we always discard the trivial eigenvector/eigenvalue,
and we know that $\mu_1 = 1$. So the denominator is $\|P_{sym}\|_F^2-1$. Equivalently, because we
are working with the shifted version $I + P_{sym}$, we can use:

$$\|I + P_{sym}\|_F^2 - N - 2\operatorname{tr}(P_{sym}) - 1.$$

As a further wrinkle, the trace term is zero only when there are no self-loops. If you are working
with a dense diffusion kernel where $w_{ii}=1$, you need to keep it.

For integer $t>1$, replace each $\mu_\ell^2$ in the fraction by $\mu_\ell^{2t}$. The denominator
becomes

$$\sum_{\ell=2}^{N}\mu_\ell^{2t}
=\operatorname{tr}\left(P_{sym}^{2t}\right)-1.$$

This looks more expensive, but we don't have to form $P_{sym}^{2t}$. There is a standard randomized
trace estimator due to [Hutchinson](https://doi.org/10.1080/03610919008812866); section 1.1 of the
paper by [Meyer and co-workers](https://arxiv.org/abs/2010.09649) gives a recent explanation.
Applied to our matrix, the recipe is:

1.  Generate a vector $z$ of length $N$, choosing each entry independently to be 1 or -1 with equal
    probability. These entries have what is called a Rademacher distribution.
2.  Multiply $z$ by $P_{sym}$ a total of $t$ times, then take the squared length of the result.
3.  Repeat with several independently generated vectors and average the results.

With $s$ random vectors, this gives

$$\operatorname{tr}\left(P_{sym}^{2t}\right)
\approx\frac{1}{s}\sum_{r=1}^{s}\left\|P_{sym}^{t}z_r\right\|^2.$$

Then subtract 1 for the trivial eigenvalue. The squared length works because $P_{sym}$ is symmetric,
so $z^{\mathsf{T}}P_{sym}^{2t}z$ equals $\|P_{sym}^t z\|^2$. Each vector costs $t$ sparse
matrix-vector products, so you can trade accuracy against the number of vectors you use.

As $t$ increases, the contributions tend to concentrate in fewer eigenvectors. If you kept the
largest magnitudes at $t=1$, that choice remains conservative for this fraction as you increase $t$.
So exploring larger diffusion times needn't require further eigenvector calculations.

All of the above should still apply to the anisotropic case, using the
$P^{\left(\alpha\right)}_{sym}$ matrix in place of $P_{sym}$ (the trivial eigenvalue remains 1).

### Truncated SVD in Practice

Is the truncated SVD approach worth doing, especially given that none of the diffusion map packages
I looked at in R or Python use SVD directly?

When comparing this approach using [irlba](https://cran.r-project.org/package=irlba) versus getting
the eigenvalues more directly via [RSpectra](https://cran.r-project.org/package=RSpectra), I didn't
notice any slowdown. However this was in the context of initializing a
[UMAP](https://github.com/lmcinnes/umap) embedding in the R package
[uwot](https://cran.r-project.org/package=uwot), and the spectral decomposition was not a noticeable
computational bottleneck in the first place. The main advantage for me would be that `uwot` already
uses `irlba` for PCA in various places, so I would be able to remove `RSpectra` as a dependency.
This is a decision that has already been reached independently in
[umappp](https://github.com/LTLA/umappp/pull/4) (another UMAP implementation, this one in C++).

I was able to find one dataset where the truncated SVD approach was slower than using `RSpectra`:
embedding a 1D line from a 3D to 2D: i.e. a dataset with 3 columns: one increasing in value from `1`
to `N`, the other two columns being all zero. `irlba` was much slower in this case. This experiment
was inspired by a [bug report in the UMAP project](https://github.com/lmcinnes/umap/issues/360), but
may not be representative of real-world data.

There could be some more advanced uses of spectral clustering where SVD is the best choice. For
example, in 2001 [Inderjit Dhillon](https://dl.acm.org/doi/10.1145/502512.502550) published a paper on
bipartite spectral graph clustering, where SVD is applied to the normalized rectangular affinity
matrix between the two sets. If the sets have $m$ and $n$ vertices, this is an $m$ by $n$ matrix,
rather than the $(m+n)$ by $(m+n)$ matrix needed for the full graph eigenvalue problem.

## Repeated Eigendirections and Other Problems

In the previous section I mentioned getting the Laplacian Eigenmap for a 1D line embedded in 3D.
Assuming you get a good converged result, the output using the first two eigenvectors is a parabola.
This is a generic problem when there is a high "aspect ratio" in a dataset, i.e. the manifold
extends much more in one direction than another. Successive eigenvectors will contain information
about the same coordinate. This problem was described as "repeated eigendirections" by [Gerber and
co-workers](https://dl.acm.org/doi/abs/10.1145/1273496.1273532) and is also discussed at length by
[Goldberg and co-workers](https://arxiv.org/abs/0806.2646).

Apart from repeated eigendirections, closely spaced eigenvalues can make individual eigenvectors
sensitive to small changes in the graph. This matters especially when the embedding keeps some of
these directions and discards others: rotating a complete retained block doesn't itself distort
an unweighted Euclidean embedding. Distortions can also occur near boundaries and holes, leaving
the resulting embedding looking twisted.

For more on this, see the discussion by [Kohli and
co-workers](https://www.jmlr.org/papers/v22/21-0131.html) and especially the references they point
to (under 'Laplacian Eigenmaps' in section 1.3). At this point I'd love to say "and here's the easy
solution that's been discovered", but that doesn't seem to be the case, so see the above papers and
the references therein for more suggested fixes.

*2 January 2025*: here's a nice [review of manifold
learning](https://doi.org/10.1146/annurev-statistics-040522-115238) which devotes an entire section
to the (surprisingly limited) literature on repeated eigenvectors. No easy fixes are presented
though.

## The Kernel PCA Connection

Kernel PCA, introduced by [Schölkopf and co-workers](https://doi.org/10.1162%2F089976698300017467),
describes quite a similar process to everything described above: you create a square affinity matrix
based on the kernel function (usually called the Kernel matrix, Gram matrix or Gramian matrix), but
instead of forming a graph Laplacian matrix from it, do SVD on the kernel matrix directly. This
relies on the "kernel trick" (I don't actually know who first coined that phrase): the value of the
kernel matrix element $w_{ij}$ can be seen as the result of transforming $\mathbf{x_i}$ and
$\mathbf{x_j}$ into some high-dimensional space (this mapping function is usually labelled as
$\Phi$), and then taking the dot product:

$$w_{ij}=k(\mathbf{x_i},\mathbf{x_j})
=\langle\Phi(\mathbf{x_i}),\Phi(\mathbf{x_j})\rangle.$$

Here $\Phi$ is the feature map, and $k$ is the positive-definite kernel: the function which gives us
the dot product without having to calculate the mapped points. This makes $W$ symmetric and positive
semidefinite. The dot products themselves can be negative, just as they can in the original space.
Hence you don't need to ever actually map the data into the high dimensional space or to even know
what $\Phi$ is.

A difference between the spectral methods and kernel PCA is the normalization of $W$. Doing PCA
requires mean-centered data and in kernel PCA this means the transformed data should also be
mean-centered. But the whole point of the kernel trick is to avoid having to actually form the
transformed data. Instead, the kernel matrix is double centered:
$W_{norm} = \left(I - \frac{1}{N}\mathbf{1} \right) W \left(I - \frac{1}{N}\mathbf{1} \right)$ where
$\mathbf{1}$ is an $N$ by $N$ matrix of all 1s. This normalization results in the row and column
means all being zero.

For more on this, [Bengio and co-workers
(PDF)](http://www.iro.umontreal.ca/~lisa/pointeurs/TR1232.pdf) have a technical report connecting
spectral clustering with kernel PCA and [Ham and
co-workers](https://dl.acm.org/doi/10.1145/1015330.1015417) describe how graph Laplacians can
themselves be considered kernels.

Standard linear PCA fits into Kernel PCA by using the "linear kernel", i.e. the dot product of the
input vectors: you can get the principal components or loadings or whatever you are looking for
whether you use the scatter/covariance matrix ($X'X$) or the Gram matrix ($XX'$), subject to some
scaling of eigenvalues here or a matrix multiplication with $X$ there. Here $X$ is the input data
matrix. If you are doing PCA, you usually need to center your input data anyway, and the
double-centering that kernel PCA does has no effect on the eigendecomposition.

So if the linear kernel is good enough for PCA, is it good choice for an affinity matrix for
spectral methods? Unfortunately not because there's nothing to stop a dot product from being
negative. That's fine for kernel PCA, as noted above, but attempting to construct a graph Laplacian
matrix from negative weights won't give a matrix with the properties we need. There seems to be a
small amount of literature on graph Laplacians with negative edge weights but it doesn't seem like
something I would get very excited about at the moment.

Also be aware that k-nearest neighbor sparsification can destroy positive semidefiniteness, so a
sparsified graph affinity matrix may no longer be suitable for kernel PCA.

## Further Reading

I try to link to official DOI URLs and the like where possible, and not post links to
copyright-busting PDFs. The tutorials and reports by von Luxburg, Horaud, and Socher are good places
to start.

- The main reference on the properties of graph Laplacians is the monograph [Spectral Graph
  Theory](http://www.math.ucsd.edu/~fan/research/revised.html) by Fan Chung. The actual amount of
  the book I have read can be rounded down to 0 though.

- The connection between t-SNE and spectral clustering is discussed in detail in [Clustering with
  t-SNE, provably](https://arxiv.org/abs/1706.02582).

- Although I haven't looked very hard, the earliest example of a mention of t-SNE with a spectral
  method I'm aware of is in [The Elastic Embedding Algorithm for Dimensionality Reduction
  (PDF)](http://faculty.ucmerced.edu/mcarreira-perpinan/papers/icml10.pdf), which draws a connection
  between t-SNE and Laplacian Eigenmaps. That paper also mentions that Diffusion Maps use normalized
  affinities (i.e. t-SNE-like normalization to probabilities), but I haven't seen this point made
  elsewhere.

- For more on the properties of the heat kernel, see e.g. [Sun and
  co-workers](https://doi.org/10.1111/j.1467-8659.2009.01515.x) or [Tsitsulin and
  co-workers](https://www.forskningsdatabasen.dk/en/catalog/2491832216).

- Von Luxburg's [A Tutorial on Spectral Clustering](https://arxiv.org/abs/0711.0189) collects a lot
  of material on the definitions of graph Laplacians and the relationship of eigenvalues and
  eigenvectors. Doesn't get into Diffusion Maps, though.

- Radu Horaud's [Graph Laplacian tutorial
  (PDF)](https://csustan.csustan.edu/~tom/Lecture-Notes/Clustering/GraphLaplacian-tutorial.pdf) also
  covers some of the ground of the Von Luxburg tutorial and expresses the relationship between the
  eigenvectors of $L_{sym}$ and $L_{rw}$.

- Another [review on spectral clustering](https://arxiv.org/abs/1901.10204) by Tremblay and Loukas.

- The [Locally Linear Embedding](https://cs.nyu.edu/~roweis/lle/) method turns out to be related to
  Laplacian Eigenmaps. A wonderfully practical tutorial, packed with R code is given in [Cosma
  Shalizi's lecture (PDF)](http://www.stat.cmu.edu/~cshalizi/350/lectures/14/lecture-14.pdf)

- The [Laplacian Eigenmap](https://doi.org/10.1162/089976603321780317) paper is quite readable and
  demonstrates the connection with LLE. It also attempts to justify the now-ubiquitous use of a
  Gaussian kernel for at least the input affinities.

- The [wikipedia page on Diffusion Maps](https://en.wikipedia.org/wiki/Diffusion_map) has one of the
  clearer statements of the algorithm which includes the diffusion parameter, but confusingly
  redefines the matrix it's called $W$ as $L$, even though it was already using $L$ for an entirely
  different purpose.

- Richard Socher's report "Manifold Learning and Dimensionality Reduction with Diffusion Maps" on
  [Diffusion Maps (PDF)](http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.162.3118) is a very
  good place to start on all this, but tragically I think there are some missing symbols in the
  description of the algorithm.

- The [first PNAS paper on Diffusion Maps](http://www.pnas.org/content/102/21/7426.long) is where
  you should go for the definitive statement on the method, but for someone with my level of
  mathematical sophistication (i.e. close to zero) it's very hard going.

- An [easier-going diffusion maps paper](https://doi.org/10.1016/j.acha.2006.04.006) by Coifman and
  Lafon, which explicitly positions diffusion maps as a generalization of ideas expressed in
  Laplacian eigenmaps.

- An [introduction to diffusion maps
  (PDF)](https://inside.mines.edu/~whereman/papers/delaPorte-Herbst-Hereman-vanderWalt-PRASA-2008.pdf)
  which has a fairly clear appendix laying out the definition and properties of $P_{sym}$.

- The [Shi and Malik](https://doi.org/10.1109/34.868688) paper on spectral clustering.

- The [Ng, Jordan and
  Weiss](https://papers.nips.cc/paper/2092-on-spectral-clustering-analysis-and-an-algorithm)
  spectral clustering paper.

- A nice [visual tool](https://dominikschmidt.xyz/spectral-clustering-exp/) for spectral clustering.

- The earliest mention (as far as I know) of the [repeated
  eigendirections](https://dl.acm.org/doi/abs/10.1145/1273496.1273532) problem as the cause for the
  distortion in Laplacian eigenmaps.

- There is a longer version of Inderjit Dhillon's paper on bipartite spectral graph clustering
  published as a technical report (TR-01-05), but it's only available from The University of Texas
  at Austin via FTP: <ftp://ftp.cs.utexas.edu/pub/techreports/tr01-05.pdf>. This is no longer a very
  web-browser friendly protocol, so have fun scrounging around for a more conveniently hosted
  version.

- [NetMF](https://arxiv.org/abs/1710.02971) and variations like
  [NetSMF](https://arxiv.org/abs/1906.11156) and [NetMF+](https://arxiv.org/abs/2110.12782) make use
  of Truncated SVD to get eigenvalues of graph Laplacians in the field of network embedding.

- A popular [review of kernel
  methods](https://projecteuclid.org/journals/annals-of-statistics/volume-36/issue-3/Kernel-methods-in-machine-learning/10.1214/009053607000000677.full).

- A more recent (as of 2021) [review on kernels in machine
  learning](https://arxiv.org/abs/2106.08443) that has references to other papers that connect
  kernel PCA to spectral methods.

- A [review of manifold learning](https://doi.org/10.1146/annurev-statistics-040522-115238) which
  connects spectral methods with dimensionality reduction in general. Note that their definition of
  diffusion maps assumes $\alpha = 1$.

## Some code

For a bit of experimentation on graph Laplacians, there is some R code at:
<https://gist.github.com/jlmelville/772060a26001d7d25d7453b0df4feff9>

Python code is at <https://gist.github.com/jlmelville/8b7655fb4803ce49e4f560d316b04a46>.

Up: [Documentation Home](https://jlmelville.github.io/smallvis/).
