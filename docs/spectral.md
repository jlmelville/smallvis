---
title: "Spectral Methods"
output:
  html_document:
    theme: cosmo
    toc: true
    toc_float:
      collapsed: false
---
Up: [Documentation Home](https://jlmelville.github.io/smallvis/).

Originally written in December 2017, with some minor corrections and typos
over the years. Some substantial new additions were added in December 2021
(these are called out in the text).

Spectral methods are only tangentially related to `smallvis`, in that a spectral
method is available for initialization (`Y_init = "laplacian"` or
`Y_init = "normlaplacian"`) and that the attractive part of most of the cost
functions approximates a spectral method. Various methods are related to each
other, but the nomenclature is not always consistent between papers, and a lot
of the more accessible material seems to be a bit unclear, so I need something
brief(ish) to remind myself without having to work it all out over and over
again every time. So it's going here.

This page is just to remind me of the definitions of matrices and the 
algorithmic procedures, not anything to do with theoretical properties. For
that, see the further reading section.

## The Usual Preliminaries

We have a matrix of input data, $X$, with $N$ observations of $K$ features. 

Now let's represent it as a fully-connected undirected graph, where each object
in the dataset is a vertex, and the strength of the connection between each
vertex is given by a weighted edge.

The graph is represented as a matrix, $W$: an $N$ x $N$ matrix of edge weights
(or "similarities" or "affinities": you may also see the matrix written as $A$),
where if $w_{ij}$ is large, that means object $i$ and $j$ are considered
similar. $W$ is sometimes also called a kernel matrix and the function that
generated the similarities as a kernel.

It doesn't particularly matter for this discussion what the kernel function is.
Let's just say it's a Gaussian function of the Euclidean distances (which it
often is):

$$w_{ij} = \exp\left(-r_{ij}^2 / \sigma\right)$$

where $r_{ij}$ is the Euclidean distance between point $i$ and $j$ and $\sigma$
is a bandwidth parameter of some kind (you can set it to one and forget it 
exists if you prefer).

This function is sometimes referred to as a radial basis function (RBF) or the
heat kernel.

What *does* matter is that $W$ should be:

* Symmetric.
* Contain all non-negative values (i.e. positive semi-definite).

If $N$ is large then it is usual to sparsify $W$ by only keeping the 
$k$-largest values in each column. This creates the k-nearest neighbor
similarity graph, which must then be resymmetrized (e.g. by adding the 
self-transpose).

If your data is naturally a graph, then you skip all of the above, you already
have the data you need to create $W$, which is now an adjacency matrix. It's
likely that in that case $W$ is a sparse matrix, but that doesn't matter for
anything that follows.

### Dealing with the Diagonal

If your $W$ matrix is derived from a graph and you don't have loops (i.e. no
edges from a node to itself), then the diagonal of $W$ will be all zeros. Kernel
matrices usually have ones on the diagonal (an obvious consequence of a Gaussian
kernel). What ends up on the diagonal doesn't affect the relationship between 
the different graph Laplacians but can affect the numerical values themselves.

*December 10 2021*: I recently stumbled upon an example of where this difference
lead to some confusion for me. The 
[sklearn's SpectralEmbedding class](https://scikit-learn.org/stable/modules/generated/sklearn.manifold.SpectralEmbedding.html)
generates an affinity matrix (i.e. the diagonal contains all 1s), but to
create the graph Laplacian derived from that affinity matrix, uses
[scipy's csgraph.laplacian function](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csgraph.laplacian.html).
This returns a symmetrized normalized graph Laplacian (see below) calculated
under the assumption that there are all zeros on the diagonal. This doesn't have
an enormous numerical effect on the output, but caused me some bewilderment.

*December 11 2021*: I just noticed that the spectral clustering method of Ng,
Jordan and Weiss explicitly calls for setting the diagonal of the affinity
matrix to zero. The Shi and Malik spectral clustering paper does not mention
this (and the weight heat map figures in their paper suggest the diagonals are
all 1s). I'm pretty sure that the diffusion map literature doesn't mention
doing anything special to the diagonal of the affinity matrix either. So I
can't give you any specific advice here. It's probably not very important.

## The Degree Matrix

The degree matrix is a diagonal matrix where each value in the diagonal is the
sum of each edge associated with a vertex. In other words, sum the rows of $W$
and put that in the diagonal of an $N$ x $N$ matrix. Or sum the columns; $W$ is,
after all, symmetric.

$$d_{ii} = \sum_{j} w_{ij}$$

## Some Graph Laplacians

Now that we have $W$ and $D$, we can create some Laplacians, using the naming
scheme given by von Luxburg:

### Unnormalized Graph Laplacian

$$L = D - W$$

Also referred to as the *combinatorial* graph Laplacian by Tremblay and Loukas.

This is the one graph Laplacian in this document which is unaffected by the
values on the diagonal.

### Symmetrized Normalized Laplacian

$$L_{sym} = D^{-1/2} L D^{-1/2} = I - D^{-1/2} W D^{-1/2}$$

also sometimes called just the *normalized* Laplacian and referred to by the 
nomenclature $L_{n}$, but see below for another normalized Laplacian, which
makes this naming ambiguous.

$D^{-1/2} W D^{-1/2}$ is sometimes called the (symmetrized) normalized adjacency
matrix.

You will often see this normalization written as:

$$L_{sym \left(ij\right)} = \frac{W_{ij}}{\sqrt{D_i D_j}}$$

### Random Walk Normalized Laplacian

$$L_{rw} = D^{-1} L = I - D^{-1} W$$

von Luxburg notes that this is *also* sometimes referred to as the normalized
Laplacian, so it's best to use the longer names von Luxburg uses to avoid
confusion. And yes, $D^{-1} W$ is *also* sometimes referred to as a normalized
adjacency matrix.

### Random Walk Transition Matrix

People really like giving $D^{-1} W$ lots of different names:

$$P = D^{-1} W$$

Also called the *Diffusion Operator* in Socher's report. This also means you
can write the Random Walk Normalized Laplacian as:

$$L_{rw} = I - P$$

$P$ is row-normalized (i.e. all rows add up to 1).

It's also pretty easy to see that:

$$L_{sym} = I - D^{1/2} P D^{-1/2}$$

The matrix $P_{sym} = D^{1/2} P D^{-1/2}$ shows up in the discussion of 
diffusion maps. The doesn't seem to be a fixed name or symbol for it, so I will
go with $P_{sym}$ here in analogy with $L_{sym}$ to indicate it's a symmetrized
version of $P$.

## Eigenvectors

As the "spectral" bit indicates, a lot of eigenvectors are going to be mentioned
below. Usually the eigenvectors are sorted according to the associated
eigenvalues.

### Smallest and Largest

When I talk about "largest" and "smallest" eigenvectors, I am referring
to the eigenvectors associated with the largest and smallest eigenvalues,
respectively. The eigenvalues of interest are all non-negative and real, so
there shouldn't be any ambiguity about whether I mean a big negative or positive
eigenvalue (none of them are negative except close to zero due to numerical 
issues).

### "Top" Eigenvectors

The "top" eigenvectors also means the eigenvectors associated with the largest
eigenvalues. This is commonly found in discussions of singular value
decomposition, especially when used as part of principal component analysis, but
can also find its way into discussions of spectral clustering and graph
laplacians. The "top" eigenvector is sometimes referred to as the "dominant"
eigenvector.

### "First" Eigenvectors

I will also borrow the nomenclature of von Luxburg which refers to the "first"
$k$ eigenvectors as being the eigenvectors associated with the $k$ *smallest*
eigenvalues. There is a bit of a clash of naming here, as you might think that
"top" and "first" eigenvectors refer to the same thing, but they don't. In fact
they have opposite meanings. Great.

### Eigenvalues

Both versions of the normalized adjacency matrix $D^{-1} W$ and
$D^{-1/2} W D^{-1/2}$, have eigenvalues that vary between -1 and 1.

For graph Laplacians, the the smallest eigenvalue is 0. For the normalized
Laplacians $L_{sym}$ and $L_{rw}$, the maximum value an eigenvalue can attain is
2. *December 9 2021*: For proof, see 
[Spectral Graph Theory](http://www.math.ucsd.edu/~fan/research/revised.html), 
Part 5 of Lemma 1.7 in Section 1.3, 'Basic facts about the spectrum of a graph'.
At the time I write this, chapter 1 is freely available at the link above.

## Connections Between Laplacian Eigenvectors

*December 8 2021* While I have made a few typo corrections and clarifications
since originally writing this, I am calling out this brief section because it
adds some relationships that may be of practical interest: it's better to choose
the Laplacian format and spectral decomposition method that meets your
robustness and speed needs, and then convert the result into what you want. You
probably aren't writing the software to do the linear algebra yourself so having
flexibility over the choice of package you use for that can be very helpful.

The eigenvalues of $L_{sym}$ and $L_{rw}$ are the same. So if you just need the
eigenvalues it doesn't matter which of those Laplacians you use. We can
therefore refer to $\lambda_i$ as the ith eigenvalue in most cases of interest
without having to worry about which of the two Laplacians we are referring to.

The eigenvectors are also related. If $v_{rw}$ is an eigenvector of $L_{rw}$
and $v_{sym}$ is an eigenvector of $L_{sym}$, then:

$$v_{rw} = D^{-1/2} v_{sym}$$

As a spoiler for later, the eigenvectors of $P$ are the same as $L_{rw}$ but
after sorting, they are in reverse order, because the eigenvalues of $P$ are
equal to $1 - \lambda$.

Further, $P$ and $P_{sym}$ have the same eigenvalues and the relationship
between eigenvectors is the same as that between $L_{rw}$ and $L_{sym}$:

$$v_{P} = D^{-1/2} v_{Psym}$$

I mainly discuss the eigenvectors of $L_{rw}$ below as these are most closely
related to Laplacian Eigenmaps and Diffusion Maps. Assume if you see reference
to a vector $v$ without any subscript that this refers to the eigenvector of 
$L_{rw}$.

### Scaling the Eigenvectors

Eigenvectors aren't defined to have any particular length. Different software 
will return the eigenvectors with different lengths, so you will need to make
a decision about their length. The Laplacian Eigenmaps literature often refers
to the smallest eigenvectors as a vector of 1s, so you could do that. The actual
length will then depend on the dimensions of your matrices. The other obvious
choice is to scale all the vectors to unit length. Forgetting to scale 
eigenvector output has been a common cause of temporary panic and wild goose 
chasing on my part when writing this document.

## Laplacian Eigenmaps

Solve the generalized eigenvalue problem:

$$Lv = \lambda D v$$

The Laplacian Eigenmap uses the smallest eigenvectors. But not the very smallest
eigenvector, $v_1$, which is constant (we can scale it to be a vector of 1s),
and corresponds to an eigenvalue of zero. So if you want to reduce to two
dimensions, use the second-smallest and third-smallest eigenvectors.

Let's do a brief bit of rearranging:

$$
Lv = \lambda D v \\
D^{-1} L v = \lambda v \\
\left(D^{-1} D - D^{-1} W \right) v = \lambda v \\
\left(I - P \right) v = \lambda v \\
L_{rw} v = \lambda v
$$

So it turns out that the standard eigenvalue problem with $L_{rw}$ will produce
the same results as the generalized eigenvalue problem with $L$ and $D$. A
non-generalized eigenvalue problem is preferable to the generalized problem, at
least in R, because generalized problems require installing the CRAN package
[`geigen`](https://cran.r-project.org/package=geigen). You could even use the
eigenvectors of $P$, although you have to bear in mind that the eigenvalues of
$P$ differ from $L_{rw}$. Compared to $L_{rw}$, when using $P$ the order of the
eigenvectors are reversed, i.e. you want those associated with the *largest*
eigenvalues, ignoring the top eigenvector (which is the constant eigenvector).
$P$ might even be preferable because it's ever-so-slightly less work to
calculate than $L_{rw}$. We'll revisit the relationship between $L_{rw}$ and $P$
when we talk about diffusion maps.

### Output

Now that you have $k$ eigenvectors, stack them columnwise to form an 
$N$ x $k$ matrix. I'll label it $Y$.

$$Y = \left[v_2 | v_3 | \dots v_k \right]$$

where, as noted above, we are not using the uninformative smallest eigenvector,
$v_1$. The rows of that matrix are the coordinates of the graph vertices in the
reduced dimension, i.e. the ith row of the 2D Laplacian Eigenmap representing
vertex i would be:

$$y_i = \left(v_{i,2}, v_{i,3} \right)$$

### The Connection with Locally Linear Embedding

The Laplacian Eigenmap paper demonstrates a connection between LE and LLE, in
that LLE is approximately computing the eigenvectors of $L^2$, which has the
same eigenvectors as $L$ (and the square of the eigenvalues).

## Spectral Clustering and Normalization

von Luxburg describes three different spectral clustering algorithms, which all
involve forming a Laplacian matrix, calculating some eigenvectors, and the
forming the reduced-dimension matrix $Y$ from column-stacking the eigenvectors.

1. Un-normalized: compute the first $k$ eigenvectors of $L$.
2. Normalized ([Shi and Malik](https://ieeexplore.ieee.org/document/868688)): 
compute the first $k$ *generalized* eigenvectors of $L$. This is just what
Laplacian Eigenmaps do, so from the above discussion we know that it is
equivalent to computing the first $k$ eigenvectors of $L_{rw}$ (hence justifying
the term "normalized").
3. Normalized ([Ng, Jordan and Weiss](https://papers.nips.cc/paper/2001/hash/801272ee79cfde7fa5960571fee36b9b-Abstract.html)):
compute the first $k$ eigenvectors of $L_{sym}$. This version requires an
additional row normalization step of the output matrix, $Y$, before you can do
clustering: normalize the rows so each row sums to one (unit norm
normalization).

After some additional theoretical discussions, von Luxburg concludes that
clustering on the un-normalized graph Laplacian has some undesirable properties,
so you definitely want to use one of the normalized Laplacians for clustering.
Of the two normalized Laplacians, the Shi-Malik approach (once cast in terms of
using $L_{rw}$) is the least effort. A slight downside to $L_{rw}$ is that
unlike $L$ and $L_{sym}$, it is not symmetric, and symmetric matrices usually
have access to slightly more methods (or more efficient methods) for solving the
eigenproblem.

Conclusion: use $L_{rw}$ (Laplacian Eigenmaps). *December 8 2021* Probably I
should have said: use the *eigenvectors* of $L_{rw}$, but given the
straightforward conversion between the eigenvectors of $L_{sym}$ and those of
$L_{rw}$, you don't need to form $L_{rw}$ directly for computational purposes.
In fact, see the 'Using Truncated SVD' section below for why you might
want to operate on a matrix related to $L_{sym}$ instead of $L_{rw}$.

## Diffusion Maps

Practically speaking, we can consider diffusion maps as a generalization of
Laplacian Eigenmaps, but solving the eigenproblem for $P$ rather than $L_{rw}$:

$$P v = \mu v$$

The eigenvectors are the same whether you use $L_{rw}$ or $P$, but the 
eigenvalues are different, so I am using $\mu$ instead of $\lambda$ to 
differentiate them from the eigenvalues associated with Laplacian Eigenmaps
and the usual spectral clustering algorithms. As it happens, the eigenvalues
are related by:

$$\mu = 1 - \lambda$$

This does have a book-keeping consequence because we order the eigenvectors by
eigenvalue. In diffusion maps, you keep the *top* eigenvectors, discarding
ignoring the very top eigenvector which is constant (and has an eigenvalue of
1).

Because of the relationship between $L_{rw}$ and $P$, it turns
out that you still end up using the same eigenvectors for diffusion maps as you
do for Laplacian Eigenmaps, just the ordering is reversed.

### Book-keeper Beware!

*December 8 2021* There is a nuance to the book-keeping I just described that
isn't spelled out directly in most discussions of diffusion maps. It might be
obvious to you dear reader, but I admit it took me four years to spot it. There
is a (potential) clash between mathematical convention and the practicalities of
creating a diffusion map via spectral decomposition of $P$.

In mathematical discussions of the properties of eigenvalues, the *absolute*
values of eigenvalues are often what's of interest. This means if you ask for
the largest eigenvalues for a matrix, you can get back the largest negative
as well as positive eigenvalues.

This is of no relevance when ordering the eigenvectors of graph Laplacians:
as mentioned in the "Eigenvalues" section above, the range of $\lambda$ for
$L_{rw}$ is 0 to 2. No negative eigenvalues can be involved (except numerically
close to zero). However, the story is different for the random walk matrices.
The range of the equivalent eigenvalues ($\mu$) is between 1 and -1.

This means that you should be careful when using spectral decomposition
libraries and check the order of returned eigenvalues. For example, both R's
built-in `eigen` function and the `eigs` function in
[RSpectra](https://cran.r-project.org/web/packages/RSpectra/index.html)
(the latter used in the diffusion map package 
[destiny](https://github.com/theislab/destiny)), return the eigenvectors ordered
by absolute value by default. The same is true for the `eigs` function in
the Python package 
[scipy](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.eigs.html).

Therefore, unless you go out of your way to specify the order in which
to return the eigenvectors of $P$, the largest eigenvectors of $P$ are not
necessarily going to be the same as the smallest eigenvectors of $L_{rw}$.

I looked at a couple of implementations of diffusion maps (for example, a [gist
by Rahul
Raj](https://gist.github.com/rahulrajpl/36a5724d0c261b915292182b1d741393), the
[pyDiffMap](https://github.com/DiffusionMapsAcademics/pyDiffMap) package, the
[megaman](https://github.com/mmp2/megaman) package) and the ordering of the
eigenvectors seems to be consistent with using $L_{rw}$, i.e. they take into
account the sign of the eigenvalues, or shift the $P$ matrix to avoid the
eigenvalues changing sign.

The diffusion map literature doesn't really get into any details about ordering
the eigenvectors (or specifically what "largest" means), but diffusion maps are
usually positioned in these papers as being a generalization or extension of
spectral methods on graph Laplacians, so following the same ordering would make
sense there. Muddying the waters a bit for me is this
[early paper on diffusion maps](https://arxiv.org/abs/math/0506090) which says
the eigenvalues of $P$ (called $M$ in the paper) "form a non-increasing sequence
of non-negative numbers", which confuses me: while if that was true it would
remove the ambiguity of ordering the eigenvalues, it doesn't seem to be true.
I must be missing something there.

Anyway, I am fairly sure that the eigenvectors you use in diffusion maps are the
same ones you use with Laplacian Eigenmaps, but let's hold onto the possibility
that I am wrong about that and you should actually order the eigenvectors based
on the largest magnitude of the eigenvalues of $P$, which would be an
interesting departure from the Laplacian Eigenmap approach.

### Scaling the Eigenvectors

However the eigenvectors are ordered, once you have them, the eigenvalues are
used to scale the eigenvectors when forming the $Y$ matrix. For example, here's
what one row of $Y$ would look like in the simplest 2D diffusion map case:

$$y_i = \left(\mu_{N-1} v_{i,N-1}, \mu_{N-2} v_{i,N-2} \right)$$

where $v_N$ is the uninteresting top eigenvector of all 1s (and $\mu_N = 1$).

Where does the diffusion come in? $P$ can also be thought of as a transition
matrix of a Markov chain: a large $p_{ij}$ means that $i$ has a high probability
of transitioning to $j$. And because you can evaluate the probabilities at time
step $t$ by creating the iterated matrix $P^{t}$, you can get a sense of the
geometry of the data at different scales by seeing how the probability changes
over time. And there's not even that much extra work to do: the eigenvectors of
the iterated matrix are the same as the original matrix $P$, and the eigenvalues
are given by $\mu^{t}$.

For a give value of $t$, the 2D diffusion map at time $t$ is therefore:

$$y_i = \left(\mu_{N-1}^{t} v_{i,N-1}, \mu_{N-2}^{t} v_{i,N-2} \right)$$

where $N$ is the total number of eigenvalues. If you want an even more
obvious connection to Laplacian Eigenmaps (and slightly clearer notation):

$$y_i = \left[\left(1 - \lambda_{2}\right)^{t} v_{i,2}, \left(1 - \lambda_{3}\right)^{t} v_{i,3}\right]$$

Because the values of $\lambda$ are in increasing order, this makes it easier
to see that the effect of $t$ is to relatively upweight the contribution of 
later eigenvectors in determining the distance between points on the map.

### Anisotropic Diffusion

So Diffusion Maps are just like Laplacian Eigenmaps but with some slight 
stretching of the axes. Big deal! Actually, there's a bit more to them than
that. You can also define different diffusion operators. At this point, sadly,
I've never actually seen the procedure written down unambiguously. Both the
Socher report and Wikipedia page seem to leave out some steps. Here's my 
attempt. I'll try and stick with the notation I've seen used elsewhere, which 
is truly unfortunate. 

*December 8 2021*: The 
[pyDiffMap](https://github.com/DiffusionMapsAcademics/pyDiffMap) package seems
to carry out a procedure very close to what follows, so this should be reasonably
reliable.

The full procedure is something like:

* Form $W$ and $D$ as usual.
* Specify $\alpha$, the anisotropic diffusion parameter, a value between 0 and 1.
* Normalize $W$ according to $\alpha$:

$$W^{\left( \alpha \right)} = D^{-\alpha} W D^{-\alpha}$$

* Form a new diagonal degree matrix, $D^{\left( \alpha \right)}$, based on the 
new kernel matrix, $W^{\left( \alpha \right)}$:

$$d^{\left( \alpha \right)}_{ii} = \sum_{j} w^{\left(\alpha\right)}_{ij}$$

* Form a new diffusion operator using the new kernel matrix and the inverse of
the new degree matrix (notation is not great here):

$$P^{\left( \alpha \right)} = D^{\left( \alpha \right)-1} W^{\left( \alpha \right)}$$

* Generate the diffusion map using $P^\left( \alpha \right)$, and scale the 
eigenvectors with the eigenvalues as described above. 

* If you want to use a symmetrized version of $P^{\left( \alpha \right)}$, then
form:

$$P^{\left( \alpha \right)}_{sym} = D^{1/2} P^{\left( \alpha \right)} D^{-1/2} $$

and remember to convert the the eigenvectors back via:

$$v_{P\left(\alpha\right)} = D^{-1/2} v_{P\left(\alpha\right)sym}$$

As far as I can tell, you could also define 
$P^{\left( \alpha \right)}_{sym} = D^{\left(\alpha\right)1/2} P^{\left( \alpha \right)} D^{\left(\alpha\right)-1/2}$,
as long as you remember to convert back by
$v_{P\left(\alpha\right)} = D^{\left(\alpha\right)-1/2} v_{P\left(\alpha\right)sym}$,
so you could just choose whichever is most computationally convenient to you.

It's particularly lamentable that you have to deal with reading both 
$D^{\left( \alpha \right)-1}$ (the inverse of the diagonal matrix 
$D^{\left( \alpha \right)})$ *and* the entirely different $D^{-\alpha}$ 
(invert the diagonal matrix $D$, then raise the resulting diagonal values to 
the power of $\alpha$).

When $\alpha = 0$, you get back the diffusion map based on the random walk-style
diffusion operator (and Laplacian Eigenmaps). For $\alpha = 1$, the diffusion
operator approximates the Laplace-Beltrami operator and for $\alpha = 0.5$, you
get Fokker-Planck diffusion.

## Using Truncated SVD

*December 8 2021*. For symmetric graph Laplacians, it should be possible to use
a truncated SVD to get the eigenvectors and eigenvalues. In terms of packages
and routines available, it's always seemed to me that there are a lot more
options for SVD than for eigenvalue problems.

For this approach to be feasible for large matrices, you still need access to
both a truncated SVD routine (i.e. one that doesn't need to calculate all the
singular vectors at once) and one that works with sparse matrices. It's likely
you won't want to code all that up from scratch. And it's probably not going to
be as fast in all cases as a dedicated eigenvalue library.

But if you find yourself in a situation where you do have access to fast SVD
routines, but not to eigenvalue problem solvers, this might be worth
considering.

### Connection Between SVD and Spectral Decomposition

For a symmetric real matrices like $W$ and $L_{sym}$, the complete set of
singular vectors and the complete set of eigenvectors are the same. Also,
the singular values and the *absolute* values of the eigenvalues coincide.

Can we use this to our advantage? To a certain extent, yes. If we are able to do
a full SVD on $L_{sym}$, we can replace the spectral decomposition. However, as
we want the smallest eigenvectors of $L_{sym}$, the equivalent singular vectors
of the SVD are the bottom singular vectors. Therefore we need to do a full SVD
and that gets more and more computationally demanding as the size of $L_{sym}$
grows.

It would be preferable to use the top singular vectors rather than the bottom
singular vectors. Then we can make use of truncated SVD approaches like
[irlba](https://cran.r-project.org/package=irlba) which run in far less time and
memory and hence scale to larger matrices.

Something like the $P$ matrix would be ideal because we want the top
eigenvalues from that matrix. A couple of minor problems to overcome though:

1. $P$ isn't symmetric, an obvious deal-breaker. However, $P_{sym}$ *is*
symmetric, has the same eigenvalues and it's easy to convert the eigenvectors.
1. The eigenvalues of $P$ and $P_{sym}$ range from -1 to 1. This is a slightly
subtler problem. Singular values are always positive, i.e. we lose track of the
sign of the eigenvalues and only get their magnitudes through SVD.
Unfortunately, this ruins the ordering of the singular vectors. For $P$, when we
say we need the largest eigenvectors, we mean the most positive. We can get
around that by taking advantage of the fact that shifting a matrix by $\alpha I$
($\alpha$ being a scalar value) shifts the eigenvalues by $\alpha$.
The eigenvalues of $\alpha I + P$ are therefore $\alpha + \mu$. Picking 
$\alpha=1$ means the eigenvalues of $I + P_{sym}$ are in the range of 0-2, which
is safe for us to run SVD on. *December 26 2021*: the 
[megaman](https://github.com/mmp2/megaman) python package forms a similar
matrix, but adds $\epsilon I$ (with $\epsilon = 2$) to ensure the smallest
eigenvalue is larger than zero and hence the matrix is positive definite rather
than positive semi-definite (see See comments in the
[spectral_embedding](https://github.com/mmp2/megaman/blob/ebf7b6f38c512be22d35ddf69c109506e62cac5b/megaman/embedding/spectral_embedding.py#L103).
Probably the reason is that a positive definite matrix is a requirement for
using the
[LOBPCG](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.lobpcg.html)
eigensolver, so this shouldn't be an issue with the SVD approach I advocate
here. We don't use the singular vectors corresponding to the zero eigenvalue so
the fact that there is no guarantee that the right and left singular vectors are
the same for that singular value in the positive semi-definite case shouldn't be
relevant.
1. Ordering the vectors is now very confusing: for graph Laplacians, we were
talking about the $k$th eigenvector as that associated with the $k$th *smallest*
eigenvalue, but when using $P$, it's actually the $k$th *largest* eigenvalue.
And now we have brought SVD into the mix where the singular vectors are ordered
by decreasing singular value: the $k$th singular vector is that associated with
the $k$th *largest* singular value.

### A Truncated SVD recipe

*December 30 2021*: Thanks to [a reddit comment on spectral
embeddings](https://www.reddit.com/r/MachineLearning/comments/rrrjrz/comment/hqirwix)
I have discovered that several network embedding methods (e.g.
[NetSMF](https://arxiv.org/abs/1906.11156)) use truncated SVD to get
eigenvectors of graph Laplacians, similar to what is described below.

Here's a procedure for calculating eigenvectors I've talked about here via
truncated SVD. As noted above, indexing can get confusing, but it *is*
consistent in the sense that if you stick with the convention of ordering
eigenvectors from smallest to largest, and ordering singular vectors from
largest to smallest, the $k$th eigenvector and $k$th singular vector are the
same (and the eigenvalues and singular values are converted easily).

1. Form $W$ and $D$.
    * If you want a diffusion map, form $W^{\left( \alpha \right)}$ and
    $D^{\left( \alpha \right) - 1}$.
1. Form the symmetric matrix $I + D^{-1/2} W D^{-1/2} = I + P_{sym}$, or
$I + P^{\left(\alpha\right)}_{sym}$ for diffusion maps.
1. Via truncated SVD find the first $k + 1$ singular vectors.
    * The top singular vectors correspond to the smallest eigenvectors of 
    $L_{sym}$.
    * To put it another way: the kth *largest* singular vector is the same as
    the kth *smallest* eigenvector. 
    * The kth *smallest* eigenvalue $\lambda_k$ and the kth *largest* singular
    value, $d_k$ are related by $\lambda_{k} = 2 - d_{k}$.
    * Just as the case with the *smallest* eigenvector, the singular vector
    associated with the *largest* singular value is constant.
1. If you want the eigenvectors of $L_{rw}$, convert in the normal way, i.e. 
$v_{rw} = D^{-1/2} v_{sym}$.
1. The diffusion map eigenvalues are $\mu_{k} = 1 - \lambda_{k} = d_{k} - 1$.

Is this actually worth considering, especially given that none of the diffusion
map packages I looked at in R or Python use SVD directly?

When comparing this approach using [irlba](https://cran.r-project.org/package=irlba)
versus getting the eigenvalues more directly via
[RSpectra](https://cran.r-project.org/package=RSpectra), I didn't
notice any slowdown. However this was in the context of initializing a 
[UMAP](https://github.com/lmcinnes/umap) embedding in the R package 
[uwot](https://cran.r-project.org/package=uwot), and the spectral decomposition
was not a noticeable computational bottleneck in the first place. The main
advantage for me would be that `uwot` already uses `irlba` for PCA in various
places, so I would be able to remove `RSpectra` as a dependency. This is a
decision that has already been reached independently in
[umappp](https://github.com/LTLA/umappp/pull/4) (another UMAP implementation,
this one in C++).

I was able to find one dataset where the truncated SVD approach was slower
than using `RSpectra`: embedding a 1D line from a 3D to 2D: i.e. a dataset with
3 columns: one increasing in value from `1` to `N`, the other two columns being
all zero. `irlba` was much slower in this case. This experiment was inspired by
a [bug report in the UMAP project](https://github.com/lmcinnes/umap/issues/360),
but may not be representative of real-world data.

There could be some more advanced uses of spectral clustering where SVD is the
best choice. For example, in 2001 Inderjit Dhillon published a paper on
[bipartite spectral graph clustering](https://dl.acm.org/doi/10.1145/502512.502550),
where SVD can be applied to a smaller matrix than would be needed if solving
a generalized eigenvalue problem directly (the difference in matrix size doesn't
seem that big, though).

## Repeated Eigendirections

In the previous section I mentioned getting the Laplacian Eigenmap for a 1D line
embedded in 3D. Assuming you get a good converged result, the output using the
first two eigenvectors is a parabola. This is a generic problem when there is a
high "aspect ratio" in a dataset, i.e. the manifold extends much more in one
direction than another. Successive eigenvectors will contain information about
the same coordinate. This problem was described as "repeated eigendirections" by
[Gerber and co-workers](https://dl.acm.org/doi/abs/10.1145/1273496.1273532) and
is also discussed at length by [Goldberg and
co-workers](https://arxiv.org/abs/0806.2646). For more on this, see the
discussion by [Kohli and
co-workers](https://www.jmlr.org/papers/v22/21-0131.html) and especially the
references they point to (under 'Laplacian Eigenmaps' in section 1.3). At this
point I'd love to say "and here's the easy solution that's been discovered", but
that doesn't seem to be the case, so see the above papers and the references
therein for more suggested fixes.

## The Kernel PCA Connection

[Kernel PCA](https://doi.org/10.1162%2F089976698300017467) describes quite a
similar process to everything described above: you create a square affinity
matrix based on the kernel function (usually called the Kernel matrix, 
Gram matrix or Gramian matrix), but instead of forming a graph Laplacian matrix
from it, do SVD on the kernel matrix directly. This relies on the "kernel trick"
(I don't actually know who first coined that phrase): the value of the kernel
matrix element $w_{ij}$ can be seen as the result of transforming $\mathbf{x_i}$
and $\mathbf{x_j}$ into some high-dimensional space (this mapping function is 
usually labelled as $\Phi$), and then taking the dot product. Hence you don't
need to ever actually map the data into the high dimensional space or to even
know what $\Phi$ is. As discussed way at the beginning of this document as long
as $W$ is symmetric and positive semi-definite, then you have a kernel that
works with the kernel trick (such functions are often referred to as a Mercer 
kernel).

A difference between the spectral methods and kernel PCA is the normalization of
$W$. Doing PCA requires mean-centered data and in kernel PCA this means the
transformed data should also be mean-centered. But the whole point of the kernel
trick is to avoid having to actually form the transformed data. Instead, the
kernel matrix is double centered: 
$W_{norm} = \left(I - \frac{1}{N}\mathbf{1} \right) W \left(I - \frac{1}{N}\mathbf{1} \right)$
where $\mathbf{1}$ is an $N$ by $N$ matrix of all 1s. This normalization results
in the row and column means all being zero.

For more on this, [Bengio and
co-workers (PDF)](http://www.iro.umontreal.ca/~lisa/pointeurs/TR1232.pdf) have a
technical report connecting spectral clustering with kernel PCA and [Ham and
co-workers](https://dl.acm.org/doi/10.1145/1015330.1015417) describe how graph
Laplacians can themselves be considered kernels.

Standard linear PCA fits into Kernel PCA by using the "linear kernel", i.e. the 
dot product of the the input vectors: you can get the principal components or
loadings or whatever you are looking for whether you use the scatter/covariance 
matrix ($W'W$) or the Gram matrix ($WW'$), subject to some scaling of eigenvalues
here or a matrix multiplication with $W$ there. If you are doing PCA, you
usually need to center your input data anyway, and the double-centering that
kernel PCA does has no effect on the eigendecomposition.

So if the linear kernel is good enough for PCA, is it good choice for an
affinity matrix for spectral methods? Unfortunately not because there's nothing
to stop a dot product from being negative. This isn't a deal breaker for kernel
PCA because a positive definite kernel doesn't actually have to produce positive
values from its inputs. The eigenvalues of the resulting kernel matrix *do* have
to be all positive, but attempting to construct a graph Laplacian matrix from
that won't give a matrix with the properties we need. There seems to be a small
amount of literature on graph Laplacians with negative edge weights but it
doesn't seem like something I would get very excited about at the moment.

## Further Reading

I try to link to official DOI URLs and the like where possible, and not post
links to copyright-busting PDFs. The tutorials and reports by von Luxburg,
Horaud, and Socher are good places to start.

* The main reference on the properties of graph Laplacians is the monograph 
[Spectral Graph Theory](http://www.math.ucsd.edu/~fan/research/revised.html)
by Fan Chung. The actual amount of the book I have read can be rounded down
to 0 though.

* The connection between t-SNE and spectral clustering is discussed in detail in
[Clustering with t-SNE, provably](https://arxiv.org/abs/1706.02582).

* Although I haven't looked very hard, the earliest example of a mention of
t-SNE with a spectral method I'm aware of is in
[The Elastic Embedding Algorithm for Dimensionality Reduction (PDF)](http://faculty.ucmerced.edu/mcarreira-perpinan/papers/icml10.pdf),
which draws a connection between t-SNE and Laplacian Eigenmaps. That paper also
mentions that Diffusion Maps use normalized affinities (i.e. t-SNE-like
normalization to probabilities), but I haven't seen this point made elsewhere.

* For more on the properties of the heat kernel, see e.g. 
[Sun and co-workers](https://doi.org/10.1111/j.1467-8659.2009.01515.x) or
[Tsitsulin and co-workers](https://www.forskningsdatabasen.dk/en/catalog/2491832216).

* Von Luxburg's 
[A Tutorial on Spectral Clustering](https://arxiv.org/abs/0711.0189)
collects a lot of material on the definitions of graph Laplacians and the
relationship of eigenvalues and eigenvectors. Doesn't get into Diffusion Maps,
though.

* Radu Horaud's 
[Graph Laplacian tutorial (PDF)](https://csustan.csustan.edu/~tom/Lecture-Notes/Clustering/GraphLaplacian-tutorial.pdf)
also covers some of the ground of the Von Luxburg tutorial and expresses the
relationship between the eigenvectors of $L_{sym}$ and $L_{rw}$.

* Another [review on spectral clustering](https://arxiv.org/abs/1901.10204) by
Tremblay and Loukas.

* The [Locally Linear Embedding](https://cs.nyu.edu/~roweis/lle/) method turns 
out to be related to Laplacian Eigenmaps. A wonderfully practical tutorial, packed with R
code is given in [Cosma Shalizi's lecture (PDF)](http://www.stat.cmu.edu/~cshalizi/350/lectures/14/lecture-14.pdf)

* The [Laplacian Eigenmap](https://doi.org/10.1162/089976603321780317) paper is
quite readable and demonstrates the connection with LLE. It also attempts to 
justify the now-ubiquitous use of a Gaussian kernel for at least the input
affinities.

* The [wikipedia page on Diffusion Maps](https://en.wikipedia.org/wiki/Diffusion_map) 
has one of the clearer statements of the algorithm which includes the diffusion 
parameter, but confusingly redefines the matrix it's called $W$ as $L$, even 
though it was already using $L$ for an entirely different purpose.

* Richard Socher's report "Manifold Learning and Dimensionality Reduction with
Diffusion Maps" on [Diffusion Maps (PDF)](http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.162.3118)
is a very good place to start on all this, but tragically I think there are
some missing symbols in the description of the algorithm.

* The [first PNAS paper on Diffusion Maps](http://www.pnas.org/content/102/21/7426.long)
is where you should go for the definitive statement on the method, but for 
someone with my level of mathematical sophistication (i.e. close to zero) it's
very hard going.

* An [easier-going diffusion maps paper](https://doi.org/10.1016/j.acha.2006.04.006)
by Coifman and Lafon, which explicitly positions diffusion maps as a
generalization of ideas expressed in Laplacian eigenmaps.

* An [introduction to diffusion maps (PDF)](https://inside.mines.edu/~whereman/papers/delaPorte-Herbst-Hereman-vanderWalt-PRASA-2008.pdf)
which has a fairly clear appendix laying out the definition and properties of
$P_{sym}$.

* The [Shi and Malik](https://doi.org/10.1109/34.868688) paper on spectral 
clustering.

* The [Ng, Jordan and Weiss](https://papers.nips.cc/paper/2092-on-spectral-clustering-analysis-and-an-algorithm) 
spectral clustering paper.

* A nice [visual tool](https://dominikschmidt.xyz/spectral-clustering-exp/) for
spectral clustering.

* The earliest mention (as far as I know) of the [repeated
eigendirections](https://dl.acm.org/doi/abs/10.1145/1273496.1273532) problem as
the cause for the distortion in Laplacian eigenmaps.

* There is a longer version of Inderjit Dhillon's paper on bipartite spectral 
graph clustering published as a technical report (TR-01-05), but it's only
available from The University of Texas at Austin via FTP:
<ftp://ftp.cs.utexas.edu/pub/techreports/tr01-05.pdf>. This is no longer a
very web-browser friendly protocol, so have fun scrounging around for a more
conveniently hosted version.

* [NetMF](https://arxiv.org/abs/1710.02971) and variations like
[NetSMF](https://arxiv.org/abs/1906.11156) and 
[NetMF+](https://arxiv.org/abs/2110.12782) make use of Truncated SVD to get
eigenvalues of graph Laplacians in the field of network embedding.

* A popular [review of kernel methods](https://projecteuclid.org/journals/annals-of-statistics/volume-36/issue-3/Kernel-methods-in-machine-learning/10.1214/009053607000000677.full).

* A more recent (as of 2021) [review on kernels in machine
learning](https://arxiv.org/abs/2106.08443) that has references to other papers
that connect kernel PCA to spectral methods.

## Some code

For a bit of experimentation on graph Laplacians, there is some R code at: 
<https://gist.github.com/jlmelville/772060a26001d7d25d7453b0df4feff9>

Python code is at 
<https://gist.github.com/jlmelville/8b7655fb4803ce49e4f560d316b04a46>.

Up: [Documentation Home](https://jlmelville.github.io/smallvis/).
