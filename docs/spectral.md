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

Spectral methods are only tangentially related to `smallvis`, in that a spectral
method is available for initialization (`Y_init = "laplacian"` or
`Y_init = "normlaplacian"`) and that the attractive part of most of the cost
functions approximates a spectral method. Various methods are related to each
other, but the nomenclature is not always consistent between papers, and a lot
of the more accessible material seems to be a bit unclear, so I need something
brief to remind myself without having to work it all out over and over again
every time. So it's going here.

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

Also, it's common for the diagonal of $W$ to be all zeros.

If your data is naturally a graph, then you skip all of the above, you already
have the data you need to create $W$.

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

### Symmetrized Normalized Laplacian

$$L_{sym} = D^{-1/2} L D^{-1/2} = I - D^{-1/2} W D^{-1/2}$$

also sometimes called just the *normalized* Laplacian and referred to by the 
nomenclature $L_{n}$, but see below for another normalized Laplacian, which
makes this naming ambiguous.

### Random Walk Normalized Laplacian

$$L_{rw} = D^{-1} L = I - D^{-1} W$$

von Luxburg notes that this is *also* sometimes referred to as the normalized
Laplacian, so it's best to use the longer names von Luxburg uses to avoid
confusion.

### Random Walk Transition Matrix

$$P = D^{-1} W$$

Also called the *Diffusion Operator* in Socher's report. This also means you
can write the Random Walk Normalized Laplacian as:

$$L_{rw} = I - P$$

$P$ is row-normalized (i.e. all rows add up to 1).

## Eigenvectors

### Smallest and Largest

As the "spectral" bit indicates, a lot of eigenvectors are going to be mentioned
below. When I talk about "largest" and "smallest" eigenvectors, I am referring
to the eigenvectors associated with the largest and smallest eigenvalues,
respectively. The eigenvalues of interest are all non-negative and real, so
there shouldn't be any ambiguity about whether I mean a big negative or positive
eigenvalue (none of them are negative except close to zero due to numerical 
issues).

### "First" Eigenvectors

I will also borrow the nomenclature of von Luxburg which refers to the "first"
$k$ eigenvectors as being the eigenvectors associated with the $k$ *smallest*
eigenvalues.

### "Top" Eigenvectors

I should mention that you will see discussion of "top" eigenvectors sometimes,
which refers to the eigenvectors associated with the largest eigenvalues in
every place I have seen it. Usually this nomenclature is seen in the context of
SVD, especially when used for PCA, rather than spectral methods. I will *not*
use this anywhere in this document.

### Eigenvalues

For Laplacians, the the smallest eigenvalue is 0. For the normalized Laplacians
$L_{sym}$ and $L_{rw}$, the maximum value an eigenvalue can attain is 2.

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
$N$ x $k$ matrix (let's call it $Y$):


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
2. Normalized (Shi and Malik): compute the first $k$ *generalized* eigenvectors
of $L$. This is just what Laplacian Eigenmaps do, so from the above discussion
we know that it is equivalent to computing the first $k$ eigenvectors of 
$L_{rw}$ (hence justifying the term "normalized").
3. Normalized (Ng, Jordan and Weiss): compute the first $k$ eigenvectors of
$L_{sym}$. This version requires an additional row normalization step of the
output matrix, $Y$, before you can do clustering: normalize each row so the
each row sums to one (unit norm normalization).

After some additional theoretical discussions, von Luxburg concludes that
clustering on the un-normalized graph Laplacian has some undesirable properties,
so you definitely want to use one of the normalized Laplacians for clustering.
Of the two normalized Laplacians, the Shi-Malik approach (once cast in terms of
using $L_{rw}$) is the least effort. A slight downside to $L_{rw}$ is that
unlike $L$ and $L_{sym}$, it is not symmetric, and symmetric matrices usually
have access to slightly more methods (or more efficient methods) for solving the
eigenproblem.

Conclusion: use $L_{rw}$ (Laplacian Eigenmaps).

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
eigenvalue. In diffusion maps, you keep the *top* eigenvectors ignoring the very
top eigenvector. Because of the relationship between $L_{rw}$ and $P$, it turns
out that you still end up using the same eigenvectors for diffusion maps as you
do for Laplacian Eigenmaps, just the ordering is reversed.

Additionally, in a diffusion map the eigenvalues are used to scale the
eigenvectors when forming the $Y$ matrix. For example, here's what one row of
$Y$ would look like in the simplest 2D diffusion map case:

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

The full procedure is something like:

* Form $W$ and $D$ as usual.
* Specify $\alpha$, the anisotropic diffusion parameter, a value between 0 and 1.
* Normalize $W$ according to $\alpha$:

$$W^{\left( \alpha \right)} = D^{-\alpha} W D^{-\alpha}$$

* Form a new diagonal degree matrix, $D^{\left( \alpha \right)}$, based on the new kernel matrix,
$W^{\left( \alpha \right)}$:

$$d^{\left( \alpha \right)}_{ii} = \sum_{j} w^{\left(\alpha\right)}_{ij}$$

* Form a new diffusion operator:

$$P^{\left( \alpha \right)} = D^{\left( \alpha \right)-1} W^{\left( \alpha \right)}$$ 

* Generate the diffusion map using $P^\left( \alpha \right)$, and scale the 
eigenvectors with the eigenvalues as described above.

It's particularly lamentable that you have to deal with reading both 
$D^{\left( \alpha \right)-1}$ (the inverse of the diagonal matrix 
$D^{\left( \alpha \right)})$ *and* the entirely different $D^{-\alpha}$ 
(invert the diagonal matrix $D$, then raise the resulting diagonal values to 
the power of $\alpha$).

When $\alpha = 0$, you get back the diffusion map based on the random walk-style
diffusion operator (and Laplacian Eigenmaps). For $\alpha = 1$, the diffusion
operator approximates the Laplace-Beltrami operator and for $\alpha = 0.5$, you
get Fokker-Planck diffusion.

## Further Reading

I try to link to official DOI URLs and the like where possible, and not post
links to copyright-busting PDFs. The tutorials and reports by von Luxburg,
Horaud, and Socher are good places to start.

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
also covers some of the ground of the Von Luxburg tutorial.

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

* Richard Socher's report on [Diffusion Maps (PDF)](http://www.socher.org/uploads/Main/DiffusionMapsSeminarReport_RichardSocher.pdf)
is a very good place to start on all this, but tragically I think there are
some missingly symbols in the description of the algorithm.

* The [first PNAS paper on Diffusion Maps](http://www.pnas.org/content/102/21/7426.long)
is where you should go for the definitive statement on the method, but for 
someone with my level of mathematical sophistication (i.e. close to zero) it's
very hard going.

* The [Shi and Malik](https://doi.org/10.1109/34.868688) paper on spectral 
clustering.

* The [Ng, Jordan and Weiss](https://papers.nips.cc/paper/2092-on-spectral-clustering-analysis-and-an-algorithm) 
spectral clustering paper.

* A nice [visual tool](https://dominikschmidt.xyz/spectral-clustering-exp/) for
spectral clustering.

## Some code

For a bit of experimentation on graph Laplacians, there is some R code at: 
<https://gist.github.com/jlmelville/772060a26001d7d25d7453b0df4feff9>

Python code is at 
<https://gist.github.com/jlmelville/8b7655fb4803ce49e4f560d316b04a46>.

Up: [Documentation Home](https://jlmelville.github.io/smallvis/).
