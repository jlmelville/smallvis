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
where if $w_{ij}$ is large, that meansobject $i$ and $j$ are considered similar.
$W$ is sometimes also called a kernel matrix and the function that generated the
similarities as a kernel.

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

### Eigenvalues

For Laplacians, the the smallest eigenvalue is 0. For the normalized Laplacians
$L_{sym}$ and $L_{rw}$, the maximum value an eigenvalue can attain in 2.

## Laplacian Eigenmaps

Solve the generalized eigenvalue problem:

$$Lv = \lambda D v$$

The Laplacian Eigenmap uses the smallest eigenvectors. But not the very smallest
eigenvector, $v_1$, which is constant (we can scale it to be a vector of 1s), and 
corresponds to an eigenvalue of zero. So if you want to reduce to two 
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
the same results as the generalized eigenvalue problem with $L$ and $D$.
A non-generalized eigenvalue problem is preferable to the generalized problem, 
at least in R, because generalized problems require installing the CRAN package
`geigen`.

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
of $L$. This is just what Laplacian Eigenamps do, so from the above discussion
we know that it is equivalent to computing the first $k$ eigenvectors of 
$L_{rw}$ (hence justifying the term "normalized").
3. Normalized (Ng, Jordan and Weiss): compute the first $k$ eigenvectors of
$L_{sym}$. This version requires an additional row normalization step of the
output matrix, $Y$, before you can do clustering: normalize each row so the
each row sums to one (unit norm normalization).

After some additional theoretical discussions, von Luxburg concludes that 
clustering on the un-normalized graph Laplacian has some undesirable properties,
so you definitely want to use one of the normalized Laplacians for clustering.
Of the two normalized Laplacians, the Shi-Malik approach (once cast in terms 
of using $L_{rw}$) is the least effort. I would a slight downside to $L_{rw}$
is that unlike $L$ and $L_{sym}$, it is not symmetric, and symmetric matrices
usually have access to slightly more methods (or more efficient methods) for 
solving the eigenproblem.

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

and unlike Laplacian Eigenmaps and spectral clustering, the eigenvalues *are*
used in generating the coordinates. Specifically, they are used to scale the
eigenvectors when forming the $Y$ matrix. For example, here's what one row of
$Y$ would look like in the simplest 2D case:

$$y_i = \left(\mu_{2} v_{i,2}, \mu_{3} v_{i,3} \right)$$

where $v_1$ is the uninteresting eigenvector of all 1s (and $\mu_1 = 1$).

Where does the diffusion come in? $P$ can also be thought of as a transition
matrix of a Markov chain: a large $p_{ij}$ means that $i$ has a high probability
oftransitioning to $j$. And because you can evaluate the probabilities at time
step $t$ by creating the iterated matrix $P^{t}$, you can get a sense of the
geometry of the data at different scales by seeing how the probability changes
over time. And there's not even that much extra work to do: the eigenvectors of 
the iterated matrix are the same as the original matrix $P$, and the eigenvalues 
are given by $\mu^{t}$. 

For a give value of $t$, the 2D diffusion map at time $t$ is therefore:

$$y_i = \left(\mu_{2}^{t} v_{i,2}, \mu_{3}^{t} v_{i,3} \right)$$

or if you want an even more obvious connection to Laplacian Eigenmaps:

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

* Form $W$ as usual.
* Specify $\alpha$, the anisotropic diffusion parameter, a value between 0 and 1.
* Normalize $W$ according to $\alpha$:

$$W^{\left( \alpha \right)} = D^{-1/\alpha} W D^{-1/\alpha}$$

* Form a new diagonal degree matrix, $D^{\left( \alpha \right)}$, based on the new kernel matrix,
$W^{\left( \alpha \right)}$:

$$d^{\left( \alpha \right)}_{ii} = \sum_{j} w^{\left(\alpha\right)}_{ij}$$

* Form a new diffusion operator:

$$P^{\left( \alpha \right)} = D^{\left( \alpha \right)-1} W^{\left( \alpha \right)}$$ 

* Generate the diffusion map using $P^\left( \alpha \right)$, and scaling the 
eigenvectors with the eigenvalues as described above.

It's particularly lamentable that you have to deal with reading both 
$D^{\left( \alpha \right)-1}$ (the inverse of the diagonal matrix 
$D^{\left( \alpha \right)})$ *and* the entirely different $D^{-1/\alpha}$ 
(invert the diagonal matrix $D$, then raise the resulting diagonal values to 
the power of $\alpha$).

When $\alpha = 0$, you get back the diffusion map based on the random walk-style
diffusion operator (and Laplacian Eigenmaps). For $\alpha = 1$, the diffusion
operator approximates the Laplace-Beltrami operator and for $\alpha = 0.5$, you
get Fokker-Planck diffusion.

## Further Reading

I try to link to official DOI URLs and the like where possible, and not post
links to copyright-busting PDFs. The tutorials and reports by Von Luxburg,
Horaud, and Socher are good places to start.

* The connection between t-SNE and spectral clustering is discussed in detail 
[Clustering with t-SNE, provably](https://arxiv.org/abs/1706.02582), but the
earliest example of a mention of t-SNE with a spectral method (Laplacian 
Eigenmaps) I've come across (not that I've looked that hard) is in
[The Elastic Embedding Algorithm for Dimensionality Reduction (PDF)](http://faculty.ucmerced.edu/mcarreira-perpinan/papers/icml10.pdf).
That paper also mentions that Diffusion Maps use normalized affinities (i.e. 
t-SNE like normalization to probabilities), but I haven't seen this point
made elsewhere.

* For more on the properties of the heat kernel, see e.g. 
[Sun and co-workers](https://doi.org/10.1111/j.1467-8659.2009.01515.x) or
[Tsitsulin and co-workers](https://www.forskningsdatabasen.dk/en/catalog/2491832216).

* Von Luxburg's 
[A Tutorial on Spectral Clustering (PDF)](http://www.kyb.mpg.de/fileadmin/user_upload/files/publications/attachments/luxburg06_TR_v2_4139%5b1%5d.pdf)
collects a lot of material on the definitions of graph Laplacians and the
relationship of eigenvalues and eigenvectors. Sadly doesn't get into Diffusion Maps,
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
has one of the clearer statements of the procedure including the diffusion 
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

* The [Ng, Jordan and Weiss (PDF)](https://papers.nips.cc/paper/2092-on-spectral-clustering-analysis-and-an-algorithm.pdf) 
spectral clustering paper.


## Some R code

For a bit of experimentation, here is some R code:

* `randw`: generate a list containing: `W`, a typical random affinity matrix: 
positive semi definite, symmetric, with zeros on the diagonal, and `D` the
degree matrix.

* `lapm`: takes the list of matrices returned by `randw` and returns its own
list containing the various Laplacian matrices discussed: `L`, the un-normalized
Laplacian; `Lsym`, the symmetrized normalized Laplacian; `Lrw`, the random walk
Laplacian; and `P` the random walk transition matrix.

* `eig`: calculates the eigenvectors and eigenvalues of an input matrix `X`, and
returns them in increasing order. If you set `norm = TRUE`, it will return the
eigenvectors after normalizing their length to 1. If you set `norm`, to a
numeric value, it will return the eigenvectors scaled such that their length is
that value. To get the lowest eigenvector to be all 1s, you should set `norm` to
the square root of the number of columns in `X`. As a short-cut you can set
`norm = "n"`. If `val1 = TRUE` then the eigenvalues are subtracted from 1 before
returning. This is useful for the comparison of the eigenvalues of $P$ and
$L_{rw}$, as described in the section on diffusion maps. The return value is the
`vectors` as a matrix with eigenvectors in each column, and `values`, with the
corresponding eigenvalues. Also, `lengths`, which gives the length of each
vector (after any scaling by `norm`).

* `geig`: calculates the generalized eigenvectors and eigenvalues for 
`A v = lambda B v`. Has the same parameters and return value structure as `eig`.
You need to install [geigen](https://cran.r-project.org/package=geigen) for this
to work.

* `reig`: uses the [RSpectra](https://cran.r-project.org/package=RSpectra) 
package to only calculate the first `k` eigenvectors, so you must provide `k`.
If you set `k` to be the full matrix it will just use `eigen` anyway.

```R
randw <- function(n = 3) {
  X <- matrix(rnorm(n * n), nrow = n)
  X <- X * X
  X <- t(X) + X
  diag(X) <- 0
  
  list(W = X, D = diag(colSums(X)))
}

lapm <- function(WD) {
  W <- WD$W
  D <- WD$D

  L <- D - W
  
  # Commented out code more closely follows the mathematical definitions, but
  # pointlessly stores and inverts the full diagonal matrix as well as carrying
  # out matrix multiplications
  # Dinv <- solve(D)
  Dinv <- 1 / diag(D)
  
  # Lsym <- Dinvs %*% L %*% Dinvs
  Dinvs <- sqrt(Dinv)
  Lsym <- Dinvs * sweep(L, 2, Dinvs, '*')

  #P <- Dinv %*% W
  P <- Dinv * WD$W
  
  # I <- diag(nrow = nrow(D), ncol = ncol(D))
  # Lrw <- I - P
  Lrw <- -P
  diag(Lrw) <- 1 - diag(Lrw)
  
  list(L = L, Lsym = Lsym, Lrw = Lrw, P = P)
}

eig <- function(X, norm = FALSE, val1 = FALSE) {
  res <- eigen(X)
  sorteig(res, norm = norm, val1 = val1)
}

geig <- function(A, B, norm = FALSE, val1 = FALSE) {
  res <- geigen::geigen(A, B)
  sorteig(res, norm = norm, val1 = val1)
}

reig <- function(X, k, norm = FALSE, val1 = FALSE) {
  res <- RSpectra::eigs(X, k = k, which = "SM")
  res$vectors <- Re(res$vectors)
  res$values <- Re(res$value)
  sorteig(res, norm = norm, val1 = val1)
}

sorteig <- function(X, norm = FALSE, val1 = FALSE) {
  vectors <- X$vectors
  values <- X$values

  if (val1) {
    values <- 1 - values
  }
  
  if ((is.logical(norm) && norm) || is.numeric(norm) || (is.character(norm) && norm == "n")) {
    if (is.logical(norm)) {
      m <- 1
    }
    else if (is.numeric(norm)) {
      m <- norm
    }
    else {
      m <- sqrt(nrow(vectors)) # make smallest eigenvector all 1s
    }
    sqrtcsums <- sqrt(colSums(vectors * vectors))
    vectors <- m * sweep(vectors, 2, sqrtcsums, "/")
  }

  vectors <- vectors[, order(values)]
  values <- sort(values)

  list(vectors = vectors, values = values, lengths = sqrt(colSums(vectors ^ 2)))
}
```

An example:

```R
WD <- randw(5)
lap <- lapm(WD)

# Laplacian Eigenmaps: these two should give the same results
# use norm = "n", because otherwise the eigenvectors can have different lengths
geig(lap$L, WD$D, norm = "n")
eig(lap$Lrw, norm = "n")
```


Up: [Documentation Home](https://jlmelville.github.io/smallvis/).
