---
title: "Spectral Methods"
output:
  html_document:
    theme: cosmo
    toc: true
    toc_float:
      collapsed: false
---

Spectral methods are only tangentially related to `smallvis`, in that a spectral
method is available for initialization (`Y_init = "laplacian"`) and that the 
attractive part of most of the cost functions approximates a spectral method.
Various methods are related to each other, but the nomenclature is not always
consistent between papers, and a lot of the more accessible material seems to
be a bit unclear, so I need something brief to remind myself without having to
work it all out over and over again every time. So it's going here.

This page is just to remind me of the definitions of matrices and the 
algorithmic procedures, not anything to do with theoretical properties. For
that, see the further reading section.

## The Usual Preliminaries

We have a matrix of input data, $X$, with $N$ observations of $K$ features. 

Now let's represent it as a fully-connected undirected graph, where each object
in the dataset is a vertex, and the strength of the connection between each
vertex is given by a weighted edge.

The graph is represented as a matrix, $W$: an $N$ x $N$ matrix of edge weights 
(or "similarities" or "affinities"), where if $w_{ij}$ is large, that means 
object $i$ and $j$ are considered similar. $W$ is sometimes also called a kernel
matrix and the function that generated the similarities as a kernel.

It doesn't particularly matter for this discussion what the kernel function is.
Let's just say it's a Gaussian function of the Euclidean distances (which it
often is):

$$w_{ij} = \exp\left(-r_{ij}^2 / \sigma\right)$$

where $r_{ij}$ is the Euclidean distance between point $i$ and $j$ and $\sigma$
is a bandwidth parameter of some kind (you can set it to one and forget it 
exists if you prefer).

What *does* matter is that $W$ should be:

* Symmetric
* Contain all non-negative values

## The Degree Matrix

The degree matrix is a diagonal matrix where each value in the diagonal is the
sum of each edge associated with a vertex. In other words, sum the rows of $W$
and put that in the diagonal of an $N$ x $N$ matrix. Or sum the columns; $W$ is,
after all, symmetric.

$$d_{ii} = \sum_{j} w_{ij}$$

## Some Graph Laplacians

Now that we have $W$ and $D$, we can create some Laplacians.

### Unnormalized Graph Laplacian

$$L = D - W$$

### Symmetrized Normalized Laplacian

$$L_{sym} = D^{-1/2} L D^{-1/2}$$

### Random Walk Normalized Laplacian

$$L_{rw} = D^{-1} L$$

### Random Walk Transition Matrix

$$P = D^{-1} W$$

Also called the "Diffusion Operator" in Socher's report.

## Eigenvectors

As the "spectral" bit indicates, a lot of eigenvectors are going to be mentioned
below. When I talk about "largest" and "smallest" eigenvectors, I am referring
to the eigenvectors associated with the largest and smallest eigenvalues,
respectively. The eigenvalues of interest are all non-negative and real, so
there shouldn't be any ambiguity about whether I mean a big negative or positive
eigenvalue (none of them are negative except close to zero due to numerical 
issues).

## Laplacian Eigenmaps

Solve the generalized eigenvalue problem:

$$Lv = \lambda_{L} D v$$

The Laplacian Eigenmap uses the smallest eigenvectors. But not the very smallest
eigenvector, which is constant. So if you want to reduce to two dimensions, use
the second-smallest and third-smallest eigenvectors.

It turns out that the following standard eigenvalue problem will produce the 
same eigenvectors:

$$Pv = \lambda_{P} v$$

except with this expression, you want the *largest* eigenvectors (except not
the very largest one), rather than the smallest eigenvectors. The eigenvalues
themselves are not the same, but are related via:

$$\lambda_{L} = 1 - \lambda_{P}$$

Not that you need them for a Laplacian Eigenmap. A non-generalized eigenvalue 
problem is preferable to the generalized problem, at least in R, because 
generalized problems require installing the CRAN package `geigen`.

### The Connection with Locally Linear Embedding

The Laplacian Eigenmap paper demonstrates a connection between LE and LLE, in
that LLE is approximately computing the eigenvectors of $L^2$, which has the
same eigenvectors as $L$.

## Spectral Clustering

Shi and Malik recommend clustering using the generalized eigenvectors of $L$, 
which is the same as clustering on the output of Laplacian Eigenmaps. Ng, Jordan
and Weiss suggest using the normalized bottom eigenvectors of $L_{sym}$, where
the normalization is with respect to ensuring the rows of the output matrix all 
sum to 1.

Von Luxburg notes that clustering on the un-normalized graph Laplacian has some
undesirable properties, so you definitely want to use a normalized Laplacian in
that case. Because the Shi and Malik version of spectral clustering can also be
cast as a standard eigenvalue problem and requires slightly less work than the
Ng, Jordan and Weiss scheme, the Shi-Malik approach is preferable.

## Diffusion Maps

Practically speaking, we can consider diffusion maps as a generalization of
Laplacian Eigenmaps. Our old friend the diffusion operator $P$ can also be
thought of as a transition matrix of a Markov chain: a large $p_{ij}$ means
that $i$ has a high probability of transitioning to $j$. And because you can
evaluate the probabilities at time step $t$ by creating the iterated matrix
$P^{t}$, you can get a sense of the geometry of the data at different 
scales by seeing how the probability changes over time. And because the 
eigenvectors of a matrix raised to a power the same as the original matrix, 
and the eigenvalues are just raised to the power themselves, there's not even 
that much extra work to do.

For the basic diffusion map, the only difference from Laplacian Eigenmaps is 
that the eigenvectors should be scaled by the eigenvalues. For example, the 2D 
diffusion map would be given by the pair of vectors:

$$(\lambda_{L,2}^{t} v_2, \lambda_{L,3}^{t} v_3)$$

where $\lambda_{L,n}$ means the $n$th largest eigenvalue as in LE, raised to the
power of $t$, according to the time step you want. But do note that the 
eigenvalues required are those of the generalized eigenproblem involving $L$,
not the eigenvalues of $P$, so if you are generating the map using $P$, the
correct scaling is:

$$\left[\left(1 - \lambda_{P,2}^{t}\right) v_2, \left(1 - \lambda_{P,3}^{t}\right) v_3 \right]$$

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

* Von Luxburg's 
[A Tutorial on Spectral Clustering (PDF)](http://www.kyb.mpg.de/fileadmin/user_upload/files/publications/attachments/luxburg06_TR_v2_4139%5b1%5d.pdf)
collects a lot of material on the definitions of graph Laplacians and the
relationship of eigenvalues and eigenvectors. Sadly doesn't get into Diffusion Maps,
though.

* Radu Horaud's 
[Graph Laplacian tutorial (PDF)](https://csustan.csustan.edu/~tom/Lecture-Notes/Clustering/GraphLaplacian-tutorial.pdf)
also covers some of the ground of the Von Luxburg tutorial.

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
