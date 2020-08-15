---
title: "Notes on NCVis"
date: "May 10 2020"
output:
  html_document:
    theme: cosmo
    toc: true
    toc_float:
      collapsed: false
---

[NCVis](https://arxiv.org/abs/2001.11411) uses Noice Contrastive Estimation
(NCE) to create a LargeVis/UMAP-like visualization method. It is interesting
theoretically because it demonstrates the optimization of what looks like an
un-normalized LargeVis-like method is actually a normalized method in disguise,
i.e. you can get to something like t-SNE using stochastic gradient descent.
Compared to LargeVis or UMAP it's more interesting feature as an implementation
is that as well as optimizing the coordinates, it also attempts to "learn" the
normalization factor.

Below I will:

1. Discuss the nuts and bolts of the implementation according to the 
[source code](https://github.com/stat-ml/ncvis).
2. Write out the gradients as they are expressed in the source code, sticking
as close as possuble to how it's done in the code.
3. Express the code gradients in a way that using the quantities and symbols
used in the paper as much as possible and is easier to compare with other
methods.
4. Use the NCE objective function given in the paper to derive the gradient.
5. Compare the gradient I derive to that in the code and show that they are
equivalent, although you have to be a bit careful with how you interpret some of
the variables, whose names and definitions are similar to, but not the same as,
notation in the paper.

## Input

The probability matrix $P$ is constructed from the $N$ high-dimensional input
data vectors $\mathbf{\hat{z}_{i}}$ using the mutual (approximate) nearest
neighbor adjacency matrix and then normalizing so the matrix sums to 1.
Specifically:

1. Create the approximate k-nearest neighbor adjacency matrix, $A$, using
[HNSW](https://github.com/nmslib/hnswlib), i.e. $a_{ij} = 1$ if $j$ is one of
the approximate k-nearest neighbors of $i$. 
2. Create the un-normalized symmetric affinity matrix, $V$ from $A$ by setting 
$v_{ij} = a_{ij} \lor a_{ji}$, i.e. set $v_{ij} = 1$ if $j$ is a neighbor of $i$
or vice versa, and to 0  otherwise. $V$ is therefore the mutual nearest neighbor
adjacency matrix.
3. Create the probability matrix $P$ by matrix normalizing in the usual way of
dividing by the grand sum of $V$: $p_{ij} = v_{ij} / \sum_{i, j} v_{ij}$.

## Initialization of the Output

Uses power iteration for a fixed number of iterations to generate the top $m$
eigenvectors of $P$ where $m$ is the output dimensionality. This is close to a
[graph Laplacian spectral approach](https://jlmelville.github.io/smallvis/spectral.html). The eigenvectors are used to initialize the output coordinates
$\mathbf{z_{i}}$.

## The Output Similarity Kernel

This is the same as UMAP:

$$
\hat{q}_{ij} = \frac{1}{1 + a d_{ij}^{2b}}
$$

where $d_{ij}^2$ is the Euclidean distance between $\mathbf{z_i}$ and
$\mathbf{z_j}$, i.e. the coordinates of $i$ and $j$ in the output space. $a$ and
$b$ are user-defined parameters.

### Normalization

Instead of calculating $Z = \sum_{i,j} \hat{q}_{ij}$, the similarities
$\hat{q}_{ij}$, are normalized to a probability, $q_{ij}$ by:

$$
q_{ij} = e^{-Q} \hat{q}_{ij}
$$

where $Q$ is a parameter which also gets optimized along with the probabilities.

## Sampling

The sampling method follows that of LargeVis and UMAP quite closely.

### Positive Sample

A positive sample is selected by uniform sampling of the positive edges of the
knn graph, which are all the non-zero entries in $P$. In the implementation,
this achieved by iterating over the graph edge list every epoch.

This gives us the pair $\left(i^+, j^+ \right)$

### Negative (Noise) Sample

A negative sample, $j^-$ is found by uniform sampling over all $N - 1$ points
where $j^- \neq i$.

This gives us the pair $\left(i^+, j^- \right)$

## SGD Optimization

This follows the expressions in the
[NCVis source code](https://github.com/stat-ml/ncvis).

Regardless of whether we are using $\left(i^+, j^+ \right)$ or 
$\left(i^+, j^- \right)$, $d_{ij}^2$ and $\hat{q}_{ij}$ are
constructed in the same way as given above. Additionally, the following quantity
is calculated:

$$
w_{ij}^{\prime} 
= 
\frac{e^{-Q^\prime} \hat{q}_{ij}}{\nu}
$$

where $Q^\prime$ is a quantity that is a lot like (but may not be exactly the 
same as) $Q$. $\nu$ is the noise ratio: the number of negative samples we take
for each positive sample. By default this is the same as for UMAP and LargeVis:
five negative samples for each positive sample. One change made in the Python
wrapper (not discussed in the paper) is to keep the total number of negative
samples over the whole optimization constant, but to redistribute $\nu$ so that
it is smaller (but not less than one) at the start of the optimization, and
larger at
the end.

At this point, we calculate the quantity $w_{ij}$. which is dependent on whether
we are optimizing for a positive or negative sample.

### Positive Sample

The positive sample is used to create the following

$$
w_{ij}^+ = \frac{1}{1 + w_{ij}^{\prime}}
$$

### Negative Sample

$$
w_{ij}^- = -\frac{1}{1 + \frac{1}{w_{ij}^{\prime}}} = -\frac{w_{ij}^{\prime}}{1 + w_{ij}^{\prime}}
$$

### Gradient for $Q^\prime$

$$
\nabla J_{Q*} = -w_{ij}^*
$$

where the $*$ subscript just means $+$ or $-$ depending on whether we are
updating using a positive or negative sample.

$Q^\prime$ is updated by gradient ascent using the $Q$-learning rate, 
$\alpha_Q$:

$$
Q^\prime \gets Q^\prime + \alpha_{Q} \nabla J_{Q*}
$$

### Gradient for $\mathbf{z_i}$

$$
k_{ij}^* = 2 w_{ij}^* \hat{q}_{ij} a b d_{ij}^{2\left(b - 1\right)}
$$

Again, the $*$ superscript just means $+$ or $-$ depending on whether we are
updating using a positive or negative sample.

The gradient is then:

$$
\nabla J_{\mathbf{z}*} = 
k_{ij}^* \left( \mathbf{z_j} - \mathbf{z_i} \right)
$$

It's important to note that the displacement is given as 
$\mathbf{z_j} - \mathbf{z_i}$, rather than the usual
$\mathbf{z_i} - \mathbf{z_j}$. This changes the sign of the gradient compared
to similar methods like UMAP and LargeVis, but you also need to bear in mind
that we are doing gradient ascent for NCVis (like LargeVis but unlike UMAP), so
there is another sign change in the update if you compare this expression to
methods that minimize a cost function by gradient descent.

$\mathbf{z_i}$ is updated by gradient ascent using the $\mathbf{z}$ learning
rate $\alpha_{\mathbf{z}}$:

$$
\mathbf{z_i} \gets 
\mathbf{z_i} + \alpha_{\mathbf{z}} 
\nabla J_{\mathbf{z}*} \left( \mathbf{z_i} \right) \\
\mathbf{z_j} \gets
\mathbf{z_j} - \alpha_{\mathbf{z}} 
\nabla J_{\mathbf{z}*} \left( \mathbf{z_i} \right)
$$

Like UMAP and LargeVis there is a clipping step to prevent too large an update,
but it is applied to the update rather than the gradient, i.e. to
$\alpha_{\mathbf{z}} \nabla J_{\mathbf{z}*}$, not just $\nabla J_{\mathbf{z}*}$
as in UMAP and LargeVis. The clipping extrema are $[-4, 4]$ as in UMAP.

## Expanding $w^*$

The expression above have been written to preserve as much as possible how they
are expressed in the implementation of NCVis, using the intermediate variables
$w^\prime$, $w^+$ and $w^-$ (in fact there's just one `w` variable that gets
re-used).

Now let's re-express the gradients, substituting for $w^*$. I will introduce
the $q$-like definition:

$$
q_{ij}^\prime = e^{-Q^\prime} \hat{q}_{ij}
$$

so that:

$$
w_{ij}^\prime = \frac{q_{ij}^\prime}{\nu}
$$

Also, we'll use $\mathbf{z_i} - \mathbf{z_j}$ for the displacement in
the definition of $\nabla J_{\mathbf{z}*}$.

$$
w_{ij}^+ = \frac{1}{1 + w_{ij}^{\prime}}
=
\frac{\nu}{q_{ij}^\prime + \nu}
$$

$$
w_{ij}^- =
-\frac{w_{ij}^{\prime}}{1 + w_{ij}^{\prime}}
=
-\frac{q_{ij}^\prime}{q_{ij}^\prime + \nu}
$$

So the gradients are:

$$
\nabla J_{Q+} = -\frac{\nu}{q_{ij}^\prime + \nu} 
\\
\nabla J_{Q-} = \frac{q_{ij}^\prime}{q_{ij}^\prime + \nu}
\\
\nabla J_{\mathbf{z}+} = 
-2 a b d_{ij}^{2\left(b - 1\right)} 
\hat{q}_{ij}
\frac{\nu}{q_{ij}^\prime + \nu} 
\left( \mathbf{z_i} - \mathbf{z_j} \right)
\\
\nabla J_{\mathbf{z}-} = 
2 a b d_{ij}^{2\left(b - 1\right)} 
\hat{q}_{ij}
\frac{q_{ij}^\prime}{q_{ij}^\prime + \nu} 
\left( \mathbf{z_i} - \mathbf{z_j} \right)
$$

## NCE Objective Function from the Paper

The paper uses $\mathbf{\theta}$ to represent all the parameters we intend to
optimize, which is the $\mathbf{z}$ vector and the scalar $Q$.

The part of the NCE objective applied to positive pairs is:

$$
J^+ = \log \frac{p_m \left( \mathbf{x_i}; \mathbf{\theta} \right)}
{p_m\left( \mathbf{x_i}; \mathbf{\theta} \right) 
+
\nu p_n\left( \mathbf{x_i} \right)}
$$

$\mathbf{x_i}$, is a positive pair $\left(i^+, j^+ \right)$. $p_m$ is the model
probability and $p_n$ is the noise probability. Those will be defined below.

The NCE objective applied to negative pairs is:

$$
J^- =
\log
\frac{
\nu p_n\left( \mathbf{y_i} \right)
}
{
p_m\left( \mathbf{y_i}; \mathbf{\theta} \right)
+
\nu p_n\left( \mathbf{y_i} \right)
}
$$

$\mathbf{y_i}$ is a noise example, i.e. the pair $\left(i^+, j^- \right)$.

The model distribution, $p_m$ is the output probability distribution, i.e.:

$$
p_m \left(i, j \right) = q_{ij}
$$

The noise distribution, $p_n$ is based on the sampling strategy that produces
the negative pairs, $\left(i^+, j^-\right)$.

The positive edges are sampled uniformly, so the probability of sampling $i^+$
is the probability of sampling any edge with $i^+$ as the first node, i.e. the
sum of the probabilities of the $i$th row (or column) of $P$, i.e. 
$\sum_k^N p_{ik}$.

The negative nodes $j^-$ are sampled uniformly with the constraint the 
$j \neq i$, so the probability of seeing any particular $j$ is
$1 / \left(N - 1\right)$.

Therefore the probability of seeing a specific noise pair is the product of
those two probabilities, i.e.

$$
p_n = \frac{1}{N - 1} \sum_{k=1}^N p_{ik}
$$

Instead of having to write $p_n\left(\mathbf{y_i}\right)$, going forward I am
going to use $p_{i}^-$ to mean $p_n$ evaluated for the negative pair
$\left(i^+, j^- \right)$, i.e.:

$$
p_{i}^- \equiv p_n\left(i^+, j^- \right)
$$

There is only the $i$ subscript because its value is entirely independent of $j$.
$p_{i}^-$ should not be confused with the high dimensional input probabilities
$p_{ij}$, although they are used in the calculation of $p_{i}^-$. In the typical
NCE notation, the $p_{ij}$ values are in fact the "data" probabilities, usually
written as $p_d$.

## SGD Form of the Gradient

The gradient of the NCE objective with respect to the output coordinates of $i$
is:

$$
\nabla J_{\mathbf{z*}} \equiv
\frac{\partial J_{ij}^*}{\partial \mathbf{z_i}} =
2 \frac{\partial J_{ij}^*}{\partial d_{ij}^2} 
\left(
\mathbf{z_i} - \mathbf{z_j}
\right)
$$

where as usual the $*$ just means it could be $+$ (for the positive pair) or $-$
(when considering the negative pair). I am going to drop the $*$ from now on
where it doesn't matter whether we are referring to the positive or negative
pairs.

We can expand the gradient expression to:

$$
\frac{\partial J_{ij}}{\partial \mathbf{z_i}} =
2 
\frac{\partial J_{ij}}{\partial q_{ij}} 
\frac{\partial q_{ij}}{\partial \hat{q}_{ij}} 
\frac{\partial \hat{q}_{ij}}{\partial d_{ij}^2} 
\left(
\mathbf{z_i} - \mathbf{z_j}
\right)
$$

We will need the following derivatives:

$$
\frac{\partial \hat{q}_{ij}}{\partial d_{ij}^2} = \frac{-a b d_{ij}^{2\left(b - 1\right)}} 
{\left(a d_{ij}^2 + 1 \right)^2} = -a b d_{ij}^{2\left( b - 1\right)} \hat{q}_{ij}^2 
$$

and:

$$
\frac{\partial q_{ij}}{\partial \hat{q}_{ij}} = 
e^{-Q}
$$

and as a consequence:

$$
\frac{\partial q_{ij}}{\partial d_{ij}^2} = 
-a b d_{ij}^{2\left( b - 1\right)} \hat{q}_{ij}^2 e^{-Q}
$$

or:

$$
\frac{\partial q_{ij}}{\partial d_{ij}^2} = 
-a b d_{ij}^{2\left( b - 1\right)} \hat{q}_{ij} q_{ij}
$$


So the gradient for $\mathbf{z_i}$ is:

$$
\nabla J_{\mathbf{z}} = 
-2 
a b d_{ij}^{2\left( b - 1\right)} 
\hat{q}_{ij} q_{ij}
\frac{\partial J_{ij}}{\partial q_{ij}} 
\left(
\mathbf{z_i} - \mathbf{z_j}
\right)
$$

we therefore need to calculate the derivative of the different attractive and
positive parts of the cost function with respect to $q_{ij}$ and plug that into
the expression above.

For $Q$ we need:

$$
\nabla J_{Q} \equiv
\frac{\partial J_{ij}}{\partial Q} =
\frac{\partial J_{ij}}{\partial q_{ij}}
\frac{\partial q_{ij}}{\partial Q}
$$

and also:

$$
\frac{\partial q_{ij}}{\partial Q} =
-e^{-Q} \hat{q}_{ij} = -q_{ij}
$$

So the "plug-in" gradient for Q is:

$$
\nabla J_{Q} =
-q_{ij} \frac{\partial J_{ij}}{\partial q_{ij}} 
$$

### SGD Update for Positive Sample

The contribution to the cost function for a positive pair sample is:

$$
J_{ij}^+ = \log \frac{q_{ij}}{q_{ij} + \nu p_{i}^-}
$$

The derivative with respect to $q_{ij}$ is:

$$
\frac{\partial J_{ij}^+}{\partial q_{ij}} =
\frac{\nu p_{i}^-}{q_{ij} + \nu p_{i}^-}
\frac{1}{q_{ij}}
$$

So the positive gradients are:

$$
\nabla J_{\mathbf{z+}} =
-2
a b d_{ij}^{2\left( b - 1\right)}
\hat{q}_{ij}
\frac{\nu p_{i}^-}{q_{ij} + \nu p_{i}^-}
\left(
\mathbf{z_i} - \mathbf{z_j}
\right) \\
\nabla J_{Q+} =
-\frac{\nu p_{i}^-}{q_{ij} + \nu p_{i}^-}
$$

### SGD Update for Negative Sample

$$
J_{ij}^- = \log \frac{\nu p_{i}^-}{q_{ij} + \nu p_{i}^-}
$$

The derivative with respect to $q_{ij}$ is:

$$
\frac{\partial J_{ij}^-}{\partial q_{ij}} =
-\frac{1}{q_{ij} + \nu p_{i}^-}
$$

So the negative gradients are:

$$
\nabla J_{\mathbf{z-}} =
2 
a b d_{ij}^{2\left( b - 1\right)}
\hat{q}_{ij}
\frac{q_{ij}}{q_{ij} + \nu p_{i}^-}
\left(
\mathbf{z_i} - \mathbf{z_j}
\right) \\
\nabla J_{Q-} =
\frac{q_{ij}}{q_{ij} + \nu p_{i}^-}
$$

## Comparing the code and paper gradients

Do the code and paper gradients match up? The paper doesn't actually write out
the gradients so we are relying on my ability to do the derivation. So I will
label that as "me" below.

Let's start with the simpler gradient for $Q$ first.

Code:

$$
\nabla J_{Q+} = -\frac{\nu}{q_{ij}^\prime + \nu} 
\\
\nabla J_{Q-} = \frac{q_{ij}^\prime}{q_{ij}^\prime + \nu}
$$

Me:

$$
\nabla J_{Q+} =
-\frac{\nu p_{i}^-}{q_{ij} + \nu p_{i}^-} \\
\nabla J_{Q-} =
\frac{q_{ij}}{q_{ij} + \nu p_{i}^-}
$$

These line up if we define:

$$
q_{ij}^\prime = \frac{q_{ij}}{p_{i}^-}
\implies q_{ij} = q_{ij}^\prime p_{i}^-
$$

and then the $p_{i}^-$ values all cancel out.

Now for the gradient for $\mathbf{z}$.

Code:

$$
\nabla J_{\mathbf{z}+} = 
-2 a b d_{ij}^{2\left(b - 1\right)} 
\hat{q}_{ij}
\frac{\nu}{q_{ij}^\prime + \nu} 
\left( \mathbf{z_i} - \mathbf{z_j} \right)
\\
\nabla J_{\mathbf{z}-} = 
2 a b d_{ij}^{2\left(b - 1\right)} 
\hat{q}_{ij}
\frac{q_{ij}^\prime}{q_{ij}^\prime + \nu} 
\left( \mathbf{z_i} - \mathbf{z_j} \right)
$$

Me:

$$
\nabla J_{\mathbf{z+}} =
-2
a b d_{ij}^{2\left( b - 1\right)}
\hat{q}_{ij}
\frac{\nu p_{i}^-}{q_{ij} + \nu p_{i}^-}
\left(
\mathbf{z_i} - \mathbf{z_j}
\right)
\\
\nabla J_{\mathbf{z-}} =
2 
a b d_{ij}^{2\left( b - 1\right)}
\hat{q}_{ij}
\frac{q_{ij}}{q_{ij} + \nu p_{i}^-}
\left(
\mathbf{z_i} - \mathbf{z_j}
\right) 
$$

Again, this all works out if we use the definition for $q_{ij}^\prime$ above.

What this implies for NCVis implementation is:

$$
q_{ij}^\prime
=
e^{-Q^\prime} \hat{q}_{ij}
=
\frac{1}{p_{i}^-} e^{-Q} \hat{q}_{ij}
$$

This means that the quantity $Q^\prime$ used in the NCVis source code is:

$$
Q^\prime \equiv Q + \log p_{i}^-
$$

This indicates that after sampling each positive edge, we should correct 
$Q^\prime$ by subtracting $\log p_{i}^-$ for the just-sampled edge and then
adding $\log p_{i}^-$ for the new positive edge. 

Without this correction, we assume that $p_{i}^-$ is a constant. Is this a good
assumption? While it is true that the number of non-zero entries in the rows of
the k-nearest neighbor adjacency matrix are constant, after symmetrization the
number of non-zero entries in the rows of the mutual nearest neighbor adjancency
matrix, $V$ (which is normalized to form $P$) is *not* constant.

In a lot of cases, the difference isn't very large, but for datasets where hubs
have formed, those hub nodes will have lots of mutual neighbors, very much 
larger than $k$. Conversely, where there are hub nodes there must be anti-hub
nodes, nodes which aren't in the nearest neighbor list of any other node. Their
mutual neighbor list is the same size as their neighbor list size, $k$. The
difference between $k$ and the hub mutual neighbor list size can be very skewed
and hence the distribution of values of $p_{i}-$ can be very skewed (I have
encountered differences of 1-2 orders of magnitude). The fact that we use the
$\log p_{i}$ may ameliorate this problem, although it could be entirely removed
by storing the $O\left(N\right)$ array of row sums and carrying out the 
correction which only adds one addition and one substraction in the edge 
sampling loop.

## Acknowledgement

Thank you Dmitry Kobak for bringing this paper to my attention and helping
correct the many mistakes I made while trying to write out the derivation.
