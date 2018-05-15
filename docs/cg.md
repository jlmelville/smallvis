---
title: "Conjugate Gradient"
date: "December 23, 2017"
output:
  html_document:
    theme: cosmo
    toc: true
    toc_float:
      collapsed: false
---

Up: [Documentation Home](https://jlmelville.github.io/smallvis/).

Another in my series of investigating the optimization of t-SNE. See also 
results for [L-BFGS](https://jlmelville.github.io/smallvis/opt.html) and
the approximate Hessian 
[Spectral Direction](https://jlmelville.github.io/smallvis/specd.html). The 
[spectral direction paper](https://arxiv.org/abs/1206.4646) also evaluates
a lot of the same methods.

A conjugate gradient method was used with 
[NeRV](http://www.jmlr.org/papers/v11/venna10a.html), without further details 
being provided. 

I didn't hold out high hopes for this, 
because the general advice in the literature is that they are at best 
competitive with quasi-Newton methods under some circumstances with their main
advantage being their simpler implementation and reduced storage requirements.
The price paid being they require stricter line search. 
The spectral direction paper tested a 
[Matlab implementation](http://learning.eng.cam.ac.uk/carl/code/minimize/),
which uses the Polak-Ribière update (without restarting) and it didn't seem
to do very well.

## Datasets

As usual, see the [Datasets](https://jlmelville.github.io/smallvis/datasets.html) page.

## Settings

In the following experiments, we'll test two types of CG update: the
Polak-Ribière update with restarting in the event of conjugacy being lost (often
called "PRP" or "PR+"), and the update created by Hager and Zhang for the
[CG_DESCENT](http://users.clas.ufl.edu/hager/papers/Software/) package. This
also has a restart rule, so I'll call it "HZ+".

For the line search, I followed the recommendation from Nocedal and Wright
(chapter 3) to use a Strong Wolfe line search with a curvature condition 
constant of `0.1`. However, in  
[the paper describing CG_DESCENT](https://doi.org/10.1145/1132973.1132979), 
Hager and Zhang use a value (in their notation, $\sigma$) of `0.9`, a looser
tolerance that Nocedal and Wright suggest is more appropriate for quasi-Newton
methods. So we'll try both.

```
# CG with a strong Wolfe line search and PR+ update
iris_cgprp <- smallvis(iris, scale = FALSE, perplexity = 40, Y_init = "spca", opt = list("cg", cg_update = "PR+", c2 = 0.1, step_tol = 1e-6, max_gr = 1000), ret_extra = c("DX", "DY"))

# CG with a strong Wolfe line search and HZ+ update
iris_cghzp <- smallvis(iris, scale = FALSE, perplexity = 40, Y_init = "spca", opt = list("cg", cg_update = "HZ+", c2 = 0.1, step_tol = 1e-6, max_gr = 1000), ret_extra = c("DX", "DY"))

# CG with a strong (but loose) Wolfe line search and PR+ update
iris_cgprp <- smallvis(iris, scale = FALSE, perplexity = 40, Y_init = "spca", opt = list("cg", cg_update = "PR+", c2 = 0.9, step_tol = 1e-6, max_gr = 1000), ret_extra = c("DX", "DY"))

# CG with a strong (but loose) Wolfe line search and HZ+ update
iris_cghzp <- smallvis(iris, scale = FALSE, perplexity = 40, Y_init = "spca", opt = list("cg", cg_update = "HZ+", c2 = 0.9, step_tol = 1e-6, max_gr = 1000), ret_extra = c("DX", "DY"))
```

## Evaluation

For each initialization, the mean neighbor preservation of the 40 nearest
neighbors, calculated using the
[quadra](https://github.com/jlmelville/quadra) package. The number reported is 
the mean average over all results and is labelled as `mnp@40` in the plots. 40
was chosen for these results to match the `perplexity`.

## Results

For each dataset, size results are shown below. On the first row, the left image
shows results using the PR+ update, and the right image from the HZ+ update. In
the middle row are the same CG settings except with the line search having a
looser tolerance as used in CG_DESCENT. On the bottom row, the left and right
hand images are those from using default `mize` L-BFGS settings with Wolfe line
search, and the default `smallvis` DBD settings respectively. These used the
same initialization, input scaling and perplexity.

### Iris

| |
:-------|:-------------------:|
|![iris CGPR+](../img/cg/iris_prp.png)|![iris CGHZ+](../img/cg/iris_hzp.png)
|![iris CGPR+ loose](../img/cg/iris_prpl.png)|![iris CGHZ+ loose](../img/cg/iris_hzpl.png)
|![iris l-BFGS](../img/opt/iris_lbfgs-5.png)|![iris DBD](../img/opt/iris_dbd.png)

### s1k

| |
:-------|:-------------------:|
|![s1k CGPR+](../img/cg/s1k_prp.png)|![s1k CGHZ+](../img/cg/s1k_hzp.png)
|![s1k CGPR+ loose](../img/cg/s1k_prpl.png)|![s1k CGHZ+ loose](../img/cg/s1k_hzpl.png)
|![s1k l-BFGS](../img/opt/s1k_lbfgs-5.png)|![s1k DBD](../img/opt/s1k_dbd.png)


### Olivetti Faces

| |
:-------|:-------------------:|
|![oli CGPR+](../img/cg/oli_prp.png)|![oli CGHZ+](../img/cg/oli_hzp.png)
|![oli CGPR+ loose](../img/cg/oli_prpl.png)|![oli CGHZ+ loose](../img/cg/oli_hzpl.png)
|![oli l-BFGS](../img/opt/oli_lbfgs-5.png)|![oli DBD](../img/opt/oli_dbd.png)

### Frey Faces

| |
:-------|:-------------------:|
|![frey CGPR+](../img/cg/frey_prp.png)|![frey CGHZ+](../img/cg/frey_hzp.png)
|![frey CGPR+ loose](../img/cg/frey_prpl.png)|![frey CGHZ+ loose](../img/cg/frey_hzpl.png)
|![frey l-BFGS](../img/opt/frey_lbfgs-5.png)|![frey DBD](../img/opt/frey_dbd.png)

### COIL-20

| |
:-------|:-------------------:|
|![coil20 CGPR+](../img/cg/coil20_prp.png)|![coil20 CGHZ+](../img/cg/coil20_hzp.png)
|![coil20 CGPR+ loose](../img/cg/coil20_prpl.png)|![coil20 CGHZ+ loose](../img/cg/coil20_hzpl.png)
|![coil20 l-BFGS](../img/opt/coil20_lbfgs-5.png)|![coil20 DBD](../img/opt/coil20_dbd.png)

### MNIST (6,000)

| |
:-------|:-------------------:|
|![mnist CGPR+](../img/cg/mnist_prp.png)|![mnist CGHZ+](../img/cg/mnist_hzp.png)
|![mnist CGPR+ loose](../img/cg/mnist_prpl.png)|![mnist CGHZ+ loose](../img/cg/mnist_hzpl.png)
|![mnist l-BFGS](../img/opt/mnist_lbfgs-5.png)|![mnist DBD](../img/opt/mnist_dbd.png)

### Fashion (6,000)

| |
:-------|:-------------------:|
|![fashion CGPR+](../img/cg/fashion_prp.png)|![fashion CGHZ+](../img/cg/fashion_hzp.png)
|![fashion CGPR+](../img/cg/fashion_prpl.png)|![fashion CGHZ+ loose](../img/cg/fashion_hzpl.png)
|![fashion l-BFGS](../img/opt/fashion_lbfgs-5.png)|![fashion DBD](../img/opt/fashion_dbd.png)

This is slightly surprising. On the smaller datasets, as usual, there's little
to choose between any of the methods. For MNIST, the CG results are less ugly
than the L-BFGS result. Not much to choose between the PR+ and HZ+ results,
either. Maybe the HZ+ results have the edge over PR+, but either seems fine. 
Neither is there that much of a difference between the "tight" Wolfe
line search results and the loose results. However, MNIST results do show
more of a cluster break up with the loose results. Also, the fashion MNIST
results do seem to show more outlying points between the (admittedly overlapping)
clusters.

Conjugate gradient seems to be no worse a choice for t-SNE than L-BFGS, although
a tight line search is recommended. But still nothing beats the simplicity of
the DBD defaults for t-SNE.

Up: [Documentation Home](https://jlmelville.github.io/smallvis/).
