---
title: "DBD versus SGD methods"
date: "December 16, 2017"
output:
  html_document:
    theme: cosmo
    toc: true
    toc_float:
      collapsed: false
---

Up: [Documentation Home](https://jlmelville.github.io/smallvis/).

The [delta-bar-delta](https://dx.doi.org/10.1016/0893-6080%2888%2990003-2) method 
used for optimizing t-SNE is an adaptive learning rate technique that
was originally suggested for neural network training. It is mentioned in the same
set of lecture notes that described the 
[RMSprop (PDF)](https://www.cs.toronto.edu/~tijmen/csc321/slides/lecture_slides_lec6.pdf)
stochastic gradient method.

Unlike stochastic gradient methods, DBD is meant for batch learning and only 
uses the sign of the gradient, whereas SGD methods need to take into account
the magnitude of the gradient too. Given what I consider to be the surprisingly
good performance of DBD compared to methods like L-BFGS (see the [optimization]([L-BFGS experiments](https://jlmelville.github.io/smallvis/opt.html) for details), it
would be interesting (to me) to see how SGD methods perform. They certainly
have an advantage of traditional unconstrained optimization methods as being
much easier to implement.

But which optimizer? A good reference is 
[Sebastian Ruder's overview](http://ruder.io/optimizing-gradient-descent/index.html).
They all seem about equally straight-forward to implement, and compared to 
traditional unconstrained optimization methods, they are a breeze.

I have no expectation that any of these should work well, but then again, I defy
anyone to have predicted that delta-bar-delta would work well for t-SNE, so 
let's just try all of them:

* [Adagrad](http://jmlr.org/papers/v12/duchi11a.html)
* [Adadelta](http://arxiv.org/abs/1212.5701) modifies Adagrad to stop the 
learning rate decaying to zero.
* [RMSProp (PDF)](http://www.cs.toronto.edu/~tijmen/csc321/slides/lecture_slides_lec6.pdf)
is similar to Adadelta in that is also modifies Adagrad to prevent zero learning rates.
* [Adam](https://arxiv.org/abs/1412.6980) attempts to estimate mean and variance
of the gradient estimates, incorporating a bias correction.
* [Adamax](https://arxiv.org/abs/1412.6980) a variant on Adam using the infinity
norm.
* [Nadam (PDF)](http://cs229.stanford.edu/proj2015/054_report.pdf), Adam 
modified to use Nesterov momentum.
* [AMSGrad](https://openreview.net/forum?id=ryQu7f-RZ) is a variant that 
corrects a flaw in Adam's convergence under some circumstances.

## Datasets

See the [Datasets](https://jlmelville.github.io/smallvis/datasets.html) page.

## Evaluation

Apart from visualizing the results, the mean neighbor preservation of the
40 closest neighbors is used to provide a rough quantification of the quality
of the result, labelled as `mnp@40` in the plots.

## Settings

Here's an example of generating the t-SNE results with Adam.

The Adam results differ from the recommended defaults given by the authors for
deep learning. I settled on these values by experimenting with the `iris`
dataset results.

```
# Adagrad
iris_adagrad <- smallvis(iris, scale = FALSE, perplexity = 40, Y_init = "spca", opt = list("adagrad", eta = 5))

# Adadelta
iris_adadelta <- smallvis(iris, scale = FALSE, perplexity = 40, Y_init = "spca", opt = list("adadelta"), tol = 0)

# RMSProp
iris_rmsprop <- smallvis(iris, scale = FALSE, perplexity = 40, Y_init = "spca", opt = list("rmsprop", eta = 0.25, rho = 0.8))

# Adam
iris_adam <- smallvis(iris, method = "tsne", perplexity = 40, max_iter = 1000, verbose = TRUE, Y_init = "spca", scale = FALSE, ret_extra = c("dx", "dy"), opt = list("adam", eta = 0.25, beta1 = 0.8))

# Nadam
iris_nadam <- smallvis(iris, scale = FALSE, perplexity = 40, Y_init = "spca", opt = list("nadam", eta = 0.25, beta1 = 0.8))

# Adamax
iris_adamax <- smallvis(iris, scale = FALSE, perplexity = 40, Y_init = "spca", opt = list("adamax", eta = 1, beta1 = 0.9))

# AMSGrad
iris_amsgrad <- smallvis(iris, scale = FALSE, perplexity = 40, Y_init = "spca", opt = list("amsgrad", eta = 0.25, beta1 = 0.8))
```

## Evaluation

For each initialization, the mean neighbor preservation of the
40 nearest neighbors, calculated using the 
[quadra](https://github.com/jlmelville/quadra) package: for each point the 40
nearest neighbors are calculated in the input and output space, and the fraction
of neighbors in common is recorded (0 means no neighbors are in common, 1 means
all the neighbors were preserved). The number reported is the mean average over
all results and is labelled as `mnp@40` in the plots. 40 was chosen for these
results to match the `perplexity`.

## Results

The DBD results are on the left, the Adam results are on the right.

### iris

|                             |                           |
:----------------------------:|:--------------------------:
![iris delta-bar-delta](../img/opt/iris_dbd.png)|![iris adam](../img/sgd/iris_adam.png)
![iris adagrad](../img/sgd/iris_adagrad.png)|![iris adadelta](../img/sgd/iris_adadelta.png)
![iris nadam](../img/sgd/iris_nadam.png)|![iris rmsprop](../img/sgd/iris_rmsprop.png)
![iris adamax](../img/sgd/iris_adamax.png)|![iris amsgrad](../img/sgd/iris_amsgrad.png)


### s1k

|                             |                           |
:----------------------------:|:--------------------------:
![s1k delta-bar-delta](../img/opt/s1k_dbd.png)|![s1k adam](../img/sgd/s1k_adam.png)
![s1k adagrad](../img/sgd/s1k_adagrad.png)|![s1k adadelta](../img/sgd/s1k_adadelta.png)
![s1k nadam](../img/sgd/s1k_nadam.png)|![s1k rmsprop](../img/sgd/s1k_rmsprop.png)
![s1k adamax](../img/sgd/s1k_adamax.png)|![s1k amsgrad](../img/sgd/s1k_amsgrad.png)


### Olivetti Faces

|                             |                           |
:----------------------------:|:--------------------------:
![oli delta-bar-delta](../img/opt/oli_dbd.png)|![oli adam](../img/sgd/oli_adam.png)
![oli adagrad](../img/sgd/oli_adagrad.png)|![oli adadelta](../img/sgd/oli_adadelta.png)
![oli nadam](../img/sgd/oli_nadam.png)|![oli rmsprop](../img/sgd/oli_rmsprop.png)
![oli adamax](../img/sgd/oli_adamax.png)|![oli amsgrad](../img/sgd/oli_amsgrad.png)

### Frey Faces

|                             |                           |
:----------------------------:|:--------------------------:
![frey delta-bar-delta](../img/opt/frey_dbd.png)|![frey adam](../img/sgd/frey_adam.png)
![frey adagrad](../img/sgd/frey_adagrad.png)|![frey adadelta](../img/sgd/frey_adadelta.png)
![frey nadam](../img/sgd/frey_nadam.png)|![frey rmsprop](../img/sgd/frey_rmsprop.png)
![frey adamax](../img/sgd/frey_adamax.png)|![frey amsgrad](../img/sgd/frey_amsgrad.png)

### COIL-20

|                             |                           |
:----------------------------:|:--------------------------:
![coil20 delta-bar-delta](../img/opt/coil20_dbd.png)|![coil20 adam](../img/sgd/coil20_adam.png)
![coil20 adagrad](../img/sgd/coil20_adagrad.png)|![coil20 adadelta](../img/sgd/coil20_adadelta.png)
![coil20 nadam](../img/sgd/coil20_nadam.png)|![coil20 rmsprop](../img/sgd/coil20_rmsprop.png)
![coil20 adamax](../img/sgd/coil20_adamax.png)|![coil20 amsgrad](../img/sgd/coil20_amsgrad.png)


### MNIST (6,000)

|                             |                           |
:----------------------------:|:--------------------------:
![mnist delta-bar-delta](../img/opt/mnist_dbd.png)|![mnist adam](../img/sgd/mnist_adam.png)
![mnist adagrad](../img/sgd/mnist_adagrad.png)|![mnist adadelta](../img/sgd/mnist_adadelta.png)
![mnist nadam](../img/sgd/mnist_nadam.png)|![mnist rmsprop](../img/sgd/mnist_rmsprop.png)
![mnist adamax](../img/sgd/mnist_adamax.png)|![mnist amsgrad](../img/sgd/mnist_amsgrad.png)


### Fashion (6,000)

|                             |                           |
:----------------------------:|:--------------------------:
![fashion delta-bar-delta](../img/opt/fashion_dbd.png)|![fashion adam](../img/sgd/fashion_adam.png)
![fashion adagrad](../img/sgd/fashion_adagrad.png)|![fashion adadelta](../img/sgd/fashion_adadelta.png)
![fashion nadam](../img/sgd/fashion_nadam.png)|![fashion rmsprop](../img/sgd/fashion_rmsprop.png)
![fashion adamax](../img/sgd/fashion_adamax.png)|![fashion amsgrad](../img/sgd/fashion_amsgrad.png)

## Conclusions

Everything except Adadelta works amazingly well. Adadelta, which has no 
controllable learning rate parameter, makes extremely small updates, so it gets
nowhere close to converging, or in fact reducing the initial error by even the
third decimal place for the larger datasets.

Everything else does very well, although not as well as Delta-Bar-Delta. On
the other hand, these methods only require one gradient evaluation per iteration
so are much faster than something like spectral directions. 

For smaller datasets, the usual story is apparent: pretty much any method gets
the job done (apart from Adadelta). For MNIST, there's a bit more difference
apparent, with more of the clusters being split up compared to the DBD results.
Nadam and Adagrad do the best job here. Looking at MNIST and Fashion, RMSProp
also seems to be lagging behind in terms of convergence. In general, the degree
of convergence wasn't as far along after 1000 iterations with any of these
alternative algorithms compared to DBD, probably because the learning rate might
be on the conservative side. On the other hand, having to tune the parameters
for each dataset and algorithm would remove much of its appeal anyway.

Up: [Documentation Home](https://jlmelville.github.io/smallvis/).
