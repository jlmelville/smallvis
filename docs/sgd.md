---
title: "DBD and SGD methods"
date: "December 16, 2017"
output:
  html_document:
    theme: cosmo
    toc: true
    toc_float:
      collapsed: false
---

## Comparison with ADAM

The delta-bar-delta results from the 
[L-BFGS experiments](https://jlmelville.github.io/smallvis/opt.html) are 
included for comparison.

```
tsne_iris_adam <- smallvis(iris, method = "tsne", perplexity = 40, max_iter = 1000, verbose = TRUE, Y_init = "spca", scale = FALSE, ret_extra = c("dx", "dy"), opt = list("adam", eta = 0.25, beta1 = 0.8))
```

### iris

|                             |                           |
:----------------------------:|:--------------------------:
![iris delta-bar-delta](../img/opt/iris_dbd.png)|![iris adam](../img/sgd/iris_adam.png)


### s1k

|                             |                           |
:----------------------------:|:--------------------------:
![s1k delta-bar-delta](../img/opt/s1k_dbd.png)|![s1k adam](../img/sgd/s1k_adam.png)


### Olivetti Faces

|                             |                           |
:----------------------------:|:--------------------------:
![oli delta-bar-delta](../img/opt/oli_dbd.png)|![oli adam](../img/sgd/oli_adam.png)

### Frey Faces

|                             |                           |
:----------------------------:|:--------------------------:
![Frey delta-bar-delta](../img/opt/frey_dbd.png)|![Frey adam](../img/sgd/frey_adam.png)


### COIL-20

|                             |                           |
:----------------------------:|:--------------------------:
![COIL-20 delta-bar-delta](../img/opt/coil20_dbd.png)|![COIL-20 adam](../img/sgd/coil20_adam.png)


### MNIST (6,000)

|                             |                           |
:----------------------------:|:--------------------------:
![MNIST delta-bar-delta](../img/opt/mnist_dbd.png)|![MNIST adam](../img/sgd/mnist_adam.png)

### Fashion (6,000)

|                             |                           |
:----------------------------:|:--------------------------:
![Fashion delta-bar-delta](../img/opt/fashion_dbd.png)|![Fashion adam](../img/sgd/fashion_adam.png)
