#' Dimensionality Reduction via Neighbor Embedding
#'
#' Carry out dimensionality reduction of a (small) dataset using one of a
#' variety of neighbor embedding methods.
#'
#' Currently supported embedding methods, which can be used as an argument
#' to the \code{method} parameter are:
#' \itemize{
#'   \item \code{"tsne"} t-Distributed Stochastic Neighbor Embedding
#'   (van der Maaten and Hinton, 2008).
#'   \item \code{"bhtsne"} Barnes-Hut t-SNE (van der Maaten, 2014).
#'   \item \code{"largevis"} the cost function of the LargeVis algorithm
#'   (Tang et al, 2016). Input affinities are calculated and symmetrized using
#'   the same perplexity calibration method as t-SNE, but are not normalized.
#'   \item \code{"umap"} the cost function the UMAP method (McInnes, 2017).
#'   Unlike LargeVis and t-SNE, UMAP uses un-normalized input weights, which
#'   are calibrated via calculating smoothed k-nearest-neighbor distances,
#'   rather than perplexity (the procedure is similar, however). The value
#'   of k is controlled by the \code{perplexity} parameter.
#'   \item \code{"tumap"} Like UMAP, except the output weight function is the
#'   t-distribution with one degree of freedom, like t-SNE and LargeVis. This
#'   simplifies the gradient calculation compared to UMAP, so it runs faster.
#'   \item \code{"mmds"} Metric Multidimensional Scaling (Borg and Groenen,
#'   2005), which reduces the "strain", i.e. the square loss between the input
#'   and output distances.
#'   \item \code{"gmmds"} Replaces the input Euclidean distances in
#'   \code{"mmds"} with geodesic distances calculated using k-nearest-neighbors,
#'   in the style of ISOMAP (Tenenbaum and co-workers, 2000). The value of
#'   k is set by the \code{perplexity} parameter.
#'   \item \code{"asne"} The original asymmetric Stochastic Neighbor Embedding
#'   method of Roweis and Hinton (2002).
#'   \item \code{"ssne"} The symmetric SNE method of Cook and co-workers (2007).
#'   \item \code{"hssne"} The heavy-tailed symmetric SNE method of Yang and
#'   co-workers (2009).
#'   \item \code{"wtsne"} The weighted t-SNE Method of Yang and co-workers
#'   (2014). The importance matrix used is the degree centrality (column sums)
#'   of the normalized input weight matrix.
#'   \item \code{"ee"} The Elastic Embedding method of Carreira-Perpinan (2010).
#'   \item \code{"tee"} A t-distributed version of Elastic Embedding (see the
#'   documentation at \url{https://jlmelville.github.io/smallvis/tee.html} for a
#'   description).
#'   \item \code{"nerv"} The Neighbor Retrieval Visualizer method of Venna
#'   and co-workers (2010). This version differs from the presentation in Venna
#'   and co-workers (2010), in that it does not use the precision values
#'   determined from the perplexity calibration in the output kernel functon,
#'   which seems to reflect its presentation in later publications, see e.g.
#'   Yang and co-workers (2014). If you really want the version from Venna and
#'   co-workers, use method \code{"bnerv"}.
#'   \item \code{"jse"} The Jensen-Shannon Embedding method of Lee and
#'   co-workers (2013).
#'   \item \code{"absne"} The alpha-beta SNE method of Narayan and co-workers
#'   (2015).
#'   \item \code{"chsne"} The chi-squared divergence version of t-SNE
#'   (Im and co-workers, 2018).
#'   \item \code{"hlsne"} The Hellinger distance divergence version of t-SNE
#'   (Im and co-workers, 2018).
#'   \item \code{"rklsne"} The reverse Kullback-Leibler divergence version of
#'   t-SNE (Im and co-workers, 2018).
#'   \item \code{"jssne"} The Jensen-Shannon divergence version of t-SNE (Im and
#'   co-workers, 2018).
#'   \item \code{"gsne"}, The global SNE (g-SNE) method of Zhou and Sharpee
#'   (2018).
#' }
#'
#' Note that only the cost function is used from these methods in the context
#' of creating the full distance matrix as part of the gradient calculation.
#' None of the approximations or other speed-ups (e.g. Barnes-Hut or approximate
#' nearest neighbors routines) are used.
#'
#' @section Method-specific options:
#'
#' Some of these methods require the use of extra parameters. To avoid
#' cluttering up the \code{smallvis} function interface, default values are
#' used. To control these, instead of passing a name to the \code{method}
#' parameter, pass a list. The first element is the name of the method you wish
#' to use. Subsequent elements must be named values specifying the parameters.
#'
#' Some parameters are available for all (nor nearly all) methods.
#'
#' \itemize{
#'    \item{\code{inp_kernel}} the input kernel function. Can be one of:
#'       \code{"gauss"} (the default), \code{"exp"} or \code{"knn"}.
#'       \code{"exp"} uses unsquared distances to calculate similarities and is
#'       not the usual literature function, but matches the implementation
#'       in \code{\link[tsne]{tsne}} (and probably doesn't matter very much).
#'       \code{"knn"} uses the symmetric k-nearest neighbors graph with
#'       element (i, j) being set to one if i is one of j's k-nearest
#'       neighbors or vice versa. Other elements are set to zero. Note that no
#'       sparsification is carried out with this kernel, so there are no
#'       memory or performance improvements to be had with this setting.
#'       \code{"skd"} uses the smooth knn distances method as used by UMAP.
#'    \item{\code{symmetrize}} the type of symmetrization, used by symmetric
#'       methods only. Can be one of:
#'       \code{"symmetric"} symmetric nearest neighbor style, by arithmetic
#'       averaging, as in t-SNE.
#'       \code{"fuzzy"} symmetrization by fuzzy set union as used in UMAP.
#'       \code{"mutual"} mutual nearest neighbor style as suggested by Schubert
#'       and Gertz (2017).
#' }
#'
#' \itemize{
#'    \item \code{"LargeVis"}
#'    \itemize{
#'    \item{\code{gamma}} Weighting term for the repulsive versus attractive
#'     forces. Default is \code{1}. The implementation by the creators of
#'     LargeVis uses a default \code{gamma = 7}, but note that this is for
#'     stochastic gradient descent with limited sampling of the repulsive
#'     contributions so it's unlikely to be a good choice with the
#'     implementation used in this package.
#'    \item{\code{gr_eps}} Epsilon used in the gradient to prevent
#'     division by zero. Default is \code{0.1}.
#'    }
#'    \item \code{"UMAP"}
#'    \itemize{
#'    \item{\code{gr_eps}} Epsilon used in the gradient to prevent
#'     division by zero. Default is \code{0.1}.
#'    \item{\code{spread}} Parameter controlling the output kernel function.
#'     Controls the length over which the output kernel decays from 1 to 0.
#'     Default is \code{1}.
#'    \item{\code{min_dist}} Parameter controlling the output kernel function.
#'     According to the UMAP documentation, controls "how tightly the embedding
#'     is allowed compress points together. Larger values ensure embedded
#'     points are more evenly distributed, while smaller values allow the
#'     algorithm to optimise more accurately with regard to local structure.
#'     Sensible values are in the range 0.001 to 0.5". Default is \code{0.001}.
#'    }
#'    \item \code{"tUMAP"}
#'    \itemize{
#'    \item{\code{gr_eps}} Epsilon used in the gradient to prevent
#'     division by zero. Default is \code{0.1}.
#'    }
#'    \item \code{"HSSNE"}
#'    \itemize{
#'    \item{\code{alpha}} Heavy-tailedness parameter, a positive numeric scalar.
#'    Increasing this value increases the amount of stretching that occurs in
#'    the output weighting function. A value of 0 performs like method
#'    \code{"SSNE"}, and a value of 1 performs like method \code{"TSNE"}.
#'    Default is \code{0.5}.
#'    }
#'    \item \code{"EE"}
#'    \itemize{
#'    \item{\code{lambda}} Weighting term for the repulsive versus attractive
#'     forces. Default is \code{100}.
#'    }
#'    \item \code{"tEE"}
#'    \itemize{
#'    \item{\code{lambda}} Weighting term for the repulsive versus attractive
#'     forces. Default is \code{0.01}.
#'    }
#'    \item \code{"NeRV"}
#'    \itemize{
#'    \item{\code{lambda}} Weighting term for relative cost of false positive
#'    versus false negatives (in terms of output distances). Should be a value
#'    between 0 and 1. Setting to 1 performs like method \code{"ASNE"}. Default
#'    is \code{0.9} based on work in Yang and co-workers (2015).
#'    }
#'    \item \code{"JSE"}
#'    \itemize{
#'    \item{\code{kappa}} Weighting term for the combination of KL divergences
#'    used in the JSE cost function. Should be a value between 0 and 1. Setting
#'    to 0 performs like method \code{"ASNE"}. Setting to 1 performs like
#'    method \code{"NeRV"} with its \code{lambda} parameter set to \code{0}.
#'    Yes, that's confusing. The default is \code{0.5}.
#'    }
#'    \item \code{"ABSNE"}
#'    \itemize{
#'    \item{\code{alpha}} Alpha value for the alpha-beta divergence. Set
#'    \code{alpha < 1} to produce more smaller, finer-grained clusters, and
#'    \code{alpha > 1} to produce fewer, larger clusters, with more emphasis
#'    on global structure. Default is \code{1.0}, to give t-SNE-like behavior.
#'    \item{\code{lambda}} Sum of alpha + beta, where beta is the beta value
#'    for the alpha-beta divergence. Set \code{lambda < 1} to increase cluster
#'    separation, and \code{lambda > 1} to decrease cluster separation.
#'    Default is \code{1.0}, to give t-SNE-like behavior.
#'    }
#'    \item \code{"gsne"}
#'    \itemize{
#'    \item{\code{lambda}} Weighting factor to put increasing emphasis on
#'    preserving global similarities. Set to \code{0} to get t-SNE (no
#'    extra emphasis on global structure), and to \code{1.0} to get equal
#'    weighting between the local and global divergences. Default is \code{1.0}.
#'    }
#' }
#'
#' The examples demonstrate how to use \code{method} with these optional
#' parameters.
#'
#' @section Output initialization:
#'
#' For initializing the output coordinates, set the \code{Y_init} parameter
#' to one of the following:
#'
#' \itemize{
#'   \item{A matrix}: which must have dimensions \code{n} by \code{k}, where
#'   \code{n} is the number of rows in \code{X}.
#'   \item{\code{"rand"}}: initialize from a Gaussian distribution with mean 0
#'   and standard deviation 1e-4, the default used by t-SNE. The standard
#'   deviation can be controlled with \code{Y_init_sdev} (see below).
#'   \item{\code{"pca"}}: use the first \code{k} scores of the
#'   PCA: columns are centered, but no scaling beyond that which is applied by
#'   the \code{scale} parameter is carried out.
#'   \item{\code{"spca"}}: uses the PCA scores and then scales each score to a
#'   standard deviation of 1e-4, similar to the method advocated by Kobak and
#'   Berens, 2018. The standard deviation can be controlled with
#'   \code{Y_init_sdev} (see below).
#'   \item{\code{"laplacian"}}: initialize from Laplacian Eigenmap (Belkin and
#'   Niyogi, 2002). The affinity matrix used as input is the result of the
#'   perplexity or smoothed k-nearest neighbor distance calibration. Either
#'   normalized or un-normalized input data will be used, with un-normalized
#'   values used if both are available (although this shouldn't make a
#'   difference). This initialization method cannot be used if the chosen
#'   embedding \code{method} does not calculate a suitable input matrix. If the
#'   input matrix is present but isn't symmetric (e.g. \code{method} choices
#'   such as \code{"jse"}, \code{"nerv"} and \code{"asne"}) then the
#'   initialization will proceed, but the results may not be reliable. If the
#'   \href{https://cran.r-project.org/package=RSpectra}{RSpectra} package is
#'   available, this will be used, otherwise the \code{\link[base]{eigen}}
#'   function is used, which is much slower. Using RSpectra is highly
#'   recommended.
#'   \item{\code{"normlaplacian"}}: initialize from the eigenvectors of the
#'   normalized Laplacian. This is a similar procedure to using Laplacian
#'   Eigenmaps, but with a slightly different Laplacian and is the same
#'   initialization procedure as carried out by UMAP (as of version 0.2.1). The
#'   same limitations apply as for \code{Y_init = "laplacian"}, including the
#'   recommendation to install the RSpectra package.
#' }
#'
#' As spectral methods, using \code{"laplacian"} or \code{"normlaplacian"} is
#' similar to turning off the repulsive interactions in many cost functions
#' related to t-SNE (see Carreira-Perpinan, 2010, and Linderman and
#' Steinerberger, 2017): it may therefore unnecessary to use the
#' \code{exaggeration_factor} setting.
#'
#' The \code{Y_init_sdev} parameter, if provided, will scale the input
#' coordinates such that the standard deviation of each dimension is the
#' provided value. The default is to do no scaling, except for
#' \code{Y_init = "spca"} and \code{Y_init = "rand"} where a scaling to a
#' standard deviation of \code{1e-4} is used, as in t-SNE initialization.
#' \code{Y_init = "spca"} is effectively an alias for \code{Y_init = "pca",
#' Y_init_sdev = 1e-4}.
#'
#' Depending on the embedding method, a particular initialization method may
#' result in initial coordinates with too small or large inter-point distances,
#' which can result in too large or small gradients, respectively. In turn this
#' will lead to difficulties in optimization. Shrinking the initial embedding by
#' rescaling can help under these circumstances. \code{Y_init= "spca"} is a
#' good default place to start.
#'
#' @section Visualization callback:
#'
#' During the optimization, the vizier package
#' (\url{https://www.github.com/jlmelville/vizier}) is used for visualization.
#' To use a custom callback, provide to the \code{epoch_callback} parameter a
#' function with the following signature:
#'
#' \code{function(Y, iter, cost = NULL)}
#'
#' where \code{Y} is the matrix of coordinates, \code{iter} is the current
#' iteration and \code{cost} is the current error value, which will be
#' \code{NULL} the first time this function is called (at iteration 0).
#' The function should have no return value, and presumably will call a plot
#' function. See the "Examples" section for the use of a custom callback.
#' Explicitly set \code{epoch_callback} to \code{NULL} or \code{FALSE} to turn
#' this off.
#'
#' @section Intrinsic dimensionality perplexity:
#'
#' Instead of using a numeric value for the \code{perplexity} argument, you may
#' also set this argument to \code{"idp"} (Intrinsic Dimensionality Perplexity),
#' and the perplexity associated with the estimate of the intrinsic
#' dimensionality of the dataset will be used. This method is only available for
#' methods that use the Gaussian input kernel. To find the IDP, candidate values
#' of perplexities in increasing powers of 2 up to 128 are used by default. To
#' use a custom set of perplexities, provide a list to the \code{perplexity}
#' argument, with the first item being \code{"idp"} and the second item being a
#' vector of perplexities, e.g. \code{perplexity = "idp"} for default behavior
#' or \code{perplexity = list("idp", c(10, 20, 30))} to use the perplexity out
#' of 10, 20, 30 that estimates the intrinsic dimensionality the best.
#' Perplexities are evaluated in the order provided, and the search will
#' terminate early if the estimate of intrinsic dimensionality is judged to be
#' getting worse. It is recommended to provide the candidate perplexities in
#' increasing order. For more details on IDP, see
#' \url{https://jlmelville.github.io/smallvis/idp.html}.
#'
#' @section Multiscale perplexities:
#'
#' Another technique to combine multiple perplexities is to use the multiscale
#' approach given by de Bodt and co-workers (2018). As with IDP, a series of
#' candidate perplexities are used, but all the affinity matrices are used to
#' create an average matrix which is used as the final probability matrix.
#' To use this method, set \code{perplexity = "multiscale"}. Default and custom
#' list of perplexities to use can be provided in the same way as with IDP.
#'
#' Note that previous work by this group described a slightly more complex
#' approach where the number of individual perplexity results are introduced
#' into the average sequentially over the course of the optimization, and
#' the output probabilities are also generated by an averaging. Although also
#' referred to as "multiscale", these variations are not implemented. Also,
#' if \code{ret_extra = TRUE} is used, extra data associated with a specific
#' perplexity (e.g. degree centrality, intrinsic dimensionality) will not
#' be returned.
#'
#' @section Alternative optimizers:
#'
#' Instead of the delta-bar-deta optimizer used in t-SNE, other optimizers can
#' be used, by passing a list to the \code{opt} parameter. The first item
#' is the name of the optimizer. Subsequent items should be named values of
#' the parameters specific to that optimizer. Optimizer and their parameters
#' are as follows:
#'
#' \itemize{
#'   \item{\code{"adagrad"}} Adagrad.
#'   \itemize{
#'     \item{\code{eta}} Initial learning rate.
#'     \item{\code{eps}} Epsilon.
#'   }
#'   \item{\code{"adadelta"}} Adadelta.
#'   \itemize{
#'     \item{\code{rho}} Momentum-like term.
#'     \item{\code{eps}} Epsilon.
#'   }
#'   \item{\code{"rmsprop"}} RMSProp.
#'   \itemize{
#'     \item{\code{eta}} Initial learning rate.
#'     \item{\code{rho}} Momentum-like term.
#'     \item{\code{eps}} Epsilon.
#'   }
#'   \item{\code{"ndbd"}} Delta-bar-delta with Nesterov momentum.
#'   \itemize{
#'     \item{\code{eta}} Initial learning rate.
#'     \item{\code{momentum}} Momentum before iteration given by
#'       \code{mom_switch_iter}.
#'     \item{\code{final_momentum}} Momentum after iteration given by
#'       \code{mom_switch_iter}.
#'     \item{\code{mom_switch_iter}} Switch from \code{momentum} to
#'       \code{final_momentum} at this iteration.
#'     \item{\code{min_gain}} Minimum step size.
#'   }
#'   Arguments to \code{smallvis} that apply to the default delta-bar-delta
#'   optimizer can also be used to set these arguments.
#'   \item{\code{"adam"}} Adam.
#'   \itemize{
#'     \item{\code{eta}} Initial learning rate.
#'     \item{\code{beta1}} Momentum-like term.
#'     \item{\code{beta2}} Momentum-like term for variance decay.
#'     \item{\code{eps}} Epsilon.
#'   }
#'   \item{\code{"adamax"}} Adamax.
#'   \itemize{
#'     \item{\code{eta}} Initial learning rate.
#'     \item{\code{beta1}} Momentum-like term.
#'     \item{\code{beta2}} Momentum-like term for variance decay.
#'     \item{\code{eps}} Epsilon.
#'   }
#'   \item{\code{"nadam"}} Adam with Nesterov momentum.
#'   \itemize{
#'     \item{\code{eta}} Initial learning rate.
#'     \item{\code{beta1}} Momentum-like term.
#'     \item{\code{beta2}} Momentum-like term for variance decay.
#'     \item{\code{eps}} Epsilon.
#'   }
#'   \item{\code{"steepd"}} Steepest descent with fixed step size.
#'   \itemize{
#'     \item{\code{eta}} Step size.
#'   }
#'   \item{\code{"mom"}} Classical momentum with fixed step size.
#'   \itemize{
#'     \item{\code{eta}} Step size.
#'     \item{\code{mu}} Momentum.
#'   }
#' }
#'
#' Note that the majority of these are intended for use in stochastic gradient
#' descent for deep learning. For \code{method = "tsne"}, you probably are better
#' to stick with the default delta-bar-delta method. For other methods, one of
#' these methods may provide greater stability.
#'
#' Alternatively, optimization methods from the mize package,
#' \url{https://cran.r-project.org/package=mize}, can also be used, with the
#' same format: provide the name of the optimizer as the first element of the
#' list, followed by named pairs of arguments. Some possible settings are:
#'
#' \itemize{
#'   \item{\code{"BFGS"}} The Broyden-Fletcher-Goldfarb-Shanno quasi-Newton
#'   method. Best used with small datasets.
#'   \item{\code{"L-BFGS"}} The limited-memory BFGS, a good default choice
#'   if the standard delta-bar-delta optimizer fails to converge.
#'   \item{\code{"SPECD"}} A non-sparse version of the spectral direction
#'   method of Vladymyrov and Carreira-Perpinan (2012). Suitable only for
#'   methods which generate a normalized input affinity matrix (like t-SNE).
#'   Due to its non-sparse forumulation, this is even more time-consuming to use
#'   than BFGS in its current state.
#' }
#'
#' Because these optimizers can carry out multiple gradient and function
#' evaluations per iteration (compared to the default which only calculates
#' one gradient calculation), you may wish to set the \code{"step_tol"}
#' (minimum step size) and/or \code{"max_gr"} (maximum number of gradient
#' evaluation) settings in the list to keep the number of calculations these
#' optimization methods use comparable to the default settings.
#'
#' See the examples for the use of alternative optimizers interface. For best
#' results, you may also consider using the \code{\link{smallvis_perpstep}}
#' function with these optimizers.
#'
#' @param X Input coordinates or distance matrix.
#' @param k Number of output dimensions for the embedding.
#' @param scale If \code{TRUE}, scale each column to zero mean and unit
#'   variance. Alternatively, you may specify one of the following strings:
#'   \code{"range"}, which range scales the matrix elements between 0 and 1;
#'   \code{"absmax"}, here the columns are mean centered and then the elements
#'   divided by absolute maximum value; \code{"scale"} does the same as using
#'   \code{TRUE}. To use the input data as-is, use \code{FALSE}, \code{NULL}
#'   or \code{"none"}.
#' @param Y_init How to initialize the output coordinates. See
#'  the 'Output initialization' section.
#' @param Y_init_sdev If non-\code{NULL}, scales each dimension of the
#'   initialized coordinates (including any user-supplied matrix) to this
#'   standard deviation. See 'Output initialization' section.
#' @param perplexity The target perplexity for parameterizing the input
#'   probabilities. For method \code{"umap"}, controls the neighborhood size
#'   for parameterizing the smoothed k-nearest neighbor distances. See also the
#'   'Intrinsic dimensionality perplexity' and 'Multiscale perplexities'
#'   sections.
#' @param max_iter Maximum number of iterations in the optimization.
#' @param pca If \code{TRUE}, apply PCA to reduce the dimensionality of
#'   \code{X} before any perplexity calibration, but after apply any scaling
#'   and filtering. The number of principal components to keep is specified by
#'   \code{initial_dims}. You may alternatively set this value to
#'   \code{"whiten"}, in which case \code{X} is also whitened, i.e. the
#'   principal components are scaled by the inverse of the square root of the
#'   equivalent eigenvalues, so that the variance of component is 1.
#' @param initial_dims If carrying out PCA or whitening, the number of
#'   principal components to keep. Must be no greater than the rank of the input
#'   or no PCA or whitening will be carried out.
#' @param method A neighbor embedding method. See "Details" and
#'   "Method-specific options" for greater control.
#' @param min_cost If the cost falls below this value, the optimization will
#'   stop early.
#' @param tol If the relative tolerance, averaged over the number of iterations
#'   between successive cost values calculations, falls below this value, the
#'   optimization will stop early. If there is another stage in the optimization
#'   to come then the next stage begins immediately. Otherwise, the optimization
#'   finishes. For \code{tol_wait} iterations after early exaggeration finished,
#'   early stopping will not occur.
#' @param g2tol If the 2-norm (i.e the magnitude) of the gradient vector falls
#'   below this value, stop early. This is off by default and the \code{tol}
#'   parameter is probably more reliable, but for some combinations of
#'   parameters, you may see the cost value plateauing and apparently the
#'   tolerance value falling below \code{tol} but termination not occurring:
#'   this means that the cost function is very slightly increasing. This is
#'   probably a numerical precision issue rather than divergence occuring, so in
#'   this case, use the gradient norm to stop early.
#' @param epoch_callback Function to call after each epoch. See the
#'   "Visualization callback" section. By default the current set of
#'   coordinates will be plotted. Set to\code{FALSE} or \code{NULL} to turn
#'   this off.
#' @param epoch After every \code{epoch} number of steps, calculates and
#'   displays the cost value and calls \code{epoch_callback}, if supplied.
#' @param momentum Initial momentum value.
#' @param final_momentum Final momentum value. If
#'   \code{late_exaggeration_factor > 1}, then during late exaggeration, the
#'   momentum is switched back to \code{momentum} from this value.
#' @param mom_switch_iter Iteration at which the momentum will switch from
#'   \code{momentum} to \code{final_momentum}. If
#'   \code{exaggeration_factor > 1}, then this should occur at some point
#'   after \code{stop_lying_iter} (default is 150 iterations after). If
#'   the early exaggeration phase stops early, this value is treated as being
#'   relative to when early exaggeration stops, to avoid wasting iterations
#'   at the lower momentum value. For example, if \code{stop_lying_iter = 100}
#'   and \code{mom_switch_iter = 250} (the defaults), but early exaggeration
#'   converges at iteration 50, the switch iteration will occur at iteration
#'   150.
#' @param eta Learning rate value, a positive number. Or set to \code{"optsne"},
#'   to use the formula suggested by Belkina and co-workers (2018) in their
#'   opt-SNE package (the size of the dataset divided by the
#'   \code{exaggeration_factor}).
#' @param min_gain Minimum gradient descent step size.
#' @param opt Optional list specifying alternative minimization method. See
#'   "Alternative optimizers" section.
#' @param exaggeration_factor Numerical value to multiply input probabilities
#'   by, during the early exaggeration phase. A value between
#'   4-12 is normal. Usually un-necessary if not using \code{Y_init = "rand"}.
#' @param stop_lying_iter Iteration at which early exaggeration is turned
#'   off.
#' @param late_exaggeration_factor Numerical value to multiply input
#'   probabilities by, during the late exaggeration phase (Linderman and
#'   co-workers, 2017). Smaller values than used with \code{exaggeration_factor}
#'   seem better here (e.g. around 1.5-2).
#' @param start_late_lying_iter Iteration at which to start late exaggeration
#'   (Linderman and co-workers, 2017). Applies \code{late_exaggeration_factor}
#'   after the specified iteration until the end of the optimization. For
#'   \code{max_iter = 1000}, \code{start_late_lying_iter = 900} is suggested.
#'   If the main optimization stops early due to the convergence tolerance being
#'   reached, the late exaggeration stage will start immediately, for the same
#'   relative number of iterations, e.g. if \code{max_iter = 1000} and
#'   \code{start_late_lying_iter = 900}, but the optimization converges at
#'   iteration 800, then late exaggeration will be started between iterations
#'   801 and 900, so that you still only get 100 iterations of late
#'   exaggeration, rather than 200 (which is unlikely to be what you wanted).
#' @param iter0_cost If \code{TRUE}, calculate the cost for the initial
#'   configuration. This while be logged to the console if \code{verbose = TRUE}
#'   or returned in the \code{itercosts} vector if \code{ret_extra = TRUE}.
#' @param ee_mon_epoch If non-\code{NULL}, then this is the number of iterations
#'   between epochs during early exaggeration, and which overrides the value
#'   provided by \code{epoch}. In this mode, the early exaggeration monitoring
#'   of the average relative rate of change in the cost as used in the opt-SNE
#'   method (Belkina and co-workers, 2018) is used and the \code{tol} parameter
#'   is ignored.
#' @param ee_mon_wait If \code{ee_mon_epoch} is non-\code{NULL}, then wait this
#'   number of iterations before monitoring the relative rate of change of the
#'   cost function during early exaggeration.
#' @param ee_mon_buffer If \code{ee_mon_epoch} is non-\code{NULL}, then ignore
#'   this number of occurences of the relative rate of change of the cost
#'   function decreasing, which would otherwise signal termination of
#'   the early exaggeration stage. This is to prevent erroneous termination of
#'   early exaggeration under conditions when the cost can fluctuate noisily.
#' @param tol_wait Wait this number of iterations during standard optimization
#'   (i.e. after early exaggeration, if any), before applying the relative
#'   tolerance early stopping criterion controlled by \code{tol}.
#' @param ret_extra If \code{TRUE}, return value is a list containing additional
#'   values associated with the embedding; otherwise just the output
#'   coordinates. You may also provide a vector of names of potentially large or
#'   expensive-to-calculate values to return, which will be returned in addition
#'   to those value which are returned when this value is \code{TRUE}. See the
#'   \code{Value} section for details.
#' @param n_threads Number of threads to use in multi-threaded code. Default is
#'   0, which means no multi-threading. Mainly affects the calculation of things
#'   like distance matrices and perplexity calibration if you set
#'   \code{use_cpp = TRUE}. Otherwise, only methods that need to calculate
#'   nearest neighbors will be affected.
#' @param use_cpp If \code{TRUE} use multi-threaded C++ code for some
#'   calculations. Default is \code{FALSE}. This won't speed up all steps and
#'   you will want to use this in conjunction with \code{n_threads}.
#' @param eps Set epsilon for avoiding division-by-zero errors. Default is
#'   \code{.Machine$double.eps}, but if you see inconsistent convergence results
#'   with optimizer that should be reducing the cost each iteration, then try
#'   setting this to a larger value, e.g. between \code{1e-3 - 1e-9}.
#' @param theta Barnes-Hut approximation accuracy. Default is \code{1.0} which
#'   more or less corresponds to the default degree of approximation used in
#'   the \code{Rtsne} package. Set to 0.0 for exact t-SNE. Applies only for
#'   \code{method = "bhtsne"}.
#' @param inp_kernel For t-SNE-like methods, determines how the input
#'   similarities are calculated. Possible values are:
#'   \itemize{
#'   \item{\code{"gaussian"}} Gaussian kernel using the target
#'   \code{perplexity}.
#'   \item{\code{"nngaussian"}} Gaussian kernel using the target
#'   \code{perplexity} but similarities beyond the nearest neighbors of each
#'   observation are set to zero. The number of nearest neighbors is three
#'   times the \code{perplexity}. For large datasets, this is faster than 
#'   \code{"gaussian"}, with only a negligible loss of accuracy.
#'   \item \code{"knn"}} k-nearest neighbors kernel. Only the \code{perplexity}
#'   nearest neighbors of each observation are considered and each are given
#'   a similarity of \code{1/perplexity}. This can be a lot faster than using
#'   \code{"nngaussian"} as fewer neighbors need to be found.
#' @param nn Method for finding nearest neighbors. Possible values are:
#'  \itemize{
#'  \item{\code{"exact"}} Exact nearest neighbors.
#'  \item{\code{"approximate"}} Approximate nearest neighbors.
#'  }
#'  Approximate nearest neighbors are faster for large datasets and usually the
#'  small loss in accuracy is not a problem. The default is \code{NULL}, where
#'  each method gets to choose the best method for its needs.
#' @param verbose If \code{TRUE}, log progress messages to the console.
#' @return If \code{ret_extra} is \code{FALSE}, the embedded output coordinates
#'   as a matrix. Otherwise, a list with the following items:
#' \itemize{
#' \item{\code{Y}} Matrix containing the embedded output coordinates.
#' \item{\code{N}} Number of objects.
#' \item{\code{origD}} Dimensionality of the input data.
#' \item{\code{scale}} Scaling applied to input data, as specified by the
#'   \code{scale} parameter.
#' \item{\code{Y_init}} Initialization type of the output coordinates, as
#'   specified by the \code{Y_init} parameter, or if a matrix was used, this will
#'   contain the string \code{"matrix"}.
#' \item{\code{iter}} Number of iterations the optimization carried out.
#' \item{\code{time_secs}} Time taken for the embedding, in seconds.
#' \item{\code{perplexity}} Target perplexity of the input probabilities, as
#'   specified by the \code{perplexity} parameter.
#' \item{\code{costs}} Embedding error associated with each observation. This is
#'   the sum of the absolute value of each component of the KL cost that the
#'   observation is associated with, so don't expect these to sum to the
#'   reported KL cost.
#' \item{\code{itercosts}} KL cost at each epoch.
#' \item{\code{opt}} List containing optimization parameters. For the default
#'   delta-bar-delta method, \code{mom_switch_iter}, \code{momentum},
#'   \code{final_momentum} and \code{eta}. Otherwise, the contents of the
#'   \code{opt} parameter argument.
#' \item{\code{exaggeration_factor}} Multiplier of the input probabilities
#'   during the early exaggeration phase.
#' \item{\code{stop_lying_iter}} Iteration at which early exaggeration is
#'   stopped, as specified by the \code{stop_lying_iter} parameter.
#' \item{\code{late_exaggeration_factor}} Multiplier of the input probabilities
#'   during the late exaggeration phase.
#' \item{\code{stop_late_lying_iter}} Iteration at which late exaggeration is
#'   stopped, as specified by the \code{stop_late_lying_iter} parameter.
#' \item{\code{pca_dims}} If PCA was carried out to reduce the initial
#'   dimensionality of the input, the number of components retained, as
#'   specified by the \code{initial_dims} parameter.
#' \item{\code{whiten_dims}} If PCA whitening was carried out to reduce the
#'   dimensionality of the input, the number of components retained, as
#'   specified by the \code{initial_dims} parameter.
#' \item{\code{G2norm}} The 2-norm of the gradient on the final iteration.
#'   A rough measure of convergence.
#' }
#' Additionally, if you set \code{ret_extra} to a vector of names, these will
#' be returned in addition to the values given above. These values are optional
#' and must be explicitly asked for, because they are either expensive to
#' calculate, take up a lot of memory, or both. The available optional values
#' are:
#' \itemize{
#' \item{\code{X}} The input data, after filtering and scaling.
#' \item{\code{P}} The input probabilities.
#' \item{\code{Q}} The output probabilities.
#' \item{\code{DX}} Input distance matrix. The same as \code{X} when the input
#'   data is already a distance matrix.
#' \item{\code{DY}} Output coordinate distance matrix.
#' }
#'
#' The following additional optional values are available only if the SNE-style
#' perplexity calibration was used as part of the embedding method:
#'
#' \itemize{
#' \item{\code{V}} Input affinities (un-normalized weights used to generate
#' \code{P}).
#' \item{\code{beta}} Precisions (inverse standard deviation) of the input
#' kernel functions. Only available if the gaussian kernel was used.
#' \item{\code{dint}} Intrinsic dimensionalities associated with the input
#' affinities. Only available if the gaussian kernel was used.
#' \item{\code{adegin}} In-degree (column sums) from the un-normalized,
#' un-symmetrized affinity matrix.
#' \item{\code{adegout}} Out-degree (row sums) from the un-normalized,
#' un-symmetrized affinity matrix.
#' \item{\code{adegc}} Degree centrality (average of in-degree and out-degree)
#' from the un-normalized, un-symmetrized affinity matrix.
#' \item{\code{pdeg}} Degree values from the probability matrix. This is always
#' the column sums of the P matrix: For row-normalized P matrices, rows sum to
#' one, so aren't informative and for symmetrized P matrices, the column sum and
#' row sums are identical.
#' }
#'
#' @examples
#' \dontrun{
#'
#' # tsne is the default. verbose = TRUE logs progress to console
#' # Also automatically uses github vizier package for plotting coordinates
#' # during optimization
#' tsne_iris <- smallvis(iris, perplexity = 50, verbose = TRUE)
#'
#' # Can use a custom epoch_callback for visualization
#' colors = grDevices::rainbow(length(unique(iris$Species)))
#' names(colors) = unique(iris$Species)
#' ecb = function(x, y) {
#'   plot(x, t = 'n')
#'   text(x, labels = iris$Species, col = colors[iris$Species])
#' }
#' tsne_iris <- smallvis(iris, epoch_callback = ecb, perplexity = 50, verbose = TRUE)
#'
#' # To turn off visualization entirely:
#' tsne_iris <- smallvis(iris, epoch_callback = FALSE, perplexity = 50, verbose = TRUE)
#'
#' # Try the LargeVis cost function, using its default parameters:
#' largevis_iris <- smallvis(iris, method = "largevis", epoch_callback = ecb,
#'                           perplexity = 50, verbose = TRUE)
#'
#' # Use a list with the method paramter for more control over LargeVis
#' largevis_iris <- smallvis(iris, method = list("largevis", gamma = 1),
#'                           epoch_callback = ecb, perplexity = 50, verbose = TRUE)
#'
#' # Use the UMAP cost function and input weights (perplexity here refers to the
#' # smoothed number of nearest neighbors)
#' umap_iris <- smallvis(iris, method = "umap", eta = 0.1,
#'                       epoch_callback = ecb, perplexity = 50, verbose = TRUE)
#'
#' # Can also control extra UMAP parameters via a list
#' umap_iris <- smallvis(iris, method = list("umap", lv_eps = 0.01), eta = 0.01,
#'                       epoch_callback = ecb, perplexity = 50, verbose = TRUE)
#'
#' # Make embedding deterministic by initializing with scaled PCA scores
#' tsne_iris_spca <- smallvis(iris, epoch_callback = ecb, perplexity = 50, Y_init = "spca",
#'                            exaggeration_factor = 4)
#'
#' # Or use Laplacian Eigenmap for initialization (no exaggeration needed)
#' tsne_iris_lap <- smallvis(iris, epoch_callback = ecb, perplexity = 50,
#'                           Y_init = "laplacian")
#'
#' # Return extra details about the embedding
#' tsne_iris_extra <- smallvis(iris, epoch_callback = ecb, perplexity = 50,
#'                             Y_init = "spca", ret_extra = TRUE)
#'
#' # Return even more details (which can be slow to calculate or take up a lot of memory)
#' tsne_iris_xextra <- smallvis(iris, epoch_callback = ecb, perplexity = 50,
#'                              Y_init = "spca", ret_extra = c("P", "Q", "X", "DX", "DY"))
#'
#' # Reduce initial dimensionality to 3 via PCA
#' # (But you would normally do this with a much larger dataset)
#' tsne_iris_pca <- smallvis(iris, epoch_callback = ecb, perplexity = 50,
#'                           pca = TRUE, initial_dims = 3)
#'
#' # Or use PCA whitening, so all columns of X have variance = 1
#' tsne_iris_whiten <- smallvis(iris, epoch_callback = ecb, perplexity = 50,
#'                              pca = "whiten", initial_dims = 3)
#'
#' # Classical momentum optimization instead of delta-bar-delta
#' umap_iris_mom <- smallvis(iris, scale = FALSE, opt = list("mom", eta = 1e-2, mu = 0.8),
#'                           method = "umap", Y_init = "spca")
#'
#' # L-BFGS optimization via the mize package
#' umap_iris_lbfgs <- smallvis(iris, scale = FALSE, opt = list("l-bfgs", c1 = 1e-4, c2 = 0.9),
#'                             method = "umap", Y_init = "spca", max_iter = 300)
#' # L-BFGS controlling maximum number of gradient evaluations and a further stopping
#' # criterion for when the step size gets too small
#' umap_iris_lbfgs <- smallvis(iris, scale = FALSE, opt = list("l-bfgs", max_gr = 1000,
#'                             step_tol = 1e-6), method = "umap", Y_init = "spca",
#'                             max_iter = 300)
#' ## Example using distance matrix
#' # Use gmmds: metric MDS using geodesic distances estimated from kNN distances
#' # k is set via the perplexity argument. Use an epoch callback to get useful color
#' # in the visualization
#' iris_gmmds <- smallvis(dist(iris[, -5]),
#'                        epoch_callback = function(Y) { plot(Y, col = iris$Species) },
#'                        method = "gmmds", max_iter = 200, perplexity = 8, eta = 1e-4,
#'                        Y_init = "spca")
#' The above should give equivalent results to this non-distance matrix version
#' iris_gmmds <- smallvis(iris, method = "gmmds", max_iter = 400, perplexity = 8, eta = 1e-4,
#'                        Y_init = "spca")
#'
#' # Let smallvis choose an appropriate perplexity for you using the
#' # Intrinsic Dimensionality Perplexity
#' tsne_iris_idp <- smallvis(iris, epoch_callback = ecb, perplexity = "idp", Y_init = "spca",
#'                            exaggeration_factor = 4)
#'
#' # Use your own IDP candidate perplexities by supplying a vector
#' tsne_iris_idp <- smallvis(iris, epoch_callback = ecb, perplexity = list("idp", 3:12),
#'                           Y_init = "spca", exaggeration_factor = 4)
#'
#' # Early Exaggeration
#' tsne_early_ex <- smallvis(iris, exaggeraton_factor = 4, stop_lying_iter = 100)
#'
#' # Also late exaggeration to improve cluster separation
#' tsne_late_ex <- smallvis(iris, exaggeraton_factor = 4, stop_lying_iter = 100,
#'                          late_exaggeration_factor = 1.5, start_late_lying_iter = 900)
#' }
#' @references
#' Belkin, M., & Niyogi, P. (2002).
#' Laplacian eigenmaps and spectral techniques for embedding and clustering.
#' In \emph{Advances in neural information processing systems}
#' (pp. 585-591).
#' \url{http://papers.nips.cc/paper/1961-laplacian-eigenmaps-and-spectral-techniques-for-embedding-and-clustering.pdf}
#'
#' Belkina, A. C., Ciccolella, C. O., Anno, R., Spidlen, J., Halpert, R., & Snyder-Cappione, J. (2018).
#' Automated optimal parameters for T-distributed stochastic neighbor embedding improve visualization and allow analysis of large datasets.
#' \emph{bioRxiv}, 451690.
#' \url{https://www.biorxiv.org/content/10.1101/451690v2.abstract}
#'
#' De Bodt, C., Mulders, D., Verleysen, M., & Lee, J. A. (2018).
#' Perplexity-free t-SNE and twice Student tt-SNE.
#' In \emph{European Symposium on Artificial Neural Networks, Computational Intelligence and Machine Learning (ESANN 2018)} (pp. 123-128).
#' \url{http://hdl.handle.net/2078.1/200844}
#'
#' Borg, I., & Groenen, P. J. (2005).
#' \emph{Modern multidimensional scaling: Theory and applications.}
#' Springer Science & Business Media.
#'
#' Carreira-Perpinan, M. A. (2010, June).
#' The Elastic Embedding Algorithm for Dimensionality Reduction.
#' In \emph{Proceedings of the 27th International Conference on Machine Learning (ICML-10)} (pp. 167-174).
#'
#' Cook, J., Sutskever, I., Mnih, A., & Hinton, G. E. (2007).
#' Visualizing similarity data with a mixture of maps.
#' In \emph{International Conference on Artificial Intelligence and Statistics} (pp. 67-74).
#'
#' Hinton, G. E., & Roweis, S. T. (2002).
#' Stochastic neighbor embedding.
#' In \emph{Advances in neural information processing systems} (pp. 833-840).
#'
#' Im, D. J., Verma, N., & Branson, K. (2018).
#' Stochastic Neighbor Embedding under f-divergences.
#' \emph{arXiv preprint} \emph{arXiv}:1811.01247.
#' \url{https://arxiv.org/abs/1811.01247}
#'
#' Kobak, D., & Berens, P. (2018).
#' The art of using t-SNE for single-cell transcriptomics.
#' \emph{bioRxiv}, 453449.
#' \url{https://doi.org/10.1101/453449}
#'
#' Lee, J. A., Renard, E., Bernard, G., Dupont, P., & Verleysen, M. (2013).
#' Type 1 and 2 mixtures of Kullback-Leibler divergences as cost functions in
#' dimensionality reduction based on similarity preservation.
#' \emph{Neurocomputing}, \emph{112}, 92-108.
#'
#' Linderman, G. C., & Steinerberger, S. (2017).
#' Clustering with t-SNE, provably.
#' \emph{arXiv preprint} \emph{arXiv}:1706.02582.
#' \url{https://arxiv.org/abs/1706.02582}
#'
#' Linderman, G. C., Rachh, M., Hoskins, J. G., Steinerberger, S., & Kluger, Y. (2017).
#' Efficient Algorithms for t-distributed Stochastic Neighborhood Embedding.
#' \emph{arXiv preprint} \emph{arXiv}:1712.09005.
#' \url{https://arxiv.org/abs/1712.09005}
#'
#' McInnes, L., & Healey, J. (2018).
#' UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction
#' \emph{arXiv preprint} \emph{arXiv}:1802.03426.
#' \url{https://arxiv.org/abs/1802.03426}
#'
#' Narayan, K. S., Punjani, A., & Abbeel, P. (2015, June).
#' Alpha-Beta Divergences Discover Micro and Macro Structures in Data.
#' In \emph{Proceedings of the 32nd International Conference on Machine Learning (ICML-14)}
#' (pp 796-804).
#' \url{http://proceedings.mlr.press/v37/narayan15.html}
#'
#' Schubert, E., & Gertz, M. (2017, October).
#' Intrinsic t-stochastic neighbor embedding for visualization and outlier detection.
#' In \emph{International Conference on Similarity Search and Applications}
#' (pp. 188-203). Springer, Cham.
#' \url{https://doi.org/10.1007/978-3-319-68474-1_13}
#'
#' Tang, J., Liu, J., Zhang, M., & Mei, Q. (2016, April).
#' Visualizing large-scale and high-dimensional data.
#' In \emph{Proceedings of the 25th International Conference on World Wide Web}
#' (pp. 287-297).
#' International World Wide Web Conferences Steering Committee.
#' \url{https://arxiv.org/abs/1602.00370}
#'
#' Tenenbaum, J. B., De Silva, V., & Langford, J. C. (2000).
#' A global geometric framework for nonlinear dimensionality reduction.
#' \emph{Science}, \emph{290}(5500), 2319-2323.
#' \url{https://dx.doi.org/10.1126/science.290.5500.2319}
#'
#' Van der Maaten, L., & Hinton, G. (2008).
#' Visualizing data using t-SNE.
#' \emph{Journal of Machine Learning Research}, \emph{9} (2579-2605).
#' \url{http://www.jmlr.org/papers/v9/vandermaaten08a.html}
#'
#' Van der Maaten, L. (2014). 
#' Accelerating t-SNE using tree-based algorithms.
#' \emph{Journal of Machine Learning Research}, \emph{15}(1) (3221-3245).
#' \url{https://jmlr.org/papers/v15/vandermaaten14a.html}
#'
#' Venna, J., Peltonen, J., Nybo, K., Aidos, H., & Kaski, S. (2010).
#' Information retrieval perspective to nonlinear dimensionality reduction for
#' data visualization.
#' \emph{Journal of Machine Learning Research}, \emph{11}, 451-490.
#'
#' Vladymyrov, M., & Carreira-Perpinan, M. A. (2012).
#' Partial-Hessian Strategies for Fast Learning of Nonlinear Embeddings.
#' In \emph{Proceedings of the 29th International Conference on Machine Learning (ICML-12)}
#' (pp. 345-352).
#'
#' Yang, Z., King, I., Xu, Z., & Oja, E. (2009).
#' Heavy-tailed symmetric stochastic neighbor embedding.
#' In \emph{Advances in neural information processing systems}
#' (pp. 2169-2177).
#'
#' Yang, Z., Peltonen, J., & Kaski, S. (2014).
#' Optimization equivalence of divergences improves neighbor embedding.
#' In \emph{Proceedings of the 31st International Conference on Machine Learning (ICML-14)}
#' (pp. 460-468).
#'
#' Yang, Z., Peltonen, J., & Kaski, S. (2015).
#' Majorization-Minimization for Manifold Embedding.
#' In \emph{Proceedings of the 18th International Conference on Artificial Intelligence and Statistics (AISTATS 2015)}
#' (pp. 1088-1097).
#'
#' Zhou, Y., & Sharpee, T. (2018).
#' Using global t-SNE to preserve inter-cluster data structure.
#' \emph{bioRxiv}, 331611.
#' \url{https://doi.org/10.1101/331611}
#'
#' @export
smallvis <- function(X,
                     k = 2,
                     scale = "absmax",
                     Y_init = "rand",
                     Y_init_sdev = NULL,
                     perplexity = 30,
                     max_iter = 1000,
                     pca = FALSE,
                     initial_dims = 50,
                     method = "tsne",
                     epoch_callback = TRUE,
                     epoch = max(1, base::round(max_iter / 10)),
                     min_cost = 0,
                     tol = 1e-7,
                     g2tol = NULL,
                     momentum = 0.5,
                     final_momentum = 0.8,
                     mom_switch_iter = stop_lying_iter + 150,
                     eta = 500,
                     min_gain = 0.01,
                     opt = list("dbd"),
                     exaggeration_factor = 1,
                     stop_lying_iter = max(1, floor(max_iter / 10)),
                     late_exaggeration_factor = 1,
                     start_late_lying_iter = max_iter - max(1, floor(max_iter / 10)),
                     iter0_cost = FALSE,
                     ee_mon_epoch = NULL,
                     ee_mon_wait = 15,
                     ee_mon_buffer = 2,
                     tol_wait = 15,
                     ret_extra = FALSE,
                     n_threads = 0,
                     use_cpp = FALSE,
                     theta = 1.0,
                     eps = .Machine$double.eps,
                     inp_kernel = NULL,
                     nn = NULL,
                     verbose = TRUE) {
  if (is.logical(epoch_callback)) {
    if (epoch_callback) {
      epoch_callback <- make_smallvis_cb(X)
    } else {
      epoch_callback <- NULL
    }
  } else if (is.function(epoch_callback)) {
    force(epoch_callback)
  }

  # The embedding method
  cost_fn <- create_cost(
    method,
    perplexity,
    eps,
    n_threads,
    use_cpp,
    theta = theta,
    inp_kernel = inp_kernel,
    nn = nn
  )
  
  if (exaggeration_factor != 1) {
    if (stop_lying_iter < 1) {
      stop("stop_lying_iter must be >= 1")
    }
  }

  if (late_exaggeration_factor != 1) {
    if (start_late_lying_iter < 1) {
      stop("start_late_lying_iter must be >= 1")
    }
  }

  if (exaggeration_factor != 1 && late_exaggeration_factor != 1) {
    if (start_late_lying_iter < stop_lying_iter) {
      stop("start_late_lying_iter must be >= stop_lying_iter")
    }
  }

  if (methods::is(pca, "character") && pca == "whiten") {
    pca <- TRUE
    whiten <- TRUE
  } else {
    whiten <- FALSE
  }
  if (pca && initial_dims < k) {
    stop(
      "Initial PCA dimensionality must be larger than desired output ",
      "dimension"
    )
  }

  start_time <- NULL
  ret_optionals <- c()
  if (methods::is(ret_extra, "character")) {
    ret_optionals <- ret_extra
    ret_extra <- TRUE
  }

  if (ret_extra) {
    start_time <- Sys.time()
  }

  if (methods::is(X, "dist")) {
    n <- attr(X, "Size")
  } else {
    if (methods::is(X, "data.frame")) {
      indexes <- which(vapply(X, is.numeric, logical(1)))
      tsmessage("Found ", length(indexes), " numeric columns")
      if (length(indexes) == 0) {
        stop("No numeric columns found")
      }
      X <- X[, indexes]
    }

    X <- scale_input(X, scale, verbose = verbose)
    X <- pca_preprocess(X, pca, whiten, initial_dims, verbose = verbose)
    n <- nrow(X)
  }

  # Check for NA
  if (any(is.na(X))) {
    stop("Input data contains NA: missing data is not allowed")
  }

  # Fail early as possible if matrix initializer is invalid
  if (methods::is(Y_init, "matrix")) {
    if (nrow(Y_init) != n || ncol(Y_init) != k) {
      stop("Y_init matrix does not match necessary configuration for X")
    }
  }

  # Optimizer
  if (opt[[1]] == "dbd" || opt[[1]] == "ndbd") {
    if (eta == "optsne") {
      eta <- n / exaggeration_factor
      tsmessage("Using opt-SNE learning rate = ", formatC(eta))
    }
    opt_list <- lmerge(
      opt,
      list(
        momentum = momentum,
        final_momentum = final_momentum,
        mom_switch_iter = mom_switch_iter,
        eta = eta,
        min_gain = min_gain,
        verbose = verbose
      )
    )
  } else {
    opt_list <- opt
  }
  opt <- opt_create(opt_list, verbose = verbose)

  # Initialize the cost function and create P
  cost_fn <- cost_init(
    cost_fn,
    X,
    max_iter = max_iter,
    verbose = verbose,
    ret_extra = ret_optionals
  )

  # Output Initialization
  if (!is.null(Y_init)) {
    if (methods::is(Y_init, "matrix")) {
      Y <- Y_init
      Y_init <- "matrix"
    } else {
      Y_init <- match.arg(
        tolower(Y_init),
        c("rand", "pca", "spca", "laplacian", "normlaplacian")
      )

      if (!is_spectral_init(Y_init)) {
        # We now handle scaling coordinates below, so spca is treated like pca
        non_spectral_init <- Y_init
        if (non_spectral_init == "spca") {
          non_spectral_init <- "pca"
        }
        Y <- init_out(non_spectral_init,
          X,
          n,
          k,
          pca_preprocessed = pca,
          verbose = verbose
        )
      } else {
        # Use normalized or unnormalized input weight for spectral initialization
        # Pretty sure it doesn't matter, given the current options
        A <- cost_fn$V
        if (is.null(A)) {
          A <- cost_fn$P
          if (!is.null(A)) {
            tsmessage("Using P for spectral initialization")
          }
        } else {
          tsmessage("Using V for spectral initialization")
        }

        if (is.null(A)) {
          if (!is.null(cost_fn$knn)) {
            tsmessage("Using knn 1 / 1 + R for spectral initialization")
            A <- 1 / (1 + cost_fn$R)
            A[cost_fn$knn == 0] <- 0
          } else if (!is.null(cost_fn$R)) {
            tsmessage("Using 1/ (1 + R) for spectral initialization")
            A <- 1 / (1 + cost_fn$R)
          } else {
            stop("No suitable input for spectral initialization")
          }
        }

        if (Y_init == "laplacian") {
          tsmessage("Initializing from Laplacian Eigenmap")
          Y <- laplacian_eigenmap(A, ndim = k)
        } else {
          tsmessage("Initializing from normalized Laplacian")
          Y <- normalized_spectral_init(A, ndim = k)
        }
      }
    }
    if (!is.null(Y_init_sdev) ||
      Y_init == "spca" || Y_init == "rand") {
      if (is.null(Y_init_sdev)) {
        Y_init_sdev <- 1e-4
      }
      tsmessage("Scaling initial coords to sdev = ", formatC(Y_init_sdev))
      Y <- shrink_coords(Y, Y_init_sdev)
    }
  }

  cost <- NA
  itercosts <- c()
  if (iter0_cost && (verbose || ret_extra)) {
    cost_eval_res <- cost_eval(cost_fn, Y)
    cost_fn <- cost_eval_res$cost
    cost <- cost_eval_res$value

    tsmessage("Iteration #0 error: ", formatC(cost))

    if (ret_extra) {
      names(cost) <- 0
      itercosts <- c(itercosts, cost)
    }
    cost_fn <- cost_clear(cost_fn)
  }

  # Display initialization
  if (!is.null(epoch_callback)) {
    do_callback(epoch_callback, Y, 0, cost, cost_fn, opt)
  }
  if (max_iter < 1) {
    return(
      ret_value(
        Y,
        ret_extra,
        method,
        X,
        scale,
        Y_init,
        iter = 0,
        cost_fn = cost_fn,
        itercosts = itercosts,
        start_time = start_time,
        optionals = ret_optionals,
        pca = ifelse(pca && !whiten, initial_dims, 0),
        whiten = ifelse(pca &&
          whiten, initial_dims, 0),
        use_cpp = use_cpp,
        n_threads = n_threads
      )
    )
  }

  opt_stages <- c()
  if (exaggeration_factor != 1) {
    opt_stages <- c(opt_stages, "early")
  }
  opt_stages <- c(opt_stages, "opt")
  if (late_exaggeration_factor != 1) {
    opt_stages <- c(opt_stages, "late")
  }
  opt_stage_idx <- 1

  opt_epoch <- epoch
  opt_tol <- tol
  old_cdiffrc <- NULL
  # Start exaggerating
  if (!is.null(cost_fn$P) && exaggeration_factor != 1) {
    cost_fn <- start_exaggerating(cost_fn, exaggeration_factor)
    if (!is.null(ee_mon_epoch)) {
      epoch <- ee_mon_epoch
      # prevent tolerance based convergence if we monitor EE
      ee_mon_tol <- 0
      tol <- ee_mon_tol
      stop_lying_iter <- max_iter
      tsmessage(
        "Applying early exaggeration factor = ",
        exaggeration_factor,
        ", epoch every ",
        epoch,
        " iterations"
      )
    }
  } else {
    stop_lying_iter <- 0
  }

  old_cost <- NULL
  tolval <- NULL
  tsmessage("Optimizing coordinates")

  # Use a while loop, so we can change max_iter inside the loop
  iter <- 0
  while (iter < max_iter) {
    iter <- iter + 1
    opt_res <- opt_step(opt, cost_fn, Y, iter)
    opt <- opt_res$opt
    cost_fn <- opt_res$cost_fn
    Y <- opt_res$Y

    if (!is.null(cost_fn$P) && iter == stop_lying_iter &&
      exaggeration_factor != 1) {
      tsmessage("Switching off exaggeration at iter ", iter)
      cost_fn <- stop_exaggerating(cost_fn, exaggeration_factor)
      opt_stage_idx <- opt_stage_idx + 1
      epoch <- opt_epoch
      tol <- opt_tol
    }

    if (!is.null(cost_fn$P) && iter == start_late_lying_iter &&
      late_exaggeration_factor != 1) {
      tsmessage(
        "Starting late exaggeration = ",
        formatC(late_exaggeration_factor),
        " at iter ",
        iter,
        " until iter ",
        max_iter
      )
      cost_fn <- start_exaggerating(cost_fn, late_exaggeration_factor)

      opt_stage_idx <- opt_stage_idx + 1
    }

    if (nnat(opt$is_terminated)) {
      tsmessage(
        "Iteration #",
        iter,
        " stopping early: optimizer reports convergence: ",
        opt$terminate$what
      )
      max_iter <- iter
      break
    }

    # Recenter after each iteration
    Y <- sweep(Y, 2, colMeans(Y))

    if ((epoch > 0 && iter %% epoch == 0) || iter == max_iter) {
      stop_early <- FALSE

      cost_eval_res <- cost_eval(cost_fn, Y, opt_res)
      cost_fn <- cost_eval_res$cost
      cost <- cost_eval_res$value

      if (!is.null(old_cost)) {
        tolval <- reltol(cost, old_cost) / epoch
      }

      if (verbose) {
        tsmessage(
          "Iteration #",
          iter,
          " error: ",
          formatC(cost),
          " ||G||2 = ",
          formatC(norm2(opt_res$G)),
          appendLF = FALSE
        )
        if (!is.null(tolval)) {
          message(" tol = ", formatC(tolval), appendLF = FALSE)
        }

        if (!is.null(opt$counts)) {
          message(" nf = ",
            opt$counts$fn,
            " ng = ",
            opt$counts$gr,
            appendLF = FALSE
          )
        }

        # special treatment for mize innards
        if (!is.null(opt$stages$gradient_descent$step_size$value)) {
          opt$stages$gradient_descent$step_size$value
          message(
            " alpha = ",
            formatC(opt$stages$gradient_descent$step_size$value),
            appendLF = FALSE
          )
        }

        if (!is.null(old_cost) && cost > old_cost) {
          message(" !", appendLF = FALSE)
        }

        if (!is.null(opt$epoch)) {
          opt <- opt$epoch(Y, iter, cost, cost_fn, opt)
        }

        message()
        utils::flush.console()
      }
      if (!is.null(epoch_callback)) {
        do_callback(epoch_callback, Y, iter, cost, cost_fn, opt)
      }

      if (ret_extra) {
        names(cost) <- iter
        itercosts <- c(itercosts, cost)
      }

      # Early stopping tests
      if (!is.null(ee_mon_epoch)) {
        cdiff <- old_cost - cost
        cdiffrc <- cdiff / old_cost
        if (length(cdiffrc) > 0 &&
          opt_stages[opt_stage_idx] == "early") {
          if (iter > ee_mon_wait &&
            length(cdiffrc) > 0 && length(old_cdiffrc) > 0 &&
            cdiffrc < old_cdiffrc) {
            if (ee_mon_buffer < 1) {
              stop_early <- TRUE
              tsmessage("Stopping early: EE relative change reached maximum")
            } else {
              ee_mon_buffer <- ee_mon_buffer - 1
            }
          }
        }
        old_cdiffrc <- cdiffrc
      }
      if (cost < min_cost) {
        stop_early <- TRUE
        tsmessage("Stopping early: cost fell below min_cost")
      }

      if (!nnat(opt$is_terminated) && !is.null(tolval) &&
        tolval < tol && cost <= old_cost &&
        (iter > stop_lying_iter + tol_wait ||
          opt_stages[opt_stage_idx] != "opt")) {
        stop_early <- TRUE
        tsmessage(
          "Stopping early: relative tolerance (",
          formatC(tol),
          ") met"
        )
      }

      # Alternative tolerance grad 2norm doesn't need cost to decrease to stop
      # (Use this for certain settings with e.g. LargeVis where numerical issues
      # can cause the cost function to increase almost negligibly)
      g2tolval <- (norm2(opt_res$G))
      if (!nnat(opt$is_terminated) &&
        !is.null(g2tol) && g2tolval < g2tol &&
        (iter > stop_lying_iter + tol_wait ||
          opt_stages[opt_stage_idx] != "opt")) {
        stop_early <- TRUE
        tsmessage(
          "Stopping early: ||G||2 tolerance (",
          formatC(g2tol),
          ") met"
        )
      }

      # Stop current stage early if we aren't making progress
      if (stop_early) {
        if (opt_stage_idx == length(opt_stages)) {
          # run out of optimization stages, so give up
          break
        }

        opt_stage <- opt_stages[opt_stage_idx]
        if (opt_stage == "early") {
          # stop early exaggeration and adjust mom_switch_iter accordingly
          # Only applies to DBD optimizer
          n_low_mom_iters <- mom_switch_iter - stop_lying_iter
          stop_lying_iter <- iter + 1
          mom_switch_iter <- stop_lying_iter + n_low_mom_iters
          opt$mom_switch_iter <- mom_switch_iter
        }

        # we know there is at least one more stage or we would have hit the
        # break earlier
        next_opt_stage <- opt_stages[opt_stage_idx + 1]

        switch(next_opt_stage,
          opt = tsmessage("Proceeding to main optimization stage"),
          late = {
            n_late_exagg_iters <- max_iter - start_late_lying_iter
            start_late_lying_iter <- iter + 1
            max_iter <- start_late_lying_iter + n_late_exagg_iters
            tsmessage("Proceeding to late exaggeration stage")
          },
          stop("BUG: unknown optimization stage '", next_opt_stage, "'")
        )
      }

      if (nnat(opt$is_terminated)) {
        break
      }

      old_cost <- cost

      # Any special custom epoch stuff
      epoch_res <- do_epoch(opt, cost_fn, iter, Y, cost)
      opt <- epoch_res$opt
      cost_fn <- epoch_res$cost
    }
  }

  if (opt_stages[opt_stage_idx] == "early") {
    cost_fn <- stop_exaggerating(cost_fn, exaggeration_factor)
  }
  if (opt_stages[opt_stage_idx] == "late") {
    cost_fn <- stop_exaggerating(cost_fn, late_exaggeration_factor)
  }

  # Recenter before output
  Y <- sweep(Y, 2, colMeans(Y))
  res <- ret_value(
    Y,
    ret_extra,
    method,
    X,
    scale,
    Y_init,
    iter,
    start_time,
    cost_fn = cost_fn,
    opt_res$G,
    perplexity,
    itercosts,
    stop_lying_iter,
    start_late_lying_iter,
    opt_list,
    opt,
    exaggeration_factor,
    late_exaggeration_factor,
    optionals = ret_optionals,
    pca = ifelse(pca && !whiten, initial_dims, 0),
    whiten = ifelse(pca && whiten, initial_dims, 0),
    use_cpp = use_cpp,
    n_threads = n_threads
  )

  res
}

create_cost <- function(method,
                        perplexity,
                        eps,
                        n_threads,
                        use_cpp,
                        theta,
                        inp_kernel,
                        nn) {
  method_names <- c(
    "tsne",
    "largevis",
    "umap",
    "tumap",
    "ntumap",
    "mmds",
    "gmmds",
    "asne",
    "ssne",
    "wtsne",
    "wssne",
    "hssne",
    "ee",
    "nerv",
    "snerv",
    "jse",
    "sjse",
    "smmds",
    "sammon",
    "tasne",
    "trmsne",
    "trsrsne",
    "tmsne",
    "arsrsne",
    "rsrjse",
    "rsrnerv",
    "btsne",
    "bssne",
    "basne",
    "btasne",
    "bnerv",
    "ballmmds",
    "knnmmds",
    "dhssne",
    "pstsne",
    "tsneu",
    "skdtsne",
    "usne",
    "cetsne",
    "tee",
    "absne",
    "chsne",
    "hlsne",
    "rklsne",
    "jssne",
    "gsne",
    "abssne",
    "bhssne",
    "bhtsne"
  )
  if (is.character(method)) {
    method <- match.arg(tolower(method), method_names)
    cost_fn <- switch(method,
      tsne = tsne(
        perplexity = perplexity,
        n_threads = n_threads,
        eps = eps,
        use_cpp = use_cpp
      ),
      bhtsne = bhtsne(
        perplexity = perplexity,
        n_threads = n_threads,
        eps = eps,
        theta = theta,
        inp_kernel = inp_kernel,
        nn = nn
      ),
      umap = umap(
        perplexity = perplexity,
        n_threads = n_threads,
        eps = eps,
        use_cpp = use_cpp
      ),
      largevis = largevis(
        perplexity = perplexity,
        n_threads = n_threads,
        eps = eps,
        use_cpp = use_cpp
      ),
      tumap = tumap(
        perplexity = perplexity,
        n_threads = n_threads,
        eps = eps,
        use_cpp = use_cpp
      ),
      ntumap = ntumap(
        perplexity = perplexity,
        n_threads = n_threads,
        eps = eps,
        use_cpp = use_cpp
      ),
      mmds = mmds(
        n_threads = n_threads,
        eps = eps,
        use_cpp = use_cpp
      ),
      gmmds = gmmds(
        k = perplexity,
        n_threads = n_threads,
        eps = eps,
        use_cpp = use_cpp
      ),
      asne = asne(
        perplexity = perplexity,
        n_threads = n_threads,
        eps = eps,
        use_cpp = use_cpp
      ),
      ssne = ssne(
        perplexity = perplexity,
        n_threads = n_threads,
        eps = eps,
        use_cpp = use_cpp
      ),
      wtsne = wtsne(
        perplexity = perplexity,
        n_threads = n_threads,
        eps = eps,
        use_cpp = use_cpp
      ),
      wssne = wssne(
        perplexity = perplexity,
        n_threads = n_threads,
        eps = eps,
        use_cpp = use_cpp
      ),
      hssne = hssne(
        perplexity = perplexity,
        n_threads = n_threads,
        eps = eps,
        use_cpp = use_cpp
      ),
      ee = ee(
        perplexity = perplexity,
        n_threads = n_threads,
        eps = eps,
        use_cpp = use_cpp
      ),
      nerv = nerv(
        perplexity = perplexity,
        n_threads = n_threads,
        eps = eps,
        use_cpp = use_cpp
      ),
      snerv = snerv(
        perplexity = perplexity,
        n_threads = n_threads,
        eps = eps,
        use_cpp = use_cpp
      ),
      jse = jse(
        perplexity = perplexity,
        n_threads = n_threads,
        eps = eps,
        use_cpp = use_cpp
      ),
      sjse = sjse(
        perplexity = perplexity,
        n_threads = n_threads,
        eps = eps,
        use_cpp = use_cpp
      ),
      smmds = smmds(
        n_threads = n_threads,
        eps = eps,
        use_cpp = use_cpp
      ),
      sammon = sammon(
        n_threads = n_threads,
        eps = eps,
        use_cpp = use_cpp
      ),
      tasne = tasne(
        perplexity = perplexity,
        n_threads = n_threads,
        eps = eps,
        use_cpp = use_cpp
      ),
      trmsne = trmsne(
        perplexity = perplexity,
        n_threads = n_threads,
        eps = eps,
        use_cpp = use_cpp
      ),
      trsrsne = trsrsne(
        perplexity = perplexity,
        n_threads = n_threads,
        eps = eps,
        use_cpp = use_cpp
      ),
      tmsne = tmsne(
        perplexity = perplexity,
        n_threads = n_threads,
        eps = eps,
        use_cpp = use_cpp
      ),
      arsrsne = arsrsne(
        perplexity = perplexity,
        n_threads = n_threads,
        eps = eps,
        use_cpp = use_cpp
      ),
      rsrjse = rsrjse(
        perplexity = perplexity,
        n_threads = n_threads,
        eps = eps,
        use_cpp = use_cpp
      ),
      rsrnerv = rsrnerv(
        perplexity = perplexity,
        n_threads = n_threads,
        eps = eps,
        use_cpp = use_cpp
      ),
      btsne = btsne(
        perplexity = perplexity,
        n_threads = n_threads,
        eps = eps,
        use_cpp = use_cpp
      ),
      bssne = bssne(
        perplexity = perplexity,
        n_threads = n_threads,
        eps = eps,
        use_cpp = use_cpp
      ),
      basne = basne(
        perplexity = perplexity,
        n_threads = n_threads,
        eps = eps,
        use_cpp = use_cpp
      ),
      btasne = btasne(
        perplexity = perplexity,
        n_threads = n_threads,
        eps = eps,
        use_cpp = use_cpp
      ),
      bnerv = bnerv(
        perplexity = perplexity,
        n_threads = n_threads,
        eps = eps,
        use_cpp = use_cpp
      ),
      ballmmds = ballmmds(
        n_threads = n_threads,
        eps = eps,
        use_cpp = use_cpp
      ),
      knnmmds = knnmmds(
        k = perplexity,
        n_threads = n_threads,
        eps = eps,
        use_cpp = use_cpp
      ),
      dhssne = dhssne(
        perplexity = perplexity,
        n_threads = n_threads,
        eps = eps,
        use_cpp = use_cpp
      ),
      tsneu = tsneu(
        perplexity = perplexity,
        n_threads = n_threads,
        eps = eps,
        use_cpp = use_cpp
      ),
      pstsne = pstsne(
        perplexity = perplexity,
        n_threads = n_threads,
        eps = eps,
        use_cpp = use_cpp
      ),
      skdtsne = skdtsne(
        perplexity = perplexity,
        n_threads = n_threads,
        eps = eps,
        use_cpp = use_cpp
      ),
      usne = usne(
        perplexity = perplexity,
        n_threads = n_threads,
        eps = eps,
        use_cpp = use_cpp
      ),
      cetsne = cetsne(
        perplexity = perplexity,
        n_threads = n_threads,
        eps = eps,
        use_cpp = use_cpp
      ),
      tee = tee(
        perplexity = perplexity,
        n_threads = n_threads,
        eps = eps,
        use_cpp = use_cpp
      ),
      absne = absne(
        perplexity = perplexity,
        n_threads = n_threads,
        eps = eps,
        use_cpp = use_cpp
      ),
      chsne = chsne(
        perplexity = perplexity,
        n_threads = n_threads,
        eps = eps,
        use_cpp = use_cpp
      ),
      hlsne = hlsne(
        perplexity = perplexity,
        n_threads = n_threads,
        eps = eps,
        use_cpp = use_cpp
      ),
      gsne = gsne(
        perplexity = perplexity,
        n_threads = n_threads,
        eps = eps,
        use_cpp = use_cpp
      ),
      rklsne = rklsne(
        perplexity = perplexity,
        n_threads = n_threads,
        eps = eps,
        use_cpp = use_cpp
      ),
      jssne = jssne(
        perplexity = perplexity,
        n_threads = n_threads,
        eps = eps,
        use_cpp = use_cpp
      ),
      abssne = abssne(
        perplexity = perplexity,
        n_threads = n_threads,
        eps = eps,
        use_cpp = use_cpp
      ),
      bhssne = bhssne(
        perplexity = perplexity,
        n_threads = n_threads,
        eps = eps,
        use_cpp = use_cpp
      ),
      stop("BUG: someone forgot to implement option: '", method, "'")
    )
  } else {
    if (is.list(method)) {
      methodlist <- method
      method_name <- methodlist[[1]]
      methodlist[[1]] <- NULL
      method <- match.arg(tolower(method_name), method_names)
      if (exists(method_name)) {
        fn <- get(method_name)
      } else {
        stop("Unknown method: '", method_name, "'")
      }
    } else if (is.function(method)) {
      fn <- method
      methodlist <- list()
    } else {
      stop("Bad method parameter type")
    }

    param_names <- names(formals(fn))
    if ("perplexity" %in% param_names) {
      methodlist$perplexity <- perplexity
    }
    if ("k" %in% param_names) {
      methodlist$k <- perplexity
    }
    cost_fn <- do.call(fn, methodlist)
  }
}


# Result Export -----------------------------------------------------------

# Prepare the return value.
# If ret_extra is TRUE, return a list with lots of extra info.
# Otherwise, Y is returned directly.
# If ret_extra is TRUE and iter > 0, then all the NULL-default parameters are
# expected to be present. If iter == 0 then the return list will contain only
# scaling and initialization information.
ret_value <- function(Y,
                      ret_extra,
                      method,
                      X,
                      scale,
                      Y_init,
                      iter,
                      start_time = NULL,
                      cost_fn = NULL,
                      G = NULL,
                      perplexity = NULL,
                      pca = 0,
                      whiten = 0,
                      itercosts = NULL,
                      stop_lying_iter = NULL,
                      start_late_lying_iter = NULL,
                      opt_input = NULL,
                      opt_res = NULL,
                      exaggeration_factor = 1,
                      late_exaggeration_factor = 1,
                      optionals = c(),
                      use_cpp = FALSE,
                      n_threads = 1) {
  attr(Y, "dimnames") <- NULL
  if (ret_extra) {
    end_time <- Sys.time()

    if (methods::is(X, "dist")) {
      N <- attr(X, "Size")
      origD <- NULL
    } else {
      N <- nrow(X)
      origD <- ncol(X)
    }

    res <- list(
      Y = Y,
      N = N,
      origD = origD,
      scale = scale,
      Y_init = Y_init,
      method = method,
      iter = iter,
      time_secs = as.numeric(end_time - start_time, units = "secs")
    )

    if (!is.null(G)) {
      res$G2norm <- norm2(G)
    }

    if (pca > 0) {
      res$pca_dims <- pca
    } else if (whiten > 0) {
      res$whiten_dims <- whiten
    }

    if (is.null(cost_fn$pcost)) {
      cost_fn <- cost_grad(cost_fn, Y)
      cost_fn <- cost_point(cost_fn, Y)
    }
    res$costs <- cost_fn$pcost

    if (!is.null(opt_input)) {
      res$opt <- opt_input
      if (!is.null(opt_res$counts)) {
        res$opt$counts <- opt_res$counts
      }
    }


    # Don't report exaggeration settings if they didn't do anything
    if (exaggeration_factor == 1 || late_exaggeration_factor == 1) {
      if (exaggeration_factor == 1) {
        exaggeration_factor <- NULL
        stop_lying_iter <- NULL
      }

      if (late_exaggeration_factor == 1) {
        late_exaggeration_factor <- NULL
        start_late_lying_iter <- NULL
      }
    } else {
      # Don't report start_late_lying_iter if we never got there
      if (start_late_lying_iter > iter) {
        start_late_lying_iter <- NULL
        late_exaggeration_factor <- NULL
      }
    }

    res <- c(
      res,
      list(
        perplexity = perplexity,
        itercosts = itercosts,
        exaggeration_factor = exaggeration_factor,
        stop_lying_iter = stop_lying_iter,
        late_exaggeration_factor = late_exaggeration_factor,
        start_late_lying_iter = start_late_lying_iter
      )
    )

    # If using the Intrinsic Dimensionality method, use the chosen perplexity
    if (!is.null(cost_fn$idp)) {
      res$perplexity <- cost_fn$idp
    }


    optionals <- tolower(unique(optionals))
    for (o in optionals) {
      exported <- NULL
      if (!is.null(cost_fn)) {
        # Could be NULL if max_iter was 0
        exported <- cost_fn$export(cost_fn, o)
      }
      if (!is.null(exported)) {
        if (nchar(o) < 3) {
          res[[toupper(o)]] <- exported
        } else {
          res[[o]] <- exported
        }
      }
      # For DX and DY, we can calculate if not done so already
      else if (o == "dx") {
        if (methods::is(X, "dist")) {
          res$DX <- X
        } else {
          res$DX <- calc_d(X, use_cpp = use_cpp, n_threads = n_threads)
        }
      } else if (o == "dy") {
        res$DY <- calc_d(Y, use_cpp = use_cpp, n_threads = n_threads)
      }

      if (o == "x") {
        res$X <- X
      }
    }

    res <- remove_nulls(res)
    res
  } else {
    Y
  }
}
