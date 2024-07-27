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
#'   0, which means no multi-threading.
#' @param eps Set epsilon for avoiding division-by-zero errors. Default is
#'   \code{.Machine$double.eps}, but if you see inconsistent convergence results
#'   with optimizer that should be reducing the cost each iteration, then try
#'   setting this to a larger value, e.g. between \code{1e-3 - 1e-9}.
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
#' colors = rainbow(length(unique(iris$Species)))
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
smallvis <- function(X, k = 2, scale = "absmax", 
                     Y_init = "rand", Y_init_sdev = NULL,
                 perplexity = 30, max_iter = 1000,
                 pca = FALSE, initial_dims = 50,
                 method = "tsne",
                 epoch_callback = TRUE,
                 epoch = max(1, base::round(max_iter / 10)),
                 min_cost = 0, tol = 1e-7, g2tol = NULL,
                 momentum = 0.5, final_momentum = 0.8, 
                 mom_switch_iter = stop_lying_iter + 150,
                 eta = 500, min_gain = 0.01,
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
                 eps = .Machine$double.eps,
                 verbose = TRUE) {

  if (is.logical(epoch_callback)) {
    if (epoch_callback) {
      epoch_callback <- make_smallvis_cb(X)
    }
    else {
      epoch_callback <- NULL
    }
  }
  else if (is.function(epoch_callback)) {
    force(epoch_callback)
  }

  # The embedding method
  method_names <- c("tsne", "largevis", "umap", "tumap",
                    "ntumap", "mmds", "gmmds",
                    "asne", "ssne", "wtsne", "wssne",
                    "hssne", "ee", "nerv", "snerv", "jse", "sjse",
                    "smmds", "sammon",
                    "tasne", "trmsne", "trsrsne", "tmsne", "arsrsne",
                    "rsrjse", "rsrnerv",
                    "btsne", "bssne", "basne", "btasne", "bnerv",
                    "ballmmds", "knnmmds",
                    "dhssne", "pstsne", "tsneu",
                    "skdtsne", "usne", "cetsne",
                    "tee", "absne", "chsne", "hlsne", "rklsne", "jssne",
                    "gsne", "abssne", "bhssne")
  if (is.character(method)) {
    method <- match.arg(tolower(method), method_names)
    cost_fn <- switch(method,
         tsne = tsne(perplexity = perplexity, n_threads = n_threads, eps = eps),
         umap = umap(perplexity = perplexity, n_threads = n_threads, eps = eps),
         largevis = largevis(perplexity = perplexity, n_threads = n_threads, eps = eps),
         tumap = tumap(perplexity = perplexity, n_threads = n_threads, eps = eps),
         ntumap = ntumap(perplexity = perplexity, n_threads = n_threads, eps = eps),
         mmds = mmds(eps = eps),
         gmmds = gmmds(k = perplexity, n_threads = n_threads, eps = eps),
         asne = asne(perplexity = perplexity, n_threads = n_threads, eps = eps),
         ssne = ssne(perplexity = perplexity, n_threads = n_threads, eps = eps),
         wtsne = wtsne(perplexity = perplexity, n_threads = n_threads, eps = eps),
         wssne = wssne(perplexity = perplexity, n_threads = n_threads, eps = eps),
         hssne = hssne(perplexity = perplexity, n_threads = n_threads, eps = eps),
         ee = ee(perplexity = perplexity, n_threads = n_threads, eps = eps),
         nerv = nerv(perplexity = perplexity, n_threads = n_threads, eps = eps),
         snerv = snerv(perplexity = perplexity, n_threads = n_threads, eps = eps),
         jse = jse(perplexity = perplexity, n_threads = n_threads, eps = eps),
         sjse = sjse(perplexity = perplexity, n_threads = n_threads, eps = eps),
         smmds = smmds(eps = eps),
         sammon = sammon(eps = eps),
         tasne = tasne(perplexity = perplexity, n_threads = n_threads, eps = eps),
         trmsne = trmsne(perplexity = perplexity, n_threads = n_threads, eps = eps),
         trsrsne = trsrsne(perplexity = perplexity, n_threads = n_threads, eps = eps),
         tmsne = tmsne(perplexity = perplexity, n_threads = n_threads, eps = eps),
         arsrsne = arsrsne(perplexity = perplexity, n_threads = n_threads, eps = eps),
         rsrjse = rsrjse(perplexity = perplexity, n_threads = n_threads, eps = eps),
         rsrnerv = rsrnerv(perplexity = perplexity, n_threads = n_threads, eps = eps),
         btsne = btsne(perplexity = perplexity, n_threads = n_threads, eps = eps),
         bssne = bssne(perplexity = perplexity, n_threads = n_threads, eps = eps),
         basne = basne(perplexity = perplexity, n_threads = n_threads, eps = eps),
         btasne = btasne(perplexity = perplexity, n_threads = n_threads, eps = eps),
         bnerv = bnerv(perplexity = perplexity, n_threads = n_threads, eps = eps),
         ballmmds = ballmmds(eps = eps),
         knnmmds = knnmmds(k = perplexity, n_threads = n_threads, eps = eps),
         dhssne = dhssne(perplexity = perplexity, n_threads = n_threads, eps = eps),
         tsneu = tsneu(perplexity = perplexity, n_threads = n_threads, eps = eps),
         pstsne = pstsne(perplexity = perplexity, n_threads = n_threads, eps = eps),
         skdtsne = skdtsne(perplexity = perplexity, n_threads = n_threads, eps = eps),
         usne = usne(perplexity = perplexity, n_threads = n_threads, eps = eps),
         cetsne = cetsne(perplexity = perplexity, n_threads = n_threads, eps = eps),
         tee = tee(perplexity = perplexity, n_threads = n_threads, eps = eps),
         absne = absne(perplexity = perplexity, n_threads = n_threads, eps = eps),
         chsne = chsne(perplexity = perplexity, n_threads = n_threads, eps = eps),
         hlsne = hlsne(perplexity = perplexity, n_threads = n_threads, eps = eps),
         gsne = gsne(perplexity = perplexity, n_threads = n_threads, eps = eps),
         rklsne = rklsne(perplexity = perplexity, n_threads = n_threads, eps = eps),
         jssne = jssne(perplexity = perplexity, n_threads = n_threads, eps = eps),      
         abssne = abssne(perplexity = perplexity, n_threads = n_threads, eps = eps),
         bhssne = bhssne(perplexity = perplexity, n_threads = n_threads, eps = eps),
         stop("BUG: someone forgot to implement option: '", method, "'")
    )
  }
  else {
    if (is.list(method)) {
      methodlist <- method
      method_name <- methodlist[[1]]
      methodlist[[1]] <- NULL
      method <- match.arg(tolower(method_name), method_names)
      if (exists(method_name)) {
        fn <- get(method_name)
      }
      else {
        stop("Unknown method: '", method_name, "'")
      }
    }
    else if (is.function(method)) {
      fn <- method
      methodlist <- list()
    }
    else {
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

  if (class(pca) == "character" && pca == "whiten") {
    pca <- TRUE
    whiten <- TRUE
  }
  else {
    whiten <- FALSE
  }
  if (pca && initial_dims < k) {
    stop("Initial PCA dimensionality must be larger than desired output ",
         "dimension")
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
  }
  else {
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
    opt_list <- lmerge(opt, list(momentum = momentum,
                       final_momentum = final_momentum,
                       mom_switch_iter = mom_switch_iter,
                       eta = eta, min_gain = min_gain,
                       verbose = verbose))
  }
  else {
    opt_list <- opt
  }
  opt <- opt_create(opt_list, verbose = verbose)

  # Initialize the cost function and create P
  cost_fn <- cost_init(cost_fn, X, max_iter = max_iter, verbose = verbose,
                       ret_extra = ret_optionals)

  # Output Initialization
  if (!is.null(Y_init)) {
    if (methods::is(Y_init, "matrix")) {
      Y <- Y_init
      Y_init <- "matrix"
    }
    else {
      Y_init <- match.arg(tolower(Y_init),
                          c("rand", "pca", "spca", "laplacian", "normlaplacian"))

      if (!is_spectral_init(Y_init)) {
        # We now handle scaling coordinates below, so spca is treated like pca
        non_spectral_init <- Y_init
        if (non_spectral_init == "spca") {
          non_spectral_init <- "pca"
        }
        Y <- init_out(non_spectral_init, X, n, k, pca_preprocessed = pca,
                      verbose = verbose)
      }
      else {
        # Use normalized or unnormalized input weight for spectral initialization
        # Pretty sure it doesn't matter, given the current options
        A <- cost_fn$V
        if (is.null(A)) {
          A <- cost_fn$P
          if (!is.null(A)) {
            tsmessage("Using P for spectral initialization")
          }        
        }
        else {
          tsmessage("Using V for spectral initialization")
        }
        
        if (is.null(A)) {
          if (!is.null(cost_fn$knn)) {
            tsmessage("Using knn 1 / 1 + R for spectral initialization")
            A <- 1 / (1 + cost_fn$R)
            A[cost_fn$knn == 0] <- 0
          }
          else if (!is.null(cost_fn$R)) {
            tsmessage("Using 1/ (1 + R) for spectral initialization")
            A <- 1 / (1 + cost_fn$R)
          }
          else {
            stop("No suitable input for spectral initialization")
          }
        }
        
        if (Y_init == "laplacian") {
          tsmessage("Initializing from Laplacian Eigenmap")
          Y <- laplacian_eigenmap(A, ndim = k)
        }
        else {
          tsmessage("Initializing from normalized Laplacian")
          Y <- normalized_spectral_init(A, ndim = k)
        }
      }
    }
    if (!is.null(Y_init_sdev) || Y_init == "spca" || Y_init == "rand") {
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
    return(ret_value(Y, ret_extra, method, X, scale, Y_init, iter = 0,
                     cost_fn = cost_fn, itercosts = itercosts,
                     start_time = start_time, optionals = ret_optionals,
                     pca = ifelse(pca && !whiten, initial_dims, 0),
                     whiten = ifelse(pca && whiten, initial_dims, 0)))
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
      tsmessage("Applying early exaggeration factor = ", exaggeration_factor, 
                ", epoch every ", epoch, " iterations")
    }
  }
  else {
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
      tsmessage("Starting late exaggeration = ",
              formatC(late_exaggeration_factor), " at iter ", iter,
              " until iter ", max_iter)
      cost_fn <- start_exaggerating(cost_fn, late_exaggeration_factor)
      
      opt_stage_idx <- opt_stage_idx + 1
    }

    if (nnat(opt$is_terminated)) {
      tsmessage("Iteration #", iter,
                " stopping early: optimizer reports convergence: ",
                opt$terminate$what)
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
        tsmessage("Iteration #", iter, " error: ",
                  formatC(cost)
                  , " ||G||2 = ", formatC(norm2(opt_res$G))
                  , appendLF = FALSE
                  )
        if (!is.null(tolval)) {
          message(" tol = ", formatC(tolval), appendLF = FALSE)
        }
        
        if (!is.null(opt$counts)) {
          message(" nf = ", opt$counts$fn, " ng = ", opt$counts$gr,
                  appendLF = FALSE)
        }

        # special treatment for mize innards
        if (!is.null(opt$stages$gradient_descent$step_size$value)) {
          opt$stages$gradient_descent$step_size$value
          message(" alpha = ", 
                  formatC(opt$stages$gradient_descent$step_size$value), 
                  appendLF = FALSE)
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
        if (length(cdiffrc) > 0 && opt_stages[opt_stage_idx] == "early") {
          if (iter > ee_mon_wait && length(cdiffrc) > 0 && length(old_cdiffrc) > 0 &&
              cdiffrc < old_cdiffrc) {
            if (ee_mon_buffer < 1) {
              stop_early <- TRUE
              tsmessage("Stopping early: EE relative change reached maximum")
            }
            else {
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
          (iter > stop_lying_iter + tol_wait 
          || opt_stages[opt_stage_idx] != "opt")) {
        stop_early <- TRUE
        tsmessage("Stopping early: relative tolerance (", formatC(tol), ") met")
      }

      # Alternative tolerance grad 2norm doesn't need cost to decrease to stop
      # (Use this for certain settings with e.g. LargeVis where numerical issues
      # can cause the cost function to increase almost negligibly)
      g2tolval <- (norm2(opt_res$G))
      if (!nnat(opt$is_terminated) && !is.null(g2tol) && g2tolval < g2tol &&
          (iter > stop_lying_iter + tol_wait 
           || opt_stages[opt_stage_idx] != "opt")) {
        stop_early <- TRUE
        tsmessage("Stopping early: ||G||2 tolerance (", formatC(g2tol), ") met")
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
  res <- ret_value(Y, ret_extra, method, X, scale, Y_init, iter, start_time,
            cost_fn = cost_fn, opt_res$G,
            perplexity, itercosts,
            stop_lying_iter, start_late_lying_iter, opt_list, opt,
            exaggeration_factor, late_exaggeration_factor,
            optionals = ret_optionals,
            pca = ifelse(pca && !whiten, initial_dims, 0),
            whiten = ifelse(pca && whiten, initial_dims, 0))

  res
}

#' Best t-SNE Result From Multiple Initializations
#'
#' Run t-SNE multiple times from a random initialization, and return the
#' embedding with the lowest cost.
#'
#' This function ignores any value of \code{Y_init} you set, and uses
#' \code{Y_init = "rand"}.
#'
#' @param nrep Number of repeats.
#' @param keep_all If \code{TRUE}, then the return value is a list of lists,
#' indexed from 1 .. \code{nrep}, with each entry the result from each
#' \code{\link{smallvis}} run. Otherwise just the result with the lowest error
#' is returned.
#' @param ... Arguments to apply to each \code{\link{smallvis}} run.
#' @return The \code{\link{smallvis}} result with the lowest final cost, or
#' if \code{keep_all} is \code{TRUE} all results as a list, indexed as 1 ..
#' \code{nrep}. If \code{ret_extra} is not \code{FALSE}, then the final costs for all
#' \code{nrep} runs are also included in the return value list as a vector
#' called \code{all_costs}. In this case, if \code{keep_all} is \code{TRUE}, then
#' \code{all_costs} appears as an extra item on all results. Additionally,
#' each result will have an extra entry \code{best_rep}, giving the index of the
#' result with the lowest cost.
#' @examples
#' \dontrun{
#' # Return best result out of five random initializations
#' tsne_iris_best <- smallvis_rep(nrep = 5, X = iris, perplexity = 50, method = "tsne",
#'                                ret_extra = TRUE)
#' # How much do the costs vary between runs?
#' range(tsne_iris_best$all_costs)
#' # Display best embedding found
#' plot(tsne_iris_best$Y)
#'
#' # Keep all results
#' # First result is in tsne_iris_rep[[1]], second in tsne_iris_rep[[2]] etc.
#' tsne_iris_rep <- smallvis_rep(nrep = 5, X = iris, perplexity = 50, method = "tsne",
#'                               ret_extra = TRUE, keep_all = TRUE)
#' # Index of result with smallest error is in special list item 'best_rep'
#' best_iris <- tsne_iris_rep[[tsne_iris_rep[[1]]$best_rep]]
#'
#' }
#' @export
smallvis_rep <- function(nrep = 10, keep_all = FALSE, ...) {
  if (nrep < 1) {
    stop("nrep must be 1 or greater")
  }
  varargs <- list(...)
  best_res <- NULL
  best_cost <- Inf
  all_costs <- c()
  # Keep requested return type for final result
  ret_extra <- varargs$ret_extra

  # always return extra so we can find the cost
  if (!should_ret_extra(ret_extra)) {
    varargs$ret_extra <- TRUE
  }

  varargs$Y_init <- "rand"
  ret <- list()

  for (i in 1:nrep) {
    # If verbose is not explicitly set to FALSE, it's TRUE by default
    if (nnat(varargs$verbose) || is.null(varargs$verbose)) {
      tsmessage("Starting embedding # ", i, " of ", nrep)
    }
    res <- do.call(smallvis, varargs)

    if (keep_all) {
      ret[[i]] <- res
    }

    final_cost <- res$itercosts[length(res$itercosts)]
    names(final_cost) <- NULL
    all_costs <- c(all_costs, final_cost)

    if (!keep_all) {
      if (final_cost < best_cost) {
        best_cost <- final_cost
        best_res <- res
      }
    }
  }

  if (keep_all) {
    # if keep_all is TRUE and we asked for extra return info
    # also add the final costs and the index of the best result to each result
    if (should_ret_extra(ret_extra)) {
      best_rep <- which.min(all_costs)
      for (i in 1:nrep) {
        ret[[i]]$all_costs <- all_costs
        ret[[i]]$best_rep <- best_rep
      }
    }
    # otherwise just store the Y-coordinates of each result
    else {
      for (i in 1:nrep) {
        ret[[i]] <- ret[[i]]$Y
      }
    }
  }
  else {
    # Only keeping one result
    if (should_ret_extra(ret_extra)) {
      # store info about other results on best result list
      best_res$all_costs <- all_costs
      ret <- best_res
    }
    else {
      ret <- best_res$Y
    }
  }
  ret
}

# If ret_extra is NULL or FALSE, we aren't returning extra info
# If ret_extra is TRUE or a vector, we are returning extra info
should_ret_extra <- function(ret_extra) {
  !is.null(ret_extra) && ((methods::is(ret_extra, "logical") && ret_extra) ||
    methods::is(ret_extra, "character"))
}

#' Dimensionality Reduction With Perplexity Stepping
#'
#' Carry out dimensionality reduction of a (small) dataset using one of a
#' variety of neighbor embedding methods, using a decreasing value of
#' perplexity to avoid bad local minima.
#'
#' This function uses ideas similar to those in the NeRV (Venna et al., 2010)
#' and JSE (Lee et al., 2013), where to avoid local minima, the initial
#' optimization steps use affinities with larger bandwidths (NeRV) or larger
#' perplexity values (JSE). This implementation uses a series of decreasing
#' perplexity values, as in JSE.
#'
#' For details on the arguments that can be passed to the dimensionality
#' reduction routine, see the help text for \code{\link{smallvis}}.
#'
#' To avoid spending too much extra time in perplexity calibrations, the extra
#' perplexities start at the power of 2 closest to, but not greater than,
#' half the dataset size (in terms of number of objects). Further calibrations
#' are then carried out halving the perplexity each time, until the perplexity
#' specified by the user is reached.
#'
#' The number of iterations spent in the larger perplexity values is specified
#' by the \code{step_iter} parameter. This determines the total number
#' of iterations, e.g. if \code{step_iter = 250} and extra optimizations
#' at a perplexity of 1024, 512, 256, 128 and 64 will be carried out, these will
#' run for 50 iterations each. To keep the number of iterations equivalent to
#' that used by a single run of \code{\link{smallvis}}, the value of
#' \code{step_iter} is subtracted from the value of \code{max_iter} before
#' the optimization at the target perplexity is carried out, e.g. if
#' \code{max_iter = 1000} and \code{step_iter = 250}, the final
#' optimization will run for 750 iterations only.
#'
#' Any value of \code{tol}, \code{exaggeration_factor} and
#' \code{stop_lying_iter} provided is used only with the final optimization.
#'
#' @param step_iter Number of iterations to carry out the perplexity
#'  stepping. Must be < the value of \code{max_iter}.
#' @param ... Arguments to be passed to \code{\link{smallvis}}. See 'Details'
#'  for information on which arguments may be modified or ignored during certain
#'  parts of the embedding.
#' @return The result of the final run of \code{\link{smallvis}} at the target
#'   perplexity.
#' @examples
#' \dontrun{
#' # t-SNE on the iris with L-BFGS optimization
#' # The 1000 max_iter is split between 250 iterations at perplexity = 64
#' # and then 750 iterations at perplexity = 40.
#' iris_lbfgs_pstep <- smallvis_perpscale(
#'   step_iter = 250, X = iris, scale = FALSE, verbose = TRUE, Y_init = "spca",
#'   ret_extra = c("DX", "DY"), perplexity = 40, max_iter = 1000, opt = list("l-bfgs"))
#' }
#' @export
#' @references
#' Venna, J., Peltonen, J., Nybo, K., Aidos, H., & Kaski, S. (2010).
#' Information retrieval perspective to nonlinear dimensionality reduction for
#' data visualization.
#' \emph{Journal of Machine Learning Research}, \emph{11}, 451-490.
#'
#' Lee, J. A., Renard, E., Bernard, G., Dupont, P., & Verleysen, M. (2013).
#' Type 1 and 2 mixtures of Kullback-Leibler divergences as cost functions in
#' dimensionality reduction based on similarity preservation.
#' \emph{Neurocomputing}, \emph{112}, 92-108.
smallvis_perpstep <- function(step_iter = 250, ...) {
  varargs <- list(...)
  max_iter <- varargs$max_iter
  if (is.null(max_iter)) {
    max_iter <- 1000
  }
  if (max_iter <= step_iter) {
    stop("max_iter must be > step_iter")
  }

  target_perplexity <- varargs$perplexity
  if (is.null(target_perplexity)) {
    target_perplexity <- 30
  }

  X <- varargs$X
  if (methods::is(X, "dist")) {
    n <- attr(X, "Size")
  }
  else {
    n <- nrow(X)
  }
  perps <- scale_perps(n = n, target_perp = target_perplexity)

  nperps <- length(perps)
  if (nperps > 0) {
    max_iter_step <- max(1, floor(step_iter / nperps))
    max_iter_target <- max(1, max_iter - step_iter)

    # Save/Modify some options between step iterations and final optimization
    ret_extra <- varargs$ret_extra
    varargs$ret_extra <- FALSE
    varargs$max_iter <- max_iter_step
    tol <- varargs$tol
    varargs$tol <- 0
    exaggeration_factor <- varargs$exaggeration_factor
    varargs$exaggeration_factor <- 1

    epoch <- varargs$epoch
    varargs$epoch <- max_iter_step

    # Loop over initial perplexities
    res <- NULL
    for (i in 1:nperps) {
      if (nnat(varargs$verbose)) {
        tsmessage("Optimizing at step perplexity ", formatC(perps[i]),
                  " for ", max_iter_step, " iterations")
      }
      varargs$perplexity <- perps[i]
      if (i > 1) {
        varargs$Y_init <- res
      }
      res <- do.call(smallvis, varargs)
    }

    varargs$Y_init <- res
    varargs$max_iter <- max_iter_target
    # Put the old arguments back before final optimization
    varargs$perplexity <- target_perplexity
    if (!is.null(ret_extra)) {
      varargs$ret_extra <- ret_extra
    }
    if (!is.null(tol)) {
      varargs$tol <- tol
    }
    if (!is.null(epoch)) {
      varargs$epoch <- epoch
    }
    if (!is.null(exaggeration_factor)) {
      varargs$exaggeration_factor <- exaggeration_factor
    }
  }

  if (nnat(varargs$verbose)) {
    tsmessage("Optimizing at target perplexity ", formatC(target_perplexity),
              " for ", max_iter_target, " iterations")
  }
  do.call(smallvis, varargs)
}

# Utility function for perplexity step
scale_perps <- function(n, target_perp) {
  max_perp <- n / 2
  max_perp <- 2 ^ floor(log(max_perp, 2))
  perp <- max_perp

  if (max_perp > target_perp) {
    perps <- c()
    while (perp > target_perp) {
      perps <- c(perps, perp)
      perp <- perp / 2
    }
  }

  perps
}


# Input Preprocess --------------------------------------------------------

# Scale X according to various strategies
scale_input <- function(X, scale, verbose = FALSE) {
  if (is.null(scale)) {
    scale <- "none"
  }
  if (is.logical(scale)) {
    if (scale) {
      scale <- "scale"
    }
    else {
      scale <- "none"
    }
  }
  scale <- match.arg(tolower(scale), c("none", "scale", "range", "absmax"))

  switch(scale,
         range = {
           tsmessage("Range scaling X")
           X <- as.matrix(X)
           X <- X - min(X)
           X <- X / max(X)
         },
         absmax = {
           tsmessage("Normalizing by abs-max")
           X <- base::scale(X, scale = FALSE)
           X <- X / abs(max(X))
         },
         scale = {
           tsmessage("Scaling to zero mean and unit variance")
           X <- Filter(stats::var, X)
           tsmessage("Kept ", ncol(X), " non-zero-variance columns")
           X <- base::scale(X, scale = TRUE)
         },
         none = {
           X <- as.matrix(X)
         }
  )
  X
}

# Reduce input dimensionality via PCA and also optionally whiten data
pca_preprocess <- function(X, pca, whiten, initial_dims, verbose = FALSE) {
  # We won't do PCA if the rank of the input is less than the requested
  # initial dimensionality
  if (pca) {
    pca <- min(nrow(X), ncol(X)) >= initial_dims
  }
  if (pca) {
    if (whiten) {
      tsmessage("Reducing initial dimensionality with PCA and ",
                "whitening to ", initial_dims, " dims")
      X <- pca_whiten(X = X, ncol = initial_dims, verbose = verbose)
    }
    else {
      tsmessage("Reducing initial dimensionality with PCA to ",
                initial_dims, " dims")
      X <- pca_scores(X = X, ncol = initial_dims, verbose = verbose)
    }
  }
  X
}

# Output Initialization ---------------------------------------------------

# Initialization of the output coordinates
init_out <- function(Y_init, X, n, ndim, pca_preprocessed, verbose = FALSE) {
  switch(Y_init,
         pca = {
           tsmessage("Initializing from PCA scores")
           pca_init(X, ndim, pca_preprocessed, verbose)
         },
         rand = {
           tsmessage("Initializing from random Gaussian")
           matrix(stats::rnorm(ndim * n, sd = 1), n)
         }
  )
}

pca_init <- function(X, ndim, pca_preprocessed, verbose = FALSE) {
  # If we've already done PCA, we can just take the first two columns
  if (pca_preprocessed) {
    X[, 1:2]
  }
  else {
    pca_scores(X, ncol = ndim, verbose = verbose)
  }
}

# Laplacian Eigenmap (Belkin & Niyogi, 2002)
# Original formulation solves the generalized eigenvalue problem of the
# unnormalized graph Laplacian: L v = lambda D v, where L = D - A
# and uses the bottom eigenvectors v that result
# (ignoring the constant eigenvector associated with the smallest eigenvalue).
#
# This is equivalent to using the top eigenvectors from the usual
# eigendecomposition of a row-normalized Laplacian P = D^-1 A: P v = lambda' v
# so we don't need to depend on an external package for generalized eigenvalues.
# Note that while the eigenvectors are the same, the eigenvalues are
# different: lambda' = 1 - lambda, but we don't use them with Laplacian
# Eigenmaps anyway.
#
# As we only need to calculate the top ndim + 1 eigenvectors (i.e. normally 3)
# it's incredibly wasteful to calculate all of them. Therefore, if the
# RSpectra library is available, we use that instead, which allows for only the
# top eigenvectors to be extracted. Otherwise, use the slower eigen routine.
# A must be symmetric and positive semi definite, but not necessarily
# normalized in any specific way.
laplacian_eigenmap <- function(A, ndim = 2, use_RSpectra = TRUE) {
  # Equivalent to: D <- diag(colSums(A)); M <- solve(D) %*% A
  # This effectively row-normalizes A: colSums is normally faster than rowSums
  # and because A is symmetric, they're equivalent
  M <- A / colSums(A)
  if (use_RSpectra && requireNamespace("RSpectra", quietly = TRUE,
                                       warn.conflicts = FALSE)) {
    tsmessage("Using RSpectra for eigenvectors")
    Re(RSpectra::eigs(M, k = ndim + 1)$vectors[, 2:(ndim + 1)])
  }
  else {
    tsmessage("Using eigen for eigenvectors")
    eigen(M, symmetric = FALSE)$vectors[, 2:(ndim + 1)]
  }
}

# Use a normalized Laplacian. The UMAP approach, taken from version 0.2.1.
normalized_spectral_init <- function(A, ndim = 2, use_RSpectra = TRUE) {
  n <- nrow(A)
  # Normalized Laplacian: clear and close to UMAP code, but very slow in R
  # I <- diag(1, nrow = n, ncol = n)
  # D <- diag(1 / sqrt(colSums(A)))
  # L <- I - D %*% A %*% D

  # A lot faster (order of magnitude when n = 1000)
  Dsq <- sqrt(colSums(A))
  L <- -t(A / Dsq) / Dsq
  diag(L) <- 1 + diag(L)

  if (use_RSpectra && requireNamespace("RSpectra", quietly = TRUE,
                                       warn.conflicts = FALSE)) {
    tsmessage("Using RSpectra for eigenvectors")
    k <- ndim + 1
    ncv <- max(2 * k + 1, floor(sqrt(n)))
    opt <- list(
      ncv = ncv,
      maxitr = 5 * n,
      tol = 1e-4
    )
    res <- RSpectra::eigs(L, k = k, which = "SM", opt = opt)
    vec_indices <- rev(order(res$values, decreasing = TRUE)[1:ndim])
    res <- Re(res$vectors[, vec_indices])
  }
  else {
    tsmessage("Using eigen for eigenvectors")
    res <- eigen(L, symmetric = FALSE)
    vec_indices <- order(res$values, decreasing = FALSE)[2:(ndim + 1)]
    res <- Re(res$vectors[, vec_indices])
  }
  res
}

is_spectral_init <- function(init) {
  tolower(init) %in% c("laplacian", "normlaplacian")
}

# Rescale embedding so that the standard deviation is the specified value.
# Default gives initialization like t-SNE, but not random.
shrink_coords <- function(X, sdev = 1e-4) {
  scale(X, scale = apply(X, 2, stats::sd) / sdev)
}


# Epoch Functions ---------------------------------------------------------

do_epoch <- function(opt, cost, iter, Y, fn_val) {
  if (!is.null(cost$epoch)) {
    res <- cost$epoch(opt, cost, iter, Y, fn_val)
    if (!is.null(res$opt)) {
      opt <- res$opt
    }
    if (!is.null(res$cost)) {
      cost <- res$cost
    }
  }

  list(
    opt = opt,
    cost = cost
  )
}

# Helper function for epoch callback, allowing user to supply callbacks with
# multiple arities.
do_callback <- function(cb, Y, iter, cost = NULL, cost_fn = NULL, opt = NULL) {
  nfs <- length(formals(cb))
  switch(nfs,
    "1" = cb(Y),
    "2" = cb(Y, iter),
    "3" = cb(Y, iter, cost),
    "4" = cb(Y, iter, cost, cost_fn),
    "5" = cb(Y, iter, cost, cost_fn, opt)
  )
}

# Create a callback for visualization
make_smallvis_cb <- function(df) {
  force(df)
  palette <- NULL
  function(Y, iter, cost = NULL) {
    if (is.null(palette)) {
      palette <- vizier:::make_palette(ncolors = nrow(Y), color_scheme = rainbow)
    }
    title <- paste0("iter: ", iter)
    if (!(is.null(cost) || is.na(cost))) {
      title <- paste0(title, " cost = ", formatC(cost))
    }
    vizier::embed_plot(Y, df, title = title, colors = palette)
  }
}

# Result Export -----------------------------------------------------------

# Prepare the return value.
# If ret_extra is TRUE, return a list with lots of extra info.
# Otherwise, Y is returned directly.
# If ret_extra is TRUE and iter > 0, then all the NULL-default parameters are
# expected to be present. If iter == 0 then the return list will contain only
# scaling and initialization information.
ret_value <- function(Y, ret_extra, method, X, scale, Y_init, iter, start_time = NULL,
                      cost_fn = NULL,
                      G = NULL,
                      perplexity = NULL, pca = 0, whiten = 0,
                      itercosts = NULL,
                      stop_lying_iter = NULL, start_late_lying_iter = NULL,
                      opt_input = NULL, opt_res = NULL,
                      exaggeration_factor = 1, late_exaggeration_factor = 1,
                      optionals = c()) {
  attr(Y, "dimnames") <- NULL
  if (ret_extra) {
    end_time <- Sys.time()

    if (methods::is(X, "dist")) {
      N <- attr(X, "Size")
      origD <- NULL
    }
    else {
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
    }
    else if (whiten > 0) {
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
    }
    else {
      # Don't report start_late_lying_iter if we never got there
      if (start_late_lying_iter > iter) {
        start_late_lying_iter <- NULL
        late_exaggeration_factor <- NULL
      }
    }

    res <- c(res, list(
      perplexity = perplexity,
      itercosts = itercosts,
      exaggeration_factor = exaggeration_factor,
      stop_lying_iter = stop_lying_iter,
      late_exaggeration_factor = late_exaggeration_factor,
      start_late_lying_iter = start_late_lying_iter
    ))

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
        }
        else {
          res[[o]] <- exported
        }
      }
      # For DX and DY, we can calculate if not done so already
      else if (o == "dx") {
        if (methods::is(X, "dist")) {
          res$DX <- X
        }
        else {
          res$DX <- sqrt(safe_dist2(X))
        }
      }
      else if (o == "dy") {
        res$DY <- sqrt(safe_dist2(Y))
      }

      if (o == "x") {
        res$X <- X
      }
    }

    res <- remove_nulls(res)
    res
  }
  else {
    Y
  }
}

# PCA ---------------------------------------------------------------------


# Calculates a matrix containing the first ncol columns of the PCA scores.
# Returns the score matrix unless ret_extra is TRUE, in which case a list
# is returned also containing the eigenvalues
pca_scores <- function(X, ncol = min(dim(X)), verbose = FALSE,
                       ret_extra = FALSE) {
  if (methods::is(X, "dist")) {
    res_mds <- stats::cmdscale(X, x.ret = TRUE, eig = TRUE, k = ncol)

    if (ret_extra || verbose) {
      lambda <- res_mds$eig
      varex <- sum(lambda[1:ncol]) / sum(lambda)
      tsmessage("PCA (using classical MDS): ", ncol, " components explained ",
                formatC(varex * 100), "% variance")
    }
    scores <- res_mds$points
  }
  else {
    X <- scale(X, center = TRUE, scale = FALSE)
    # do SVD on X directly rather than forming covariance matrix
    s <- svd(X, nu = ncol, nv = 0)
    D <- diag(c(s$d[1:ncol]))
    if (verbose || ret_extra) {
      # calculate eigenvalues of covariance matrix from singular values
      lambda <- (s$d ^ 2) / (nrow(X) - 1)
      varex <- sum(lambda[1:ncol]) / sum(lambda)
      tsmessage("PCA: ", ncol, " components explained ", formatC(varex * 100),
                "% variance")
    }
    scores <- s$u %*% D
  }

  if (ret_extra) {
    list(
      scores = scores,
      lambda = lambda[1:ncol]
    )
  }
  else {
    scores
  }
}

# Whiten the data by PCA. This both reduces the dimensionality, but also
# scales the scores by the inverse square root of the equivalent eigenvalue
# so that the variance of each column is 1.
pca_whiten <- function(X, ncol = min(dim(X)), eps = 1e-5, verbose = FALSE) {
  pca <- pca_scores(X, ncol = ncol, verbose = verbose, ret_extra = TRUE)
  sweep(pca$scores, 2, sqrt(pca$lambda + eps), "/")
}


# Perplexity Calibration --------------------------------------------------

# Calculates the input affinities from X, such that each normalized row of the
# affinity matrix has the specified perplexity (within the supplied tolerance).
# Returns a list containing the affinities, beta values and intrinsic
# dimensionalities.
# NB set default kernel to "exp" to get results closer to the tsne package.
# This differs from the procedure in the t-SNE paper by exponentially weighting
# the distances, rather than the squared distances.
x2aff <- function(X, perplexity = 15, tol = 1e-5, kernel = "gauss",
                  verbose = FALSE, guesses = NULL) {
  x_is_dist <- methods::is(X, "dist")
  if (x_is_dist) {
    D <- X
    n <- attr(D, "Size")

    D <- as.matrix(D)
    if (kernel == "gauss") {
      D <- D * D
    }
  }
  else {
    XX <- rowSums(X * X)
    n <- nrow(X)
  }

  nperps <- length(perplexity)
  if (nperps > 1 && nperps != n) {
    stop("Must provide one perplexity per point")
  }

  if (!is.null(guesses) && length(guesses) != n) {
    stop("Initial guess vector must match number of observations in X")
  }

  W <- matrix(0, n, n)
  intd <- rep(0, n)
  if (!is.null(guesses)) {
    beta <- guesses
  }
  else {
    beta <- rep(1, n)
  }
  if (nperps == 1) {
    logU <- log(perplexity)
  }
  else {
    perps <- perplexity
  }
  bad_perp <- 0

  for (i in 1:n) {
    if (nperps > 1) {
      perplexity <- perps[i]
      logU <- log(perplexity)
    }
    betamin <- -Inf
    betamax <- Inf

    if (x_is_dist) {
      Di <- D[i, -i]
    }
    else {
      Di <- (XX[i] + XX - 2 * colSums(tcrossprod(X[i, ], X)))[-i]
      Di[Di < 0] <- 0
      if (kernel == "exp") {
        Di <- sqrt(Di)
      }
    }

    # If we haven't been provided with guesses, then try the initialization used
    # for all points in ELKI according to Schubert & Gertz in "Intrinsic
    # t-Stochastic Neighbor Embedding for Visualization and Outlier Detection: A
    # Remedy Against the Curse of Dimensionality?" Using the last optimized beta
    # seems to be better most of the time based on my testing though, so we'll
    # only use it for the first point.
    if (is.null(guesses) && i == 1) {
      beta[1] <- 0.5 * perplexity / mean(Di)
    }

    sres <- shannon(Di, beta[i])
    H <- sres$H
    Wi <- sres$W
    sumWi <- sres$Z

    Hdiff <- H - logU
    tries <- 0
    while (abs(Hdiff) > tol && tries < 50) {
      if (Hdiff > 0) {
        betamin <- beta[i]
        if (is.infinite(betamax)) {
          beta[i] <- beta[i] * 2
        } else {
          beta[i] <- (beta[i] + betamax) / 2
        }
      } else {
        betamax <- beta[i]
        if (is.infinite(betamin)) {
          beta[i] <- beta[i] / 2
        } else {
          beta[i] <- (beta[i] + betamin) / 2
        }
      }

      sres <- shannon(Di, beta[i])
      H <- sres$H
      Wi <- sres$W
      sumWi <- sres$Z

      Hdiff <- H - logU
      tries <- tries + 1
    }
    if (abs(Hdiff) > tol) {
      bad_perp <- bad_perp + 1
      knn_idx <- order(Di, decreasing = FALSE)[1:max(floor(perplexity), 1)]
      knn_idx[knn_idx >= i] <- knn_idx[knn_idx >= i]  + 1
      Wi <- rep(0, length(Di))
      Wi[knn_idx] <- 1

      intd[i] <- 0
    }
    else {
      # if we didn't supply estimates for beta manually, then initialize guess for
      # next point with optimized beta for this point: doesn't save many
      # iterations, but why not?
      if (is.null(guesses) && i < n) {
        beta[i + 1] <- beta[i]
      }
      intd[i] <- intd_x2aff(Di, beta[i], Wi, sumWi, H)
    }
    W[i, -i] <- Wi
  }
  sigma <- sqrt(1 / beta)

  if (bad_perp > 0) {
    tsmessage("Warning: ", bad_perp, " perplexity calibrations failed!")
    warning(bad_perp, " perplexity calibrations failed")
  }

  if (verbose) {
    summarize(sigma, "sigma summary", verbose = verbose)
    summarize(intd, "Dint", verbose = verbose)
  }
  list(W = W, beta = beta, dint = intd)
}

# Calculates affinites based on exponential weighting of D2 with beta
# and returns a list containing:
# W, the affinities; Z, the sum of the affinities; H, the Shannon entropy
# This routine relies specifically on input weights being = exp(-beta * D)
# and calculates the Shannon entropy as log(Z) + beta * sum(W * D) / Z
# where Z is the sum of W.
shannon <- function(D2, beta) {
  W <- exp(-D2 * beta)
  Z <- sum(W)

  if (Z == 0) {
    H <- 0
  }
  else {
    H <- log(Z) + beta * sum(D2 * W) / Z
  }
  list(
    W = W,
    Z = Z,
    H = H
  )
}

x2aff_sigma <- function(X, sigma = 1e-3, verbose = FALSE) {
  x_is_dist <- methods::is(X, "dist")
  if (x_is_dist) {
    D <- X

    D <- as.matrix(D)
    D <- D * D
  }
  else {
    D <- safe_dist2(X)
  }
  beta <- 1 / (sigma * sigma)
  sres <- shannon(D, beta)
  W <- sres$W
  diag(W) <- 0
  
  list(W = W, beta = beta)
}

# Create a symmetrized distance matrix based on the k-nearest neighbors
# Non-neighbor distances are set to Inf
knn_dist <- function(X, k, n_threads, verbose) {
  if (methods::is(X, "dist")) {
    # If it's already a distance matrix, find k-smallest distance per column
    # (ignoring self-distances of zero) and set everything larger to Inf
    # (potentially more than k finite distances in the event of ties, not
    # going to worry about that)
    D <- as.matrix(X)
    n <- nrow(D)
    if (k > n - 1) {
      stop("k must be not be > n - 1")
    }
    kdists <- Rfast::colnth(D, rep(k + 1, n))
    for (i in 1:n) {
      Di <- D[, i]
      Di[Di > kdists[i]] <- Inf
      D[, i] <- Di
    }
  }
  else {
    # Find the k-nearest indexes and distances of X, and set the corresponding
    # distance matrix elements
    n <- nrow(X)
    if (k > n - 1) {
      stop("k must be not be > n - 1")
    }
    tsmessage("Finding ", k + 1, " nearest neighbors")
    knn <- rnndescent::brute_force_knn(X, k = k + 1, n_threads = n_threads)
    knn$idx <- knn$idx[, 2:(k + 1)]
    knn$dist <- knn$dist[, 2:(k + 1)]
    
    D <- matrix(Inf, nrow = n, ncol = n)
    diag(D) <- 0
    for (i in 1:n) {
      D[i, knn$idx[i, ]] <- knn$dist[i, ]
    }
  }

  # symmetrize
  pmin(D, t(D))
}

# Create the knn graph: D[i, j] = 1 if j is one of i's k-nearest neighbors.
# i is NOT considered a neighbor of itself.
# No symmetrization is carried out.
# Used by knnmmds
knn_graph <- function(X, k, n_threads, verbose) {
  if (methods::is(X, "dist")) {
    D <- as.matrix(X)
    n <- nrow(D)
    if (k > n - 1) {
      stop("k must be not be > n - 1")
    }
    kdists <- Rfast::colnth(D, rep(k + 1, n))
    for (i in 1:n) {
      Di <- D[, i]
      Di[Di > kdists[i]] <- 0
      D[, i] <- 1
    }
  }
  else {
    # Find the k-nearest indexes and distances of X, and set the corresponding
    # distance matrix elements
    n <- nrow(X)
    if (k > n - 1) {
      stop("k must be not be > n - 1")
    }
    
    tsmessage("Finding ", k + 1, " nearest neighbors")
    knn <- rnndescent::brute_force_knn(X, k = k + 1, n_threads = n_threads)
    knn$idx <- knn$idx[, 2:(k + 1)]
    
    D <- matrix(0, nrow = n, ncol = n)
    for (i in 1:n) {
      D[i, knn$idx[i, ]] <- 1
    }
  }
  D
}

# Given data X and k nearest neighbors, return a geodisic distance matrix
# Disconnections are treated by using the Euclidean distance.
geodesic <- function(X, k, fill = TRUE, n_threads = 0, verbose = FALSE) {
  tsmessage("Calculating geodesic distances with k = ", k)

  # The hard work is done by Rfast's implementation of Floyd's algorithm
  R <- knn_dist(X, k, n_threads = n_threads, verbose = verbose)
  G <- Rfast::floyd(R)
  if (any(is.infinite(G)) && fill) {
    tsmessage("k = ", k, " resulted in disconnections: filling with Euclidean distances")
    if (methods::is(X, "dist")) {
      R <- as.matrix(X)
    }
    else {
      R <- sqrt(safe_dist2(X))
    }
    G[is.infinite(G)] <- R[is.infinite(G)]
  }
  G
}

# Multiscale perplexities: P is an average over the results of multiple 
# perplexities
# as described by de Bodt et al in
# Perplexity-free t-SNE and twice Student tt-SNE (2018)
msp <- function(X, perplexities = NULL, tol = 1e-5,
                symmetrize = "symmetric", row_normalize = TRUE,
                normalize = TRUE,
                verbose = FALSE, guesses = NULL) {
  if (methods::is(X, "dist")) {
    n <- attr(X, "Size")
  }
  else {
    n <- nrow(X)
  }

  if (is.null(perplexities)) {
    perplexities <- idp_perps(n)
  }
  tsmessage("Calculating multi-scale P with perplexities from ",
            formatC(perplexities[1]), " to ", formatC(last(perplexities)))

  res <- NULL
  for (perplexity in perplexities) {
    tsmessage("Commencing calibration for perplexity = ",
              format_perps(perplexity))
    x2a_res <- x2aff(X = X, perplexity = perplexity, tol = tol, kernel = "gauss",
                     verbose = verbose, guesses = guesses)
    P <- x2a_res$W
    Q <- scale_affinities(P, 
                          symmetrize = "symmetric", row_normalize = TRUE,
                          normalize = TRUE)
    if (is.null(res)) {
      res$P <- Q
    }
    else {
      res$P <- res$P + Q
    }
  }

  if (length(perplexities) > 1) {
    res$P <- res$P / length(perplexities)
  }
  
  if (is.logical(row_normalize)) {
    tsmessage("Effective perplexity of multiscale P approx = ", 
              formatC(stats::median(perpp(res$P))))
  }

  res
}

# Use the Intrinsic Dimensionality Perplexity (IDP)
# Scan through the provided perplexities and use the result which maximizes
# the mean correlation dimension (which is an estimate for the intrinsic
# dimensionality). Stops at the first maxmimum found.
idp <- function(X, perplexities = NULL, tol = 1e-5,
                verbose = FALSE, guesses = NULL) {
  if (methods::is(X, "dist")) {
    n <- attr(X, "Size")
  }
  else {
    n <- nrow(X)
  }
  
  if (is.null(perplexities)) {
    perplexities <- idp_perps(n)
  }
  if (verbose) {
    tsmessage("Searching for intrinsic dimensionality with perplexities from ",
              formatC(perplexities[1]), " to ", formatC(last(perplexities)))
  }
  
  corr_dim_max <- -Inf
  idp <- 0
  idp_res <- NULL
  for (perplexity in perplexities) {
    if (verbose) {
      tsmessage("Commencing calibration for perplexity = ",
                format_perps(perplexity))
    }
    x2a_res <- x2aff(X = X, perplexity = perplexity, tol = tol, kernel = "gauss",
                     verbose = verbose, guesses = guesses)
    corr_dim <- mean(x2a_res$dint)
    if (corr_dim <= corr_dim_max) {
      break
    }
    else {
      corr_dim_max <- corr_dim
      idp <- perplexity
      idp_res <- x2a_res
    }
  }
  if (idp <= 0) {
    stop("Unable to find an IDP: all correlation dimensions were -ve")
  }
  if (verbose) {
    tsmessage("Found IDP at perplexity = ", formatC(idp),
              " intrinsic dimensionality = ", formatC(corr_dim_max))
  }
  
  idp_res$idp <- idp
  idp_res
}

# Come up with a set of candidate perplexities for finding the Intrinsic
# Dimensionality Perplexity. Use powers of 2 up to around half the data set
# size, or a perplexity of 128, whichever is smaller. Tries to provide a balance
# of coverage of useful perplexities vs time consumption.
idp_perps <- function(n) {
  max_u <- min(128, max(2, ceiling(n / 2)))
  max_uexp <- floor(log2(max_u))
  min_uexp <- min(2, max_uexp)
  2 ^ (min_uexp:max_uexp)
}

# Is the perplexity argument a string or a list with the first element is a
# string?
perp_method <- function(perplexity) {
  method <- ""
  if (is.character(perplexity) || is.list(perplexity)) {
    if (is.list(perplexity)) {
      method <- perplexity[[1]]
    }
    else {
      method <- perplexity
    }
  }
  tolower(method)
}

# If the user provided a list like ("idp", c(10, 20, 30)), extract the list
# of numbers as the candidate list of perplexities. Otherwise, return NULL
user_idp_perps <- function(perplexity) {
  perplexities <- NULL
  if (is.list(perplexity) && length(perplexity) == 2) {
    perplexities <- perplexity[[2]]
  }
  perplexities
}

# Utility Functions -------------------------------------------------------

# Create matrix of squared Euclidean distances
# For low dimension, X %*% t(X) seems to a bit faster than tcrossprod(X)
# Small -ve distances are possible
dist2 <- function(X) {
  D2 <- rowSums(X * X)
  D2 + sweep(X %*% t(X) * -2, 2, t(D2), `+`)
}

# Squared Euclidean distances, ensuring no small -ve distances can occur
safe_dist2 <- function(X) {
  D2 <- dist2(X)
  D2[D2 < 0] <- 0
  D2
}

# 2-norm of a vector or matrix
norm2 <- function(X) {
  sqrt(sum(X * X))
}

# Simple time stamp
stime <- function() {
  format(Sys.time(), "%T")
}

# message with a time stamp
tsmessage <- function(..., domain = NULL, appendLF = TRUE, force = FALSE,
                      time_stamp = TRUE) {
  verbose <- get0("verbose", envir = sys.parent())
  
  if (force || (!is.null(verbose) && verbose)) {
    msg <- ""
    if (time_stamp) {
      msg <- paste0(stime(), " ")
    }
    message(msg, ..., domain = domain, appendLF = appendLF)
    utils::flush.console()
  }
}

# merge lists, where anything non-NULL in l is kept
# e.g.
# all(unlist(update_list(list(a = 1, b = 2), list(a = 10, c = 3))) ==
#     unlist(list(a = 1, b = 2, c = 3)))
lmerge <- function(l, l2) {
  for (name in names(l2)) {
    if (is.null(l[[name]])) {
      l[[name]] <- l2[[name]]
    }
  }
  l
}

# replaces the contents of l with the named arguments
# e.g.
# all(lreplace(c(a = 1, b = 2), a = 10, c = 3) ==
#     c(a = 10, b = 2, c = 3))
lreplace <- function(l, ...) {
  varargs <- list(...)
  for (i in names(varargs)) {
    l[[i]] <- varargs[[i]]
  }
  l
}

# relative tolerance between x and y
reltol <- function(x, y) {
  abs(x - y) / min(abs(x), abs(y))
}

# Check if a value is non-null and true
nnat <- function(x) {
  !is.null(x) && is.logical(x) && x
}

# log vector information
summarize <- function(X, msg = "", verbose = FALSE) {
  summary_X <- summary(X, digits = max(3, getOption("digits") - 3))
  tsmessage(msg, ": ", paste(names(summary_X), ":", summary_X, "|",
                             collapse = ""))
}

# Format perplexity as a string. Could be a scalar or a vector. In the latter
# case, just list the first two values and then ellipses
format_perps <- function(perplexity) {
  if (length(perplexity) > 1) {
    paste0(formatC(perplexity[1]), ", ",
           formatC(perplexity[2]), "...")
  }
  else {
    formatC(perplexity)
  }
}

# last item of a vector
last <- function(x) {
  x[length(x)]
}

# remove NULL items from a list
remove_nulls <- function(l) {
  l[!vapply(l, is.null, logical(1))]
}

# UMAP  -------------------------------------------------------------------

# Fits a kernel for the output distances of the form w = 1 / (1 + a dsq ^ b)
# where dsq is the squared Euclidean distance.
# Standard t-SNE function is a = 1, b = 1.
# Default UMAP values are a = 1.929, b = 0.7915.
find_ab_params <- function(spread = 1, min_dist = 0.001) {
  xv <- seq(from = 0, to = spread * 3, length.out = 300)
  yv <- rep(0, length(xv))
  yv[xv < min_dist] <- 1
  yv[xv >= min_dist] <- exp(-(xv[xv >= min_dist] - min_dist) / spread)
  result <- try({
    stats::nls(yv ~ 1 / (1 + a * xv ^ (2 * b)),
             start = list(a = 1, b = 1))$m$getPars()
    }, silent = TRUE)
    if (class(result) == "try-error") {
      stop("Can't find a, b for provided spread/min_dist values")
    }
  result
}

# The UMAP equivalent of perplexity calibration in x2aff. k is continuous rather
# than integral and so is analogous to perplexity.
# Some differences:
# 1. The target value is the log2 of k, not the Shannon entropy associated
# with the desired perplexity.
# 2. Input weights are exponential, rather than Gaussian, with respect to the
# distances. The distances are also centered with respect to the smoothed
# distance to the nearest (non-zero distance) neighbor. A non-integral
# 'local_connectivity' value can result in this shortest distance between an
# interpolated value between two distances.
# 3. The weights are not normalized. Their raw sum is compared to the target
# value.
# 4. Distances beyond the k-nearest neighbors are not used in the calibration.
# The equivalent weights are set to 0.
# 5. Weights associated with distances shorter than the smoothed nearest
# neighbor distance are clamped to 1.
# This code has been converted from the original Python and may not be very
# idiomatic (or vectorizable).
# tol is SMOOTH_K_TOLERANCE in the Python code.
smooth_knn_distances <-
  function(X,
           k,
           n_iter = 64,
           local_connectivity = 1.0,
           bandwidth = 1.0,
           tol = 1e-5,
           min_k_dist_scale = 1e-3,
           cardinality = log2(k),
           n_threads = 0,
           verbose = FALSE) {

    tsmessage("Commencing smooth kNN distance calibration for k = ", formatC(k))

    if (methods::is(X, "dist")) {
      X <- as.matrix(X)
      nn_idx <- t(apply(X, 2, order))[, 1:k]
      nn_dist <- matrix(0, nrow = nrow(X), ncol = k)
      for (i in 1:nrow(nn_idx)) {
        nn_dist[i, ] <- X[i, nn_idx[i, ]]
      }
    }
    else {
      tsmessage("Finding ", k + 1, " nearest neighbors")
      knn <- rnndescent::brute_force_knn(X, k = k, n_threads = n_threads)
      knn$idx <- knn$idx[, 2:k]
      knn$dist <- knn$dist[, 2:k]

      nn_idx <- matrix(nrow = nrow(X), ncol = k)
      nn_idx[, 1] <- 1:nrow(nn_idx)
      nn_idx[, 2:ncol(nn_idx)] <- knn$idx
      nn_dist <- matrix(0, nrow = nrow(X), ncol = k)
      nn_dist[, 2:ncol(nn_dist)] <- knn$dist
    }

    n <- nrow(nn_dist)
    target <- cardinality * bandwidth
    rho <- rep(0, n)
    sigma <- rep(0, n)
    P <- matrix(0, nrow = n, ncol = n)
    mean_distances <- NULL

    for (i in 1:n) {
      lo <- 0.0
      hi <- Inf
      mid <- 1.0

      ith_distances <- nn_dist[i, ]
      non_zero_dists <- ith_distances[ith_distances > 0.0]
      if (length(non_zero_dists) >= local_connectivity) {
        index <- floor(local_connectivity)
        interpolation <- local_connectivity - index
        if (index > 0) {
          if (interpolation <= tol) {
            rho[i] <- non_zero_dists[index]
          }
          else {
            rho[i] <- non_zero_dists[index] + interpolation *
              (non_zero_dists[index + 1] - non_zero_dists[index])
          }
        }
        else {
          rho[i] <- interpolation * non_zero_dists[1]
        }
      } else if (length(non_zero_dists) > 0) {
        rho[i] <- max(non_zero_dists)
      }
      else {
        rho[i] <- 0.0
      }

      for (iter in 1:n_iter) {
        psum <- 0.0
        for (j in 2:ncol(nn_dist)) {
          dist <- max(0, (nn_dist[i, j] - rho[i]))
          psum <- psum + exp(-(dist / mid))
        }
        val <- psum

        if (abs(val - target) < tol) {
          break
        }

        if (val > target) {
          hi <- mid
          mid <- (lo + hi) / 2.0
        }
        else {
          lo <- mid
          if (is.infinite(hi)) {
            mid <- mid * 2
          }
          else {
            mid <- (lo + hi) / 2.0
          }
        }
      }
      sigma[i] <- mid

      if (rho[i] > 0.0) {
        sigma[i] <- max(sigma[i], min_k_dist_scale * mean(ith_distances))
      }
      else {
        if (is.null(mean_distances)) {
          mean_distances <- mean(nn_dist)
        }
        sigma[i] <- max(sigma[i], min_k_dist_scale * mean_distances)
      }

      prow <- exp(-(nn_dist[i, ] - rho[i]) / (sigma[i] * bandwidth))
      prow[nn_dist[i, ] - rho[i] <= 0] <- 1
      P[i, nn_idx[i, ]] <- prow
    }
    diag(P) <- 0

    if (verbose) {
      summarize(sigma, "sigma summary", verbose = verbose)
    }
    list(sigma = sigma, rho = rho, P = P)
  }


# set_op_mix_ratio = between 0 and 1 mixes in fuzzy set intersection
# set to 0 for intersection only
fuzzy_set_union <- function(X, set_op_mix_ratio = 1) {
  XX <- X * t(X)
  set_op_mix_ratio * (X + t(X) - XX) + (1 - set_op_mix_ratio) * XX
}

init_ab <- function(cost, spread = 1, min_dist = 0.001, verbose = FALSE) {
  ab_params <- find_ab_params(spread = spread, min_dist = min_dist)
  a <- ab_params[1]
  b <- ab_params[2]
  if (verbose) {
    message("Umap curve parameters = ", formatC(a), ", ", formatC(b))
  }
  cost$a <- a
  cost$b <- b
  cost
}
