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

  list(opt = opt, cost = cost)
}

# Helper function for epoch callback, allowing user to supply callbacks with
# multiple arities.
do_callback <- function(cb,
                        Y,
                        iter,
                        cost = NULL,
                        cost_fn = NULL,
                        opt = NULL) {
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
      palette <- vizier:::color_helper(df, color_scheme = grDevices::rainbow)$palette
    }
    title <- paste0("iter: ", iter)
    if (!(is.null(cost) || is.na(cost))) {
      title <- paste0(title, " cost = ", formatC(cost))
    }
    vizier::embed_plot(Y, df, title = title, color_scheme = palette)
  }
}
