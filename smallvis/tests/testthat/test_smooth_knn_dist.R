test_that("smooth distances support eight neighbors", {
  res <- smooth_knn_distances(iris10, k = 8)
  expect_equal(
    res$sigma,
    c(
      0.2567215,
      0.22098923,
      0.08285332,
      0.09981823,
      0.28608322,
      0.17873764,
      0.15968704,
      0.17134094,
      0.25434113,
      0.19572449
    )
  )
  expect_equal(
    res$rho,
    c(
      0.14142136,
      0.17320508,
      0.24494897,
      0.24494897,
      0.14142136,
      0.6164414,
      0.26457513,
      0.17320508,
      0.3,
      0.17320508
    ),
    tolerance = 1e-6
  )
})

test_that("smooth distances support four neighbors", {
  res <- smooth_knn_distances(iris10, k = 4)
  expect_equal(
    res$sigma,
    c(
      0.17993927,
      0.20488739,
      0.0493803,
      0.09060478,
      0.24940491,
      0.00390625,
      0.15367126,
      0.13551712,
      0.24542618,
      0.20633698
    ),
    tolerance = 1e-5
  )
  expect_equal(
    res$rho,
    c(
      0.14142136,
      0.17320508,
      0.24494897,
      0.24494897,
      0.14142136,
      0.6164414,
      0.26457513,
      0.17320508,
      0.3,
      0.17320508
    ),
    tolerance = 1e-6
  )
})

test_that("distance input and fuzzy set operations match UMAP references", {
  # distance matrix
  res <- smooth_knn_distances(dist(iris10), k = 4)
  expect_equal(
    res$sigma,
    c(
      0.17993927,
      0.20488739,
      0.0493803,
      0.09060478,
      0.24940491,
      0.00390625,
      0.15367126,
      0.13551712,
      0.24542618,
      0.20633698
    )
  )
  expect_equal(
    res$rho,
    c(
      0.14142136,
      0.17320508,
      0.24494897,
      0.24494897,
      0.14142136,
      0.6164414,
      0.26457513,
      0.17320508,
      0.3,
      0.17320508
    ),
    tolerance = 1e-6
  )

  # numbers come from inserting a print statement into the appropriate bit of UMAP
  # source code
  P <- matrix(0, nrow = 10, ncol = 10)
  P[1, 5] <- 1
  P[1, 8] <- 0.838084924053286
  P[1, 10] <- 0.1619080617920443

  P[2, 3] <- 0.5385624888941846
  P[2, 4] <- 0.46144714836635475
  P[2, 10] <- 1.0

  P[3, 2] <- 0.32796848436705417
  P[3, 4] <- 1.0
  P[3, 7] <- 0.6720321472619579

  P[4, 3] <- 1.0
  P[4, 9] <- 0.5446591215696793
  P[4, 10] <- 0.4553449084802392

  P[5, 1] <- 1.0
  P[5, 8] <- 0.7192644685443685
  P[5, 7] <- 0.28072806955243307

  P[6, 1] <- 1.0
  P[6, 5] <- 1.0
  P[6, 8] <- 5.128685589309561e-10

  P[7, 3] <- 1.0
  P[7, 4] <- 0.6462531113433841
  P[7, 8] <- 0.35375192447644943

  P[8, 1] <- 1.0
  P[8, 5] <- 0.6894084285873151
  P[8, 10] <- 0.3105906239414311

  P[9, 4] <- 1.0
  P[9, 3] <- 0.5748250840738938
  P[9, 2] <- 0.4251747827928667

  P[10, 2] <- 1.0
  P[10, 3] <- 0.4999980790544043
  P[10, 4] <- 0.49999807905434984

  expect_equal(res$P, P)

  # numbers from python fuzzy_simplicial_set
  Ps <- matrix(0, nrow = 10, ncol = 10)
  Ps[1, 5] <- 1.0
  Ps[1, 6] <- 1.0
  Ps[1, 8] <- 1.0
  Ps[1, 10] <- 0.1619080617920443

  Ps[2, 3] <- 0.6898994500416646
  Ps[2, 4] <- 0.46144714836635475
  Ps[2, 9] <- 0.4251747827928667
  Ps[2, 10] <- 1.0

  Ps[3, 2] <- 0.6898994500416646
  Ps[3, 4] <- 1.0
  Ps[3, 7] <- 1.0
  Ps[3, 9] <- 0.5748250840738938
  Ps[3, 10] <- 0.4999980790544043

  Ps[4, 2] <- 0.46144714836635475
  Ps[4, 3] <- 1.0
  Ps[4, 7] <- 0.6462531113433841
  Ps[4, 9] <- 1.0
  Ps[4, 10] <- 0.7276714079872907

  Ps[5, 1] <- 1.0
  Ps[5, 6] <- 1.0
  Ps[5, 7] <- 0.28072806955243307
  Ps[5, 8] <- 0.9128059101338203

  Ps[6, 1] <- 1.0
  Ps[6, 5] <- 1.0
  Ps[6, 8] <- 5.128685589309561e-10

  Ps[7, 3] <- 1.0
  Ps[7, 4] <- 0.6462531113433841
  Ps[7, 5] <- 0.28072806955243307
  Ps[7, 8] <- 0.35375192447644943

  Ps[8, 1] <- 1.0
  Ps[8, 5] <- 0.9128059101338203
  Ps[8, 6] <- 5.128685589309561e-10
  Ps[8, 7] <- 0.35375192447644943
  Ps[8, 10] <- 0.3105906239414311

  Ps[9, 2] <- 0.4251747827928667
  Ps[9, 3] <- 0.5748250840738938
  Ps[9, 4] <- 1.0

  Ps[10, 1] <- 0.1619080617920443
  Ps[10, 2] <- 1.0
  Ps[10, 3] <- 0.4999980790544043
  Ps[10, 4] <- 0.7276714079872907
  Ps[10, 8] <- 0.3105906239414311

  expect_equal(fuzzy_set_union(res$P), Ps)

  # mix intersection with union
  Psm <- matrix(0, nrow = 10, ncol = 10)
  Psm[1, 5] <- 1.0
  Psm[1, 6] <- 0.5
  Psm[1, 8] <- 0.9190424620266431
  Psm[1, 10] <- 0.08095403089602216

  Psm[2, 3] <- 0.43326548663061937
  Psm[2, 4] <- 0.23072357418317738
  Psm[2, 9] <- 0.21258739139643335
  Psm[2, 10] <- 1.0

  Psm[3, 2] <- 0.43326548663061937
  Psm[3, 4] <- 1.0
  Psm[3, 7] <- 0.836016073630979
  Psm[3, 9] <- 0.2874125420369469
  Psm[3, 10] <- 0.24999903952720215

  Psm[4, 2] <- 0.23072357418317738
  Psm[4, 3] <- 1.0
  Psm[4, 7] <- 0.3231265556716921
  Psm[4, 9] <- 0.7723295607848397
  Psm[4, 10] <- 0.4776714937672945

  Psm[5, 1] <- 1.0
  Psm[5, 6] <- 0.5
  Psm[5, 7] <- 0.14036403477621653
  Psm[5, 8] <- 0.7043364485658419

  Psm[6, 1] <- 0.5
  Psm[6, 5] <- 0.5
  Psm[6, 8] <- 2.5643427946547806e-10

  Psm[7, 3] <- 0.836016073630979
  Psm[7, 4] <- 0.3231265556716921
  Psm[7, 5] <- 0.14036403477621653
  Psm[7, 8] <- 0.17687596223822472

  Psm[8, 1] <- 0.9190424620266431
  Psm[8, 5] <- 0.7043364485658419
  Psm[8, 6] <- 2.5643427946547806e-10
  Psm[8, 7] <- 0.17687596223822472
  Psm[8, 10] <- 0.15529531197071555

  Psm[9, 2] <- 0.21258739139643335
  Psm[9, 3] <- 0.2874125420369469
  Psm[9, 4] <- 0.7723295607848397

  Psm[10, 1] <- 0.08095403089602216
  Psm[10, 2] <- 1.0
  Psm[10, 3] <- 0.24999903952720215
  Psm[10, 4] <- 0.4776714937672945
  Psm[10, 8] <- 0.15529531197071555

  expect_equal(fuzzy_set_union(res$P, set_op_mix_ratio = 0.5), Psm)

  Psi <- matrix(0, nrow = 10, ncol = 10)
  Psi[1, 5] <- 1.0
  Psi[1, 8] <- 0.838084924053286

  Psi[2, 3] <- 0.17663152321957415
  Psi[2, 10] <- 1.0

  Psi[3, 2] <- 0.17663152321957415
  Psi[3, 4] <- 1.0
  Psi[3, 7] <- 0.6720321472619579

  Psi[4, 3] <- 1.0
  Psi[4, 9] <- 0.5446591215696793
  Psi[4, 10] <- 0.22767157954729833

  Psi[5, 1] <- 1.0
  Psi[5, 8] <- 0.49586698699786347

  Psi[7, 3] <- 0.6720321472619579

  Psi[8, 1] <- 0.838084924053286
  Psi[8, 5] <- 0.49586698699786347

  Psi[9, 4] <- 0.5446591215696793

  Psi[10, 2] <- 1.0
  Psi[10, 4] <- 0.22767157954729833

  expect_equal(fuzzy_set_union(res$P, set_op_mix_ratio = 0), Psi)
})

test_that("fractional local connectivity matches the UMAP reference", {
  res <- smooth_knn_distances(iris10, k = 4, local_connectivity = 1.5)
  expect_equal(
    res$sigma,
    c(
      0.13976288,
      0.11270523,
      0.03319359,
      0.05047989,
      0.17582703,
      0.00390625,
      0.1002903,
      0.09278679,
      0.14465714,
      0.10316849
    ),
    tolerance = 1e-5
  )
  expect_equal(
    res$rho,
    c(
      0.15731322,
      0.23660254,
      0.25476205,
      0.27247449,
      0.18251408,
      0.6164414,
      0.29811881,
      0.19840594,
      0.36794495,
      0.24471642
    ),
    tolerance = 1e-6
  )
  Pl <- matrix(0, nrow = 10, ncol = 10)

  Pl[1, 5] <- 1.0
  Pl[1, 8] <- 0.8925204178669899
  Pl[1, 10] <- 0.10748447876830035

  Pl[2, 3] <- 0.569778945418246
  Pl[2, 4] <- 0.4302286057469135
  Pl[2, 10] <- 1.0

  Pl[3, 2] <- 0.255929829694113
  Pl[3, 4] <- 1.0
  Pl[3, 7] <- 0.7440613820135572

  Pl[4, 3] <- 1.0
  Pl[4, 9] <- 0.5796812988186044
  Pl[4, 10] <- 0.42031672929651537

  Pl[5, 1] <- 1.0
  Pl[5, 8] <- 0.7915905048023626
  Pl[5, 7] <- 0.20840631064514223

  Pl[6, 1] <- 1.0
  Pl[6, 5] <- 1.0
  Pl[6, 8] <- 5.128685589309561e-10

  Pl[7, 3] <- 1.0
  Pl[7, 4] <- 0.7157203287723206
  Pl[7, 8] <- 0.28427839448271264

  Pl[8, 1] <- 1.0
  Pl[8, 5] <- 0.762159337934477
  Pl[8, 10] <- 0.2378396564365533

  Pl[9, 2] <- 0.3748106659309435
  Pl[9, 3] <- 0.6251919975207915
  Pl[9, 4] <- 1.0

  Pl[10, 2] <- 1.0
  Pl[10, 3] <- 0.4999980790544043
  Pl[10, 4] <- 0.49999807905429533

  expect_equal(res$P, Pl, tolerance = 1e-5)
})

test_that("duplicated inputs retain k minus one non-self UMAP neighbors", {
  X <- rbind(matrix(0, 6, 2), cbind(1:10, (1:10)^2))
  # Inspect the directed graph before fuzzy symmetrization can increase degree.
  for (input in list(X, stats::dist(X))) {
    set.seed(42)
    result <- smooth_knn_distances(input, k = 4, n_threads = 1)
    expect_equal(diag(result$P), rep(0, 16))
    expect_true(all(is.finite(result$P)))
    expect_equal(rowSums(result$P[1:6, ] > 0), rep(3, 6))
    expect_true(all(rowSums(result$P > 0) <= 3))
  }
})

test_that("neighbor extraction removes self by identity before truncation", {
  # A controlled tied pool tests membership independently of backend tie choices.
  # These private graph builders own the directed support, hidden by public
  # symmetrization or shortest-path completion.
  candidates <- rbind(
    c(2L, 1L, 3L),
    c(3L, 4L, 2L),
    c(1L, 2L, 4L),
    c(4L, 1L, 2L)
  )
  local_mocked_bindings(get_nn = function(X, k, ...) {
    expect_equal(k, 3)
    list(idx = candidates, dist = matrix(0, 4, 3))
  })
  expected <- matrix(0, 4, 4)
  expected[1, c(2, 3)] <- 1
  expected[2, c(3, 4)] <- 1
  expected[3, c(1, 2)] <- 1
  expected[4, c(1, 2)] <- 1
  X <- matrix(0, 4, 2)
  expect_equal(smooth_knn_distances(X, k = 3)$P, expected)
  graph <- matrix(Inf, 4, 4)
  diag(graph) <- 0
  graph[expected > 0] <- 0
  expect_equal(
    knn_dist(X, 2, n_threads = 1, verbose = FALSE),
    pmin(graph, t(graph))
  )
})

test_that("fractional UMAP k retains its continuous cardinality target", {
  X <- cbind(sin(1:12), cos(1:12))
  for (input in list(X, stats::dist(X))) {
    result <- smooth_knn_distances(input, k = 4.5, n_threads = 1)
    expect_equal(diag(result$P), rep(0, 12))
    expect_true(all(rowSums(result$P > 0) <= 3))
    expect_equal(rowSums(result$P), rep(log2(4.5), 12), tolerance = 1e-4)
  }
})

test_that("public UMAP neighborhoods enforce the reserved-self count contract", {
  X <- cbind(sin(1:12), cos(1:12))
  for (input in list(X, stats::dist(X))) {
    run <- function(perplexity) {
      smallvis(
        input,
        method = "umap",
        perplexity = perplexity,
        scale = FALSE,
        Y_init = cbind(sin(1:12), cos(1:12)),
        max_iter = 0,
        ret_extra = "P",
        epoch_callback = FALSE,
        verbose = FALSE
      )
    }
    for (perplexity in c(1, 13)) {
      expect_error(run(perplexity), "UMAP k.*between 2")
    }
    for (perplexity in c(2, 4.5, 12)) {
      result <- run(perplexity)
      expect_equal(diag(result$P), rep(0, 12))
      expect_true(all(is.finite(result$P)))
      expect_true(all(rowSums(result$P) > 0))
      expect_true(all(is.finite(result$costs)))
    }
  }
})
