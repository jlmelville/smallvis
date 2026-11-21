#include "vendor/zufft.hpp"
#include <Rcpp.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <vector>

using namespace Rcpp;

namespace {
// Keep only the point data read by the existing spread and gather loops.
struct DirectGrid {
  enum class Status { ready, coincident, too_large };
  std::vector<int> starts;
  std::vector<double> weights;
  std::array<std::vector<double>, 2> xy;
  std::array<int, 2> side{};
  double width = 0;

  static int cell(double position, int boxes) {
    return std::min(static_cast<int>(position), boxes - 1);
  }

  Status prepare(const NumericMatrix &points, double interval,
                 const NumericVector &fast_boxes, const zufft_api &api,
                 double cell_limit) {
    const int n = points.nrow();
    if (points.ncol() != 2 || n < 1 ||
        n > std::numeric_limits<int>::max() / 9 || !std::isfinite(interval) ||
        interval <= 0)
      Rcpp::stop("Invalid FFT grid input");
    double minimum[2], extent[2];
    for (int axis = 0; axis < 2; ++axis) {
      const double *column = points.begin() + axis * n;
      double low = column[0], high = column[0];
      for (int i = 0; i < n; ++i) {
        if (!std::isfinite(column[i]))
          Rcpp::stop("FFT coordinates must be finite");
        low = std::min(low, column[i]);
        high = std::max(high, column[i]);
      }
      minimum[axis] = low;
      extent[axis] = high - low;
    }
    const double span = std::max(extent[0], extent[1]);
    if (span == 0)
      return Status::coincident;
    const double requested =
        std::max(fast_boxes[0], std::floor(span / interval));
    if (requested > cell_limit)
      return Status::too_large;
    int boxes = static_cast<int>(requested);
    const auto fast =
        std::lower_bound(fast_boxes.begin(), fast_boxes.end(), boxes);
    if (fast != fast_boxes.end())
      boxes = static_cast<int>(*fast);
    width = span / boxes;
    std::array<int, 2> nodes{};
    // Each axis starts at its own minimum in box zero. Its maximum therefore
    // determines the occupied extent without constructing a point-to-cell
    // vector.
    for (int axis = 0; axis < 2; ++axis) {
      nodes[axis] = 3 * (1 + cell((extent[axis] / span) * boxes, boxes));
    }
    if (4.0 * nodes[0] * nodes[1] > cell_limit)
      return Status::too_large;
    for (int axis = 0; axis < 2; ++axis) {
      side[axis] = static_cast<int>(api.good_size(2 * nodes[axis], 1));
    }
    if (side[0] <= 0 || side[1] <= 0 ||
        static_cast<double>(side[0]) * side[1] > cell_limit)
      return Status::too_large;
    starts.resize(n);
    weights.resize(9 * n);
    for (auto &column : xy)
      column.resize(n);
    for (int i = 0; i < n; ++i) {
      double basis[2][3];
      int start[2];
      for (int axis = 0; axis < 2; ++axis) {
        const double local = points(i, axis) - minimum[axis];
        const double position = (local / span) * boxes;
        const int box = cell(position, boxes);
        const double u = 3 * (position - box) - 1.5;
        basis[axis][0] = u * (u - 1) / 2;
        basis[axis][1] = 1 - u * u;
        basis[axis][2] = u * (u + 1) / 2;
        start[axis] = 3 * box;
        xy[axis][i] = local - extent[axis] / 2;
      }
      starts[i] = start[0] + side[0] * start[1];
      for (int b = 0; b < 3; ++b)
        for (int a = 0; a < 3; ++a)
          weights[9 * i + 3 * b + a] = basis[0][a] * basis[1][b];
    }
    return Status::ready;
  }
};

void fill_kernel(double *out, int nx, int ny, double width, bool squared) {
  const double step = width / 3;
  std::vector<double> dx(nx);
  for (int x = 0; x < nx; ++x) {
    double v = (2 * x >= nx ? x - nx : x) * step;
    dx[x] = v * v;
  }
  for (int y = 0; y < ny; ++y) {
    double v = (2 * y >= ny ? y - ny : y) * step;
    double dy = v * v;
    for (int x = 0; x < nx; ++x) {
      double w = 1 / (1 + (dx[x] + dy));
      out[y * nx + x] = squared ? w * w : w;
    }
  }
}

void spread(double *out, int nx, int ny, const DirectGrid &grid, int field) {
  std::fill_n(out, nx * ny, 0.0);
  for (size_t i = 0; i < grid.starts.size(); ++i) {
    const double q = field == 0 ? 1 : grid.xy[field - 1][i];
    for (int b = 0; b < 3; ++b) {
      const int offset = grid.starts[i] + b * nx;
      for (int a = 0; a < 3; ++a)
        out[offset + a] += grid.weights[9 * i + 3 * b + a] * q;
    }
  }
}

void gather(const double *in, double *out, int nx, int ny,
            const DirectGrid &grid) {
  for (size_t i = 0; i < grid.starts.size(); ++i) {
    double value = 0;
    for (int b = 0; b < 3; ++b) {
      const int offset = grid.starts[i] + b * nx;
      for (int a = 0; a < 3; ++a)
        value += grid.weights[9 * i + 3 * b + a] * in[offset + a];
    }
    out[i] = value / (nx * ny);
  }
}

struct Workspace {
  DirectGrid grid;
  std::vector<double> real, kernel, source, potential;
  // Keep the largest grid seen during this FFT phase. A new shape always gets
  // fresh layouts, weights and spectra, including when it fits the existing
  // storage.
  static void ensure(std::vector<double> &buffer, size_t size) {
    if (buffer.size() < size)
      std::vector<double>(size).swap(buffer);
  }
};
} // namespace

// [[Rcpp::export(rng = false)]]
SEXP zufft_workspace_cpp() { return XPtr<Workspace>(new Workspace(), true); }

// [[Rcpp::export(rng = false)]]
SEXP zufft_tsne_sums_cpp(SEXP workspace, const NumericMatrix &points,
                         double interval, const NumericVector &fast_boxes,
                         int workers, bool mass_only, bool reuse,
                         bool hybrid = false) {
  XPtr<Workspace> ws(workspace);
  // Resolve the optional provider on the R thread, after requireNamespace().
  const zufft_api *api = zufft_api_get();
  if (!api)
    stop("Incompatible zufft C interface (expected version 1)");
  if (workers < 1 || fast_boxes.size() == 0)
    stop("Invalid FFT settings");
  auto &grid = ws->grid;
  const int n = points.nrow();
  constexpr double max_cells = 4000000;
  const double limit = hybrid ? std::min(max_cells, 14.0 * n) : max_cells;
  const auto status = grid.prepare(points, interval, fast_boxes, *api, limit);
  if (status == DirectGrid::Status::too_large) {
    if (!hybrid)
      stop("FFT padded grid exceeds four million cells. Increase grid_interval "
           "or use Barnes-Hut.");
    // Free the FFT buffers before the caller allocates the mutual tree.
    if (!ws->real.empty())
      *ws = Workspace{};
    return R_NilValue;
  }
  if (status == DirectGrid::Status::coincident) {
    List result = List::create(_["Z"] = static_cast<double>(n) * (n - 1));
    if (!mass_only)
      result["force"] = NumericMatrix(n, 2);
    return result;
  }

  const int nx = grid.side[0], ny = grid.side[1];
  const size_t shape[] = {static_cast<size_t>(nx), static_cast<size_t>(ny)};
  const size_t half[] = {shape[0] / 2 + 1, shape[1]};
  const size_t spectrum_size = 2 * half[0] * half[1];
  const zufft_layout real_layout = {2, shape, nullptr},
                     packed_layout = {2, half, nullptr};
  // Pack the contiguous x axis. All transforms are unnormalized.
  const int axes[] = {1, 0};
  Workspace::ensure(ws->real, shape[0] * shape[1]);
  Workspace::ensure(ws->kernel, spectrum_size);
  Workspace::ensure(ws->source, spectrum_size);
  Workspace::ensure(ws->potential, n);
  const bool save_unit = reuse && !mass_only;
  auto check = [&](int status) {
    if (status != ZUFFT_OK)
      stop("zufft: %s", api->message());
  };
  auto forward = [&](double *out) {
    check(api->r2c(&real_layout, ws->real.data(), &packed_layout, out, 2, axes,
                   1, 1, workers));
  };

  // First gather the Cauchy mass. The squared kernel then acts on unit,
  // x and y charges for the two components of the unnormalized force.
  List result;
  NumericMatrix force = no_init_matrix(mass_only ? 0 : n, 2);
  for (int col = 0; col < (mass_only ? 1 : 4); ++col) {
    if (col < 2) {
      fill_kernel(ws->real.data(), nx, ny, grid.width, col == 1);
      forward(ws->kernel.data());
    }
    if (!(save_unit && col == 1)) {
      spread(ws->real.data(), nx, ny, grid, std::max(0, col - 1));
      forward(ws->source.data());
    }
    // The mass kernel is used once. Overwrite it with the first product
    // so the unit spectrum survives in source for the squared kernel.
    double *product =
        save_unit && col == 0 ? ws->kernel.data() : ws->source.data();
    const double *source_data = ws->source.data();
    const double *kernel_data = ws->kernel.data();
    for (size_t i = 0; i < spectrum_size; i += 2) {
      const double a = source_data[i], b = source_data[i + 1];
      const double c = kernel_data[i], d = kernel_data[i + 1];
      product[i] = a * c - b * d;
      product[i + 1] = a * d + b * c;
    }
    // The product spectrum can be overwritten. Invert y in place, then
    // unpack x into the real grid for gathering.
    check(api->c2c(&packed_layout, product, &packed_layout, product, 1, axes, 0,
                   1, workers));
    check(api->c2r(&packed_layout, product, &real_layout, ws->real.data(), 1,
                   axes + 1, 0, 1, workers));
    gather(ws->real.data(), ws->potential.data(), nx, ny, grid);
    if (col == 0) {
      long double mass = 0;
      for (int i = 0; i < n; ++i)
        mass += ws->potential[i];
      const double Z = static_cast<double>(mass) - n;
      if (!(Z > 0) || !std::isfinite(Z))
        stop("FFT could not estimate a positive normalization. Reduce "
             "grid_interval.");
      result["Z"] = Z;
    } else if (col == 1) {
      for (int axis = 0; axis < 2; ++axis)
        for (int i = 0; i < n; ++i)
          force(i, axis) = grid.xy[axis][i] * ws->potential[i];
    } else {
      for (int i = 0; i < n; ++i)
        force(i, col - 2) -= ws->potential[i];
    }
  }
  if (!mass_only)
    result["force"] = force;
  return result;
}
