#include "threads.h"
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <memory>
#include <utility>
#include <vector>

using Rcpp::_;
using Rcpp::List;
using Rcpp::NumericMatrix;

namespace {
struct Point {
  double x, y;
  int row;
};

// Local affine vector field around the cell's centroid.
struct Field {
  double x = 0, y = 0, xx = 0, xy = 0, yy = 0;

  void inherit(const Field &p, double dx, double dy) {
    x += p.x + p.xx * dx + p.xy * dy;
    y += p.y + p.xy * dx + p.yy * dy;
    xx += p.xx;
    xy += p.xy;
    yy += p.yy;
  }
};

struct Cell {
  int lo, hi, left = -1, right = -1;
  double x = 0, y = 0, radius = 0, xx = 0, xy = 0, yy = 0;

  int size() const { return hi - lo; }
};

struct Force {
  double x = 0, y = 0;
};

struct Sums {
  Sums(size_t n, size_t cells) : forces(n), fields(cells) {}
  std::vector<Force> forces;
  std::vector<Field> fields;
  double mass = 0;
  void merge(const Sums &other) {
    mass += other.mass;
    for (size_t i = 0; i < forces.size(); ++i) {
      forces[i].x += other.forces[i].x;
      forces[i].y += other.forces[i].y;
    }
    for (size_t i = 0; i < fields.size(); ++i) {
      Field &f = fields[i];
      const Field &g = other.fields[i];
      f.x += g.x;
      f.y += g.y;
      f.xx += g.xx;
      f.xy += g.xy;
      f.yy += g.yy;
    }
  }
};

using Pair = std::pair<int, int>;

class MutualTree {
public:
  explicit MutualTree(const NumericMatrix &Y) {
    points.reserve(Y.nrow());
    for (int i = 0; i < Y.nrow(); ++i) {
      const double x = Y(i, 0), y = Y(i, 1);
      if (!std::isfinite(x) || !std::isfinite(y))
        Rcpp::stop("Mutual tree coordinates must be finite");
      points.push_back({x, y, i});
    }
    build(0, points.size());
  }

  void refit(const NumericMatrix &Y) {
    if (static_cast<size_t>(Y.nrow()) != points.size())
      Rcpp::stop("Mutual tree refit requires the same number of points");
    for (Point &p : points) {
      p.x = Y(p.row, 0);
      p.y = Y(p.row, 1);
      if (!std::isfinite(p.x) || !std::isfinite(p.y))
        Rcpp::stop("Mutual tree coordinates must be finite");
    }
    for (Cell &cell : cells)
      update_geometry(cell);
  }

  template <bool Gradient> List evaluate(double theta, int threads) {
    theta2 = theta * theta;
    const size_t n = Gradient ? points.size() : 0;
    const size_t nc = Gradient ? cells.size() : 0;
    Sums total(n, nc);
    const size_t capacity =
        std::min<size_t>((points.size() + 255) / 256,
                         std::max(1u, std::thread::hardware_concurrency()));
    const size_t workers = smallvis::worker_count(capacity, threads);
    if (workers == 1) {
      visit<Gradient>(0, 0, total);
    } else {
      // Stop the ordinary walk at a shallow frontier. Accepted ancestors
      // and direct leaf pairs above it have already contributed exactly once.
      std::vector<Pair> tasks;
      int depth = 0;
      for (size_t width = 1; width < 32 * workers; width *= 2)
        ++depth;
      visit<Gradient>(0, 0, total, &tasks, depth);
      const size_t count = std::min(workers, tasks.size());
      std::vector<Sums> partial;
      for (size_t i = 0; i < count; ++i)
        partial.emplace_back(n, nc);
      smallvis::parallel_for(
          count, count, [&](size_t start, size_t end, size_t) {
            for (size_t w = start; w < end; ++w)
              for (size_t t = w; t < tasks.size(); t += count)
                visit<Gradient>(tasks[t].first, tasks[t].second, partial[w]);
          });
      for (const auto &part : partial)
        total.merge(part);
    }
    if (!Gradient)
      return List::create(_["Z"] = 2 * total.mass);
    NumericMatrix force(points.size(), 2);
    distribute(0, total, force);
    return List::create(_["Z"] = 2 * total.mass, _["force"] = force);
  }

private:
  std::vector<Point> points;
  std::vector<Cell> cells;
  double theta2 = 0;

  void update_geometry(Cell &cell) const {
    cell.x = cell.y = cell.xx = cell.xy = cell.yy = cell.radius = 0;
    for (int i = cell.lo; i < cell.hi; ++i) {
      cell.x += points[i].x;
      cell.y += points[i].y;
    }
    cell.x /= cell.size();
    cell.y /= cell.size();
    for (int i = cell.lo; i < cell.hi; ++i) {
      double x = points[i].x - cell.x, y = points[i].y - cell.y;
      cell.radius = std::max(cell.radius, x * x + y * y);
      cell.xx += x * x;
      cell.xy += x * y;
      cell.yy += y * y;
    }
    cell.radius = std::sqrt(cell.radius);
    cell.xx /= cell.size();
    cell.xy /= cell.size();
    cell.yy /= cell.size();
  }

  int build(int lo, int hi) {
    const int id = cells.size();
    Cell cell{lo, hi};
    update_geometry(cell);
    cells.push_back(cell);
    if (cell.size() > 12) {
      int mid = lo + (hi - lo) / 2;
      bool along_x = cell.xx >= cell.yy;
      std::nth_element(points.begin() + lo, points.begin() + mid,
                       points.begin() + hi,
                       [along_x](const Point &a, const Point &b) {
                         return along_x ? a.x < b.x : a.y < b.y;
                       });
      int left = build(lo, mid), right = build(mid, hi);
      cells[id].left = left;
      cells[id].right = right;
    }
    return id;
  }

  // Second-order cell-average force and mass, first-order variation within
  // each target. Source and target offsets have zero mean by construction.
  template <bool Gradient>
  void summarize(int ia, int ib, double x, double y, double w, Sums &s) const {
    const Cell &a = cells[ia], &b = cells[ib];
    double xx = a.xx + b.xx, xy = a.xy + b.xy, yy = a.yy + b.yy;
    double cx = xx * x + xy * y, cy = xy * x + yy * y, trace = xx + yy;
    double radial = x * cx + y * cy, w2 = w * w, w3 = w2 * w;
    s.mass += static_cast<double>(a.size()) * b.size() *
              (w - w2 * trace + 4 * w3 * radial);
    if (!Gradient)
      return;
    double factor = w2 - 2 * w3 * trace + 12 * w3 * w * radial;
    double fx = factor * x - 4 * w3 * cx, fy = factor * y - 4 * w3 * cy;
    auto add = [&](Field &f, double mass, double sign) {
      f.x += sign * mass * fx;
      f.y += sign * mass * fy;
      f.xx += mass * (w2 - 4 * w3 * x * x);
      f.xy -= mass * 4 * w3 * x * y;
      f.yy += mass * (w2 - 4 * w3 * y * y);
    };
    add(s.fields[ia], b.size(), 1);
    add(s.fields[ib], a.size(), -1);
  }

  template <bool Gradient>
  void visit(int ia, int ib, Sums &s, std::vector<Pair> *frontier = nullptr,
             int depth = 0) const {
    if (frontier && depth == 0) {
      frontier->emplace_back(ia, ib);
      return;
    }
    const Cell &a = cells[ia], &b = cells[ib];
    bool same = ia == ib;
    if (!same) {
      double x = a.x - b.x, y = a.y - b.y, scale = 1 + x * x + y * y;
      double radius = a.radius + b.radius;
      double radius2 = radius * radius;
      // The Cauchy kernel varies on scale sqrt(1 + squared distance).
      if (radius2 < theta2 * scale) {
        summarize<Gradient>(ia, ib, x, y, 1 / scale, s);
        return;
      }
    }
    if (a.left < 0 && b.left < 0) {
      double mass = 0;
      for (int i = a.lo; i < a.hi; ++i) {
        const Point &p = points[i];
        for (int j = same ? i + 1 : b.lo; j < b.hi; ++j) {
          const Point &q = points[j];
          double x = p.x - q.x, y = p.y - q.y, w = 1 / (1 + x * x + y * y);
          mass += w;
          if (Gradient) {
            double fx = x * w * w, fy = y * w * w;
            s.forces[i].x += fx;
            s.forces[i].y += fy;
            s.forces[j].x -= fx;
            s.forces[j].y -= fy;
          }
        }
      }
      s.mass += mass;
    } else if (same) {
      visit<Gradient>(a.left, a.left, s, frontier, depth - 1);
      visit<Gradient>(a.left, a.right, s, frontier, depth - 1);
      visit<Gradient>(a.right, a.right, s, frontier, depth - 1);
    } else if (b.left < 0 || (a.left >= 0 && a.radius >= b.radius)) {
      visit<Gradient>(a.left, ib, s, frontier, depth - 1);
      visit<Gradient>(a.right, ib, s, frontier, depth - 1);
    } else {
      visit<Gradient>(ia, b.left, s, frontier, depth - 1);
      visit<Gradient>(ia, b.right, s, frontier, depth - 1);
    }
  }

  void distribute(int idx, Sums &s, NumericMatrix &force) const {
    const Cell &a = cells[idx];
    if (a.left >= 0) {
      for (int id : {a.left, a.right}) {
        const Cell &child = cells[id];
        s.fields[id].inherit(s.fields[idx], child.x - a.x, child.y - a.y);
        distribute(id, s, force);
      }
    } else {
      for (int i = a.lo; i < a.hi; ++i) {
        const Point &p = points[i];
        double x = p.x - a.x, y = p.y - a.y;
        force(p.row, 0) = s.forces[i].x + s.fields[idx].x +
                          s.fields[idx].xx * x + s.fields[idx].xy * y;
        force(p.row, 1) = s.forces[i].y + s.fields[idx].y +
                          s.fields[idx].xy * x + s.fields[idx].yy * y;
      }
    }
  }
};

struct MutualWorkspace {
  explicit MutualWorkspace(int interval) : interval(interval) {}
  const int interval;
  int age = 0;
  std::unique_ptr<MutualTree> tree;

  List gradient(const NumericMatrix &Y, double theta, int threads) {
    if (age == 0) {
      tree.reset();
      tree = std::make_unique<MutualTree>(Y);
    } else {
      tree->refit(Y);
    }
    List result = tree->evaluate<true>(theta, threads);
    age = (age + 1) % interval;
    return result;
  }
};
} // namespace

// [[Rcpp::export(rng = false)]]
SEXP mutual_workspace_cpp(int rebuild_interval) {
  if (rebuild_interval < 1)
    Rcpp::stop("rebuild_interval must be a positive integer");
  return Rcpp::XPtr<MutualWorkspace>(new MutualWorkspace(rebuild_interval),
                                     true, Rf_install("smallvis_mutual"));
}

// [[Rcpp::export(rng = false)]]
List mutual_cauchy_cpp(const NumericMatrix &Y, double theta, bool mass_only,
                       int n_threads = 0, SEXP workspace = R_NilValue) {
  if (n_threads < 0)
    Rcpp::stop("n_threads must be nonnegative");
  if (Y.ncol() != 2 || Y.nrow() < 2)
    Rcpp::stop("Mutual tree requires at least two points in two dimensions");
  if (!std::isfinite(theta) || theta < 0)
    Rcpp::stop("Mutual tree theta must be finite and nonnegative");
  if (!mass_only && workspace != R_NilValue) {
    if (TYPEOF(workspace) != EXTPTRSXP ||
        R_ExternalPtrTag(workspace) != Rf_install("smallvis_mutual"))
      Rcpp::stop("Invalid mutual tree workspace");
    Rcpp::XPtr<MutualWorkspace> state(workspace);
    return state->gradient(Y, theta, n_threads);
  }
  MutualTree tree(Y);
  return mass_only ? tree.evaluate<false>(theta, n_threads)
                   : tree.evaluate<true>(theta, n_threads);
}
