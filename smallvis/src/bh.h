// BSD 3-Clause License
//
// Copyright (c) 2024, James Melville
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice, this
//   list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
//
// * Neither the name of the copyright holder nor the names of its
//   contributors may be used to endorse or promote products derived from
//   this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

// This is a rough translation of the Barnes-Hut t-SNE cython code from the
// OpenTSNE package originally written by Pavlin Poličar. The original code can
// found at https://github.com/pavlin-policar/openTSNE

#ifndef SMALLVIS_BH_H
#define SMALLVIS_BH_H

#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <numeric>
#include <vector>

#include "threads.h"

namespace smallvis {

class QuadTree {
public:
  // Node structure
  struct Node {
    // Center of mass of points within this node
    double center_of_mass_x = 0.0;
    double center_of_mass_y = 0.0;
    std::size_t num_points = 0;
    Node *children[4] = {nullptr, nullptr, nullptr, nullptr};

    Node(double center_x, double center_y, double length) noexcept
        : center_x(center_x), center_y(center_y), length(length),
          length2(length * length) {}

    ~Node() noexcept {
      for (Node *child : children) {
        delete child;
      }
    }

    bool is_close(double point_x, double point_y, double eps) const noexcept {
      return std::abs(center_of_mass_x - point_x) < eps &&
             std::abs(center_of_mass_y - point_y) < eps;
    }

    bool is_skippable(double point_x, double point_y,
                      double eps) const noexcept {
      return num_points == 0 || (is_leaf && is_close(point_x, point_y, eps));
    }

    // Barnes-Hut criterion rearranged for squared distances and squared theta
    // to avoid square roots and divisions
    bool is_summary(double d2, double theta2) const noexcept {
      return is_leaf || length2 < d2 * theta2;
    }

    void add_point(double point_x, double point_y, double eps) {
      if ((is_leaf && num_points == 0) || is_close(point_x, point_y, eps)) {
        update_center_of_mass(point_x, point_y);
        return;
      }

      if (is_leaf) {
        split();
        std::size_t qid = find_quadrant(center_of_mass_x, center_of_mass_y);
        children[qid]->add_point(center_of_mass_x, center_of_mass_y, eps);
      }

      update_center_of_mass(point_x, point_y);
      std::size_t qid = find_quadrant(point_x, point_y);
      children[qid]->add_point(point_x, point_y, eps);
    }

  private:
    // Center of the node
    double center_x;
    double center_y;
    // (Longest) side length of the node
    double length;
    // length squared used in Barnes Hut criterion
    double length2;
    bool is_leaf = true;

    std::size_t find_quadrant(double point_x, double point_y) const noexcept {
      // Determine which quadrant the point is in -- 0: SW, 1: SE, 2: NW, 3: NE
      if (point_x <= center_x) {
        return point_y <= center_y ? 0 : 2;
      }
      return point_y <= center_y ? 1 : 3;
    }

    void update_center_of_mass(double point_x, double point_y) noexcept {
      num_points += 1;
      center_of_mass_x += (point_x - center_of_mass_x) / num_points;
      center_of_mass_y += (point_y - center_of_mass_y) / num_points;
    }

    void split() {
      is_leaf = false;
      double new_len = length * 0.5;
      double half_len = new_len * 0.5;

      children[0] = new Node(center_x - half_len, center_y - half_len, new_len);
      children[1] = new Node(center_x + half_len, center_y - half_len, new_len);
      children[2] = new Node(center_x - half_len, center_y + half_len, new_len);
      children[3] = new Node(center_x + half_len, center_y + half_len, new_len);
    }
  };
  std::unique_ptr<Node> root;

  QuadTree(const std::vector<double> &data, double eps = EPS) {
    double coords_min_x = std::numeric_limits<double>::max();
    double coords_max_x = std::numeric_limits<double>::lowest();
    double coords_min_y = std::numeric_limits<double>::max();
    double coords_max_y = std::numeric_limits<double>::lowest();

    const std::size_t data_len = data.size();
    for (std::size_t i = 0; i < data_len; i += 2) {
      double x = data[i];
      double y = data[i + 1];

      coords_min_x = std::min(coords_min_x, x);
      coords_max_x = std::max(coords_max_x, x);
      coords_min_y = std::min(coords_min_y, y);
      coords_max_y = std::max(coords_max_y, y);
    }

    double center_x = (coords_max_x + coords_min_x) * 0.5;
    double center_y = (coords_max_y + coords_min_y) * 0.5;
    double length =
        std::max(coords_max_x - coords_min_x, coords_max_y - coords_min_y);

    root = std::make_unique<Node>(center_x, center_y, length);
    for (std::size_t i = 0; i < data_len; i += 2) {
      root->add_point(data[i], data[i + 1], eps);
    }
  }

private:
  constexpr static double EPS = std::numeric_limits<double>::epsilon();
};

void _bh_tsne_negative_gradient_single(const QuadTree::Node &node,
                                       double point_x, double point_y,
                                       double theta2, double eps,
                                       double &gradient_x, double &gradient_y,
                                       double &Zi) {
  // Ensure we do not process empty nodes or simple self-interactions
  if (node.is_skippable(point_x, point_y, eps)) {
    return;
  }

  // Compute the squared Euclidean distance between the point and the center of
  // mass
  double diff_x = node.center_of_mass_x - point_x;
  double diff_y = node.center_of_mass_y - point_y;
  double d2_ij = eps + diff_x * diff_x + diff_y * diff_y;

  // Check if we can use this node as a summary
  if (node.is_summary(d2_ij, theta2)) {
    double w_ij = 1.0 / (1.0 + d2_ij);

    Zi += node.num_points * w_ij;

    // un-normalized -ve grad W*W*dY
    w_ij *= w_ij;
    gradient_x += node.num_points * w_ij * diff_x;
    gradient_y += node.num_points * w_ij * diff_y;

    return;
  }

  // Recursively apply Barnes-Hut to the children
  for (const auto &child : node.children) {
    _bh_tsne_negative_gradient_single(*child, point_x, point_y, theta2, eps,
                                      gradient_x, gradient_y, Zi);
  }
}

// Function to estimate the negative gradient using the Barnes-Hut approximation
void bh_tsne_negative_gradient(const QuadTree &tree,
                               const std::vector<double> &embedding,
                               double theta2, double eps, std::size_t n_threads,
                               std::vector<double> &gradient) {
  std::size_t embedding_len = embedding.size();
  std::vector<double> Zi(std::max(static_cast<std::size_t>(1), n_threads), 0.0);

  const QuadTree::Node *root = tree.root.get();

  // Function to calculate the negative gradient for a single point
  auto worker = [&](std::size_t start, std::size_t end, std::size_t thread_id) {
    for (std::size_t i = start * 2; i < end * 2; i += 2) {
      _bh_tsne_negative_gradient_single(*root, embedding[i], embedding[i + 1],
                                        theta2, eps, gradient[i],
                                        gradient[i + 1], Zi[thread_id]);
    }
  };
  parallel_for(embedding_len / 2, n_threads, worker);

  double Z = std::accumulate(Zi.begin(), Zi.end(), eps);
  // Normalize the gradient
  for (std::size_t i = 0; i < embedding_len; i += 2) {
    gradient[i] /= Z;
    gradient[i + 1] /= Z;
  }
}

void bh_tsne_positive_gradient(const std::vector<std::size_t> &indices,
                               const std::vector<std::size_t> &indptr,
                               const std::vector<double> &P_data,
                               const std::vector<double> &embedding,
                               std::size_t n_threads,
                               std::vector<double> &gradient) {

  auto worker = [&](std::size_t start, std::size_t end) {
    for (std::size_t i = start, i2 = start * 2; i < end; ++i, i2 += 2) {
      for (std::size_t k = indptr[i]; k < indptr[i + 1]; ++k) {
        std::size_t j2 = indices[k] * 2;

        // Compute the direction of the points' attraction and the squared
        // Euclidean distance
        double diff_x = embedding[i2] - embedding[j2];
        double diff_y = embedding[i2 + 1] - embedding[j2 + 1];
        double w_ij_p_ij =
            P_data[k] / (1.0 + diff_x * diff_x + diff_y * diff_y);

        // Compute F_{attr} of point `j` on point `i`
        // W x P x dY
        gradient[i2] += w_ij_p_ij * diff_x;
        gradient[i2 + 1] += w_ij_p_ij * diff_y;
      }
    }
  };
  parallel_for(gradient.size() / 2, n_threads, worker);
}

void _bh_Zi(const QuadTree::Node &node, double point_x, double point_y,
            double theta2, double eps, double &Zi) {
  // Ensure we do not process empty nodes or simple self-interactions
  if (node.is_skippable(point_x, point_y, eps)) {
    return;
  }

  // Compute the squared Euclidean distance between the point and the center of
  // mass
  double diff_x = node.center_of_mass_x - point_x;
  double diff_y = node.center_of_mass_y - point_y;
  double d2_ij = eps + diff_x * diff_x + diff_y * diff_y;

  // Check if we can use this node as a summary
  if (node.is_summary(d2_ij, theta2)) {
    Zi += node.num_points / (1.0 + d2_ij);
    return;
  }

  // Recursively apply Barnes-Hut to the children
  for (const auto &child : node.children) {
    _bh_Zi(*child, point_x, point_y, theta2, eps, Zi);
  }
}

double bh_Z(const QuadTree &tree, const std::vector<double> &embedding,
            double theta2, double eps, std::size_t n_threads) {
  std::vector<double> Zi(std::max(static_cast<std::size_t>(1), n_threads), 0.0);

  const QuadTree::Node *root = tree.root.get();

  // Function to calculate the negative gradient for a single point
  auto worker = [&](std::size_t start, std::size_t end, std::size_t thread_id) {
    const std::size_t end2 = end * 2;
    for (std::size_t i = start * 2; i < end2; i += 2) {
      _bh_Zi(*root, embedding[i], embedding[i + 1], theta2, eps, Zi[thread_id]);
    }
  };
  parallel_for(embedding.size() / 2, n_threads, worker);

  return std::accumulate(Zi.begin(), Zi.end(), eps);
}

double bh_plogq(const std::vector<std::size_t> &indices,
                const std::vector<std::size_t> &indptr,
                const std::vector<double> &P_data,
                const std::vector<double> &embedding, double theta, double eps,
                std::size_t n_threads) {
  QuadTree tree(embedding, eps);
  double mlogZ =
      -std::log(bh_Z(tree, embedding, theta * theta, eps, n_threads));

  std::vector<double> partial_plogq(
      std::max(static_cast<std::size_t>(1), n_threads));
  auto worker = [&](std::size_t start, std::size_t end, std::size_t thread_id) {
    double window_plogq = 0.0;
    for (std::size_t i = start, i2 = start * 2; i < end; ++i, i2 += 2) {
      for (std::size_t k = indptr[i]; k < indptr[i + 1]; ++k) {
        std::size_t j2 = indices[k] * 2;

        // p * log(q) = p * log(w/Z) = p * (log w - log Z) =
        // p * (log(1/(1+d2)) - log Z) = p * (-log(1+d2) - log Z)
        double diff_x = embedding[i2] - embedding[j2];
        double diff_y = embedding[i2 + 1] - embedding[j2 + 1];

        window_plogq +=
            P_data[k] *
            (mlogZ - std::log(1.0 + diff_x * diff_x + diff_y * diff_y));
      }
    }
    partial_plogq[thread_id] += window_plogq;
  };
  parallel_for(embedding.size() / 2, n_threads, worker);

  return std::accumulate(partial_plogq.begin(), partial_plogq.end(), 0.0);
}

std::vector<double> bh_tsne_gradient(const std::vector<std::size_t> &indices,
                                     const std::vector<std::size_t> &indptr,
                                     const std::vector<double> &P_data,
                                     const std::vector<double> &embedding,
                                     double theta, double eps,
                                     std::size_t n_threads) {

  QuadTree tree(embedding, eps);
  std::vector<double> gradient(embedding.size(), 0.0);

  // don't re-arrange the order of the following two function calls!
  bh_tsne_negative_gradient(tree, embedding, theta * theta, eps, n_threads,
                            gradient);
  bh_tsne_positive_gradient(indices, indptr, P_data, embedding, n_threads,
                            gradient);
  return gradient;
}
} // namespace smallvis

#endif // SMALLVIS_BH_H
