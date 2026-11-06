// BSD 3-Clause License
//
// Copyright (c) 2018, Pavlin Poličar
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

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <numeric>
#include <stdexcept>
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
    std::unique_ptr<Node> children[4];

    Node(double center_x, double center_y, double length) noexcept
        : center_x(center_x), center_y(center_y), length(length),
          length2(length * length) {}

    // Never summarize a subtree containing the query, even at a large theta.
    bool is_summary(double d2, double theta2,
                    bool contains_query) const noexcept {
      return is_leaf || (!contains_query && length2 < d2 * theta2);
    }

    void add_point(double point_x, double point_y, std::size_t count = 1) {
      if (num_points == 0) {
        center_of_mass_x = point_x;
        center_of_mass_y = point_y;
        num_points = count;
        return;
      }

      if (is_leaf) {
        if (center_of_mass_x == point_x && center_of_mass_y == point_y) {
          num_points += count;
          return;
        }
        split();
        const auto qid = find_quadrant(center_of_mass_x, center_of_mass_y);
        children[qid]->add_point(center_of_mass_x, center_of_mass_y,
                                 num_points);
      }

      update_center_of_mass(point_x, point_y, count);
      const auto qid = find_quadrant(point_x, point_y);
      children[qid]->add_point(point_x, point_y, count);
    }

    std::size_t find_quadrant(double point_x, double point_y) const noexcept {
      // 0: SW, 1: SE, 2: NW, 3: NE
      if (point_x <= center_x) {
        return point_y <= center_y ? 0 : 2;
      }
      return point_y <= center_y ? 1 : 3;
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

    void update_center_of_mass(double point_x, double point_y,
                               std::size_t count) noexcept {
      num_points += count;
      center_of_mass_x += (point_x - center_of_mass_x) * count / num_points;
      center_of_mass_y += (point_y - center_of_mass_y) * count / num_points;
    }

    void split() {
      double new_len = length * 0.5;
      double half_len = new_len * 0.5;

      children[0] = std::make_unique<Node>(center_x - half_len,
                                           center_y - half_len, new_len);
      children[1] = std::make_unique<Node>(center_x + half_len,
                                           center_y - half_len, new_len);
      children[2] = std::make_unique<Node>(center_x - half_len,
                                           center_y + half_len, new_len);
      children[3] = std::make_unique<Node>(center_x + half_len,
                                           center_y + half_len, new_len);
      is_leaf = false;
    }
  };
  std::unique_ptr<Node> root;

  explicit QuadTree(const std::vector<double> &data) {
    double coords_min_x = std::numeric_limits<double>::max();
    double coords_max_x = std::numeric_limits<double>::lowest();
    double coords_min_y = std::numeric_limits<double>::max();
    double coords_max_y = std::numeric_limits<double>::lowest();

    const std::size_t data_len = data.size();
    for (std::size_t i = 0; i < data_len; i += 2) {
      double x = data[i];
      double y = data[i + 1];
      if (!std::isfinite(x) || !std::isfinite(y)) {
        throw std::invalid_argument("Barnes-Hut coordinates must be finite");
      }

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
      root->add_point(data[i], data[i + 1]);
    }
  }
};
} // namespace smallvis

#endif // SMALLVIS_BH_H
