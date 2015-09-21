#ifndef BG_NAIVE_H
#define BG_NAIVE_H

#include <vector>

#include "steiner/utils/point.hpp"
#include "steiner/utils/bottleneck_graph/bottleneck_graph.hpp"

typedef Utils::Point Point;

/**
 * @class BottleneckGraphNaive
 *
 * Implements a naive version of the Bottleneck Graph.
 * Each time a set of points are merged, the entire graph will be recomputed.
 */
class BottleneckGraphNaive : public BottleneckGraph {
public:
  BottleneckGraphNaive(const std::vector<Point> &points);
  ~BottleneckGraphNaive();

  /**
   * Merges all of the points in the given vector, and recomputes the Bottleneck Graph
   */
  void mergePoints(const std::vector<unsigned int> &points);
protected:
private:
  /**
   * Merges the two points indexed by i and j
   */
  void _mergePoints(const unsigned int i, const unsigned int j);
};

#endif // BG_NAIVE_H
