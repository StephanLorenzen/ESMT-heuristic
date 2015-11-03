#ifndef BG_NAIVE_H
#define BG_NAIVE_H

#include <vector>

#include "steiner/graph.hpp"
#include "steiner/utils/point.hpp"
#include "steiner/utils/bottleneck_graph/bg_simple.hpp"

typedef Utils::Point Point;

/**
 * @class BottleneckGraphNaive
 *
 * Implements a naive version of the Bottleneck Graph.
 * Each time a set of points are merged, the entire graph will be recomputed.
 */
class BottleneckGraphNaive : public BottleneckGraphSimple {
public:
  /**
   * Constructor.
   * Constructs the initial BottleneckGraph
   *
   * @param g   The graph to construct the Bottleneck Graph for.
   */
  BottleneckGraphNaive(Graph &g);
  
  /**
   * Destructor
   */
  ~BottleneckGraphNaive();

  /**
   * Contracts all of the points in the given vector, and recomputes the Bottleneck Graph
   *
   * @param points  The points (indicies) to be contracted
   */
  void contract(const std::vector<unsigned int> &points);

protected:
private:  
  /**
   * Recomputes all bottleneck distances.
   */
  void _recompute();
};

#endif // BG_NAIVE_H
