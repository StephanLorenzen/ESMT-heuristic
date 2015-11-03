#ifndef BG_LAZY_H
#define BG_LAZY_H

#include <vector>

#include "steiner/graph.hpp"
#include "steiner/utils/point.hpp"
#include "steiner/utils/bottleneck_graph/bg_simple.hpp"

typedef Utils::Point Point;

/**
 * @class BottleneckGraphLazy
 *
 * Implements a lazy version of the Bottleneck Graph.
 * Bottleneck distances are calculated (and stored) when needed.
 */
class BottleneckGraphLazy : public BottleneckGraphSimple {
public:
  /**
   * Constructor.
   * Constructs the initial BottleneckGraph
   *
   * @param g   The graph to construct the Bottleneck Graph for.
   */
  BottleneckGraphLazy(Graph &g);

  /**
   * Destructor
   */
  ~BottleneckGraphLazy();

  /**
   * Contracts all of the points in the given vector, and recomputes the Bottleneck Graph
   *
   * @param points  The points (indicies) to be contracted
   */
  void contract(const std::vector<unsigned int> &points);
protected:
private:
  /**
   * Get the index in this->edges of the bottleneck edge for the two given points
   *
   * @param i  The first point
   * @param j  The second point
   */
  unsigned int _getEdgeIndex(const unsigned int i, const unsigned int j);
  
  /**
   * Recomputes the bottleneck distance between i and j.
   */
  void _recompute(const unsigned int i, const unsigned int j);
};

#endif // BG_LAZY_H
