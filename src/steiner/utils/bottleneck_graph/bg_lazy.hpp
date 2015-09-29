#ifndef BG_LAZY_H
#define BG_LAZY_H

#include <vector>

#include "steiner/graph.hpp"
#include "steiner/utils/point.hpp"
#include "steiner/utils/bottleneck_graph/bottleneck_graph.hpp"

typedef Utils::Point Point;

/**
 * @class BottleneckGraphLazy
 *
 * Implements a lazy version of the Bottleneck Graph.
 * Bottleneck distances are calculated (and stored) when needed.
 */
class BottleneckGraphLazy : public BottleneckGraph {
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
   * Merges all of the points in the given vector, and recomputes the Bottleneck Graph
   *
   * @param points  The points (indicies) to be merged
   */
  void mergePoints(const std::vector<unsigned int> &points);
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
  
  /**
   * Traverses the tree to find the Bottleneck distance for a single point p.
   */
  void _traverse(const unsigned int p, const unsigned int cur,
		 const unsigned int prevEdge, const unsigned int mEdge);
  
  /**
   * Computes the key into this->bdist, for an edge (i,j), i < j
   */
  unsigned long _key(const unsigned int i, const unsigned int j);

  /** Stores the bottleneck distances */
  std::unordered_map<unsigned long, unsigned int> bdist;
};

#endif // BG_LAZY_H
