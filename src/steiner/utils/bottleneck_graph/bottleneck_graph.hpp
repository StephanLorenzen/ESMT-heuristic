#ifndef BOTTLENECK_GRAPH_H
#define BOTTLENECK_GRAPH_H

#include <vector>

#include "steiner/utils/point.hpp"

typedef Utils::Point Point;

/**
 * @class BottleneckGraph
 *
 * Implements the Bottleneck Graph.
 */
class BottleneckGraph {
public:
  /**
   * Constructor.
   * Constructs the initial BottleneckGraph
   */
  BottleneckGraph(Graph &g);
  
  /**
   * Destructor
   */
  virtual ~BottleneckGraph();
  
  /**
   * Gets the Bottleneck distance between points indexed by i and j
   */
  double distance(const unsigned int i, const unsigned int j);
  
  /**
   * Gets the Bottleneck Minimum Spanning Tree length for a set of points
   */
  double getBMSTLength(std::vector<unsigned int> &points);
  
protected:
  /**
   * Recomputes all bottleneck distances.
   */
  void _recompute();
  
  /**
   * Recomputes the bottleneck distances for all given points (indices).
   */
  void _recompute(std::vector<unsigned int> &points);

  /**
   * Initial traverse of the tree to find the Bottleneck distance for a single point.
   */
  void _traverse(const unsigned int p, const unsigned int cur,
		 const unsigned int prevEdge, const unsigned int mEdge);

  struct _edge {
    double dist;
    unsigned int p, q;
  };

  std::vector<Point> &points;
  unsigned int** bdist;
  std::vector< std::vector<unsigned int> > connections;
  std::vector<_edge> edges;
private:
};

#endif /* BOTTLENECK_GRAPH_H */
