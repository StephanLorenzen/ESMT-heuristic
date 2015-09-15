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
  BottleneckGraph(const std::vector<Point> &points);
  
  /**
   * Destructor
   */
  ~BottleneckGraph();
  
  /**
   * Gets the Bottleneck distance between points indexed by i and j
   */
  float distance(const unsigned int i, const unsigned int j);
  
  /**
   * Merges all of the points in the given vector, and updates the Bottleneck Graph
   */
  void mergePoints(const std::vector<Point> &points);
  
protected:
private:
  /**
   * Merges the two points indexed by i and j, and updates the Bottleneck Graph
   */
  unsigned int _mergePoints(const unsigned int i, const unsigned int j);
  
  struct _edge {
    float dist;
    unsigned int p, q;
  };

  std::vector<Point> points;
  std::vector< std::vector<unsigned int> > connections;
  std::vector<_edge> edges;
};

#endif /* BOTTLENECK_GRAPH_H */
