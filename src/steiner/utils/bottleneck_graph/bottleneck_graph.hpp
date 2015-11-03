#ifndef BOTTLENECK_GRAPH_H
#define BOTTLENECK_GRAPH_H

#include <vector>

#include "steiner/utils/point.hpp"

typedef Utils::Point Point;

/**
 * @class BottleneckGraph
 *
 * Implements the abstract Bottleneck Graph.
 */
class BottleneckGraph {
public:
  /**
   * Constructor.
   * Constructs the initial BottleneckGraph
   *
   * @param g The graph to construct the BottleneckGraph for
   */
  BottleneckGraph(Graph &g);
  
  /**
   * Destructor
   */
  virtual ~BottleneckGraph();
  
  /**
   * Gets the Bottleneck distance between points indexed by i and j
   *
   * @param i  Index of the first point
   * @param j  Index of the second point
   * @return   The Bottleneck distance between i and j
   */
  virtual double distance(const unsigned int i, const unsigned int j) = 0;
  
  /**
   * Gets the Bottleneck Minimum Spanning Tree length for a set of points
   *
   * @param points  The set of points (indices) to calculate the BMST length for
   * @return        The BMST length 
   */
  double getBMSTLength(std::vector<unsigned int> &points);
  
  /**
   * Sets the Bottleneck Minimum Spanning Tree length of the given graph
   *
   * @param g  The graph to be updated with new BMST length
   */
  void setBMSTLength(Graph &g);

  /**
   * Contracts the given points
   *
   * @param points The set of points (indices) to be contracted
   */
  virtual void contract(const std::vector<unsigned int> &points) = 0;
  
protected:
  std::vector<Point> &points;
private:
};

#endif /* BOTTLENECK_GRAPH_H */
