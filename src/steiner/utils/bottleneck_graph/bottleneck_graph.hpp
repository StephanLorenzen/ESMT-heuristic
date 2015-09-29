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
   * @param g   The graph to construct the Bottleneck Graph for.
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
  double distance(const unsigned int i, const unsigned int j);
  
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
   * Merges the given points
   *
   * @param points The set of points (indices) to be merged
   */
  virtual void mergePoints(const std::vector<unsigned int> &points);
  
protected:
  /** struct used for storing edges */
  struct _edge {
    double dist;
    unsigned int p, q;
  };

  /**
   * Get the index in this->edges of the bottleneck edge for the two given points
   *
   * @param i  The first point
   * @param j  The second point
   */
  virtual unsigned int _getEdgeIndex(const unsigned int i, const unsigned int j) = 0;
  
  /**
   * Remove the bottleneck edge between the two points
   *
   * @param i  The first point
   * @param j  The second point
   */
  virtual void _removeEdge(const unsigned int i, const unsigned int j);
  
  /**
   * Merge the two points
   *
   * @param i  The first point
   * @param j  The second point
   */
  virtual void _mergePoints(const unsigned int i, const unsigned int j);
  
  /** Stores a reference to the actual points */
  std::vector<Point> &points;
  /** Stores the connecting edges (index into this->edges) for each point */
  std::vector< std::vector<unsigned int> > connections;
  /** Stores info (struct _edge) on each edge */
  std::vector<_edge> edges;
private:
};

#endif /* BOTTLENECK_GRAPH_H */
