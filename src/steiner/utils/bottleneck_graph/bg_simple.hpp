#ifndef BG_SIMPLE_H
#define BG_SIMPLE_H

#include <vector>
#include <unordered_map>

#include "steiner/graph.hpp"
#include "steiner/utils/point.hpp"
#include "steiner/utils/bottleneck_graph/bottleneck_graph.hpp"

typedef Utils::Point Point;

/**
 * @class BottleneckGraphSimple
 *
 * Implements an abstract version of the Bottleneck Graph, with a simple structure for storing
 * the MST.
 */
class BottleneckGraphSimple : public BottleneckGraph {
public:
  /**
   * Constructor.
   * Constructs the initial BottleneckGraph
   *
   * @param g   The graph to construct the Bottleneck Graph for.
   */
  BottleneckGraphSimple(Graph &g);

  /**
   * Destructor
   */
   ~BottleneckGraphSimple();

  /**
   * Gets the Bottleneck distance between points indexed by i and j
   *
   * @param i  Index of the first point
   * @param j  Index of the second point
   * @return   The Bottleneck distance between i and j
   */
  double distance(const unsigned int i, const unsigned int j);

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
  virtual unsigned int _getEdgeIndex(const unsigned int i, const unsigned int j);
  
  /**
   * Remove the bottleneck edge between the two points
   *
   * @param i  The first point
   * @param j  The second point
   */
  void _removeEdge(const unsigned int i, const unsigned int j);
  
  /**
   * Merge the two points
   *
   * @param i  The first point
   * @param j  The second point
   */
  void _mergePoints(const unsigned int i, const unsigned int j);

  /**
   * Traverses the tree to find the Bottleneck distance for a single point p.
   */
  void _traverse(const unsigned int p, const unsigned int cur,
		 const unsigned int prevEdge, const unsigned int mEdge);
  
  /**
   * Computes the key into this->bdist, for an edge (i,j), i < j
   */
  unsigned long _key(const unsigned int i, const unsigned int j);
  
  /** Stores the connecting edges (index into this->edges) for each point */
  std::vector< std::vector<unsigned int> > connections;
  /** Stores info (struct _edge) on each edge */
  std::vector<_edge> edges;
  /** Stores the bottleneck distances */
  std::unordered_map<unsigned long, unsigned int> bdist;

private:
};

#endif // BG_SIMPLE_H
