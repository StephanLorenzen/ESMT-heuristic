#ifndef BG_SLEATOR_H
#define BG_SLEATOR_H

#include <vector>

#include "steiner/graph.hpp"
#include "steiner/utils/point.hpp"
#include "steiner/utils/bottleneck_graph/bottleneck_graph.hpp"

typedef Utils::Point Point;

/**
 * @class BottleneckGraphSleator
 *
 * Implements a lazy version of the Bottleneck Graph.
 * Bottleneck distances are calculated (and stored) when needed.
 */
class BottleneckGraphSleator : public BottleneckGraph {
public:
  /**
   * Constructor.
   * Constructs the initial BottleneckGraph
   *
   * @param g   The graph to construct the Bottleneck Graph for.
   */
  BottleneckGraphSleator(Graph &g);

  /**
   * Destructor
   */
  ~BottleneckGraphSleator();
  
  /**
   * Gets the Bottleneck distance between points indexed by i and j
   *
   * @param i  Index of the first point
   * @param j  Index of the second point
   * @return   The Bottleneck distance between i and j
   */
  double distance(const unsigned int i, const unsigned int j);

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

  //////////////////////////////////////////////////////
  // Functions needed for Sleator/Trajan implementation
  
  // Structures needed
  struct Vertex {
    Vertex *bparent;
    Vertex *bleft;
    Vertex *bright;
    Vertex *bhead;
    Vertex *btail;
    bool external, reversed;
    double netcost, netmin;

    Vertex *dparent;
    double dcost;
  };
  typedef Vertex Path;

  // Static tree operations
  Vertex *parent(Vertex *v);
  Vertex *root(Vertex *v);
  double cost(Vertex *v);
  double mincost(Vertex *v);
  void update(Vertex *v, double x);
  // Dynamic tree operations
  Vertex *link(Vertex *v, Vertex *u, double x);
  Vertex *cut(Vertex *v);
  void evert(Vertex *v);
  
  /////////////////////////////
  // Operations on paths
  
  // Struct for storing result of a split
  struct SplitResult {
    Vertex *q;
    Vertex *p;
    double x, y;
  };

  // Path operations
  Path *path(Vertex *v);
  Vertex *head(Path *p);
  Vertex *tail(Path *p);
  Vertex *before(Vertex *v);
  //bool before_rec(const unsigned int i, int &r, bool &st);
  Vertex *after(Vertex *v);
  //bool after_rec(const unsigned int i, int &r, bool &st);
  double pcost(Vertex *v);
  //bool pcost_rec(unsigned int i, int &r, bool &st, double &g, double &grossmin);
  double pmincost(Path *p);
  //double pmincost_rec(unsigned int u, bool &r);
  void pupdate(Path *p, double x);
  void reverse(Path *p);
  Path *concatenate(Path *p, Path *q, double x);
  void split(Path *p, SplitResult &res);
  
  // Splice and expose
  Path *splice(Path *p);
  Path *expose(Vertex *v);
  
  std::vector<Vertex> _vertices;
};

#endif // BG_SLEATOR_H
