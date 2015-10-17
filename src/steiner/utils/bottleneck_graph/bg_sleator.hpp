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
    //double netcost, netmin;

    double cost;
    double maxcost;

    Vertex *dparent;
    double dcost;

    int height;
  };
  typedef Vertex Path;

  void pathDecompose(std::vector< std::vector<unsigned int> > conns,
		     unsigned int cur, unsigned int prev, Path *path);

  // Static tree operations
  Vertex *parent(Vertex *v);
  Vertex *root(Vertex *v);
  double cost(Vertex *v);
  Vertex *maxcost(Vertex *v);
  //Vertex *mincost(Vertex *v);
  //void update(Vertex *v, double x);
  // Dynamic tree operations
  void link(Vertex *v, Vertex *w, double x);
  double cut(Vertex *v);
  void evert(Vertex *v);
  
  /////////////////////////////
  // Operations on paths
  
  // Struct for storing result of a split
  struct SplitResult {
    Path *q;
    Path *r;
    double x, y;
  };

  // Recursive result struct
  struct RRS {
    Vertex *res;
    Vertex *path;
    bool reversed;
    //double grossmin;
    double cost;
    double maxcost;
  };

  // Path operations
  Path *path(Vertex *v);
  Vertex *head(Path *p);
  Vertex *tail(Path *p);
  Vertex *before(Vertex *v);
  void before_rec(Vertex *v, RRS &c);
  Vertex *after(Vertex *v);
  void after_rec(Vertex *v, RRS &c);
  double pcost(Vertex *v);
  void pcost_rec(Vertex *v, RRS &c);
  Vertex *pmaxcost(Path *p);
  //Vertex *pmincost(Path *p);
  //void pupdate(Path *p, double x);
  void reverse(Path *p);
  Path *concatenate(Path *p, Path *q, double x);
  void split(Vertex *v, SplitResult &res);

  // Splice and expose
  Path *splice(Path *p);
  Path *expose(Vertex *v);
  
  // AVL tree operations
  struct DestructResult {
    Vertex *v;
    Vertex *w;
    double x;
  };

  Vertex *construct(Vertex *v, Vertex *w, double x);
  void destruct(Vertex *u, DestructResult &c);
  void rotateleft(Vertex *v);
  void rotateright(Vertex *v);
  void balance(Vertex *v);
  void tsplit(Vertex *r, Vertex *v, std::vector<bool> &path, int i, DestructResult &dr);
  void treepath(Vertex *v, std::vector<bool> &path, RRS &c, bool before);
    
  std::vector<Vertex> _vertices;
  std::vector<Point> &_pointsRef;
};

#endif // BG_SLEATOR_H
