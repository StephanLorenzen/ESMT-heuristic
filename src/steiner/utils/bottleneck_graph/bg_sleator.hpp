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

  void info();
  
protected:
private:

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
    
    double cost;
    double maxcost;
    
    Vertex *dparent;
    double dcost;
    
    unsigned int height;
    
    int idx;
  };
  typedef Vertex Path;

  void treewalk(Vertex *v);
  
  void _cleanUp(Vertex *v);
  
  // Static tree operations
  Vertex *parent(Vertex *v);
  Vertex *root(Vertex *v);
  double cost(Vertex *v);
  Vertex *maxcost(Vertex *v);
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
    bool reversed;
    bool resrev;
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
  Vertex *pmaxcost(Path *p);
  void reverse(Path *p);
  Path *concatenate(Path *p, Path *q, double x);
  void split(Vertex *v, SplitResult &res);

  // Splice and expose
  Path *splice(Path *p);
  Path *expose(Vertex *v);
  
  // AVL tree operations
  struct DestructResult {
    Vertex *left;
    Vertex *right;
    double x;
  };

  Vertex *construct(Vertex *v, Vertex *w, double x);
  void destruct(Vertex *u, DestructResult &c);
  void rotateleft(Vertex *v);
  void rotateright(Vertex *v);
  void unreverse(Vertex *v);
  void balance(Vertex *v);
  int balance_factor(Vertex *v);
  bool is_balanced(Vertex *v);
  void tsplit(Vertex *v, bool right, DestructResult &dr);
  void treepath(Vertex *v, Vertex **r, bool &reversed, bool before);
    
  std::vector<Vertex> _vertices;
  std::vector<Point> &_pointsRef;

  /* TEMP. DEBUG FUNCTIONS */
  void print_tree_rec(Vertex *v, std::string s) {
    if(s.length() >= 10) {
      std::cout << "more children..." << std::endl;
      return;
    }
    std::cout << s << v->idx << "["<<v->reversed<<"]["<<v->height<<"]("<<
      v->cost<<"/"<<v->maxcost<<")" << std::endl;
    if(v->bleft)
      this->print_tree_rec(v->bleft,s+" ");
    if(v->bright)
      this->print_tree_rec(v->bright,s+" ");
  }
  void print_tree(Vertex *v) {
    this->print_tree_rec(v, "");
  }
};

#endif // BG_SLEATOR_H
