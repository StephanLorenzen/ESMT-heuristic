#include <vector>

#include "steiner/graph.hpp"
#include "steiner/utils/utils.hpp"
#include "steiner/utils/point.hpp"
#include "steiner/utils/bottleneck_graph/bottleneck_graph.hpp"

typedef Utils::Point Point;
typedef Utils::Edge  Edge;

/*
 * Implementation of BottleneckGraph constructor
 */
BottleneckGraph::BottleneckGraph(Graph &g)
  : points(g.getPointsRef()), connections(g.n()), edges(g.n()-1)
{
  // Compute MST
  std::vector<Edge> edges;
  Utils::MSTKruskalEdges(g, edges);
  
  // Setup connection lists
  for(unsigned int i=0; i<edges.size(); i++) {
    int i0 = edges[i].i0, i1 = edges[i].i1;
    this->connections[i0].push_back(i);
    this->connections[i1].push_back(i);
    this->edges[i].p = i0 < i1 ? i0 : i1;
    this->edges[i].q = i0 < i1 ? i1 : i0;
    this->edges[i].dist = edges[i].length;
  }
  // Setup bottleneck distance array
  bdist = new unsigned int*[points.size()];
  for(unsigned int i=0; i<this->points.size(); i++)
    bdist[i] = new unsigned int[points.size()];
 
  this->_recompute();
}

/*
 * Implementation of BottleneckGraph destructor
 */
BottleneckGraph::~BottleneckGraph() {
  for(unsigned int i=0; i<this->points.size(); i++)
    delete bdist[i];
  delete bdist;
}

/*
 * Implementation of BottleneckGraph::distance(...)
 */
double BottleneckGraph::distance(const unsigned int i, const unsigned int j) {
  unsigned int a = i, b = j;
  if(i > j) {
    a = j;
    b = i;
  }
  return this->edges[bdist[a][b]].dist;
}

/*
 * Implementation of BottleneckGraph::getBMSTLength(...)
 */
double BottleneckGraph::getBMSTLength(std::vector<unsigned int> &pointIdxs) {
  unsigned int i, j;
  Graph g(pointIdxs, this->points);
  for(i = 0; i < g.n(); i++) {
    for(j = i+1; j < g.n(); j++)
      g.getEdges().push_back(Edge(i, j, this->distance(g.pidx(i),g.pidx(j))));
  }
  Utils::MSTKruskalMod(g);
  return g.getLength(false);
}

/*
 * Implementation of BottleneckGraph::mergePoints(...)
 */
//void BottleneckGraph::mergePoints(const std::vector<Point> &points) {
  
//}

/*
 * Implementation of BottleneckGraph::_mergePoints(...)
 */
//unsigned int BottleneckGraph::_mergePoints(const unsigned int i, const unsigned int j) {
  
  
//  return 0;
//}

/*
 * Implementation of BottleneckGraph::_recompute()
 */
void BottleneckGraph::_recompute() {
  // Calculate bottleneck tree
  for(unsigned int i=0; i<this->points.size(); i++)
    this->_traverse(i,i,0,0);
}

/*
 * Implementation of BottleneckGraph::_recompute(...)
 */
void BottleneckGraph::_recompute(std::vector<unsigned int> &points) {
  // Calculate bottleneck tree
  for(unsigned int i=0; i<points.size(); i++)
    this->_traverse(points[i],points[i],0,0);
}

/*
 * Implementation of BottleneckGraph::_traverse(...)
 */
void BottleneckGraph::_traverse(const unsigned int p, const unsigned int cur,
	       const unsigned int prevEdge, const unsigned int mEdge) {
  double carry = this->edges[mEdge].dist;
  if(p < cur)
    // Update bottleneck dist for this point
    bdist[p][cur] = mEdge;
  
  std::vector<unsigned int> &conn = this->connections[cur];
  for(unsigned int i = 0; i<conn.size(); i++) {
    unsigned int nextEdge = conn[i];
    if(p == cur || nextEdge != prevEdge) {
      // Get edge
      _edge &e = this->edges[nextEdge];
      unsigned int next = cur == e.p ? e.q : e.p;
      this->_traverse(p, next, nextEdge, e.dist > carry ? nextEdge : mEdge);
    }
  }
}
