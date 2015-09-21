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
BottleneckGraph::BottleneckGraph(const std::vector<Point> &points) : connections(points.size()), edges(points.size()-1) {
  // Compute MST
  this->points = points;
  Graph g(this->points);
  Graph mst = Utils::MSTKruskal(g);
  
  // Setup connection lists
  std::vector<Edge> edges = mst.getEdges();
  for(unsigned int i=0; i<edges.size(); i++) {
    int i0 = edges[i].i0, i1 = edges[i].i1;
    this->connections[i0].push_back(i1);
    this->connections[i1].push_back(i0);
    this->edges[i].p = i0 < i1 ? i0 : i1;
    this->edges[i].q = i0 < i1 ? i1 : i0;
    this->edges[i].dist = edges[i].length;
  }

  // Setup bottleneck distance array
  bdist = new unsigned int*[points.size()];
  for(unsigned int i=0; i<this->points.size(); i++)
    bdist[i] = new unsigned int[points.size()];
 
  std::vector<unsigned int> ps;
  for(unsigned int i=0; i<this->points.size(); i++)
    ps.push_back(i);
  this->_recompute(ps);
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
  return bdist[a][b];
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
  //std::sort(points);
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
  for(unsigned int nextEdge = 0; nextEdge<conn.size(); nextEdge++) {
    if(p == cur || nextEdge != prevEdge) {
      // Get edge
      _edge &e = this->edges[nextEdge];
      unsigned int next = cur == e.p ? e.q : e.p;
      this->_traverse(p, next, nextEdge, e.dist > carry ? nextEdge : mEdge);
    }
  }
}
