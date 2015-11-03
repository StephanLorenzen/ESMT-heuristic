#include <vector>
#include <algorithm>

#include "steiner/graph.hpp"
#include "steiner/utils/utils.hpp"
#include "steiner/utils/point.hpp"
#include "steiner/utils/bottleneck_graph/bg_simple.hpp"

typedef Utils::Point Point;
typedef Utils::Edge  Edge;

/*
 * Implementation of BottleneckGraphSimple constructor
 */
BottleneckGraphSimple::BottleneckGraphSimple(Graph &g)
  : BottleneckGraph(g), connections(g.n()), edges(g.n()-1)
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
}

/*
 * Implementation of BottleneckGraphSimple destructor
 */
BottleneckGraphSimple::~BottleneckGraphSimple() { }

/**
 * Implementation of BottleneckGraphSimple::distance(...)
 */
double BottleneckGraphSimple::distance(const unsigned int i, const unsigned int j) {
  return this->edges[this->_getEdgeIndex(i,j)].dist;
}

/**
 * Implementation of BottleneckGraphSimple::contract(...)
 */
void BottleneckGraphSimple::contract(const std::vector<unsigned int> &points) {
  if(points.size() < 2)
    return;
  
  // Remove the bottleneck edges
  for(unsigned int i=0; i<points.size(); i++)
    for(unsigned int j=i+1; j<points.size(); j++)
      this->_removeEdge(points[i],points[j]);
  
  // Merge the points together
  for(unsigned int i=1; i<points.size(); i++)
    this->_mergePoints(points[i-1], points[i]);
}

/**
 * Implementation of BottleneckGraphSimple::_getEdgeIndex(...)
 */
unsigned int BottleneckGraphSimple::_getEdgeIndex(const unsigned int i, const unsigned int j) {
  unsigned int a = i, b = j;
  if(i > j) {
    a = j;
    b = i;
  }
  unsigned long key = this->_key(a,b);
  return this->bdist[key];
}

/**
 * Implementation of BottleneckGraphSimple::_removeEdge(...)
 */
void BottleneckGraphSimple::_removeEdge(const unsigned int i, const unsigned int j) {
  unsigned int e = this->_getEdgeIndex(i,j);
  
  // Might be removed already
  _edge &edge = this->edges[e];
  if(edge.dist == 0.0)
    // Already removed or merge edge
    return;
  edge.dist = 0.0;

  // Remove from cons lists
  std::vector<unsigned int> &connsP = this->connections[edge.p];
  std::vector<unsigned int> &connsQ = this->connections[edge.q];
  connsP.erase(std::remove(connsP.begin(), connsP.end(), e), connsP.end());
  connsQ.erase(std::remove(connsQ.begin(), connsQ.end(), e), connsQ.end());
}

/**
 * Implementation of BottleneckGraphSimple::_mergePoints(...)
 */
void BottleneckGraphSimple::_mergePoints(const unsigned int i, const unsigned int j) {
  // To merge two points, we simply add a zero edge between them.
  BottleneckGraphSimple::_edge nedge;
  nedge.dist = 0.0;
  nedge.p    = i;
  nedge.q    = j;
  this->edges.push_back(nedge);
  this->connections[i].push_back(this->edges.size()-1);
  this->connections[j].push_back(this->edges.size()-1);
}

/*
 * Implementation of BottleneckGraphSimple::_traverse(...)
 */
void BottleneckGraphSimple::_traverse(const unsigned int p, const unsigned int cur,
				      const unsigned int prevEdge, const unsigned int mEdge) {
  double carry = this->edges[mEdge].dist;
  if(p < cur)
    // Update bottleneck dist for this point
    this->bdist[this->_key(p,cur)] = mEdge;
  
  std::vector<unsigned int> &conn = this->connections[cur];
  for(unsigned int i = 0; i<conn.size(); i++) {
    unsigned int nextEdge = conn[i];
    if(p == cur || nextEdge != prevEdge) {
      // Get edge
      BottleneckGraphSimple::_edge &e = this->edges[nextEdge];
      unsigned int next = cur == e.p ? e.q : e.p;
      this->_traverse(p, next, nextEdge, e.dist > carry ? nextEdge : mEdge);
    }
  }  
}

/*
 * Implementation of BottleneckGraphSimple::_key(...)
 */
unsigned long BottleneckGraphSimple::_key(const unsigned int i, const unsigned int j) {
  return (((long)i) << 32) + j;
}
