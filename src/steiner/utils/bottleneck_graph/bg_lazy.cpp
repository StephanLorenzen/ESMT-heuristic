#include <vector>
#include <unordered_map>

#include "steiner/graph.hpp"

#include "steiner/utils/bottleneck_graph/bg_lazy.hpp"

/*
 * Implementation of BottleneckGraphLazy constructor
 */
BottleneckGraphLazy::BottleneckGraphLazy(Graph &g)
  : BottleneckGraph(g)
{
  for(unsigned int i = 0; i < this->points.size(); i++)
    for(unsigned int j = i+1; j < this->points.size(); j++)
      this->_getEdgeIndex(i,j);
}

/*
 * Implementation of BottleneckGraphLazy destructor
 */
BottleneckGraphLazy::~BottleneckGraphLazy() { }

/**
 * Implementation of BottleneckGraphLazy::mergePoints(...)
 */
void BottleneckGraphLazy::mergePoints(const std::vector<unsigned int> &points) {
  BottleneckGraph::mergePoints(points);
  this->bdist.clear();
}

/**
 * Implementation of BottleneckGraphLazy::_getEdgeIndex(...)
 */
unsigned int BottleneckGraphLazy::_getEdgeIndex(const unsigned int i, const unsigned int j) {
  unsigned int a = i, b = j;
  if(i > j) {
    a = j;
    b = i;
  }
  unsigned long key = this->_key(a,b);
  if(this->bdist.find(key) == this->bdist.end())
    this->_recompute(a,b);
  return this->bdist[key];
}

/*
 * Implementation of BottleneckGraphLazy::_recompute(...)
 */
void BottleneckGraphLazy::_recompute(const unsigned int i, const unsigned int j) {
  // Calculate distance from i to everybody else, including j
  this->_traverse(i,i,0,0);
}

/*
 * Implementation of BottleneckGraphLazy::_traverse(...)
 */
void BottleneckGraphLazy::_traverse(const unsigned int p, const unsigned int cur,
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
      BottleneckGraph::_edge &e = this->edges[nextEdge];
      unsigned int next = cur == e.p ? e.q : e.p;
      this->_traverse(p, next, nextEdge, e.dist > carry ? nextEdge : mEdge);
    }
  }  
}

/*
 * Implementation of BottleneckGraphLazy::_key(...)
 */
unsigned long BottleneckGraphLazy::_key(const unsigned int i, const unsigned int j) {
  return (((long)i) << 32) + j;
}
