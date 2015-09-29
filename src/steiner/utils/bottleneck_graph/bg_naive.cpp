#include <vector>

#include "steiner/utils/point.hpp"
#include "steiner/utils/bottleneck_graph/bg_naive.hpp"

/*
 * Implementation of BottleneckGraphNaive constructor
 */
BottleneckGraphNaive::BottleneckGraphNaive(Graph &g)
  : BottleneckGraph(g) {
  // Setup bottleneck distance array
  bdist = new int*[points.size()];
  for(unsigned int i=0; i<this->points.size(); i++)
    bdist[i] = new int[points.size()];

  // Initial computation
  this->_recompute();
}

/*
 * Implementation of BottleneckGraphNaive destructor
 */
BottleneckGraphNaive::~BottleneckGraphNaive() {
  for(unsigned int i=0; i<this->points.size(); i++)
    delete bdist[i];
  delete bdist;
}

/*
 * Implementation of BottleneckGraphNaive::mergePoints(...)
 */
void BottleneckGraphNaive::mergePoints(const std::vector<unsigned int> &points) {
  BottleneckGraph::mergePoints(points);
  // Do the recompute
  this->_recompute();
}

/**
 * Implementation of BottleneckGraphNaive::_getEdgeIndex(...)
 */
unsigned int BottleneckGraphNaive::_getEdgeIndex(const unsigned int i, const unsigned int j) {
  unsigned int a = i, b = j;
  if(i > j) {
    a = j;
    b = i;
  }
  // Assumed to always be valid
  return this->bdist[a][b];
}

/*
 * Implementation of BottleneckGraph::_recompute()
 */
void BottleneckGraphNaive::_recompute() {
  // Calculate bottleneck tree
  for(unsigned int i=0; i<this->points.size(); i++)
    this->_traverse(i,i,0,0);
}

/*
 * Implementation of BottleneckGraph::_traverse(...)
 */
void BottleneckGraphNaive::_traverse(const unsigned int p, const unsigned int cur,
				     const unsigned int prevEdge, const unsigned int mEdge) {
  double carry = this->edges[mEdge].dist;
  if(p < cur)
    // Update bottleneck dist for this point
    this->bdist[p][cur] = mEdge;
  
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
