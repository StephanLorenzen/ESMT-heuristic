#include <vector>

#include "steiner/utils/point.hpp"
#include "steiner/utils/bottleneck_graph/bg_naive.hpp"

/*
 * Implementation of BottleneckGraphNaive constructor
 */
BottleneckGraphNaive::BottleneckGraphNaive(const std::vector<Point> &points)
  : BottleneckGraph(points) { }

/*
 * Implementation of BottleneckGraphNaive destructor
 */
BottleneckGraphNaive::~BottleneckGraphNaive() { }

/*
 * Implementation of BottleneckGraphNaive::mergePoints(...)
 */
void BottleneckGraphNaive::mergePoints(const std::vector<unsigned int> &points) {
  for(unsigned int i=1; i<points.size(); i++)
    this->_mergePoints(points[i-1], points[i]);
  this->_recompute();
}

/*
 * Implementation of BottleneckGraphNaive::_mergePoints(...)
 */
void BottleneckGraphNaive::_mergePoints(const unsigned int i, const unsigned int j) {
  unsigned int p = i, q = j;
  if(i > j) {
    p = j;
    q = i;
  }
  // To merge two points, we simply add a zero edge between them.
  BottleneckGraph::_edge nedge;
  nedge.dist = 0.0;
  nedge.p    = p;
  nedge.q    = q;
  this->edges.push_back(nedge);
  this->connections[p].push_back(q);
  this->connections[q].push_back(p);
  
  // Then we remove their bottleneck edge.
  unsigned int bde = this->bdist[p][q];
  BottleneckGraph::_edge &bedge = this->edges[bde];
  std::vector<unsigned int> newConnsP, newConnsQ;
  std::vector<unsigned int> oldConnsP = this->connections[bedge.p];
  std::vector<unsigned int> oldConnsQ = this->connections[bedge.q];
  for(unsigned int i=0; i<oldConnsP.size(); i++)
    if(oldConnsP[i]!=bedge.q)
      newConnsP.push_back(oldConnsP[i]);
  for(unsigned int i=0; i<oldConnsQ.size(); i++)
    if(oldConnsQ[i]!=bedge.p)
      newConnsQ.push_back(oldConnsQ[i]);
  this->connections[bedge.p] = newConnsP;
  this->connections[bedge.q] = newConnsQ;
}

