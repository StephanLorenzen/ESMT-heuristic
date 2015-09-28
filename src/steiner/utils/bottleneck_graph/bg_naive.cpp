#include <vector>

#include "steiner/utils/point.hpp"
#include "steiner/utils/bottleneck_graph/bg_naive.hpp"

/*
 * Implementation of BottleneckGraphNaive constructor
 */
BottleneckGraphNaive::BottleneckGraphNaive(Graph &g)
  : BottleneckGraph(g) {
  /*for(int i = 0; i < points.size(); i++) {
    std::cout << i << ":";
    for(int j = i+1; j < points.size(); j++) {
      std::cout << " " << this->bdist[i][j];
    }
    std::cout << std::endl;
    }*/
}

/*
 * Implementation of BottleneckGraphNaive destructor
 */
BottleneckGraphNaive::~BottleneckGraphNaive() { }

/*
 * Implementation of BottleneckGraphNaive::mergePoints(...)
 */
void BottleneckGraphNaive::mergePoints(const std::vector<unsigned int> &points) {
  if(points.size() < 2)
    return;
  
  for(unsigned int i=0; i<points.size(); i++)
    for(unsigned int j=i+1; j<points.size(); j++)
      this->_removeEdge(points[i],points[j]);
  
  for(unsigned int i=1; i<points.size(); i++)
    this->_mergePoints(points[i-1], points[i]);
  //std::cout << "Recompute now..." << std::endl;
  
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
  //std::cout << "merge: "<<p << "," << q << std::endl;
  // To merge two points, we simply add a zero edge between them.
  BottleneckGraph::_edge nedge;
  nedge.dist = 0.0;
  nedge.p    = p;
  nedge.q    = q;
  this->edges.push_back(nedge);
  //std::cout << " added: "<<this->edges.size()-1 << std::endl;
  this->connections[p].push_back(this->edges.size()-1);
  this->connections[q].push_back(this->edges.size()-1);
}

/*
 * Implementation of BottleneckGraphNaive::_removeEdge(...)
 */
void BottleneckGraphNaive::_removeEdge(const unsigned int i, const unsigned int j) {
  //std::cout << "remove: " << i << "," << j << std::endl;
  unsigned int p = i, q = j;
  if(i > j) {
    p = j;
    q = i;
  }
  // Remove the bottleneck edge.
  unsigned int bde = this->bdist[p][q];
  //std::cout << "edge: "<<bde << std::endl;
  BottleneckGraph::_edge &bedge = this->edges[bde];
  bedge.dist = -1;
  std::vector<unsigned int> newConnsP, newConnsQ;
  std::vector<unsigned int> oldConnsP = this->connections[bedge.p];
  std::vector<unsigned int> oldConnsQ = this->connections[bedge.q];
  for(unsigned int i=0; i<oldConnsP.size(); i++)
    if(oldConnsP[i]!=bde)
      newConnsP.push_back(oldConnsP[i]);
  for(unsigned int i=0; i<oldConnsQ.size(); i++)
    if(oldConnsQ[i]!=bde)
      newConnsQ.push_back(oldConnsQ[i]);
  this->connections[bedge.p] = newConnsP;
  this->connections[bedge.q] = newConnsQ;

}
