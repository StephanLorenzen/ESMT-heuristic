#include <vector>
#include <algorithm>

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
}

/*
 * Implementation of BottleneckGraph destructor
 */
BottleneckGraph::~BottleneckGraph() { }

/**
 * Implementation of BottleneckGraph::distance(...)
 */
double BottleneckGraph::distance(const unsigned int i, const unsigned int j) {
  return this->edges[this->_getEdgeIndex(i,j)].dist;
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

/**
 * Implementation of BottleneckGraph::setBMSTLength(...)
 */
void BottleneckGraph::setBMSTLength(Graph &g) {
  g.setBMSTLength(this->getBMSTLength(g.getPoints()));
}

/**
 * Implementation of BottleneckGraph::mergePoints(...)
 */
void BottleneckGraph::mergePoints(const std::vector<unsigned int> &points) {
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
 * Implementation of BottleneckGraph::_removeEdge(...)
 */
void BottleneckGraph::_removeEdge(const unsigned int i, const unsigned int j) {
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
 * Implementation of BottleneckGraph::_mergePoints(...)
 */
void BottleneckGraph::_mergePoints(const unsigned int i, const unsigned int j) {
  // To merge two points, we simply add a zero edge between them.
  BottleneckGraph::_edge nedge;
  nedge.dist = 0.0;
  nedge.p    = i;
  nedge.q    = j;
  this->edges.push_back(nedge);
  this->connections[i].push_back(this->edges.size()-1);
  this->connections[j].push_back(this->edges.size()-1);
}
