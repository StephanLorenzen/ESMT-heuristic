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
BottleneckGraph::BottleneckGraph(Graph &g) : points(g.getPointsRef()) { }

/*
 * Implementation of BottleneckGraph destructor
 */
BottleneckGraph::~BottleneckGraph() { }

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
