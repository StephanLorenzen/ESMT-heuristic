#include <vector>

#include "steiner/graph.hpp"
#include "steiner/utils/utils.hpp"
#include "steiner/utils/point.hpp"
#include "steiner/utils/bottleneck_graph.hpp"

typedef Utils::Point Point;

/*
 * Implementation of BottleneckGraph constructor
 */
BottleneckGraph::BottleneckGraph(const std::vector<Point> &points) {
    
}

/*
 * Implementation of BottleneckGraph destructor
 */
BottleneckGraph::~BottleneckGraph() { }

/*
 * Implementation of BottleneckGraph::distance(...)
 */
float BottleneckGraph::distance(const unsigned int i, const unsigned int j) {
  return 0.0f;
}

/*
 * Implementation of BottleneckGraph::mergePoints(...)
 */
void BottleneckGraph::mergePoints(const std::vector<Point> &points) {
  
}

/*
 * Implementation of BottleneckGraph::_mergePoints(...)
 */
unsigned int _mergePoints(const unsigned int i, const unsigned int j) {
  return 0;
}
