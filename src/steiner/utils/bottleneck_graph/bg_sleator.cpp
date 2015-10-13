#include <vector>
#include <algorithm>

#include "steiner/graph.hpp"
#include "steiner/utils/utils.hpp"
#include "steiner/utils/point.hpp"
#include "steiner/utils/bottleneck_graph/bg_sleator.hpp"

typedef Utils::Point Point;
typedef Utils::Edge  Edge;

/*
 * Implementation of BottleneckGraphSleator constructor
 */
BottleneckGraphSleator::BottleneckGraphSleator(Graph &g) : BottleneckGraph(g) { }

/*
 * Implementation of BottleneckGraphSleator destructor
 */
BottleneckGraphSleator::~BottleneckGraphSleator() { }

/**
 * Implementation of BottleneckGraphSleator::distance(...)
 */
double BottleneckGraphSleator::distance(const unsigned int i, const unsigned int j) {
  return 0.0;
}

/**
 * Implementation of BottleneckGraphSleator::mergePoints(...)
 */
void BottleneckGraphSleator::mergePoints(const std::vector<unsigned int> &points) {
  return;
}

/**
 * Implementation of BottleneckGraphSleator::_getEdgeIndex(...)
 */
unsigned int BottleneckGraphSleator::_getEdgeIndex(const unsigned int i, const unsigned int j) {
  return 0;
}
