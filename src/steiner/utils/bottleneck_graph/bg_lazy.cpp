#include <vector>
#include <unordered_map>

#include "steiner/graph.hpp"

#include "steiner/utils/bottleneck_graph/bg_lazy.hpp"

/*
 * Implementation of BottleneckGraphLazy constructor
 */
BottleneckGraphLazy::BottleneckGraphLazy(Graph &g) : BottleneckGraphSimple(g) { }

/*
 * Implementation of BottleneckGraphLazy destructor
 */
BottleneckGraphLazy::~BottleneckGraphLazy() { }

/**
 * Implementation of BottleneckGraphLazy::mergePoints(...)
 */
void BottleneckGraphLazy::mergePoints(const std::vector<unsigned int> &points) {
  BottleneckGraphSimple::mergePoints(points);
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
