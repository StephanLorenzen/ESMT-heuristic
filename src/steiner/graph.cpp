#include <vector>
#include <iostream>
#include <assert.h>

#include <steiner/graph.hpp>
#include <steiner/utils/utils.hpp>
#include <steiner/utils/point.hpp>

Graph::Graph(std::vector<Utils::Point> &pointsRef)
  : pointsRef(&pointsRef)
{
  // Assume all points, fill pidxs
  for(unsigned int i=0; i<pointsRef.size(); i++)
    this->pidxs.push_back(i);
  this->mst_length  = -1;
  this->bmst_length = -1;
}

Graph::Graph(std::vector<Pidx> &points, std::vector<Utils::Point> &pointsRef)
  : pointsRef(&pointsRef), pidxs(points)
{
  this->mst_length  = -1;
  this->bmst_length = -1;
}

Graph::Graph(std::vector<Pidx> &points, std::vector<Utils::Point> &pointsRef,
	     std::vector<Utils::Edge> &edges)
  : pointsRef(&pointsRef), pidxs(points), edges(edges) {
  this->mst_length  = -1;
  this->bmst_length = -1;
}

Graph::~Graph() {}

/*
 * Getter for the points
 */
std::vector<Pidx> &Graph::getPoints() {
  return this->pidxs;
}
std::vector<Utils::Point> &Graph::getPointsRef() {
  return *this->pointsRef;
}

/*
 * Getter for the edges
 */
std::vector<Utils::Edge> &Graph::getEdges() {
  return this->edges;
}

Point &Graph::getPoint(unsigned int index) const {
  assert(index < this->n());
  return (*this->pointsRef)[this->pidx(index)];
}

Pidx Graph::pidx(unsigned int index) const {
  return this->pidxs[index];
}
unsigned int Graph::n() const {
  return this->pidxs.size();
}
double Graph::getMSTLength() const {
  /*  if(this->mst_length > 0)
    return this->mst_length;
  
  Graph mst = Utils::MSTKruskal(*this);
  this->mst_length = mst.getLength();
  */
  return this->mst_length;
}

double Graph::getBMSTLength() const {
  return this->bmst_length;
}

void Graph::setMSTLength(double l) {
  this->mst_length = l;
}
void Graph::setBMSTLength(double l) {
  this->bmst_length = l;
}

void Graph::calculateMSTLength() {

}
void Graph::calculateBMSTLength() {

}

double Graph::getLength(bool use_euclidean) const {
  std::vector<Edge>::const_iterator it;
  double result = 0.0;
  for(it = this->edges.begin(); it != this->edges.end(); it++) {
    result += (use_euclidean
	       ? Utils::length(this->getPoint(it->i0), this->getPoint(it->i1))
	       : it->length);
  }
  return result;
}

unsigned int Graph::dimension() const {
  if(this->pointsRef->size() < 1)
    return 0;
  return (*this->pointsRef)[0].dim();
}


Graph &Graph::operator=(const Graph &other) {
  if(this != &other) {
    this->pointsRef   = other.pointsRef;
    this->pidxs       = other.pidxs;
    this->edges       = other.edges;
    this->mst_length  = other.mst_length;
    this->bmst_length = other.bmst_length;
  }
  return *this;
}

std::ostream& operator<<(std::ostream& os, Graph &g) {
  unsigned int i, n = g.n();
  std::string indent("  ");
  std::vector<Utils::Edge> &edges = g.getEdges();
  os << "# Terminals" << std::endl;
  for(i = 0; i < n; i++)
    os << indent << i << " " << g.getPoint(i) << std::endl;
  os << std::endl << "# Edges" << std::endl;
  for(i = 0; i < edges.size(); i++) {
    os << indent << "(" << edges[i].i0
       << " " << edges[i].i1 << ")" << std::endl;
  }
  os << std::endl << "# |MST|: " << g.getMSTLength() << std::endl;
  return os;
}
