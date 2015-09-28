#include <vector>
#include <iostream>

#include "steiner/steiner_tree.hpp"
#include "steiner/utils/point.hpp"

SteinerTree::SteinerTree(std::vector<Utils::Point> &pointsRef)
  : Graph(pointsRef)
{
  this->ratio      = -1.0;
  this->b_ratio    = -1.0;
  this->smt_length = -1.0;
}

SteinerTree::SteinerTree(std::vector<Pidx> &terminals, std::vector<Utils::Point> &pointsRef)
  : Graph(terminals, pointsRef)
{
  this->ratio      = -1.0;
  this->b_ratio    = -1.0;
  this->smt_length = -1.0;
}

SteinerTree::SteinerTree(std::vector<Pidx> &terminals,
			 std::vector<Utils::Point> &pointsRef,
			 std::vector<Utils::Edge> &edges)
  : Graph(terminals, pointsRef, edges)
{
  this->ratio      = -1.0;
  this->b_ratio    = -1.0;
  this->smt_length = -1.0;
}

SteinerTree::~SteinerTree() {}

double SteinerTree::getSMTLength() const {
  return this->smt_length;
}

void SteinerTree::setSMTLength(double l) {
  this->smt_length = l;
}

void SteinerTree::computeRatios() {
  this->ratio   = this->smt_length / this->getMSTLength();
  this->b_ratio = this->smt_length / this->getBMSTLength();
}

double SteinerTree::getSteinerRatio() const {
  return this->ratio;
}

double SteinerTree::getBSteinerRatio() const {
  return this->b_ratio;
}

void SteinerTree::setSteinerRatio(double l) {
  this->ratio = l;
}

unsigned int SteinerTree::m() const {
  return Graph::n() + this->steiner_points.size();
}

unsigned int SteinerTree::s() const {
  return this->steiner_points.size();
}

std::vector<Point> &SteinerTree::getSteinerPoints() {
  return this->steiner_points;
}

Point &SteinerTree::getPoint(unsigned int index) {
  assert(index < this->m());
  if(index >= this->n())
    return this->getSteinerPoint(index - this->n());
  else
    return Graph::getPoint(index);
}

Point &SteinerTree::getSteinerPoint(unsigned int index) {
  assert(index < this->s());
  return this->steiner_points[index];
}

double SteinerTree::getLength(bool use_euclidean) {
  std::vector<Edge>::const_iterator it;
  double result = 0.0;
  for(it = this->edges.begin(); it != this->edges.end(); it++) {
    result += (use_euclidean
	       ? Utils::length(this->getPoint(it->i0), this->getPoint(it->i1))
	       : it->length);
  }
  return result;
}


std::ostream& operator<<(std::ostream& os, SteinerTree &st) {
  unsigned int i;
  std::string indent("  ");
  std::vector<Utils::Edge> &edges = st.getEdges();
  
  os << "# Terminals" << std::endl;
  for(i = 0; i < st.n(); i++)
    os << indent << i << " " << st.getPoint(i) << std::endl;
  
  os << std::endl << "# Steiner points" << std::endl;
  for(; i < st.m(); i++)
    os << indent << i << " " << st.getPoint(i) << std::endl;
  
  os << std::endl << "# Edges" << std::endl;
  for(i = 0; i < edges.size(); i++)
    os << indent << "(" << edges[i].i0
       << " " << edges[i].i1 << ")" << std::endl;
  
  os << std::endl << "# |MST|: " << st.getMSTLength() << std::endl
     << "# |SMT|: " << st.getSMTLength() << std::endl
     << "# Steiner ratio: " << st.getSteinerRatio() << std::endl;
  return os;
}

SteinerTree &SteinerTree::operator=(const SteinerTree &other) {
  if(this != &other) {
    Graph::operator=(other);
    this->steiner_points = other.steiner_points;
    this->smt_length     = other.smt_length;
    this->ratio          = other.ratio;
    this->b_ratio        = other.b_ratio;
  }
  return *this;
}
