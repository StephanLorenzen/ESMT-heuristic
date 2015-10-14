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

// Static tree operations
BottleneckGraphSleator::Vertex *BottleneckGraphSleator::parent(Vertex *v) {
  return v == this->tail(this->path(v)) ? v->dparent : this->after(v);
}
BottleneckGraphSleator::Vertex *BottleneckGraphSleator::root(Vertex *v) {
  return this->tail(this->expose(v));
}
double BottleneckGraphSleator::cost(Vertex *v) {
  return v == this->tail(this->path(v)) ? v->dcost : this->pcost(v);
}
double BottleneckGraphSleator::mincost(Vertex *v) {
  return this->pmincost(this->expose(v));
}
void BottleneckGraphSleator::update(Vertex *v, double x) {
  this->pupdate(this->expose(v), x);
}
// Dynamic tree operations
void BottleneckGraphSleator::link(Vertex *v, Vertex *w, double x) {
  this->concatenate(this->path(v), this->expose(w), x);
}
double BottleneckGraphSleator::cut(Vertex *v) {
  this->expose(v);
  SplitResult s;
  this->split(v, s);
  v->dparent = NULL;
  return s.y;
}
void BottleneckGraphSleator::evert(Vertex *v) {
  this->reverse(this->expose(v));
  v->dparent = NULL;
}

// Path operations
BottleneckGraphSleator::Path *BottleneckGraphSleator::path(Vertex *v) {
  while(v->bparent)
    v = v->bparent;
  return v;
}
BottleneckGraphSleator::Vertex *BottleneckGraphSleator::head(Path *p) {
  return p->reversed ? p->btail : p->bhead;
}
BottleneckGraphSleator::Vertex *BottleneckGraphSleator::tail(Path *p) {
  return p->reversed ? p->bhead : p->btail;
}
BottleneckGraphSleator::Vertex *BottleneckGraphSleator::before(Vertex *v) {
  RRS c;
  c.reversed = false;
  this->before_rec(v, c);
  Vertex *w = c.res;
  Vertex *u;
  u = w->bparent->bleft == w ? w->bparent->bright : w->bparent->bleft;
  bool rev = (c.reversed != u->reversed); 
  if(u->external)
    return u;
  else
    return rev ? u->bhead : u->btail;
}
void BottleneckGraphSleator::before_rec(Vertex *v, RRS &c) {
  if(!v)
    return;
  Vertex *par = v->bparent;
  this->before_rec(par, c);
  if(par) {
    if(c.reversed) {
      if(v == par->bleft) {
	c.res = v;
      }
    }
    else if(v == par->bright) {
      c.res = v;
    }
  }
  c.reversed = (c.reversed != v->reversed);
}
//bool before_rec(const unsigned int i, int &r, bool &st);
BottleneckGraphSleator::Vertex *BottleneckGraphSleator::after(Vertex *v) {
  return NULL;
}
//bool after_rec(const unsigned int i, int &r, bool &st);
double BottleneckGraphSleator::pcost(Vertex *v) {
  return 0.0;
}
//bool pcost_rec(unsigned int i, int &r, bool &st, double &g, double &grossmin);
double BottleneckGraphSleator::pmincost(Path *p) {
  return 0.0;
}
//double pmincost_rec(unsigned int u, bool &r);
void BottleneckGraphSleator::pupdate(Path *p, double x) {

}
void BottleneckGraphSleator::reverse(Path *p) {

}
BottleneckGraphSleator::Path *BottleneckGraphSleator::concatenate(Path *p, Path *q, double x) {
  return NULL;
}
void BottleneckGraphSleator::split(Path *p, SplitResult &res) {

}
  
// Splice and expose
BottleneckGraphSleator::Path *BottleneckGraphSleator::splice(Path *p) {
  Vertex *v = this->tail(p)->dparent;
  SplitResult s;
  this->split(v, s);
  if(s.q) {
    this->tail(s.q)->dparent = v;
    this->tail(s.q)->dcost   = s.x;
  }
  p = this->concatenate(p, this->path(v), this->tail(p)->dcost);
  return s.r ? this->concatenate(p, s.r, s.y) : p;
}
BottleneckGraphSleator::Path *BottleneckGraphSleator::expose(Vertex *v) {
  Path *p;
  SplitResult s;
  this->split(v, s);
  if(s.q) {
    this->tail(s.q)->dparent = v;
    this->tail(s.q)->dcost   = s.x;
  }
  p = s.r ? this->concatenate(this->path(v), s.r, s.y) : this->path(v);
  while(this->tail(p)->dparent) p = this->splice(p);
  return p;
}
  
