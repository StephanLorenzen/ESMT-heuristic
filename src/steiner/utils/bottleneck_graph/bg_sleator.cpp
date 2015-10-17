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
BottleneckGraphSleator::BottleneckGraphSleator(Graph &g)
  : BottleneckGraph(g), _vertices(g.n()), _pointsRef(g.getPointsRef()) {
  for(unsigned int i = 0; i < g.n(); i++) {
    this->_vertices[i].bparent  = NULL;
    this->_vertices[i].bleft    = NULL;
    this->_vertices[i].bright   = NULL;
    this->_vertices[i].bhead    = &this->_vertices[i];
    this->_vertices[i].btail    = &this->_vertices[i];
    this->_vertices[i].external = true;
    this->_vertices[i].reversed = false;
    this->_vertices[i].height   = 1;
  }
  Utils::MSTKruskalMod(g);
  // Traverse tree
  std::vector< std::vector<unsigned int> > conns(g.n());
  std::vector<Edge> &edges = g.getEdges();
  for(unsigned int i = 0; i < edges.size(); i++) {
    Edge e = edges[i];
    conns[e.i0].push_back(e.i1);
    conns[e.i1].push_back(e.i0);
  }
  this->pathDecompose(conns, 0, 0, NULL);
}
void BottleneckGraphSleator::pathDecompose(std::vector< std::vector<unsigned int> > conns,
					   unsigned int cur, unsigned int prev, Path *path) {
  if(path)
    path = this->construct(path, &this->_vertices[cur], Utils::length(this->_pointsRef[cur], this->_pointsRef[prev]));
  else
    path = &this->_vertices[cur];
  for(unsigned int i = 0; i < conns[cur].size(); i++) {
    if(conns[cur][i] != prev) {
      //std::cout << i << std::endl;
      this->pathDecompose(conns, conns[cur][i], cur, path);
      path = NULL;
    }
  }
}

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
/*BottleneckGraphSleator::Vertex *BottleneckGraphSleator::mincost(Vertex *v) {
  return this->pmincost(this->expose(v));
  }*/
BottleneckGraphSleator::Vertex *BottleneckGraphSleator::maxcost(Vertex *v) {
  return this->pmaxcost(this->expose(v));
}
/*void BottleneckGraphSleator::update(Vertex *v, double x) {
  this->pupdate(this->expose(v), x);
  }*/
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
  Vertex *u = w->bparent->bleft == w ? w->bparent->bright : w->bparent->bleft;
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
BottleneckGraphSleator::Vertex *BottleneckGraphSleator::after(Vertex *v) {
  RRS c;
  c.reversed = false;
  this->after_rec(v, c);
  Vertex *w = c.res;
  Vertex *u = w->bparent->bleft == w ? w->bparent->bright : w->bparent->bleft;
  bool rev = (c.reversed != u->reversed); 
  if(u->external)
    return u;
  else
    return rev ? u->bhead : u->btail;
}
void BottleneckGraphSleator::after_rec(Vertex *v, RRS &c) {
  if(!v)
    return;
  Vertex *par = v->bparent;
  this->before_rec(par, c);
  if(par) {
    if(c.reversed) {
      if(v == par->bright) {
	c.res = v;
      }
    }
    else if(v == par->bleft) {
      c.res = v;
    }
  }
  c.reversed = (c.reversed != v->reversed);
}
/*double BottleneckGraphSleator::pcost(Vertex *v) {
  RRS c;
  this->pcost_rec(v, c);
  Vertex *w = c.res;
  return w->bparent->netcost + (c.grossmin-w->netmin);
}
void BottleneckGraphSleator::pcost_rec(Vertex *v, RRS &c) {
  if(!v)
    return;
  Vertex *par = v->bparent;
  this->pcost_rec(par, c);
  if(par) {
    if(c.reversed) {
      if(v == par->bright) {
	c.res = v;
      }
    }
    else if(v == par->bleft) {
      c.res = v;
    }
  }
  c.grossmin += v->netmin;
  c.reversed = (c.reversed != v->reversed);
}*/
double BottleneckGraphSleator::pcost(Vertex *v) {
  RRS c;
  this->pcost_rec(v, c);
  Vertex *w = c.res;
  return w->bparent->cost;
}
void BottleneckGraphSleator::pcost_rec(Vertex *v, RRS &c) {
  if(!v)
    return;
  Vertex *par = v->bparent;
  this->pcost_rec(par, c);
  if(par) {
    if(c.reversed) {
      if(v == par->bright) {
	c.res = v;
      }
    }
    else if(v == par->bleft) {
      c.res = v;
    }
  }
  c.reversed = (c.reversed != v->reversed);
}
BottleneckGraphSleator::Vertex *BottleneckGraphSleator::pmaxcost(Path *p) {
  bool rev = false;
  Vertex *u = p;
  Vertex *r;
  Vertex *l;
  while(true) {
    rev = (rev != u->reversed);
    if(rev) { r = u->bleft;  l = u->bright; }
    else    { r = u->bright; l = u->bleft;  }
    if(!r->external && r->maxcost-r->cost == 0)
      u = r;
    else if(r->maxcost-u->cost > 0)
      u = l;
    else
      break;
  }
  // Return rightmost external descendant of l
  if(l->external) return l;
  else if(rev) return l->bhead;
  else return l->btail;
}
/*BottleneckGraphSleator::Vertex *BottleneckGraphSleator::pmincost(Path *p) {
  bool rev = false;
  Vertex *u = p;
  Vertex *r;
  Vertex *l;
  while(true) {
    rev = (rev != u->reversed);
    if(rev) { r = u->bleft;  l = u->bright; }
    else    { r = u->bright; l = u->bleft;  }
    if(!r->external && r->netcost == 0)
      u = r;
    else if(u->netcost > 0)
      u = l;
    else
      break;
  }
  // Return rightmost external descendant of l
  if(l->external) return l;
  else if(rev) return l->bhead;
  else return l->btail;
  }*/
/*void BottleneckGraphSleator::pupdate(Path *p, double x) {
  p->netmin += x;
  }*/
void BottleneckGraphSleator::reverse(Path *p) {
  p->reversed = !p->reversed;
}
BottleneckGraphSleator::Path *BottleneckGraphSleator::concatenate(Path *p, Path *q, double x) {
  return this->construct(this->tail(p), this->head(q), x);
}

void BottleneckGraphSleator::treepath(Vertex *v, std::vector<bool> &path, RRS &c, bool before) {
  if(!v)
    return;
  Vertex *par = v->bparent;
  this->treepath(par, path, c, before);
  if(par) {
    if(c.reversed) {
      if((before && v == par->bleft) || (!before && v == par->bright)) {
	c.res = v;
      }
    }
    else if((!before && v == par->bleft) || (before && v == par->bright)) {
      c.res = v;
    }
    path.push_back(v==par->bright);
  }
  else
    c.path = v;
  c.reversed = (c.reversed != v->reversed);  
}
void BottleneckGraphSleator::tsplit(Vertex *r, Vertex *v, std::vector<bool> &path, int i, DestructResult &dr) {
  assert(r);
  if(r==v) {
    this->destruct(v, dr);
  }
  else if(path[i]) {
    // Right
    this->tsplit(r->bright, v, path, i+1, dr);
    dr.v = this->construct(r->bleft, dr.v, r->cost);
    delete r;
  }
  else {
    // Left
    this->tsplit(r->bleft, v, path, i+1, dr);
    dr.w = this->construct(dr.w, r->bright, r->cost);
    delete r;
  }
}
void BottleneckGraphSleator::split(Vertex *v, SplitResult &res) {
  std::vector<bool> path;
  RRS c;
  DestructResult dr1, dr2;
  this->treepath(v, path, c, true);
  this->tsplit(c.path, c.res, path, 0, dr1);
  this->treepath(v, path, c, false);
  this->tsplit(c.path, c.res, path, 0, dr2);
  
  res.q = dr1.v;
  res.x = dr1.x;
  res.r = dr2.v == v ? dr2.w : dr2.v;
  res.y = dr2.x;
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
  
// AVL tree operations
BottleneckGraphSleator::Vertex *BottleneckGraphSleator::construct(Vertex *v, Vertex *w, double x) {
  Vertex *r = new Vertex;
  r->bparent  = NULL;
  r->bleft    = v;
  v->bparent  = r;
  r->bright   = w;
  w->bparent  = r;
  r->bhead    = this->head(v);
  r->btail    = this->tail(w);
  r->external = false;
  r->reversed = false;
  //double gmin = min(v->netmin, w->netmin);
  //r->netmin   = gmin;
  //r->netcost  = x - gmin;
  r->cost     = x;
  r->maxcost  = v->maxcost > w->maxcost ? v->maxcost : w->maxcost;;
  
  r->dparent  = NULL;
  r->dcost    = 0.0;

  r->height   = (v->height > w->height ? v->height : w->height) + 1;

  this->balance(r);
  return r;
}
void BottleneckGraphSleator::destruct(Vertex *u, DestructResult &c) {
  c.v = u->bleft;
  c.w = u->bright;
  c.x = u->cost;
  delete u;
  // Make sure subtrees are balanced (they should always be)
  this->balance(c.v);
  this->balance(c.w);
}
void BottleneckGraphSleator::rotateleft(Vertex *v) {
  Vertex *rc = v->bright;
  if(!rc)
    return;
  Vertex *p = v->bleft;
  Vertex *q = rc->bleft;
  Vertex *r = rc->bright;
  // Update v
  v->bleft   = rc;
  v->bright  = r;
  // Update p
  if(p) p->bparent = rc;
  // Update r
  if(r) r->bparent = v;
  // Update rc
  rc->bleft  = p;
  rc->bright = q;

  // Update heights
  rc->height = (p->height > q->height ? p->height : q->height)+1;
  v->height  = (rc->height > r->height ? rc->height : r->height)+1;
}
void BottleneckGraphSleator::rotateright(Vertex *v) {
  Vertex *lc = v->bleft;
  if(!lc)
    return;
  Vertex *p = lc->bleft;
  Vertex *q = lc->bright;
  Vertex *r = v->bright;
  // Update v
  v->bleft   = p;
  v->bright  = lc;
  // Update p
  if(p) p->bparent = v;
  // Update q
  if(r) r->bparent = lc;
  // Update rc
  lc->bleft  = q;
  lc->bright = r;

  // Update heights
  lc->height = (q->height > r->height ? q->height : r->height)+1;
  v->height  = (lc->height > p->height ? lc->height : p->height)+1;
}
void BottleneckGraphSleator::balance(Vertex *v) {
  int lh = v->bleft ? v->bleft->height : 0;
  int rh = v->bright ? v->bright->height : 0;
  int balance = lh-rh;
  if(balance <= 1 && balance >= -1)
    return;
  // Do rotate
  //      v
  //     / \        .
  //    /   \       .
  //   l     r
  //  / \   / \     .
  // a   b c   d
  if(lh > rh) {
    int la = v->bleft->bleft ? v->bleft->bleft->height : 0;
    int lb = v->bleft->bright ? v->bleft->bright->height : 0;
    balance = la-lb;
    if(balance == -1) {
      // Must rotate b up
      this->rotateleft(v->bleft);
    }
    this->rotateright(v);
    // Balance new b-right
    this->balance(v->bright);
  }
  else { // lh < rh
    int lc = v->bright->bleft ? v->bright->bleft->height : 0;
    int ld = v->bright->bright ? v->bright->bright->height : 0;
    balance = lc-ld;
    if(balance == 1) {
      // Must rotate b up
      this->rotateright(v->bright);
    }
    this->rotateleft(v);
    // Balance new b-right
    this->balance(v->bleft);
  }
}
