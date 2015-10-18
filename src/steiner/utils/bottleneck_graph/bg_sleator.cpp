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
    this->_vertices[i].idx      = i;
    this->_vertices[i].cost     = 0;
    this->_vertices[i].maxcost  = 0;
  }
  Utils::MSTKruskalMod(g, true);
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
    path = this->concatenate(&this->_vertices[cur], path, Utils::length(this->_pointsRef[cur], this->_pointsRef[prev]));
  else {
    path = &this->_vertices[cur];
    if(cur != 0) {
      path->dparent = &this->_vertices[prev];
      path->dcost   = Utils::length(this->_pointsRef[cur], this->_pointsRef[prev]);
    }
  }
  for(unsigned int i = 0; i < conns[cur].size(); i++) {
    if(conns[cur][i] != prev) {
      this->pathDecompose(conns, conns[cur][i], cur, path);
      path = NULL;
    }
  }
}

/*
 * Implementation of BottleneckGraphSleator destructor
 */
BottleneckGraphSleator::~BottleneckGraphSleator() {
  // Assume all connected
  for(unsigned int i = 0; i < this->_vertices.size(); i++)
    if(this->_vertices[i].bparent)
      this->_cleanUp(this->path(&this->_vertices[i]));
  
}

void BottleneckGraphSleator::_cleanUp(Vertex *v) {
  if(v->external)
    v->bparent = NULL;
  else {
    this->_cleanUp(v->bleft);
    this->_cleanUp(v->bright);
    delete v;
  }
}

/**
 * Implementation of BottleneckGraphSleator::distance(...)
 */
double BottleneckGraphSleator::distance(const unsigned int i, const unsigned int j) {
  if(i==j)
    return 0.0;
  this->evert(&this->_vertices[i]);
  Vertex *v = this->maxcost(&this->_vertices[j]);
  return v ? this->cost(v) : -1.0;
}

/**
 * Implementation of BottleneckGraphSleator::mergePoints(...)
 */
void BottleneckGraphSleator::mergePoints(const std::vector<unsigned int> &points) {
  if(points.size() < 2)
    return;
  
  // Remove the bottleneck edges
  for(unsigned int i=0; i<this->_vertices.size(); i++) {
    this->evert(&this->_vertices[i]);
    for(unsigned int j=i+1; j<this->_vertices.size(); j++) {
      this->cut(&this->_vertices[j]);
    }
  }
  
  for(unsigned int i=1; i<this->_vertices.size(); i++)
    this->link(&this->_vertices[i-1], &this->_vertices[i], 0.0);
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
BottleneckGraphSleator::Vertex *BottleneckGraphSleator::maxcost(Vertex *v) {
  return this->pmaxcost(this->expose(v));
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
  c.res      = NULL;
  this->before_rec(v, c);
  Vertex *w = c.res;
  if(!w)
    return NULL;
  Vertex *u = w->bparent->bleft == w ? w->bparent->bright : w->bparent->bleft;
  bool rev = (c.resrev != u->reversed); 
  if(u->external)
    return u;
  else
    return rev ? u->bhead : u->btail;
}
void BottleneckGraphSleator::before_rec(Vertex *v, RRS &c) {
  Vertex *par = v->bparent;
  if(!par)
    // Root
    return;
  this->before_rec(par, c);
  c.reversed = (c.reversed != par->reversed);
  if(c.reversed) {
    if(v == par->bleft) {
      c.res = v;
      c.resrev = c.reversed;
    }
  }
  else if(v == par->bright) {
    c.res = v;
    c.resrev = c.reversed;
  }
}
BottleneckGraphSleator::Vertex *BottleneckGraphSleator::after(Vertex *v) {
  RRS c;
  c.reversed = false;
  c.res = NULL;
  this->after_rec(v, c);
  Vertex *w = c.res;
  if(!w)
    return NULL;
  Vertex *u = w->bparent->bleft == w ? w->bparent->bright : w->bparent->bleft;
  bool rev = (c.resrev != u->reversed); 
  if(u->external)
    return u;
  else
    return rev ? u->btail : u->bhead;
}
void BottleneckGraphSleator::after_rec(Vertex *v, RRS &c) {
  Vertex *par = v->bparent;
  if(!par)
    // Root
    return;
  this->after_rec(par, c);
  c.reversed = (c.reversed != par->reversed);
  if(c.reversed) {
    if(v == par->bright) {
      c.res = v;
      c.resrev = c.reversed;
    }
  }
  else if(v == par->bleft) {
    c.res = v;
    c.resrev = c.reversed;
  }
}
double BottleneckGraphSleator::pcost(Vertex *v) {
  RRS c;
  c.res = NULL;
  c.reversed = false;
  this->after_rec(v, c);
  Vertex *w = c.res->bparent;
  return w->cost;
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
    if(u->external || (r->external && l->external) || abs(u->maxcost-u->cost) <= 0.0001)
      // Edge has been found
      break;
    else
      u = (r->maxcost > l->maxcost && !r->external ? r : l);
  }
  if(u->external)
    return NULL;
  // Return rightmost external descendant of l
  if(l->external) return l;
  else if(rev) return l->bhead;
  else return l->btail;
}
void BottleneckGraphSleator::reverse(Path *p) {
  p->reversed = !p->reversed;
}
BottleneckGraphSleator::Path *BottleneckGraphSleator::concatenate(Path *p, Path *q, double x) {
  if(p->dparent == q)
    p->dparent = NULL;
  if(q->dparent == p)
    q->dparent = NULL;
  return this->construct(p, q, x);
}
void BottleneckGraphSleator::treepath(Vertex *v, Vertex **r, bool &reversed, bool before) {
  if(!v) {
    r = NULL;
    return;
  }
  Vertex *par = v->bparent;
  this->treepath(par, r, reversed, before);
  if(par) {
    if(reversed) {
      if((before && v == par->bleft) || (!before && v == par->bright)) {
	*r = par;
      }
    }
    else if((!before && v == par->bleft) || (before && v == par->bright)) {
      *r = par;
    }
  }
  reversed = (reversed != v->reversed);
}
void BottleneckGraphSleator::tsplit(Vertex *v, bool right, DestructResult &dr) {
  Vertex *par = v->bparent;
  bool pright = (par && par->bright == v);
  if(!dr.left)
    this->destruct(v, dr);
  else {
    if(right)
      dr.left = this->construct(v->bleft, dr.left, v->cost);
    else
      dr.right = this->construct(dr.right, v->bright, v->cost);
    dr.left->reversed  = (dr.left->reversed != v->reversed);
    dr.right->reversed = (dr.right->reversed != v->reversed);
    if(v->reversed) {
      Vertex *tmp = dr.left;
      dr.left     = dr.right;
      dr.right    = tmp;
    }
    delete v;
  }
  if(par)
    this->tsplit(par, pright, dr);
}
void BottleneckGraphSleator::split(Vertex *v, SplitResult &res) {
  Vertex *r = NULL;
  bool rev = false;
  DestructResult dr1, dr2;
  dr1.left = NULL;
  dr1.right = NULL;
  dr2.left = NULL;
  dr2.right = NULL;
  this->treepath(v, &r, rev, true);
  if(r)
    this->tsplit(r, false, dr1);
  
  r = NULL;
  rev = false;
  this->treepath(v, &r, rev, false);
  if(r)
    this->tsplit(r, false, dr2);
  
  res.q = dr1.left;
  res.x = dr1.x;
  res.r = dr2.right;
  res.y = dr2.x;
}
  
// Splice and expose
BottleneckGraphSleator::Path *BottleneckGraphSleator::splice(Path *p) {
  Vertex *v = this->tail(p)->dparent;
  assert(v);
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
  r->cost     = x;
  r->maxcost  = v->maxcost > w->maxcost ? v->maxcost : w->maxcost;
  r->maxcost  = r->maxcost > x ? r->maxcost : x;
  r->idx      = this->head(v)->idx*10+this->tail(w)->idx;
  
  r->dparent  = NULL;
  r->dcost    = 0.0;

  r->height   = (v->height > w->height ? v->height : w->height) + 1;
  this->balance(r);
  return r;
}
void BottleneckGraphSleator::destruct(Vertex *u, DestructResult &c) {
  c.left = u->bleft;
  c.left->bparent = NULL;
  c.left->reversed = (c.left->reversed != u->reversed);
  c.right = u->bright;
  c.right->bparent = NULL;
  c.right->reversed = (c.right->reversed != u->reversed);
  c.x = u->cost;
  if(u->reversed) {
    Vertex *tmp = c.left;
    c.left = c.right;
    c.right = tmp;
  }
  delete u;
  // Make sure subtrees are balanced (they should always be)
  this->balance(c.left);
  this->balance(c.right);
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
