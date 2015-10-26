#include <vector>
#include <iostream>
#include <assert.h>
#include <cstdlib>
#include <cmath>

#include "steiner/graph.hpp"
#include "steiner/steiner_tree.hpp"
#include "steiner/heuristics/concat.hpp"
#include "steiner/utils/point.hpp"
#include "steiner/utils/utils.hpp"

typedef Utils::Point Point;
typedef Utils::Edge  Edge;

/**
 * Constructor
 */
IterativeConcat::IterativeConcat(int dim, bool doCleanUp) : SubgraphHeuristic(), Iterative(dim, CNMAX), do_clean_up(doCleanUp) {}

/**
 * Destructor
 */
IterativeConcat::~IterativeConcat() {}

/*
 * Finds the steiner points and edges.
 */
void IterativeConcat::findSteinerPoints(SteinerTree &subgraph) {
  unsigned int i, j, n;
  
  n = subgraph.n();
  std::vector<Edge>  &edges  = subgraph.getEdges();
  edges.clear();
  
  this->dim = subgraph.dimension();
  
  // Check that the input is valid
  assert(this->dim < this->max_dim);
  assert(n > 2);
  assert(n < CNMAX);
  
  this->N = n;
  
  // Copy points (+ SP dummies)
  this->P = std::vector<Point>();
  this->P.reserve(2*this->N-2);
  for(i = 0; i < n; i++)
    this->P.push_back(subgraph.getPoint(i));
  for(i = 0; i < this->N-2; i++)
    this->P.push_back(Point(this->dim));
  
  double l, r;

  // Init topology vector
  for(i = 0; i < this->N - 3; i++) {
    this->topo_vec[i][0] = 0;
    this->topo_vec[i][1] = 0;
  }

  // Build tree with three terminals
  this->buildTreeConcat(0);
  this->S = 1;
  
  if(this->N == 3) {
    // Optimise and return
    l = this->length();
    r = this->error(); 
    do {
      this->optimise(this->treshold*r/this->N);
      l = this->length();
      r = this->error();
    } while(r>l*this->treshold);
  }
  else {
    for(i = 3; i < this->N; i++) {
      // Add points iteratively
      int sp      = this->N;
      double dist = this->dist(i, sp);
      double tmp;
    
      for(j = this->N+1; j < this->N+this->S; j++) {
	tmp = this->dist(i, j);
	if(tmp < dist) {
	  dist = tmp;
	  sp = j;
	}
      }
      // sp is now closest SP.
      this->topo_vec[i-3][0] = sp-this->N;
    
      double best_l       = -1;
      int    best_e       = 0;
      double low_treshold = 50*this->treshold;
      for(j = 0; j < 3; j++) {
	this->topo_vec[i-3][1] = j;
	this->buildTreeConcat(i-2);
	
	this->S++;
	
	// Optimise
	l = this->length();
	r = this->error();  
	do {
	  this->optimise(low_treshold*r/this->N);
	  l = this->length();
	  r = this->error();
	} while(r>l*low_treshold);
      
	this->S--;
	
	if(best_l < 0 || l < best_l) {
	  best_l = l;
	  best_e = j;
	}
      }
      this->topo_vec[i-3][1] = best_e;
      this->buildTreeConcat(i-2);
      this->S++;
    
      // Optimise
      l = this->length();
      r = this->error();  
      do {
	this->optimise(this->treshold*r/this->N);
	l = this->length();
	r = this->error();
      } while(r>l*this->treshold);
      // Done
    }
  }
  this->cleanUp(subgraph.getSteinerPoints(), edges, this->do_clean_up);
  
  // Update length
  subgraph.setSMTLength(subgraph.getLength());
}

void IterativeConcat::insertTerminal(SteinerTree &fst, unsigned int pi, double mst_length) {
  unsigned int i, j, m, n, s, i0, i1, is0, is1;
  double tmp, l = -1.0;

  this->dim = fst.dimension();
  
  n = fst.n();
  s = fst.s();
  m = n+s;
  std::vector<Edge> &edges = fst.getEdges();
  
  // Allocate space for new points
  this->P = std::vector<Point>();
  this->P.reserve(m+2);
  for(i = 0; i < n; i++)
    this->P.push_back(fst.getPoint(i));
  
  // Check for full tree
  assert(m == n*2-2);
  this->N = i+1;
  this->S = i-2;
  // New terminal - add to fst and this->P
  fst.getPoints().push_back(pi);
  Point &p = fst.getPoint(fst.n()-1);
  this->P.push_back(p);
  
  // Current Steiner points
  for(; i < m; i++)
    this->P.push_back(fst.getPoint(i));
  
  // Add dummy for last
  this->P.push_back(Point(this->dim));
  
  for(i = 0; i < this->S; i++)
    for(j = 0; j < 3; j++)
      this->adj[i][j] = -1;
  
  for(i = 0; i < edges.size(); i++) {
    i0 = edges[i].i0;
    i1 = edges[i].i1;
    // Increment sp indices
    i0 = i0 >= this->N-1 ? i0+1 : i0;
    i1 = i1 >= this->N-1 ? i1+1 : i1;
    if(i0 >= this->N) {
      is0 = i0 - this->N;
      for(j = 0; j < 3; j++)
	if(this->adj[is0][j] < 0) {
	  this->adj[is0][j] = i1;
	  break;
	}
    }
    if(i1 >= this->N) {
      is1 = i1 - this->N;
      for(j = 0; j < 3; j++)
	if(this->adj[is1][j] < 0) {
	  this->adj[is1][j] = i0;
	  break;
	}
    }
  }
  // Adj set up.
  // Search for closest sp.
  unsigned int sp = 0;
  for(i = 0; i < this->S; i++) {
    tmp = this->dist(i+this->N, this->N-1);
    if(l < 0.0 || tmp < l) {
      l = tmp;
      sp = i;
    }
  }
  // sp is now closest. Try to insert
  l = -1.0;
  n = 0;
  for(i = 0; i < 3; i++) {
    tmp = this->insertAndRemovePoint(this->N-1, sp, i);
    if(l < 0 || tmp < l) {
      l = tmp;
      n = i;
    }
  }
  this->insertPoint(this->N-1, sp, n);
  this->cleanUp(fst.getSteinerPoints(), edges, this->do_clean_up);
  
  // Compute ratio
  fst.setMSTLength(mst_length);
  fst.setSMTLength(fst.getLength());
}

double IterativeConcat::insertPoint(unsigned int i,
				    unsigned int sp,
				    unsigned int j) {
  unsigned int d, k, s, op, opi;
  double l, r;
  s = sp+this->N;
  op = this->adj[sp][j];
  for(d = 0; d < this->dim; d++) {
    this->P[this->N + this->S][d] = (this->P[i][d]
				     + this->P[s][d]
				     + this->P[op][d]) / 3.0 + 0.001 * Utils::frand();
  }
  this->adj[this->S][0] = i;
  this->adj[this->S][1] = s;
  this->adj[this->S][2] = op;
  // Update sp
  this->adj[sp][j] = this->N+this->S;
  // Update op if needed
  if(op >= this->N) {
    opi = op-this->N;
    for(k = 0; k < 3; k++) {
      if(this->adj[opi][k] == (int)s)
	break;
    }
    this->adj[opi][k] = this->N+this->S;
  }
  this->S++;
  
  // Optimise
  l = this->length();
  r = this->error();
  do {
    this->optimise(this->treshold*r/this->N);
    l = this->length();
    r = this->error();
  } while(r>l*this->treshold);
  return l;
}

double IterativeConcat::insertAndRemovePoint(unsigned int i,
					     unsigned int sp,
					     unsigned int j) {
  unsigned int k, prev, previ = 0;
  prev = this->adj[sp][j];
  if(prev >= this->N) {
    for(k = 0; k < 3; k++)
      if(this->adj[prev-this->N][k] == (int)(sp+this->N))
	break;
    previ = k;
  }
  double result = this->insertPoint(i, sp, j);
  this->S--;
  for(k = 0; k < 3; k++)
    this->adj[this->S][k] = -1;
  this->adj[sp][j] = prev;
  if(prev >= this->N)
    this->adj[prev-this->N][previ] = sp+this->N;
  
  return result;
}

void IterativeConcat::buildTreeConcat(unsigned int k) {
  unsigned int i, g, d;
  assert(this->N > 2);

  // First, build a 3-terminal tree.
  int si = this->N; // Steiner index
  
  for(d = 0; d < this->dim; d++) {
    this->P[si][d] = (this->P[0][d]
		      + this->P[1][d]
		      + this->P[2][d]) / 3.0 + 0.001 * Utils::frand();
  }
  
  // Set adjs
  this->adj[0][0] = 0;
  this->adj[0][1] = 1;
  this->adj[0][2] = 2;
  
  for(i = 0; i < k; i++) {
    // Best sp;
    int spi = this->topo_vec[i][0];
    int sp  = spi+this->N;
    // Best edge
    int e  = this->topo_vec[i][1];
    int nspi = i+1;
    int nsp  = nspi+this->N;
    
    // Split at SP sp, #edge e
    int point = this->adj[spi][e];
    // Set new adj
    this->adj[nspi][0] = sp;
    this->adj[nspi][1] = point;
    this->adj[nspi][2] = i+3;
    
    // And for the other
    this->adj[spi][e] = nsp;

    point -= this->N;
    if(point >= 0) {
      for(e = 0; e < 3; e++)
	if(this->adj[point][e] == sp)
	  break;
      assert(e < 3);
      this->adj[point][e] = nsp;
    }
    
    for(d = 0; d < this->dim; d++) {
      this->P[nsp][d] = 0.0;
      for(g = 0; g < 3; g++) {
	this->P[nsp][d] += this->P[this->adj[nspi][g]][d];
      }
      this->P[nsp][d] = this->P[nsp][d] / 3.0 + 0.001 * Utils::frand();
    }
  }
}
