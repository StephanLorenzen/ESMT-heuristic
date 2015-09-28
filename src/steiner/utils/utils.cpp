#include <vector>
#include <list>
#include <algorithm>
#include <iostream>

#include "steiner/graph.hpp"
#include "steiner/utils/utils.hpp"
#include "steiner/utils/disjoint_set.hpp"
#include "steiner/utils/point.hpp"

typedef Utils::Edge         Edge;
typedef Utils::Point        Point;

/*
 * Compare function for the std::sort in Utils::MSTKruskal
 */
bool compareLength(const Edge &e1, const Edge &e2) {
  return e1.length < e2.length;
}

/*
 * Get edges for a full graph with n nodes
 */
void getEdges(Graph &graph) {
  std::vector<Edge> &edges = graph.getEdges();
  edges.clear();
  for(unsigned int i=0; i<graph.n(); i++)
    for(unsigned int j=i+1; j<graph.n(); j++)
      edges.push_back(Edge(i,j,Utils::length(graph.getPoint(i),graph.getPoint(j))));
}

/*
 * Implementation of Kruskal's MST algorithm
 *
 * Runs in O(n lg n)
 */
Graph Utils::MSTKruskal(Graph &graph) {
  Graph res(graph.getPoints(), graph.getPointsRef(), graph.getEdges());
  Utils::MSTKruskalMod(res);
  return res;
}

Graph Utils::MSTKruskal(std::vector<Point> &points) {
  std::vector<Edge> edges;
  for(unsigned int i=0; i<points.size(); i++)
    for(unsigned int j=i+1; j<points.size(); j++)
      edges.push_back(Edge(i,j,Utils::length(points[i],points[j])));
  Graph g(points);
  // Set the new edges
  g.getEdges() = edges;
  // Compute
  return Utils::MSTKruskal(g);
}

void Utils::MSTKruskalMod(Graph &graph, bool add_edges) {
  // Find resulting edges
  std::vector<Edge> res;
  if(add_edges)
    getEdges(graph);
  Utils::MSTKruskalEdges(graph, res);
  
  // Update edges
  graph.getEdges() = res;
}

void Utils::MSTKruskalEdges(Graph &graph, std::vector<Edge> &res) {
  // Iterator
  std::vector<Edge>::iterator eit;

  // Edges
  std::vector<Edge> &edges = graph.getEdges();
  
  // List of disjoint sets
  std::vector< Utils::DisjointSet<unsigned int> > sets;
  
  // Create the disjoint sets
  for(unsigned int i = 0; i < graph.n(); i++) {
    Utils::DisjointSet<unsigned int> ds(i);
    sets.push_back(ds);
  }
  // Sort the list of edges according to length
  std::sort(edges.begin(), edges.end(), compareLength);
  
  // Iterate through the edges
  for(eit = edges.begin(); eit != edges.end(); ++eit) {
    // Get the two sets from the list of sets
    Utils::DisjointSet<unsigned int> &ds1 = sets[eit->i0];
    Utils::DisjointSet<unsigned int> &ds2 = sets[eit->i1];

    // Check if we may safely concatenate
    if(ds1.findSet() != ds2.findSet()) {
      // Add edge to the resulting list and union
      res.push_back(*eit);
      ds1.setUnion(ds2);
    }
  }
}

double Utils::MSTLengthKruskal(std::vector<Point> &points) {
  Graph g(points);
  return MSTLengthKruskal(g, true);
}

double Utils::MSTLengthKruskal(Graph &graph, bool add_edges) {
  Graph g = graph;
  // Compute
  Utils::MSTKruskalMod(g, add_edges);
  return g.getLength();
}

double Utils::frand() {
  return (double)rand() / RAND_MAX;
}

void validateRec(std::vector<Edge> &edges, std::vector<int> &pf,
		 std::vector<int> &ef, int cur, int prev) {
  unsigned int i;
  pf[cur]++;
  for(i = 0; i < edges.size(); i++) {
    if(edges[i].i0 == cur && edges[i].i1 != prev) {
      ef[i]++;
      if(ef[i] > 1)
	return;
      validateRec(edges, pf, ef, edges[i].i1, edges[i].i0);
    }
    else if(edges[i].i1 == cur && edges[i].i0 != prev) {
      ef[i]++;
      if(ef[i] > 1)
	return;
      validateRec(edges, pf, ef, edges[i].i0, edges[i].i1);
    }
  }
}

bool Utils::validate(SteinerTree &st) {
  unsigned int i, p = st.n(), sp = st.s();
  std::vector<Edge>  &edges  = st.getEdges();
  // No of edges check
  if(st.m() != edges.size()+1) {
    std::cerr << "Wrong number of edges!!!" << std::endl;
    return false;
  }
  
  if(sp > p-2) {
    std::cerr << "Too many Steiner points: " << sp << " vs. " << p
	      << " terminals." << std::endl;
    return false;
  }
  
  std::vector<int> pf(st.m(), 0);
  std::vector<int> ef(edges.size(), 0);
  validateRec(edges, pf, ef, 0, -1);

  for(i = 0; i < pf.size(); i++) {
    if(pf[i] == 0) {
      std::cerr << "Point not reached: " << i << " " << st.getPoint(i) << std::endl;
      return false;
    }
    else if(pf[i] > 1) {
      std::cerr << "Point reached more than once: "
		<< i << " " << st.getPoint(i) << std::endl;
      return false;
    }
  }
  for(i = 0; i < ef.size(); i++) {
    if(!ef[i]) {
      std::cerr << "Edge not reached: " << i << " (" << edges[i].i0
		<< " : " << edges[i].i1 << ")" << std::endl;
      return false;
    }
    else if(ef[i] > 1) {
      std::cerr << "Edge reached more than once: " << i << " (" << edges[i].i0
		<< " : " << edges[i].i1 << ")" << std::endl;
      return false;
    }
  }

  std::vector<Point> points;
  for(i = 0; i < p; i++)
    points.push_back(st.getPoint(i));
  double mst_length = Utils::MSTLengthKruskal(points);
  if(st.getMSTLength()-mst_length>0.0001) {
    std::cerr << "MST-length is not correct!" << std::endl
	      << "  Correct |MST|: " << mst_length << std::endl
	      << "  Actual |MST|:  " << st.getMSTLength() << std::endl;
    return false;
  }
  
  return true;
}
