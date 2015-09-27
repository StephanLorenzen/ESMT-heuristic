#ifndef GRAPH_H
#define GRAPH_H

#include <vector>

#include "steiner/utils/point.hpp"

typedef unsigned int   Pidx;
typedef Utils::Point   Point;
typedef Utils::Edge    Edge;

class Graph {
public:

  /**
   * Default constructor.
   */
  Graph(std::vector<Utils::Point> &pointsRef);
  
  /**
   * Constructor. Creates a graph with no edges.
   *
   * @param points   The points for this graph.
   */
  Graph(std::vector<Pidx> &points, std::vector<Utils::Point> &pointsRef);

  /**
   * Constructor.
   *
   * @param points   The points for this graph.
   * @param edges    The edges for this graph.
   */
  Graph(std::vector<Pidx> &points, std::vector<Utils::Point> &pointsRef,
	std::vector<Utils::Edge> &edges);
  
  /**
   * Destructor.
   */
  ~Graph();

  /**
   * Getter for the points.
   *
   * @return  A copy of the set of points for this graph.
   */
  std::vector<Pidx> &getPoints();
  std::vector<Utils::Point> &getPointsRef();
  
  /**
   * Getter for the edges.
   *
   * @return  A copy of the set of edges for this graph.
   */
  std::vector<Utils::Edge> &getEdges();
  
  Utils::Point &getPoint(unsigned int index) const;
  
  Pidx pidx(unsigned int index) const;

  unsigned int n() const;

  /**
   * Calculates the length of the MST for this graph.
   *
   * @return   The length of the MST for this graph.
   */
  double getMSTLength() const;

  /**
   * Calculates the length of the MST for this graph.
   *
   * @return   The length of the MST for this graph.
   */
  double getBMSTLength() const;
  
  /**
   * Setter for the MST length.
   *
   * Used to set the MST length to avoid overhead of
   * calculating it twice.
   *
   * @param l   The new MST length.
   */
  void setMSTLength(double l);

  void setBMSTLength(double l);

  void calculateMSTLength();
  void calculateBMSTLength();
  
  /**
   * Calculates the length of this graph.
   *
   * @return   The length of this graph.
   */
  double getLength() const;

  /**
   * Returns the dimension of this graph
   *
   * @return    The dimension of this graph
   */
  unsigned int dimension() const;

  Graph &operator=(const Graph&);
  
protected:
  
  /** The points of this graph. */
  std::vector<Utils::Point> *pointsRef;
  std::vector<Pidx> pidxs;
  /** The edges of this graph. */
  std::vector<Utils::Edge> edges;
  /** The length of the MST of this graph. */
  double mst_length;
  double bmst_length;
private:
};

/**
 * Prints the Graph
 *
 * @param os     The out-stream
 * @param st     The Graph to print.
 *
 * @return       A reference to the stream.
 */
std::ostream& operator<<(std::ostream& os, Graph &g);


#endif // GRAPH_H
