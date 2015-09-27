#ifndef STEINER_TREE_H
#define STEINER_TREE_H

#include <vector>
#include <assert.h>

#include "steiner/graph.hpp"
#include "steiner/utils/point.hpp"

/**
 * @class SteinerTree
 *
 * Implements a SteinerTree and related functions.
 */
class SteinerTree : public Graph {
public:

  /**
   * Default constructor.
   */
  SteinerTree(std::vector<Utils::Point> &pointsRef);

  /**
   * Constructor. Creates a SteinerTree with no edges
   *
   * @param terminals  The terminals for this SteinerTree.
   *                   This list contains the Steiner points too,
   *                   with the isSteiner flag set to true.
   */
  SteinerTree(std::vector<Pidx> &terminals, std::vector<Utils::Point> &pointsRef);
  
  /**
   * Constructor.
   *
   * @param terminals  The terminals for this SteinerTree.
   *                   This list contains the Steiner points too,
   *                   with the isSteiner flag set to true.
   * @param edges      The list of edges for this graph.
   */
  SteinerTree(std::vector<Pidx> &terminals,
	      std::vector<Utils::Point> &pointsRef,
	      std::vector<Utils::Edge> &edges);

  /**
   * Destructor.
   */
  virtual ~SteinerTree();

  /**
   * Gets the length of the SMT (= the length of this tree).
   *
   * @return  The length of this SMT.
   */
  double getSMTLength() const;

  /**
   * Setter for the SMT length.
   *
   * Used to set the SMT length to avoid overhead of
   * calculating it twice.
   *
   * @param l   The new SMT length.
   */
  void setSMTLength(double l);

  void calculateSMTLength();

  /**
   * Getter for the Steiner ratio.
   *
   * Computes |SMT|/|MST|.
   *
   * @return  The Steiner ratio.
   */
  double getSteinerRatio() const;
  double getBSteinerRatio() const;

  /**
   * Set the Steiner ratio for this ST.
   *
   * Used to avoid overhead by calculating the
   * ratio twice.
   *
   * @param l   The Steiner ratio of this ST.
   */
  void setSteinerRatio(double l);
  
  std::vector<Point> &getSteinerPoints();

  unsigned int m() const;
  unsigned int s() const;

  Utils::Point &getSteinerPoint(unsigned int index);
  Utils::Point &getPoint(unsigned int index);

  /**
   * Calculates the length of this graph.
   *
   * @return   The length of this graph.
   */
  double getLength();

  SteinerTree &operator=(const SteinerTree&);
  
protected:
  std::vector<Utils::Point> steiner_points;
private:
  /** The SMT length for this graph. */
  double smt_length;
  /** The Steiner ratio for this graph. */
  double ratio;
  double b_ratio;
};

/**
 * Prints the Steiner tree
 *
 * @param os     The out-stream
 * @param st     The Steiner tree to print.
 *
 * @return       A reference to the stream.
 */
std::ostream& operator<<(std::ostream& os, SteinerTree &st);

#endif // STEINER_TREE_H
