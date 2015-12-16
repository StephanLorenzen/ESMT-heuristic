#ifndef ITERATIVE_CONCAT_H
#define ITERATIVE_CONCAT_H

#include <vector>

#include "steiner/graph.hpp"
#include "steiner/heuristics/subgraph_heuristic.hpp"
#include "steiner/iterative.hpp"
#include "steiner/utils/point.hpp"

#define CNMAX 25

class IterativeConcat : public SubgraphHeuristic, public Iterative {
public:
  /**
   * Constructor
   */
  IterativeConcat(int dim, bool doCleanUp = false, bool doSortInput = false);
  
  /**
   * Destructor
   */
  ~IterativeConcat();

  /**
   * Finds the Steiner points and edges using a concatenation heuristic
   * for selecting the topology.
   */
  void findSteinerPoints(SteinerTree &subgraph);
  
  /**
   * Inserts a terminal in a full Steiner tree.
   *
   * @param fst         A pointer to the full Steiner tree.
   * @param pi          The index of the terminal to be inserted.
   * @param mst_length  The mst_length of the new tree.
   *
   * @return      The (allocated) resulting SteinerTree
   */
  void insertTerminal(SteinerTree &fst, unsigned int pi, double mst_length);

  /**
   * Setter for the doCleanUp flag
   *
   * @param doCleanUp  the new value of doCleanUp
   */
  void setDoCleanUp(bool doCleanUp);

  /**
   * Setter for the doSortInput flag
   *
   * @param doSortInput  the new value of doSortInput
   */
  void setDoSortInput(bool doSortInput);
  
protected:
private:
  void buildTreeConcat(unsigned int k);

  double insertPoint(unsigned int i,
		     unsigned int sp,
		     unsigned int j);
  double insertAndRemovePoint(unsigned int i,
			      unsigned int sp,
			      unsigned int j);
  
  int topo_vec[CNMAX][2];
  bool do_clean_up;
  bool do_sort_input;
};

#endif // ITERATIVE_CONCAT_H
