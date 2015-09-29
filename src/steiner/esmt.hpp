#ifndef ESMT_H
#define ESMT_H

#include <vector>
#include <unordered_map>

#include "steiner/heuristics/subgraph_heuristic.hpp"
#include "steiner/heuristics/concat.hpp"
#include "steiner/heuristics/smith.hpp"
#include "steiner/iterative.hpp"
#include "steiner/utils/point.hpp"
#include "steiner/utils/delaunay.hpp"
#include "steiner/utils/heap.hpp"
#include "steiner/utils/disjoint_set.hpp"
#include "steiner/utils/bottleneck_graph/bottleneck_graph.hpp"

// Bottleneck Graph - allowed values
#define BOTTLENECK_GRAPH_NONE   0
#define BOTTLENECK_GRAPH_NAIVE  1
#define BOTTLENECK_GRAPH_LAZY   2

/**
 * @class ESMT
 * @authors Andreas Olsen and Stephan Lorenzen
 *
 * A structure for describing an euclidean Steiner minimal tree in two or
 * three dimensions.
 */
class ESMT : public SteinerTree, public Iterative {
  
public:
  
  /**
   * Constructs an approximated eucliedean Steiner minimal tree, using
   * the algorithm described in our report.
   * Uses Smith's iterative algorithm for subgraphs and concates tetrahedrons.
   *
   * @param points          The given point set
   */
  ESMT(std::vector<Utils::Point> &points);

  /**
   * Constructs an approximated eucliedean Steiner minimal tree, using
   * the algorithm described in our report.
   *
   * @param points          The given point set
   * @param sh              A pointer to a subgraph heuristic.
   * @param concat_subgraphs   Determines wether subgraph concatenation should be
   *                           done.
   * @param post_optimise      Determines wether post optimisation using Smith's
   *                           should be applied
   * @param special_concat     Redo concatenation, adding non-covered faces.
   * @param verbose            If true, stats will be printed.
   */
  ESMT(std::vector<Utils::Point> &points, SubgraphHeuristic *sh,
       bool concat_subgraphs, bool post_optimise, bool special_concat,
       unsigned int use_bg = BOTTLENECK_GRAPH_NONE, bool verbose = false);


  /**
   * Constructs an approximated eucliedean Steiner minimal tree, using
   * the algorithm described in our report. This constructor takes a Delaunay
   * triangulation as input, to skip this step in the algorithm.
   *
   * @param del             The Delaunay data structure. The algorithm will
   *                        use the dimension of this datastructure.
   * @param sh              A pointer to a subgraph heuristic.
   * @param concat_subgraphs   Determines wether subgraph concatenation should be
   *                           done.
   * @param post_optimise      Determines wether post optimisation using Smith's
   *                           should be applied
   * @param special_concat     Redo concatenation, adding non-covered faces.
   * @param verbose            If true, stats will be printed.
   */
  ESMT(Utils::Delaunay &del, SubgraphHeuristic *sh,
       bool concat_subgraphs, bool post_optimise, bool special_concat,
       unsigned int use_bg = BOTTLENECK_GRAPH_NONE, bool verbose = false);

  /**
   * Destructor
   */
  ~ESMT();

  /**
   * A structure for storing statistics for later retrival.
   *
   * Set compile-flag ESMT_COLLECT_STATS = 1 to collect stats.
   */
  struct Stats {
    /** Constructor */
    Stats();

    // Delaunay
    unsigned int no_of_simplices;
    // Covered faces
    std::vector<unsigned int> covered_faces;
    // All faces
    std::vector<unsigned int> faces;
    // Sub-trees in queue
    unsigned int sub_trees_in_queue;
    // Added sub-trees
    std::vector<unsigned int> added_sub_trees;
    unsigned int add_sub_trees_total;
    // SPs before post-optimisation
    unsigned int no_of_sp;
    // SPs after post-optimisation
    unsigned int no_of_sp_post_optimisation;
    // Overlapping (very close) SPs after post-optimisation
    unsigned int no_of_sp_overlapping;
  };

  /**
   * Getter for the statistics struct.
   *
   * @return  The statistics object.
   */
  ESMT::Stats *getStats();
  
  /**
   * Getter for the list of covered faces.
   *
   * @return  The list of covered faces.
   */
  std::vector<Graph> &getCoveredFaces();

  /**
   * Getter for the list of generated SMTs.
   *
   * @return  The list of generated SMTs.
   */
  std::vector<SteinerTree> &getComponents();

protected:
private:
  
  //////////////////////////////////////////////////////////
  // Private functions

  /**
   * Performs the ESMT algorithm.
   *
   * @param points             The given point set
   * @param sh                 A pointer to a subgraph heuristic.
   * @param concat_subgraphs   Determines wether subgraph concatenation should be
   *                           done.
   * @param post_optimise      Determines wether post optimisation using Smith's
   *                           should be applied
   * @param special_concat     Redo concatenation, adding non-covered faces.
   * @param verbose            If true, stats will be printed.
   */
  void findESMT(std::vector<Utils::Point> &points,
		SubgraphHeuristic *sh,
		bool concat_subgraphs,
		bool post_optimise,
		bool special_concat,
		unsigned int use_bg,
		bool verbose = false);

  /**
   * Performs the ESMT algorithm.
   *
   * @param del             A Delaunay triangulation structure
   * @param sh              A pointer to a subgraph heuristic.
   * @param concat_subgraphs   Determines wether subgraph concatenation should be
   *                           done.
   * @param post_optimise      Determines wether post optimisation using Smith's
   *                           should be applied
   * @param special_concat     Redo concatenation, adding non-covered faces.
   * @param verbose            If true, stats will be printed.
   */
  void findESMT(Utils::Delaunay &del,
		SubgraphHeuristic *sh,
		bool concat_subgraphs,
		bool post_optimise,
		bool special_concat,
		unsigned int use_bg,
		bool verbose = false);
  
  /**
   * Performs the concatenation step of the algorithm.
   * Assumes this->queue to contain all sub-trees in order by ratio.
   *
   * @param verbose       If true, stats will be printed
   */
  void doConcatenate(bool verbose = false);

  /**
   * Performs the concatenation step of the algorithm, using the given bottleneck graph
   * for MST lengths.
   * Assumes that this->queue contains all sub-trees in order by b-ratio.
   * 
   * @param verbose       If true, stats will be printed
   */
  void doConcatenateWithBottleneck(bool verbose = false);
  
  /**
   * Performs the concatenation step of the algorithm repeatedly, using the redo technic.
   * Assumes that this->queue contains all sub-trees in order by ratio.
   * 
   * @param verbose       If true, stats will be printed
   */
  void doConcatenateWithRedo(bool verbose = false);

  /**
   * Checks if the given SteinerTree can be added
   */
  bool concatCheck(SteinerTree &st, std::vector< Utils::DisjointSet<unsigned int> > &sets, bool *flags);
  
  /**
   * Adds the given SteinerTree to the final solution.
   * NOTE!!! Should only be called after concatCheck(...).
   */
  void concatAdd(SteinerTree &st, std::vector< Utils::DisjointSet<unsigned int> > &sets);
  
  /**
   * Performs after optimisation
   *
   * Adds extra SPs in order to create a full topology. Then optimises with
   * Smith, and removes SP coincident with terminals.
   *
   */
  void postOptimisation();
  
  /**
   * The findFaces procedure finds all covered (MST) faces from the Delaunay
   * triangulation, which contains the given point.
   * These faces will be placed in components.
   * The procedure will leave out faces, which cross simplex boundaries and
   * which contains points with a lower index.
   *
   * @param handles      A list of point handles from the Delaunay triangulation.
   *                     These contain simplex indices for each point.
   * @param components   A list to add the faces to.
   * @param connections  Lists the connections for each point in the MST.
   * @param point        The point to be contained in all faces.
   */
  void findCoveredFaces(std::vector< Utils::Delaunay::PointHandle > &handles,
			std::vector< std::vector<int> > &connections);

  /**
   * The recursive part of the findFaces procedure.
   *
   * @param handles      A list of point handles from the Delaunay triangulation.
   *                     These contain simplex indices for each point.
   * @param components   A list to add the faces to.
   * @param connections  Lists the connections for each point in the MST.
   * @param map          An array of size |N|, needed for the procedure.
   * @param point        Mapmax is the current maximum value in the map.
   * @param bound        No points with index < bound will be added.
   * @param prevSet      The index of the previous set.
   */
  void findCoveredFacesRec(std::vector< Utils::Delaunay::PointHandle > &handles,
			   std::vector< std::vector<int> > &connections,
			   std::vector< int > &currentSimplices,
			   Graph prevSet,
			   std::unordered_map<unsigned int, unsigned int> map,
			   unsigned int &mapmax,
			   bool *flag,
			   unsigned int cur,
			   unsigned int prev,
			   unsigned int bound);
  
  void findAllFaces(Utils::Delaunay &del);

  /**
   * Builds sausages recursively.
   *
   * @param prevSet              The set we are currently building.
   * @param prevTree             The previously constructed SteinerTree (contained in
   *                             a SubST struct).
   * @param added                Number of added points to the initial simplex.
   * @param next_simplex_index   Neighbouring index of next simplex for cur_simplex.
   * @param cur_simplex          Index of current simplex in simplices.
   * @param org_simplex          Index of original simplex in simplices.
   * @param c                    The number of needed connections in MST for next
   *                             point.
   */
  /*void buildSausage(std::vector<int> &prevSet,
		    ESMT::SubST &prevTree,
		    unsigned int added,
		    unsigned int next_simplex_index,
		    int cur_simplex,
		    int org_simplex,
		    unsigned int c);
  */
  /**
   * Builds sausages recursively in the other direction.
   *
   * @param prevSet              The set we are currently building.
   * @param prevTree             The previously constructed SteinerTree (contained in
   *                             a SubST struct).
   * @param prev_added           Number of added points to the initial simplex
   *                             in the first build step (buildSausage).
   * @param added                Number of added points in reverse.
   * @param next_simplex_index   Neighbouring index of next simplex for cur_simplex.
   * @param cur_simplex          Index of current simplex in simplices.
   */
  /*void buildSausageReverse(std::vector<int> &prevSet,
			   ESMT::SubST &prevTree,
			   unsigned int prev_added,
			   unsigned int added,
			   unsigned int next_simplex_index,
			   int cur_simplex);
  */
  /**
   * Checks if the edge between i0 and i1 is in the MST.
   *
   * @param i0  The first end point of the edge.
   * @param i1  The second end point of the edge.
   *
   * @return    True, if i0-i1 is in the MST, otherwise false.
   */
  bool isInMST(int i0, int i1);

  void setBLength(SteinerTree &st);

  ////////////////////////////////////////////////
  // Static functions
  
  /**
   * Compares two SubST trees with respect to the Steiner ratio.
   *
   * @param st1   The first SubST.
   * @param st2   The second SubST.
   *
   * @return      True if st1.ratio > st2.ratio
   */
  static bool compareSteinerRatio(const SteinerTree &st1, const SteinerTree &st2);

  static bool compareSteinerBRatio(const SteinerTree &st1, const SteinerTree &st2);

  /**
   * Compares two edges with respect to the length.
   *
   * @param e1   The first edge.
   * @param e2   The second edge.
   *
   * @return     True, if e1.length > e2.length
   */
  static bool compareLength(Utils::Edge e1, Utils::Edge e2);

  /////////////////////////////////////////////////
  // Attributes

  /** Simplices from Delaunay */
  std::vector<Utils::Delaunay::Simplex> simplices;
  /** Flags - set to true if simplex is covered */
  std::vector<bool> is_covered_simplex;
  
  /** Components to create SMTs for */
  std::vector<Graph> components;
  std::unordered_map<unsigned int, unsigned int> simplex_id;
  
  /** List of SMTs for concatenation */
  std::vector<SteinerTree> smts;

  /** Priority queue for concatenation */
  Utils::Heap<SteinerTree> *queue;

  /** Map showing if the edge i1,i2 is in the MST,
      where i1, i2 is indexes in this->points. */
  std::unordered_map<unsigned long, bool> in_MST;
  
  /** Stats object */
  Stats stats;
  
  /** Iterative concat used for sausage construction */
  IterativeConcat iterCon;

  /** Iterative Smith used for recomputing FSTs */
  IterativeSmith iterSmith;

  /** BottleneckGraph object */
  BottleneckGraph *bgraph;
};

#endif // ESMT_H
