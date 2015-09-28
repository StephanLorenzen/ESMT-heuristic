#include <vector>
#include <algorithm>
#include <unordered_map>
#include <assert.h>

#include "steiner/graph.hpp"
#include "steiner/steiner_tree.hpp"
#include "steiner/esmt.hpp"
#include "steiner/iterative.hpp"
#include "steiner/heuristics/subgraph_heuristic.hpp"
#include "steiner/heuristics/concat.hpp"
#include "steiner/heuristics/smith.hpp"
#include "steiner/heuristics/steiner_finder.hpp"
#include "steiner/utils/point.hpp"
#include "steiner/utils/utils.hpp"
#include "steiner/utils/fermat.hpp"
#include "steiner/utils/disjoint_set.hpp"
//#include "steiner/utils/priority_queue.hpp"
#include "steiner/utils/delaunay.hpp"
#include "steiner/utils/bottleneck_graph/bg_naive.hpp"

#define ESMT_COLLECT_STATS  0

/*
 * Utils type definitions
 */
typedef Utils::Point                      Point;
typedef Utils::Edge                       Edge;
typedef Utils::DisjointSet<unsigned int>  DisjointSet;
typedef Utils::Delaunay                   Delaunay;
typedef Utils::Delaunay::Simplex          Simplex;
typedef Utils::Delaunay::PointHandle      PointHandle;
typedef Utils::PriorityQueue<SteinerTree> Queue;

/*
 * Constructor
 */
ESMT::ESMT(std::vector<Point> &points)
  : SteinerTree(points), Iterative(points[0].dim(), points.size()),
    iterCon(points[0].dim()), iterSmith(points[0].dim())
{
  SubgraphHeuristic *sh = new IterativeSmith(points[0].dim());
  this->findESMT(points, sh, true, true, 0, false);
  delete sh;
}

/*
 * Constructor
 */
ESMT::ESMT(std::vector<Point> &points, SubgraphHeuristic *sh,
	   bool concat_subgraphs, bool post_optimise, bool special_concat,
	   unsigned int use_bg, bool verbose)
  : SteinerTree(points), Iterative(points[0].dim(), points.size()),
    iterCon(points[0].dim()), iterSmith(points[0].dim())
{
  if(verbose)
    std::cout << "*** ESMT Heuristic ***" << std::endl;
  this->findESMT(points, sh, concat_subgraphs, post_optimise, special_concat, use_bg, verbose);
  if(verbose)
    std::cout << "**********************" << std::endl;
}

/*
 * Constructor
 */
ESMT::ESMT(Utils::Delaunay &del, SubgraphHeuristic *sh,
	   bool concat_subgraphs, bool post_optimise, bool special_concat,
	   unsigned int use_bg, bool verbose)
  : SteinerTree(del.getPointsRef()), Iterative(del.dimension(), del.getPoints().size()),
    iterCon(del.dimension()), iterSmith(del.dimension())
{
  if(verbose)
    std::cout << "*** ESMT Heuristic ***" << std::endl;

  this->findESMT(del, sh, concat_subgraphs, post_optimise, special_concat, use_bg, verbose);
  if(verbose)
    std::cout << "**********************" << std::endl;
}

ESMT::Stats *ESMT::getStats() {
  return &this->stats;
}

std::vector<SteinerTree> &ESMT::getComponents() {
  return this->smts;
}

/*
 * Implementation of findESMT(...)
 */
void ESMT::findESMT(std::vector<Point> &points,
		    SubgraphHeuristic *sh,
		    bool concat_subgraphs,
		    bool post_optimise,
		    bool special_concat,
		    unsigned int use_bg,
		    bool verbose) {
  // Do Delaunay
  Delaunay del(points);
  this->findESMT(del, sh, concat_subgraphs, post_optimise, special_concat, use_bg, verbose);
}

/*
 * Implementation of findESMT(...) with Delaunay as argument
 */
void ESMT::findESMT(Delaunay &del, 
		    SubgraphHeuristic *sh,
		    bool concat_subgraphs,
		    bool post_optimise,
		    bool special_concat,
		    unsigned int use_bg,
		    bool verbose) {
  unsigned int i,bits;
  this->bgraph = NULL;
  
  assert(sh);
  this->dim = del.dimension();
  this->N = this->n();
  
  if(use_bg) {
    // Don't do special concat if use_bg
    special_concat = false;
    // Check that there are not too many points:
    //  dim = 2 -> 64 / 2 ~= 32 bits available
    //  dim = 3 -> 64 / 3 ~= 21 bits available
    //  dim = 4 -> 64 / 4 ~= 16 bits available
    //  dim = 5 -> 64 / 5 ~= 12 bits available
    //  dim = 6 -> 64 / 6 ~= 10 bits available
    bits = 64 / dim;
    if(this->N > (unsigned long)(2L << bits)) {
      std::cerr << "Error in ESMT using B-Graph: too many points for this dimension." << std::endl;
      exit(1);
    }
    
    this->bgraph = new BottleneckGraphNaive(del);
  }
  this->edges              = del.getEdges();
  this->simplices          = del.getSimplices();
  this->is_covered_simplex = std::vector<bool>(this->simplices.size(), false);
  std::vector<PointHandle> &handles = del.getPointHandles();
  
#if(ESMT_COLLECT_STATS)
  this->stats.no_of_simplices = this->simplices.size();
#endif

  if(verbose) {
    std::cout << "Delaunay done!" << std::endl
	      << "  Simplices:  " << this->simplices.size() << std::endl;
  }
  
  // Get the mst.
  Utils::MSTKruskalMod(*this);
  // Mark entries in is_MST
  std::vector<Edge>::iterator eit;
  for(eit = this->edges.begin(); eit != this->edges.end(); ++eit)
    in_MST[eit->key()] = true;
  
  this->setMSTLength(this->getLength());
  
  if(verbose) {
    std::cout << "MST done!" << std::endl
	      << "  Length:  " << this->getMSTLength() << std::endl;
  }
  
  // Get components to concatenate
  this->components.reserve(this->n());
  
  //////////////////
  // Preprocessing
  
  // 1. Build the connections list and sort all handles
  std::vector< std::vector<int> > connections(this->n());
  for(i = 0; i < this->n(); i++)
    std::sort(handles[i].simplices.begin(),handles[i].simplices.end());
  unsigned int i0, i1;
  for(i = 0; i < this->edges.size(); i++) {
    i0 = this->edges[i].i0;
    i1 = this->edges[i].i1;
    connections[i0].push_back(i1);
    connections[i1].push_back(i0);
  }
  
  // 2. Generate faces
  //if(false && use_bg) {
    // We have to find all faces
    /*std::unordered_map<unsigned long,bool> flags;
    std::vector<Simplex>::iterator sit;
    for(sit = this->simplices->begin(); sit != this->simplices->end(); sit++) {
      Subset simp;
      Subset ss;
      ss.map.push_back(0);
      for(i = 0; i < sit->n; i++) {
	ss.map[0] = sit->map[i];
	this->findAllFaces(*sit, ss, flags);
	// Add simplex (simp) outside, since there is not room for key in flag.
	// (key = d * (64/d-bit) numbers).
	simp.map.push_back(sit->map[i]);
      }
      this->components.push_back(simp);
    }
    for(i = 0; i < this->edges.size(); i++) {
      // Add all MST edges
      Subset ed;
      Edge &e = this->edges[i];
      ed.map.push_back(e.i0);
      ed.map.push_back(e.i1);
      this->components.push_back(ed);
      }*/
  if(use_bg)
    // Find all faces
    this->findAllFaces(del);
  else
    // Only find covered faces
    this->findCoveredFaces(handles, connections);
  
  if(verbose) {
    std::cout << "Preprocessing done!" << std::endl
	      << "  Number of covered faces: " << this->components.size() << std::endl;
#if(ESMT_COLLECT_STATS)
    for(i = 2; i < this->stats.covered_faces.size(); i++)
      std::cout << "   [" << i << "]: "
		<< this->stats.covered_faces[i] << std::endl;
#endif
  }
  
  // We now have a list of all faces (or only those covered by the MST, if not use_bg).
  // Begin finding small SMTs.
  // Allocate vector now. Approx each element should be added.
  this->smts.reserve(this->components.size());
  
  std::vector< Graph >::iterator sit;
  for(sit = this->components.begin(); sit != this->components.end(); sit++) {
    if(sit->n() == 2) {
      std::vector<Edge> edges;
      edges.push_back(Edge(0,1));
      SteinerTree st(sit->getPoints(), sit->getPointsRef(), edges);
      double l = st.getLength();
      st.setMSTLength(l);
      st.setSMTLength(l);
      st.setBMSTLength(l);
      st.computeRatios();
      this->smts.push_back(st);
    }
    else {
      if(use_bg)
	// First compute MST
	Utils::MSTKruskalMod(*sit, true);
      SteinerTree st(sit->getPoints(), sit->getPointsRef(), sit->getEdges());
      st.setMSTLength(st.getLength());
      
      if(sit->n() == 3)
	Utils::getFermatSMT(st);
      else
	sh->findSteinerPoints(st);
      
      st.setSMTLength(st.getLength());
      st.computeRatios();
      
      if(st.getSteinerRatio() >= 1.0)
	continue;

      if(use_bg)
	st.setBMSTLength(this->bgraph->getBMSTLength(st.getPoints()));
      st.computeRatios();
      
      this->smts.push_back(st);
      // TODO disable sauages for now
      /*if(false && sit->map.size() == this->dim+1 && concat_subgraphs) {
	int simplex = sit->simplex_index;
	std::vector<int> map;
	for(i = 0; i < (*this->simplices)[simplex].n; i++)
	  map.push_back((*this->simplices)[simplex].map[i]);
	for(i = 0; i < map.size(); i++) {
	  // Swap first element with i
	  int tmp = map[0];
	  map[0]  = map[i];
	  map[i]  = tmp;
	  this->buildSausage(map, subst, 0, i, simplex, simplex, 1);
	}
	}*/
    }
  }
  
#if(ESMT_COLLECT_STATS)
  this->stats.sub_trees_in_queue = this->smts.size();
#endif

  if(verbose) {
    std::cout << "Sub-trees found!" << std::endl
	      << "  Number of sub-trees in queue: " << this->smts.size() << std::endl;
#if(ESMT_COLLECT_STATS)
    std::cout << "  Number of covered sausages: " << std::endl;
    for(i = this->dim+2; i < this->stats.covered_faces.size(); i++)
      std::cout << "   [" << i << "]: "
		<< this->stats.covered_faces[i] << std::endl;
#endif
  }
  
  // TODO - what is wrong here??? this->queue = Queue(this->smts, use_bg ? ESMT::compareSteinerBRatio : ESMT::compareSteinerRatio);
  // Sort priority queue of sub-graphs according to ratio
  if(use_bg)
    std::sort(this->smts.begin(), this->smts.end(), ESMT::compareSteinerBRatio);
  else
    std::sort(this->smts.begin(), this->smts.end(), ESMT::compareSteinerRatio);
  
  // Sort priority queue of sub-graphs according to ratio
  //if(use_bg)
  //  std::sort(this->smts.begin(), this->smts.end(), ESMT::compareSteinerBRatio);
  //else
  //std::sort(this->smts.begin(), this->smts.end(), ESMT::compareSteinerRatio);
  
  /*if(special_concat) {
    // Starting concat
    this->doConcatenate(false);
    std::vector<SubST> non_covered;
    // Preprocessing: find SMTs of all non_covered
    std::vector< Simplex >::iterator spit;
    for(i = 0, spit = this->simplices->begin();
	spit != this->simplices->end(); i++, spit++) {
      if(this->is_covered_simplex[i])
	// Simplex is covered. Continue
	continue;
      // All is of size = d+1 (since only simplices)
      Graph submst;
      for(i = 0; i < spit->n; i++) {
	submst.getPoints().push_back(this->points[spit->map[i]]);
	for(j = i+1; j < spit->n; j++) {
	  int a = spit->map[i];
	  int b = spit->map[j];
	  submst.getEdgesPtr()->push_back(Edge(i,j,Utils::length(this->points[a],this->points[b])));
	}
      }

      submst = Utils::MSTKruskal(submst);
      submst.setMSTLength(submst.getLength());

      SteinerTree *st;
      if(spit->n == 3)
	st = Utils::getFermatSMT(submst);
      else
	st = sh->findSteinerPoints(submst);
      if(st->getSteinerRatio() >= 1.0) {
	delete st;
	continue;
      }
      SubST subst(st,spit->n,spit->map);
      non_covered.push_back(subst);
    }

    std::sort(non_covered.begin(), non_covered.end(), ESMT::compareSteinerRatio);

    // Now, redo concatenation
    std::vector<SubST> pre_add;
    double best_ratio = this->getSteinerRatio();
    unsigned int k = 1000; // Maximum of max(1000,non_covered.size()) tries
    if(verbose)
      std::cout << "Ratio after initial concatenation: "
		<< best_ratio << std::endl;
    for(i = 0; i < k && i < non_covered.size(); i++) {
      pre_add.push_back(non_covered[i]);
      // Now try to concatenate again.
      this->points.erase(this->points.begin()+this->N, this->points.end());
      this->doConcatenate(pre_add, false);
      this->setSMTLength(-1);
      this->setSteinerRatio(-1);
      if(this->getSteinerRatio() < best_ratio) {
	if(verbose)
	  std::cout << "Better ratio achieved: " << this->getSteinerRatio()
		    << std::endl;
	best_ratio = this->getSteinerRatio();
      }
      else {
	// Remove this FST from pre_add again.
	pre_add.pop_back();
      }
    }

#if(ESMT_COLLECT_STATS)
    this->stats.added_sub_trees.clear();
    this->stats.add_sub_trees_total = 0;
#endif

    this->points.erase(this->points.begin()+this->N, this->points.end());
    this->doConcatenate(pre_add, verbose);
    this->setSMTLength(-1);
    this->setSteinerRatio(-1);
  }
  else if(use_bg) {
    this->doConcatenateB(verbose);
  }
  else {}*/
  if(use_bg)
    this->doConcatenateWithBottleneck(verbose);
  else
    this->doConcatenate(verbose);
    
  this->setSMTLength(this->getLength());
  this->computeRatios();

#if(ESMT_COLLECT_STATS)
  this->stats.no_of_sp = this->s;
#endif

  /*if(post_optimise) {
    if(verbose)
      std::cout << "Starting post optimisation." << std::endl
		<< "  Current |SP| = " << this->points.size()-this->N << std::endl;
    // Setup for Smiths and add extra SPs. Then optimise, do clean-up and return.
    this->postOptimisation();

#if(ESMT_COLLECT_STATS)
    // Find overlapping SPs
    double l = this->length();
    l /= this->N+this->S;
    for(i = 0; i < this->N-2; i++) {
      for(j = 0; j < 3; j++) {
	unsigned int op = this->adj[i][j];
	if(op >= this->N && op-this->N > i) {
	  if(this->EL[i][j] < 0.0001*l)
	    this->stats.no_of_sp_overlapping++;
	}
      }
    }
#endif

    // Convert to ESMT structure
    this->points.erase(this->points.begin()+this->N, this->points.end());
    this->edges.clear();

    this->cleanUp(&this->points, &this->edges);

#if(ESMT_COLLECT_STATS)
    this->stats.no_of_sp_post_optimisation = this->points.size()-this->N;
#endif

    if(verbose) {
      std::cout << "Post optimisation done!" << std::endl
		<< "  |SP| = " << this->points.size()-this->N << std::endl;
#if(ESMT_COLLECT_STATS)
      std::cout << "  Number of short SP-SP edges: "
		<< this->stats.no_of_sp_overlapping << std::endl;
#endif
    }
    }*/
  
  if(use_bg)
    delete this->bgraph;
}

/*
 * Destructor
 */
ESMT::~ESMT() {
}

// Private functions

void ESMT::doConcatenate(bool verbose) {
  unsigned int i, c;
  // Concatenation process.
  // This is really just MSTKruskal.

  // Clear edges from Delaunay
  this->edges.clear();
  
  if(verbose)
    std::cout << "Starting concatenation process." << std::endl;  

  // First, create Disjoint sets. Each set is build over vertex indicies
  std::vector< DisjointSet > sets;
  for(i = 0; i < this->n(); i++)
    sets.push_back(DisjointSet(i));

  // Flags used to determine if two sets are disjoint
  bool *flags = new bool[this->N];
  for(i = 0; i < this->N; i++)
    flags[i] = false;

  c = 1;
  i = 0;
  while(true) {
    SteinerTree st = this->smts[i];//this->queue.next();
    if(this->checkAndAdd(st,sets,flags)) {
      c += st.n()-1;
      if(c >= this->n())
	break;
    }
    i++;
  }
  
  delete flags;
}

void ESMT::doConcatenateWithBottleneck(bool verbose) {
  unsigned int i, c;
  // Concatenation process.
  // This is really just MSTKruskal.

  // Clear edges from Delaunay
  this->edges.clear();
  
  if(verbose)
    std::cout << "Starting concatenation process." << std::endl;  
  // First, create Disjoint sets. Each set is build over vertex indicies
  std::vector< DisjointSet > sets;
  for(i = 0; i < this->n(); i++)
    sets.push_back(DisjointSet(i));

  // Flags used to determine if two sets are disjoint
  bool *flags = new bool[this->N];
  for(i = 0; i < this->N; i++)
    flags[i] = false;
  
  c = 1;
  unsigned idx = 0;
  //this->queue.rebuild();
  while(true) {
    SteinerTree st = this->smts[idx];//this->queue.next();
    if(this->checkAndAdd(st,sets,flags)) {
      c += st.n()-1;
      if(c >= this->n())
	break;
      // Recompute the queue
      this->bgraph->mergePoints(st.getPoints());
      for(i = idx; i < this->smts.size(); i++) {
	SteinerTree &st = this->smts[i];
	st.setBMSTLength(this->bgraph->getBMSTLength(st.getPoints()));
	st.computeRatios();
      }
      idx++;
      std::sort(this->smts.begin()+idx, this->smts.end(), ESMT::compareSteinerBRatio);
    }
    else
      idx++;
  }
  delete flags;
}

void ESMT::doConcatenateWithRedo(bool verbose) {

}

bool ESMT::checkAndAdd(SteinerTree &st, std::vector<DisjointSet> &sets, bool *flags) {
  unsigned int i;
  //std::cout << st << std::endl;
  std::vector<unsigned int> flags_set;
  for(i = 0; i < st.n(); i++) {
    // id is the vertex index in the full graph
    unsigned int id = sets[st.pidx(i)].findSet()->getValue();
    if(flags[id]) {
      // Not good, has been encountered before = same set
      // Unset all flags and return
      for(i = 0; i < flags_set.size(); i++)
	flags[flags_set[i]] = false;
      return false;
    }
    // Set flag
    flags[id] = true;
    flags_set.push_back(id);
  }
  // No conflicts, this graph may be inserted without problem
  
  // Add points
  std::vector<unsigned int> indexes = st.getPoints();
  for(i = 0; i < st.s(); i++) {
    this->steiner_points.push_back(st.getSteinerPoint(i));
    indexes.push_back(this->m()-1);
  }
  // Add edges
  std::vector<Edge> edges = st.getEdges();
  std::vector<Edge>::iterator eit;
  for(eit = edges.begin(); eit != edges.end(); eit++) {
    Edge e(indexes[eit->i0],indexes[eit->i1]);
    this->edges.push_back(e);
  }
  // Now union all
  // All points in flag_set should be unioned
  for(i = 1; i < flags_set.size(); i++)
    sets[flags_set[0]].setUnion(sets[flags_set[i]]);

#if(ESMT_COLLECT_STATS)
  while(this->stats.added_sub_trees.size() < st.n()+1)
    this->stats.added_sub_trees.push_back(0);
  this->stats.added_sub_trees[st.n()]++;
  this->stats.add_sub_trees_total++;
#endif

  // Unset all flags and return
  for(i = 0; i < flags_set.size(); i++)
    flags[flags_set[i]] = false;
  return true;
}

//void ESMT::doConcatenateB(bool verbose) {
  /*unsigned int i, j, k, c, index;
  // Concatenation process.
  // This is really just MSTKruskal.

  //std::cout << "Components: " << this->smts.size() << std::endl;
  //std::cout << "Simplx: " << this->simplices->size() << std::endl;
  
  // Clear edges from Delaunay
  this->edges.clear();

  if(verbose) {
    std::cout << "Starting B-concatenation process." << std::endl;
  }  

  c = 0;
  // First, create Disjoint sets. Each set is build over vertex indicies
  std::vector< DisjointSet > sets;
  for(i = 0; i < this->points.size(); i++) {
    sets.push_back(DisjointSet(i));
  }
  
  // Flags used to determine if two sets are disjoint
  bool *flags = new bool[this->N];
  for(i = 0; i < this->N; i++)
    flags[i] = false;
  
  std::vector<SubST> tsmts = this->smts;
  bool done = false;
  while(!done) {
    for(j = 0; j < tsmts.size(); j++) {
      //std::cout << j << ":"<<tsmts.size() << std::endl;
      SubST &sit = tsmts[j];
      std::vector<int> flags_set;
      // Check for conflicts
      for(i = 0; i < sit.n; i++) {
	// id is the vertex index in the full graph
	int id = sets[sit.map[i]].findSet()->getValue();
	if(flags[id])
	  // Not good, has been encountered before = same set
	  break;
	// Set flag
	flags[id] = true;
	flags_set.push_back(id);
      }
      // Reset flags
      for(k = 0; k < flags_set.size(); k++) {
	flags[flags_set[k]] = false;
      }
      if(i == sit.n) {
	// No conflicts, this graph may be inserted without problem
	SteinerTree *st = sit.st;
  */
	/*std::cout << "added: ";
	for(int uu =0 ; uu < sit.n; uu++)
	  std::cout << sit.map[uu] << " ";
	std::cout << "smt="<<sit.st->getSMTLength()<< ", bmst="<<sit.bmst_length<<std::endl;
	*/
  /*	std::vector<Point> *points = st->getPointsPtr();
	std::vector<Edge>  *edges  = st->getEdgesPtr();
	std::vector<Point>::iterator pit;
	std::vector<Edge>::iterator  eit;
	unsigned int no_of_terminals = points->size();

	// Add points
	index = 0;
	std::vector<int> indexes;
	for(pit = points->begin(); pit != points->end(); pit++) {
	  if(pit->isSteiner()) {
	    no_of_terminals--;
	    this->points.push_back(*pit);
	    indexes.push_back(this->points.size()-1);
	  }
	  else
	    indexes.push_back(sit.map[index++]);
	}
	// Add edges
	for(eit = edges->begin(); eit != edges->end(); eit++) {
	  Edge e(indexes[eit->i0],indexes[eit->i1]);
	  this->edges.push_back(e);
	}
	// Now union all
	// All points in flag_set should be unioned
	for(i = 1; i < flags_set.size(); i++) {
	  sets[flags_set[0]].setUnion(sets[flags_set[i]]);
	}

#if(ESMT_COLLECT_STATS)
	while(this->stats.added_sub_trees.size() < no_of_terminals+1)
	  this->stats.added_sub_trees.push_back(0);
	this->stats.added_sub_trees[no_of_terminals]++;
	this->stats.add_sub_trees_total++;
#endif
	
	// Increment counter
	c += no_of_terminals-1;
	break;
      }
    }
    // Check if we have added N-1 mst-edges
    if(c == this->N-1) {
      // Connected graph
      break;
    }
    
    // We have to update the priority queue
    std::vector<SubST> newQueue(tsmts.begin()+j+1, tsmts.end());
    // Update B-graph
    SubST &sit = tsmts[j];
    std::vector<unsigned int> idxs;
    for(k = 0; k < sit.n; k++)
      idxs.push_back(sit.map[k]);
    this->bgraph->mergePoints(idxs);
    for(k = 0; k < newQueue.size(); k++)
      this->setBLength(newQueue[k]);
    std::sort(newQueue.begin(), newQueue.end(), ESMT::compareSteinerBRatio);
    tsmts = newQueue;
  }

  if(verbose) {
    std::cout << "Concatenation done." << std::endl;
#if(ESMT_COLLECT_STATS)
    std::cout << "  Sub-trees added: " << std::endl;
    for(i = 2; i < this->stats.added_sub_trees.size(); i++)
      std::cout << "   [" << i << "]: "
		<< this->stats.added_sub_trees[i] << std::endl;
    std::cout << "   Total: " << this->stats.add_sub_trees_total << std::endl;
#endif
  }

  delete flags;*/
//}

/* Implementation of postOptimisation() */
void ESMT::postOptimisation() {
  /*unsigned int i, j, k, v;

  this->S = this->points.size() - this->N;
  this->P.clear();
  for(i = 0; i < this->N+this->S; i++)
    this->P.push_back(this->points[i]);
  for(i = this->N+this->S; i < 2*this->N-2; i++)
    this->P.push_back(Point(this->dim));

  for(i = 0; i < this->N; i++)
    for(j = 0; j < 3; j++)
      this->adj[i][j] = -1;

  std::vector< std::vector<int> > tadj(this->N);

  // Fill adj and tadj array
  unsigned int i0, i1, is0, is1;
  for(i = 0; i < this->edges.size(); i++) {
    i0 = this->edges[i].i0;
    i1 = this->edges[i].i1;
    is0 = i0-this->N;
    is1 = i1-this->N;
    if(i0 >= this->N)
      for(j = 0; j < 3; j++) {
	if(this->adj[is0][j] < 0) {
	  this->adj[is0][j] = i1;
	  break;
	}
      }
    else
      tadj[i0].push_back(i1);
    if(i1 >= this->N)
      for(j = 0; j < 3; j++) {
	if(this->adj[is1][j] < 0) {
	  this->adj[is1][j] = i0;
	  break;
	}
      }
    else
      tadj[i1].push_back(i0);
  }

  // Adj and Tadj are now up to date.
  // Now do the inserts of new SPs.
  for(i = 0; i < this->N; i++) {
    // Insert SPs if valens > 1
    unsigned int val;
    while((val = tadj[i].size()) > 1) {
      // This terminal has valens > 1
      float smallest_angle = Utils::angle(this->P[tadj[i][0]], this->P[i], this->P[tadj[i][1]]);
      float tmp_angle;
      int   a, b;
      a = 0;
      b = 1;
      for(j = 0; j < val; j++) {
	for(k = j+1; k < val; k++) {
	  tmp_angle = Utils::angle(this->P[tadj[i][j]], this->P[i], this->P[tadj[i][k]]);
	  if(tmp_angle > smallest_angle) {
	    smallest_angle = tmp_angle;
	    a = j;
	    b = k;
	  }
	}
      }

      // Best angle is now aib. Insert point.
      this->P[this->N+this->S] = (this->P[tadj[i][a]] + this->P[i] + this->P[tadj[i][b]]) * (1.0 / 3.0);

      // Update edges
      this->adj[this->S][0] = tadj[i][a];
      this->adj[this->S][1] = tadj[i][b];
      this->adj[this->S][2] = i;

      if(tadj[i][a] >= (int)this->N) {
	int ai = tadj[i][a]-this->N;
	for(j = 0; j < 3; j++)
	  if(this->adj[ai][j] == (int)i) {
	    this->adj[ai][j] = this->N+this->S;
	    break;
	  }
      }
      else {
	v = tadj[tadj[i][a]].size();
	for(j = 0; j < v; j++)
	  if(tadj[tadj[i][a]][j] == (int)i) {
	    tadj[tadj[i][a]][j] = this->N+this->S;
	    break;
	  }
      }
      if(tadj[i][b] >= (int)this->N) {
	int bi = tadj[i][b]-this->N;
	for(j = 0; j < 3; j++)
	  if(this->adj[bi][j] == (int)i) {
	    this->adj[bi][j] = this->N+this->S;
	    break;
	  }
      }
      else {
	v = tadj[tadj[i][b]].size();
	for(j = 0; j < v; j++)
	  if(tadj[tadj[i][b]][j] == (int)i) {
	    tadj[tadj[i][b]][j] = this->N+this->S;
	    break;
	  }
      }

      // Update a
      tadj[i][a] = this->N+this->S;
      // Erase b
      tadj[i].erase(tadj[i].begin()+b);
      this->S++;
    }
  }

  // Now all terminals should have val = 1.
  // Optimise
  double l = this->length();
  double r = this->error(); 
  do {
    this->optimise(0.0001*r/this->N);
    l = this->length();
    r = this->error();
    } while(r>l*0.0001);*/
}

bool ESMT::compareSteinerRatio(const SteinerTree &st1, const SteinerTree &st2) {
  if(st1.n() == 2 && st2.n() == 2)
    return st1.getMSTLength() < st2.getMSTLength();
  return st1.getSteinerRatio() < st2.getSteinerRatio();
}

bool ESMT::compareSteinerBRatio(const SteinerTree &st1, const SteinerTree &st2) {
  if(st1.n() == 2 && st2.n() == 2)
    return st1.getMSTLength() < st2.getMSTLength();
  return st1.getBSteinerRatio() < st2.getBSteinerRatio();
}

bool ESMT::compareLength(Edge e1, Edge e2) {
  return e1.length < e2.length;
}

bool ESMT::isInMST(int i0, int i1) {
  return (this->in_MST.find(Edge::key(i0,i1)) != this->in_MST.end());
}

void ESMT::setBLength(SteinerTree &s) {
  unsigned int i, j;
  std::vector<Edge> edges;
  for(i = 0; i < s.n(); i++) {
    for(j = 0; j < i; j++) {
      edges.push_back(Edge(i,j,this->bgraph->distance(s.pidx(i),s.pidx(j))));
    }
  }
  Graph g(s.getPoints(), s.getPointsRef(), edges);
  Utils::MSTKruskalMod(g, false);
  double bmst_length = 0.0;
  for(i = 0; i < edges.size(); i++) {
    bmst_length += edges[i].length;
  }
  s.setBMSTLength(bmst_length);
}
//void ESMT::buildMST(Subset &subset, Graph &g) {
  /*std::vector<Point> *points = g.getPointsPtr();
  std::vector<Edge>  *edges  = g.getEdgesPtr();
  points->clear();
  edges->clear();

  unsigned int i, j, i0, i1;
  for(i = 0; i < subset.map.size(); i++) {
    points->push_back(this->points[subset.map[i]]);
    for(j = i+1; j < subset.map.size(); j++) {
      i0 = subset.map[i];
      i1 = subset.map[j];
      if(this->isInMST(i0, i1))
	edges->push_back(Edge(i, j,
			      Utils::length(this->points[i0], this->points[i1])));
    }
  }
  g.setMSTLength(g.getLength());*/
//}

void ESMT::findCoveredFaces(std::vector< PointHandle > &handles,
			    std::vector< std::vector<int> > &connections) {

  unsigned int i, p;
  std::unordered_map<unsigned int, unsigned int> map;
  bool *flag = new bool[this->n()];
  for(i = 0; i < this->n(); i++)
    flag[i] = false;
  for(p = 0; p < this->n(); p++) {
    unsigned int mapmax = 0;
    std::vector<Pidx> pidxs;
    pidxs.push_back(p);
    Graph sps(pidxs, this->getPointsRef());
    map[p] = mapmax++;
    flag[p] = true;
    this->findCoveredFacesRec(handles, connections, handles[p].simplices,
			      sps, map, mapmax, flag, p, p, p);
    flag[p] = false;
    map.clear();
  }
  delete flag;
}

void ESMT::findCoveredFacesRec(std::vector< PointHandle > &handles,
			       std::vector< std::vector<int> > &connections,
			       std::vector< int > &currentSimplices,
			       Graph prevSet,
			       std::unordered_map< unsigned int, unsigned int > map,
			       unsigned int &mapmax,
			       bool *flag,
			       unsigned int cur,
			       unsigned int prev,
			       unsigned int bound) {
  unsigned int c, i, j, rank, next;
  std::vector<int> intersection;
  
  if(cur != bound) {
    // Check, if we may add this point.
    // Assume point handle lists and currentSimplices to be sorted.
    std::set_intersection(currentSimplices.begin(),currentSimplices.end(),
			  handles[cur].simplices.begin(), 
			  handles[cur].simplices.end(),
			  std::back_inserter(intersection));
    if(intersection.size() == 0)
      // We cannot add this point :(
      return;
    // Add this point, and add the new set
    std::vector<Pidx> &points = prevSet.getPoints();
    std::vector<Edge> &edges  = prevSet.getEdges();
    points.push_back(cur);
    edges.push_back(Edge(prev, points.size()-1));
    this->components.push_back(prevSet);
    c = this->components.size()-1;

    if(intersection.size() == 1 && prevSet.n() == this->dim+1)
      // This is a simplex. Add index to simplex list
      this->simplex_id[c] = intersection[0];

#if(ESMT_COLLECT_STATS)
    //unsigned int size = components[prevSet].map.size();
    //while(this->stats.covered_faces.size() < size+1)
    //  this->stats.covered_faces.push_back(0);
    //this->stats.covered_faces[size]++;
#endif
  }
  else
    intersection = currentSimplices;
  
  // For each connection, try to add new point.
  rank = map[prevSet.pidx(prevSet.n()-1)];
  for(j = 0; j < prevSet.n(); j++) {
    cur = prevSet.pidx(j);
    //std::cout << "j="<<j<<" cur="<<cur << std::endl;
    for(i = 0; i < connections[cur].size(); i++) {
      next = connections[cur][i];
      bool notFound = (map.find(next) == map.end());
      //std::cout << " n = "<<next<<", nfound = "<< notFound << std::endl;
      if(next > bound && (notFound || map[next] > rank) && !flag[next]) {
	if(notFound)
	  map[next] = mapmax++;
	flag[next] = true;
	this->findCoveredFacesRec(handles, connections, intersection, prevSet, map, mapmax,
				  flag, next, j, bound);
	flag[next] = false;
      }
    }
  }
}

void ESMT::findAllFaces(Delaunay &del) {
  unsigned int i;
  std::vector<Graph> &faces = del.findFaces();
  
  this->components.clear();
  this->components.reserve(faces.size());
  for(i=0; i<faces.size(); i++)
    if(faces[i].n() > 2)
      this->components.push_back(faces[i]);
  // Get only MST edges
  for(i=0; i<this->edges.size(); i++) {
    std::vector<Pidx> pidxs;
    pidxs.push_back(this->edges[i].i0);
    pidxs.push_back(this->edges[i].i1);
    this->components.push_back(Graph(pidxs, this->getPointsRef()));
  }
}

/*void ESMT::buildSausage(std::vector<int> &prevSet,
			SubST &prevTree,
			unsigned int added,
			unsigned int next_simplex_index,
			int cur_simplex,
			int org_simplex,
			unsigned int c) {
  unsigned int i, j, z;
  double d_mst_length = 0.0;

  Simplex *cur  = &(*this->simplices)[cur_simplex];
  int next_vertex = cur->nextVertex[next_simplex_index];
  int next_simplex = cur->nextSimplex[next_simplex_index];

  // Check for neighbour?
  if(next_vertex < 0)
    return;

  Simplex *next = &(*this->simplices)[next_simplex];

  // Detect if we are connected - must have at least c connections
  for(i = 0; i < cur->n; i++)
    if(this->isInMST(cur->map[i], next_vertex)) {
      c--;
      d_mst_length += Utils::length(this->points[cur->map[i]], this->points[next_vertex]);
    }
  if(c > 0)
    // Stop - but we may have to look in other direction here!
    return;
  // Check if we are adding a duplet.
  for(i = 0; i < prevSet.size(); i++)
    if(prevSet[i] == next_vertex)
      return;
  // Check if the sharing facet has d edges in MST
  c = 0;
  for(i = 0; i < cur->n; i++)
    if(i != next_simplex_index && this->isInMST(cur->map[i],cur->map[next_simplex_index]))
      c++;
  if(c == 1 && next_simplex < cur_simplex)
    return;
  std::vector<int> curSet = prevSet;
  curSet.push_back(next_vertex);
  // Is connected. Calculate ST and add.
  SteinerTree *st = this->iterCon.insertTerminal(prevTree.st, this->points[next_vertex], prevTree.st->getMSTLength()+d_mst_length);
  if(st->getSteinerRatio() >= 1.0) {
    delete st;
    return;
  }

#if(ESMT_COLLECT_STATS)
  unsigned int size = prevTree.n+1;
  while(this->stats.covered_faces.size() < size+1)
    this->stats.covered_faces.push_back(0);
  this->stats.covered_faces[size]++;
#endif

  SubST subst(st, prevTree.n, prevTree.map, next_vertex);
  this->smts.push_back(subst);  
  added += 1;
  // We build from the largest index, so we have to choose the
  // simplex on the face with the last dim points from cur-set.
  // However, at first, we have different directions, with one
  // less for each added simplex. Thus, we check (opposite) points from
  // size-(this->dim+1) to size-(this->dim+1-max(added,dim)), where
  c = added >= this->dim ? 1 : this->dim - added;
  z = curSet.size()-(this->dim+1);
  for(i = z; i < z+c; i++) {
    // Swap first element with i
    int tmp   = curSet[z];
    curSet[z] = curSet[i];
    curSet[i] = tmp;
    // Get index (j) for next simplex
    for(j = 0; j < next->n; j++)
      if(next->map[j] == curSet[z])
	break;
    // Now, next vertex is in next->nextVertex[j]. Recursive call now
    // if next simplex is less than bound.
    this->buildSausage(curSet, subst, added, j, next_simplex, org_simplex, 1);
  }
  // We have to start the concatenation in the other direction.
  // prevSet[0..(added-1)] contains the set indices, which must be respected.
  // If(added==dim+1), direction is uniqly determined, otherwise, we must make
  // A permutation of prevSet[added..dim]. First exit will be prevSet[dim].
  int start = curSet[0];
  // Do not build if neighbour point is <= start
  for(i = this->dim; i >= added; i--) {
    if(curSet[i] <= start)
      continue;
    int tmp           = curSet[this->dim];
    curSet[this->dim] = curSet[i];
    curSet[i]         = tmp;
    // Get index (j) for next simplex
    for(j = 0; j < (*this->simplices)[org_simplex].n; j++)
      if((*this->simplices)[org_simplex].map[j] == curSet[this->dim])
	break;
    // j is now index of the next connection
    this->buildSausageReverse(curSet, subst, added, 0, j, org_simplex);
    }
}*/

/*void ESMT::buildSausageReverse(std::vector<int> &prevSet,
			       SubST &prevTree,
			       unsigned int prev_added,
			       unsigned int added,
			       unsigned int next_simplex_index,
			       int cur_simplex) {
  unsigned int i, j, z, c;
  double d_mst_length = 0.0;
  Simplex *cur  = &(*this->simplices)[cur_simplex];
  int next_vertex = cur->nextVertex[next_simplex_index];
  int next_simplex = cur->nextSimplex[next_simplex_index];
  // Check for neighbour?
  if(next_vertex < 0)
    return;

  Simplex *next = &(*this->simplices)[next_simplex];

  // Detect if we are connected - must have at least 1 connection
  for(i = 0; i < cur->n; i++)
    if(this->isInMST(cur->map[i], next_vertex)) {
      d_mst_length = Utils::length(this->points[cur->map[i]], this->points[next_vertex]);
      break;
    }
  if(i >= cur->n)
    // Not connected
    return;
  // Check if we are adding a duplet.
  for(i = 0; i < prevSet.size(); i++)
    if(prevSet[i] == next_vertex)
      return;
  // Check if the sharing facet has d edges in MST
  c = 0;
  for(i = 0; i < cur->n; i++)
    if(i != next_simplex_index && this->isInMST(cur->map[i],cur->map[next_simplex_index]))
      c++;
  if(c == 1 && next_simplex < cur_simplex)
    return;
  // Is connected. Calculate ST and add.
  std::vector<int> curSet;
  // Push front because of reverse build
  curSet.reserve(prevSet.size()+1);
  curSet.push_back(next_vertex);
  curSet.insert(curSet.end(), prevSet.begin(), prevSet.end());
  SteinerTree *st = this->iterCon.insertTerminal(prevTree.st, this->points[next_vertex], prevTree.st->getMSTLength()+d_mst_length);

  if(st->getSteinerRatio() >= 1.0) {
    delete st;
    return;
  }

#if(ESMT_COLLECT_STATS)
  unsigned int size = prevTree.n+1;
  while(this->stats.covered_faces.size() < size+1)
    this->stats.covered_faces.push_back(0);
  this->stats.covered_faces[size]++;
#endif

  SubST subst(st, prevTree.n, prevTree.map, next_vertex);
  this->smts.push_back(subst);
  added += 1;
  // We now build in reverse. Given curSet of the form:
  // curSet[0..(1st index)..((d+1)th index)..n]
  // = curSet[0..added..(added+this->dim+1)..n]
  //   c = added+prev_added is in the middle.
  // Consider all between c and (added+this->dim+1) for next direction.
  c = added+prev_added > this->dim ? this->dim : added+prev_added;
  z = this->dim;
  for(i = this->dim; i >= c; i--) {
    // Swap first element with i
    int tmp   = curSet[z];
    curSet[z] = curSet[i];
    curSet[i] = tmp;
    // Get index (j) for next simplex
    for(j = 0; j < next->n; j++)
      if(next->map[j] == curSet[z])
	break;
    // Now, next vertex is in next->nextVertex[j]. Recursive call now
    // if next simplex is less than bound.
    this->buildSausageReverse(curSet, subst, prev_added, added, j, next_simplex);
    }
}*/

/* Stats constructor */
ESMT::Stats::Stats() {
  this->no_of_simplices            = 0;
  this->sub_trees_in_queue         = 0;
  this->add_sub_trees_total        = 0;
  this->no_of_sp                   = 0;
  this->no_of_sp_post_optimisation = 0;
  this->no_of_sp_overlapping       = 0;
}
