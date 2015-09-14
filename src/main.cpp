/**
 * Test of the gdl_window
 */
#include <string>
#include <vector>
#include <iostream>
#include <stdlib.h>

#include "libclap/clap.hpp"
#include "test/test.hpp"
#include "steiner/utils/delaunay.hpp"
#include "steiner/utils/utils.hpp"
#include "steiner/utils/point.hpp"
#include "steiner/utils/point_set_generator.hpp"
#include "steiner/graph.hpp"
#include "steiner/esmt.hpp"
#include "steiner/steiner_tree.hpp"
#include "steiner/heuristics/subgraph_heuristic.hpp"
#include "steiner/heuristics/steiner_finder.hpp"
#include "steiner/heuristics/smith.hpp"
#include "steiner/heuristics/concat.hpp"

#define MAX_DIM     12

typedef Utils::Point      Point;
typedef Utils::Edge       Edge;
typedef Utils::Generator  Generator;

int set_type_from_name(std::string &settype) {
  int t = -1;
  if(settype.compare(0, 20, "random") == 0)
    t = TEST_POINT_SET_RANDOM_D;
  else if(settype.compare(0, 20, "simplex") == 0)
    t = TEST_POINT_SET_SIMPLEX;
  else if(settype.compare(0, 20, "sausage") == 0)
    t = TEST_POINT_SET_SAUSAGE;
  else if(settype.compare(0, 20, "grid") == 0)
    t = TEST_POINT_SET_GRID;
  else if(settype.compare(0, 20, "delaunay") == 0)
    t = TEST_POINT_SET_DELAUNAY_SIMPLICES;
  return t;
}

int main(int argc, char *argv[]) {
  const std::string clap_config =
    "OPTIONS:                                                              \n"
    "-os   --onlysgh                                                       \n"
    "-ct   --counttest  br dimension:i num_points:i seed:i                 \n"
    "-tt   --testtopo   br                                                 \n"
    "-v    --verbose                                                       \n"
    "-npo  --nopostopt                                                     \n"
    "-nsc  --nosubcon                                                      \n"
    "-rdc  --redocon                                                       \n"
    "-alg  --alg           name:s                                          \n"
    "-s    --seed          value:i                                         \n"
    "-pt   --printtree                                                     \n"
    "-val  --validate                                                      \n"
    "-nd   --nodelaunay                                                    \n"
    "-st   --stats                                                         \n"
    "-out  --outfile       path:s                                          \n"
    "-lt   --looptime      sec:i                                           \n"
    "-b    --batch                                                         \n"
    "-g    --generate      set_type:s dimension:i num_points:i             \n"
    "-gn   --generateN     set_type:s dimension:i num_points:i num_tests:i \n"
    "-if   --infile        file:s set_name:s                               \n"
    "-ifa  --infileall     file:s                                          \n"
    "PARAMETERS:                                                           \n"
    "                                                                      \n"
    "file:s set_name:s                                                     \n"
    "set_type:s dimension:i num_points:i                                   \n";
    
  CLAP c(clap_config, argc, argv);
  
  if(c.is_set("testtopo")) {
    std::cout << "Running topo test!" << std::endl;
    std::string name("Smith iterative");
    std::string resFile("/home/stephan/Desktop/result.dat");
    std::string path("../data/protein3D/W1.stp");
    std::string sname("W1");
    std::vector<Point> points1 = Generator::loadFromFile(path, sname);
    path = "../data/protein3D/W3.stp";
    sname = "W3";
    std::vector<Point> points2 = Generator::loadFromFile(path, sname);
    
    ESMT e1(points1);
    ESMT e2(points2);
    
    Test test;
    std::cout << test.testTopology(e1,e2) << std::endl;
    
    return 0;
  }
  else if(c.is_set("counttest")) {
    int d = c.get_int_param("counttest",0),
      n = c.get_int_param("counttest",1),
      seed = c.get_int_param("counttest",2);
    Test test;
    test.doTestESMTSpecial(d, n, seed);
    return 0;
  }
  else {
    bool verbose = c.is_set("v");
    bool print   = c.is_set("pt");
    int n, t = -1, dim, numtests;
    std::string algname, outfile, filename, settype, setname;
    SubgraphHeuristic *sh;
    SubgraphHeuristic *shsp = NULL;
    algname = c.is_set("alg") ? c.get_string_param("alg") : "NO";
    if(algname.compare(0, 10, "NO") == 0)
      sh = new IterativeSmith(MAX_DIM);
    else if(algname.compare(0, 10, "RNO") == 0)
      sh = new IterativeConcat(MAX_DIM);
    else if(algname.compare(0, 10, "SP") == 0) {
      shsp = new IterativeConcat(MAX_DIM);
      sh = new SteinerFinder(shsp);
    }
    else {
      c.error_usage("Unknown sub-graph algorithm");
      return 1;
    }

    // Heuristic parameters
    bool no_subgraph_concat = c.is_set("nsc"),
      no_post_optimise = c.is_set("npo"),
      redo_concat = c.is_set("rdc");
    int seed = c.is_set("s") ? c.get_int_param("s") : time(NULL);
    
    int p = c.get_chosen_pattern();
    if(p == 0 && !c.is_set("batch")) {
      delete sh;
      delete shsp;
      c.error_usage("Please use --batch when no inputs given.");
    }
    else if(p == 0) {      
      // Check that some set is added
      if(!c.is_set("g") && !c.is_set("gn") && !c.is_set("if") && !c.is_set("ifa")) {
	delete sh;
	delete shsp;
	c.error_usage("No input set specified");
      }
    
      Test test;
      test.setSubgraphHeuristicOne(sh, algname);
      test.doConcatSubgraphs(!no_subgraph_concat);
      test.doPostOptimise(!no_post_optimise);
      test.doUseSpecialConcat(redo_concat);
      test.inclDelaunay(!c.is_set("nd"));
      test.setLoopTime(c.is_set("lt") ? c.get_int_param("lt") : 0);
      test.setSeed(seed);
      test.doCollectStats(c.is_set("st"));
      outfile = c.is_set("out") ? c.get_string_param("out") : "";
      
      unsigned int i, num_g = c.is_set("g"), num_gn = c.is_set("gn"),
	num_if = c.is_set("if"), num_ifa = c.is_set("ifa");
      
      numtests = 1;
      for(i = 0; i < num_g; i++) {
	settype  = c.get_string_param("g", 0, i);
	dim      = c.get_int_param("g", 1, i);
	n        = c.get_int_param("g", 2, i);
	t        = set_type_from_name(settype);
	if(t < 0)
	  c.error_usage("Unknown point set type");
	test.addTest(t, dim, n, numtests);
      }
      for(i = 0; i < num_gn; i++) {
	settype  = c.get_string_param("gn", 0, i);
	dim      = c.get_int_param("gn", 1, i);
	n        = c.get_int_param("gn", 2, i);
	numtests = c.get_int_param("gn", 3, i);
	t        = set_type_from_name(settype);
	if(t < 0)
	  c.error_usage("Unknown point set type");
	test.addTest(t, dim, n, numtests);
      }
      for(i = 0; i < num_if; i++) {
	filename = c.get_string_param("if", 0, i);
	setname  = c.get_string_param("if", 1, i);
	test.addTest(filename, setname);
      }
      for(i = 0; i < num_ifa; i++) {
	filename = c.get_string_param("ifa", 0, i);
	test.addTest(filename);
      }
      
      if(c.is_set("os"))
	test.doTestSubgraphAlgorithm(verbose);
      else
	test.doTestESMT(verbose);
      if(outfile.size() > 0)
	test.createDatFile(outfile);
    }
    else {
      std::vector<Point> points;
      if(p == 1) {
	filename = c.get_string_param("file");
	setname  = c.get_string_param("set_name");
	points = Generator::loadFromFile(filename,setname);
      }
      else if(p == 2) {
	settype = c.get_string_param("set_type");
	dim     = c.get_int_param("dimension");
	n       = c.get_int_param("num_points");
	Generator::setSeed(seed);
	if(settype.compare(0, 20, "random") == 0)
	  points = Generator::randomFloatPoints(Point(dim, 100), Point(dim, -100), n);
	else if(settype.compare(0, 20, "simplex") == 0)
	  points = Generator::simplex(dim);
	else if(settype.compare(0, 20, "sausage") == 0)
	  points = Generator::sausage(dim, n);
	else if(settype.compare(0, 20, "grid") == 0)
	  points = Generator::gridFromSide(dim, n);
	else
	  c.error_usage("Unknown point set type");
      }
      
      ESMT esmt(points, sh, !no_subgraph_concat,
		!no_post_optimise, redo_concat, verbose);
      
      if(c.is_set("val")) {
	if(Utils::validate(esmt))
	  std::cout << "Validate OK!" << std::endl;
	else
	  std::cout << "Validate ERROR!" << std::endl;
      }
      
      if(print)
	std::cout << "## RESULT ##" << std::endl << esmt << std::endl;
      else {
	std::cout << "Done!" << std::endl
		  << "  |MST| = " << esmt.getMSTLength() << std::endl
		  << "  |SMT| = " << esmt.getSMTLength() << std::endl
		  << "  Ratio = " << esmt.getSteinerRatio() << std::endl;
      }
    }
    
    delete sh;
    delete shsp;
    
    return 0;
  }
}


/*
void usage() {
  std::cout << "Usage:" << std::endl
	    << " esmt [options] <points>" << std::endl
	    << " test esmt [options] <points>" << std::endl
	    << " test sgh [options] <points>" << std::endl
	    << " counttest <dimension> <no of points> <seed>" << std::endl
	    << std::endl
	    << " options:" << std::endl
	    << "   -v           Verbose" << std::endl
	    << "   -npo         Disable post optimisation" << std::endl
	    << "   -nsc         Disable sub-graph concatenation" << std::endl
	    << "   -sct         Redo concatenation, adding other, non-covered FSTs" << std::endl
	    << "   -alg <name>  Sub-graph algorithm:" << std::endl
	    << "                  NO  = Numerical optimisation" << std::endl
	    << "                  RNO = Restricted numerical optimisation" << std::endl
	    << "                  SP  = Simplex partitioning" << std::endl
	    << "                  Default NO" << std::endl
	    << "   -s           Seed for random generation of point sets." << std::endl
	    << "  # esmt only" << std::endl
	    << "   -pt          Print resulting tree" << std::endl
	    << "   -val         Validate the resulting tree" << std::endl
	    << "  # test only" << std::endl
	    << "   -nd          Do not include Delaunay tesselation in time measurement." << std::endl
	    << "   -cs          Collect extra statistics (no of simplices, etc.)." << std::endl
	    << "   -out <path>  Output file for test results (CSV format)." << std::endl
	    << "   -lt <sec>    Loop time" << std::endl
	    << std::endl
	    << " points:" << std::endl
	    << "   -in <file name> <set name>          Read set <set name>" << std::endl
	    << "                                       from file <file name>" << std::endl
	    << "   -g <set name> <dim> <no of points>  Generate a point set with"<<std::endl
	    << "                                       the given number of points in"<<std::endl
	    << "                                       the given dimension."<<std::endl
	    << "                                       Possible set names:" << std::endl
	    << "                                           random, sausage, simplex, grid" << std::endl
	    << " points (test only):" << std::endl
	    << "   -ina <file name>              Read all sets from file <file name>" << std::endl
	    << "   -gn <set name> <dim> <no of points> <no of tests>" << std::endl
	    << "                                 Generate the given number of" << std::endl
	    << "                                 point set with the given"<<std::endl
	    << "                                 number of points in the given"<<std::endl
	    << "                                 dimension."<<std::endl;
  exit(1);
}
*/
