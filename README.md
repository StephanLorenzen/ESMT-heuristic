# Heuristic for the Euclidean Steiner Minimal Tree in any dimension

![](esmt-3d.png)

This heuristic is an extended version of the heuristic presented in the paper

* [Euclidean Steiner Tree Heuristic in d-Space](http://dimacs11.cs.princeton.edu/workshop/OlsenLorenzenFonsecaWinter.pdf). A.E. Olsen, S.S. Lorenzen, R. Fonseca and P. Winter.

The code for the original heuristic can be found [here](https://github.com/RasmusFonseca/ESMT-heuristic).  

The heuristic finds solutions to the Euclidean Steiner Minimal Tree (ESTP) problem. In ESTP, we seek a network of minimal total edge length spanning a set of n terminal points while allowing for the insertion of additional points (Steiner points) to decrease the overall length of the network. 

In this extended version of the heuristic we explore the possible gains from using a Bottleneck Graph to determine good candidate sub-Steiner trees in the concatenation part of the heuristic. For more info about the extended heuristic, please see: 

* [Steiner Tree Heuristic in the Euclidean d-Space Using Bottleneck Distances]().

# Compiling

```
$ cd src
$ make
```

The executable depends on having the qdelaunay executable from [qhull](http://www.qhull.org) in the systems PATH. Theres an easy-to-follow explanation [in the wiki](http://github.com/RasmusFonseca/ESMT-heuristic/wiki/qdelaunay) for the original heuristic.

# Usage

```
Usage: 
 esmt-heuristic [Options] <file> <name>
 esmt-heuristic [Options] <type> <d> <n>

file is a STP file and name is the name of a set in that file.
type is a set type (random, grid, sausage), d is the dimension and n is the number of
points to generate.

Options:
 -h    --help               Prints this help
 -v    --verbose            Print info during execution.
 -npo  --nopostopt          Disable fine tuning.
 -nsc  --nosubcon           Disable subgraph concatenation/sausages.
 -rdc  --redocon            Enable concatenation redo.
 -ubd  --usebdist     type  Use bottleneck distances. type must be 1 (table computation),
                            2 (lazy computation) or 3 (dyn. trees).
 -k    --facemaxsize  type  Set the maximum face size to consider for concatenation.
                            Only applicable when using -ubd.
 -alg  --alg          name  Set the subgraph heuristic to use. Must be NO, RNO or SP.
 -s    --seed         s     Use seed s for generating random point sets.
 -pt   --printtree          Print tree when done.
 -val  --validate           Validate tree when done.
 -st   --stats              Collect extra stats.
```

The important option is -ubd, which enables bottleneck distances. Recommended value is 3.
Use the -k option to reduce the number of faces considered when using bottleneck distances in higher dimensions.

# Examples
Solve one of the ESTEIN-3D library instances using bottleneck instances (dyn. tree implementation):
```
$ cd src
$ ./esmt-heuristic -ubd 3 ../data/eSteiner-3d/estein250.stp estein250-03
Done!
  |MST| = 27.1691
  |SMT| = 25.7987
  Ratio = 0.949558
```

Time the computation for 1000 random points in 4 dimensions when using the restricted numerical optimisation (RNO) heuristic for finding small FSTs:
```
$ cd src
$ time ./esmt-heuristic -ubd 3 -alg RNO random 4 1000
Done!
  |MST| = 25611.4
  |SMT| = 23952.3
  Ratio = 0.935222

real	0m9.913s
user	0m9.880s
sys	0m0.037s
```
