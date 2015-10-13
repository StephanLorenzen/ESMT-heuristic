# Heuristic for the Euclidean Steiner Minimal Tree in any dimension

This heuristic is an extended version of the heuristic presented in the paper

* [Euclidean Steiner Tree Heuristic in d-Space](http://dimacs11.cs.princeton.edu/workshop/OlsenLorenzenFonsecaWinter.pdf). A.E. Olsen, S.S. Lorenzen, R. Fonseca and P. Winter.

The code for the original heuristic can be found [here](https://github.com/RasmusFonseca/ESMT-heuristic).  

The heuristic finds solutions to the Euclidean Steiner Minimal Tree (ESTP) problem. In ESTP, we seek a network of minimal total edge length spanning a set of n terminal points while allowing for the insertion of additional points (Steiner points) to decrease the overall length of the network. 

In this extended version of the heuristic we explore the possible gains from using a Bottleneck Graph to determine good candidate sub-Steiner trees in the concatenation part of the heuristic. For more info about the extended heuristic, please see: 

* [Improved Steiner Tree Heuristic Using Bottleneck Graph]().

# Compiling

```
$ cd src
$ make
```

The executable depends on having the qdelaunay executable from [qhull](http://www.qhull.org) in the systems PATH. Theres an easy-to-follow explanation [in the wiki](http://github.com/RasmusFonseca/ESMT-heuristic/wiki/qdelaunay) for the original heuristic.