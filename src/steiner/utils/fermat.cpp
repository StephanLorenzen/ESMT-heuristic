#include <vector>
#include <cmath>
#include <assert.h>

#include "steiner/graph.hpp"
#include "steiner/steiner_tree.hpp"
#include "steiner/utils/utils.hpp"
#include "steiner/utils/point.hpp"
#include "steiner/utils/fermat.hpp"

#define PI     3.14159265
#define SQRT32 0.866025403784439

typedef Utils::Point        Point;

void calculateEqPoint(Point &A, Point &B, Point &C, Point &res) {
  unsigned int i;
  double cos_A, sin_A, sin_O, ac, l, ext;

  // cos(BAC)
  cos_A = Utils::angle(B, A, C);
  // cos(CAO)
  cos_A = -cos_A;
  // sin(CAO)
  sin_A = std::sqrt(1-cos_A*cos_A);
  // sin(AOC) = sin(pi-CAO-ACO) = sin(pi-1/3*pi-ACO)
  //          = sin(2/3*pi-ACO)
  //          = sin(2/3*pi)cos(ACO) - cos(2/3*pi)sin(ACO)
  //          = sqrt(3)/2 cos(ACO) + 0.5 sin(ACO)
  sin_O   = SQRT32*cos_A + 0.5*sin_A;
  
  ac = Utils::length(A,C);
  ext = (ac / sin_O) * SQRT32;
  
  l = (1.0 / Utils::length(A,B));
  for(i = 0; i < res.dim(); i++)
    res[i] = A[i] + (A[i]-B[i])*l*ext;

  l = (1.0 / Utils::length(C,res));
  for(i = 0; i < res.dim(); i++)
    res[i] = C[i] + (res[i]-C[i])*l*ac;
}

void calculateFermatPoint(Point &corner, Point &p1, Point &p2, Point &res) {
  unsigned int dim = corner.dim(), i = 0, j = 0;
  Point e1(dim), e2(dim);
  calculateEqPoint(corner, p1, p2, e1);
  calculateEqPoint(corner, p2, p1, e2);

  // Calculate intersection of e1-p1 and e2-p2
  double d[2][2];
  double p[2];
  
  double det = 0.0;
  i = 0;
  j = 1;
  while(det == 0.0) {
    assert(i < corner.dim()-1);
    d[0][0] = e2[i] - p2[i];
    d[0][1] = p1[i] - e1[i];
    d[1][0] = e2[j] - p2[j];
    d[1][1] = p1[j] - e1[j];
    
    p[0]    = p1[i] - p2[i];
    p[1]    = p1[j] - p2[j];
    det = d[0][0]*d[1][1] - d[0][1]*d[1][0];
    j++;
    if(j == corner.dim()) {
      i++;
      j = i+1;
    }
  }
  
  double t = (d[1][1]*p[0] - d[0][1]*p[1]) / det;
  for(i = 0; i < dim; i++)
    res[i] = t * (e2[i] - p2[i]) + p2[i];
}

int Utils::getFermatPoint(Point &A, Point &B, Point &C, Point &res) {
  double cos_A, cos_B, cos_C;  
  //unsigned int dim = A.dim();

  // Special cases
  cos_A = Utils::angle(B, A, C);
  if(cos_A < -0.5) {
    res = A;
    return 0;
  }
  cos_B = Utils::angle(A, B, C);
  if(cos_B < -0.5) {
    res = B;
    return 1;
  }
  cos_C = Utils::angle(A, C, B);
  if(cos_C < -0.5) {
    res = C;
    return 2;
  }

  // Special case
  if(cos_A*1.001 <= -0.5 && cos_B*1.001 <= -0.5 && cos_C*1.001 <= -0.5) {
    res = (A + B + C) * (1.0 / 3.0);
    return -1;
  }
  
  if(cos_A < cos_B && cos_A < cos_C) {
    // A is biggest angle
    calculateFermatPoint(A, B, C, res);
  }
  else if(cos_B < cos_A && cos_B < cos_C) {
    // B is biggest angle
    calculateFermatPoint(B, A, C, res);
  }
  else {
    // C is biggest angle
    calculateFermatPoint(C, A, B, res);
  }
  return -1;
}

void Utils::getFermatSMT(SteinerTree &triangle, bool add_sp) {
  assert(triangle.n()==3);
  
  std::vector<Edge> &edges = triangle.getEdges();
  edges.clear();
 
  Point sp(triangle.dimension());
  Point &p0 = triangle.getPoint(0);
  Point &p1 = triangle.getPoint(1);
  Point &p2 = triangle.getPoint(2);
  int add = Utils::getFermatPoint(p0, p1, p2, sp);
  
  if(add == -1 || add_sp) {
    edges.push_back(Edge(0,3, Utils::length(p0, sp)));
    edges.push_back(Edge(1,3, Utils::length(p1, sp)));
    edges.push_back(Edge(2,3, Utils::length(p2, sp)));
    sp.setSteiner();
    triangle.getSteinerPoints().push_back(sp);
  }
  else {
    unsigned int pi1 = (add+1) % 3, pi2 = (add+2) % 3;
    edges.push_back(Edge(add,pi1, Utils::length(triangle.getPoint(pi1), triangle.getPoint(add))));
    edges.push_back(Edge(add,pi2, Utils::length(triangle.getPoint(pi2), triangle.getPoint(add))));
  }
  triangle.setSMTLength(triangle.getLength());
}
