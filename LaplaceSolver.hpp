#ifndef _LAPLACESOLVER_HPP
#define _LAPLACESOLVER_HPP

#include <fstream>
#include <string>
#include <vector>
#include "Mesh.hpp"
#include "GetPot.hpp"



struct CSR_pattern
{
  std::vector<long int> I;
  std::vector<long int> J;
  size_t n_zero;
  CSR_pattern()
  :I(0,0),
  J(0,0),
  n_zero(0)
  {}
  void clear()
  {
    I.clear();
    J.clear();
    n_zero=0;
  }
};

struct CSR_matrix
{
  CSR_pattern _pattern;
  std::vector<double> K;
  CSR_matrix()
  :K(0,0)
  {}
  void clear(bool clearpattern=0)
  {
    K.clear();
    if(clearpattern)
    {
      _pattern.clear();
    }
    else
    {
      K.resize(_pattern.n_zero,0);
    }
  }
};


class LaplaceSolver
{

 public:
  LaplaceSolver();
  LaplaceSolver(const Mesh * _mesh);
  LaplaceSolver(const GetPot & dfile, const Mesh * _mesh);
  ~LaplaceSolver();
  void setMesh(const Mesh * _mesh);
  // set functions
  inline void set_abs_toll(double toll){ if(toll>0.0){abs_toll=toll;} };
  inline void set_rel_toll(double toll){ if(toll>0.0){rel_toll=toll;} };
  inline void set_max_it(int maxit){ if(maxit>0){itr_max=maxit;} };
  inline void set_Krilov_dim(long int Kdim){ if(Kdim>0){dimKrilovSp=Kdim;} };
  inline void set_verbosity(short int verb) { if(verb>=0){ verbose=verb;} };
  
 private: 
  bool _consistentState;
  const Mesh *  _ptrmesh;
  double abs_toll;
  double rel_toll;
  int itr_max;
  long int dimKrilovSp;
  short int verbose;

  

};
#endif
