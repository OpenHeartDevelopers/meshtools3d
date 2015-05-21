#ifndef _LAPLACESOLVER_HPP
#define _LAPLACESOLVER_HPP

#include <fstream>
#include <string>
#include <vector>
#include<map>
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

 typedef std::map<long int, double> BCContainerType;
 typedef BCContainerType::iterator BCContainerTypeIterator;
 public:
  LaplaceSolver();
  LaplaceSolver(const Mesh * _mesh);
  LaplaceSolver(const GetPot & dfile, const Mesh * _mesh);
  void solve();
  ~LaplaceSolver();
  void setMesh(const Mesh * _mesh);
  // set functions
  inline void setBCValue(long int node, double value){DirichletBC.insert(std::pair<long int, double>(node,value));};
  void setBCValue(std::set<long int> region, double value);
  void setBCValue(BCContainerType BCS);
  inline void set_abs_toll(double toll){ if(toll>0.0){abs_toll=toll;} };
  inline void set_rel_toll(double toll){ if(toll>0.0){rel_toll=toll;} };
  inline void set_max_it(int maxit){ if(maxit>0){itr_max=maxit;} };
  inline void set_Krilov_dim(long int Kdim){ if(Kdim>0){dimKrilovSp=Kdim;} };
  inline void set_verbosity(short int verb) {  verbose=verb; };
  inline std::vector<double> & sol(){return _sol;};
  inline double & sol(size_t iP){return (_sol[iP]);};
  inline const std::vector<double> & sol() const {return _sol;};
  inline const double & sol(size_t iP) const {return _sol[iP];};
  void writeSolution(std::string filename);
  std::vector<double> ElementTetraGradient(size_t iEl) const;
 private: 
  //functions
  void eval_pattern();
  void evaldphi0();
  void matrixAssembly(bool build_pattern=1);
  std::vector<double> localStiff(size_t iTet, double k=1.0);
  short int RMIndex(short int irow, short int jcol, short int rank=4) const;
  
  //variables
  bool _consistentState;
  const Mesh *  _ptrmesh;
  double abs_toll;
  double rel_toll;
  int itr_max;
  long int dimKrilovSp;
  short int verbose;
  CSR_matrix _Matrix;
  std::vector<std::vector<double> > dphi0;
  std::vector<double> _sol;
  std::vector<double> _RHS;
  BCContainerType DirichletBC;
  

};
#endif
