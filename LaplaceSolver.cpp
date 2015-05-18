#include "LaplaceSolver.hpp"
#include<iostream>
#include "mgmres.hpp"

LaplaceSolver::LaplaceSolver()
:_consistentState(false),
_ptrmesh(NULL),
abs_toll(1e-6),
rel_toll(1e-6),
itr_max(100),
dimKrilovSp(10),
verbose(0)
{}

LaplaceSolver::LaplaceSolver(const Mesh *  _mesh)
:_consistentState(false),
_ptrmesh(NULL),
abs_toll(1e-6),
rel_toll(1e-6),
itr_max(100),
dimKrilovSp(10),
verbose(0)
{
  this->setMesh(_mesh);
}

void LaplaceSolver::setMesh(const Mesh *  _mesh)
{
  if(_mesh)
  {
    _ptrmesh=_mesh;
    _consistentState=true;
  }
}

LaplaceSolver::LaplaceSolver(const GetPot & dfile, const Mesh * _mesh)
:_consistentState(false),
_ptrmesh(NULL),
abs_toll(1e-6),
rel_toll(1e-6),
itr_max(100),
dimKrilovSp(10),
verbose(0)
{
  this->setMesh(_mesh);
  abs_toll    = dfile("laplacesolver/abs_toll",1.e-6);
  rel_toll    = dfile("laplacesolver/rel_toll",1.e-6);
  itr_max     = dfile("laplacesolver/itr_max",100);    
  dimKrilovSp = dfile("laplacesolver/dimKrilovSp",100);
  verbose     = dfile("laplacesolver/verbose",0);
  if(dimKrilovSp>static_cast<long int>(_ptrmesh->nPt()))
  {
    dimKrilovSp=_ptrmesh->nPt();
    std::cout<<"Warning: Krilov space dimension was bigger than matrix size; changed to "<<_ptrmesh->nPt()<<std::endl;
  }

}



LaplaceSolver::~LaplaceSolver()
{
  if(_ptrmesh)
  {
    _ptrmesh=NULL;
  }
}

