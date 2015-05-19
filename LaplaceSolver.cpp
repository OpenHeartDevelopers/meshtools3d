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
    this->evaldphi0();
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


void LaplaceSolver::eval_pattern()
{
  _Matrix._pattern.I.resize(1+_ptrmesh->nPt());
  std::vector<std::set<long int> > labelPerNode;
  labelPerNode.resize(_ptrmesh->nPt());
  // I create a vector of set
  // in each entry the set gives the connectivity
  // i.e. the label of neighb. nodes

  for(size_t iTet=0; iTet<_ptrmesh->nTet(); iTet++)  
  {
    for(short int iPt=0; iPt<4; iPt++)
    {
        for(short int jPt=iPt; jPt<4; jPt++)
        {
          size_t iNode=_ptrmesh->Tet(iTet).vertex[iPt];
          size_t jNode=_ptrmesh->Tet(iTet).vertex[jPt];
          labelPerNode[iNode].insert(jNode);
          labelPerNode[jNode].insert(iNode);
        }
    }
  }
  // I evaluate the number of non_zero elements into the pattern
  _Matrix._pattern.n_zero=0;
  for(size_t iPt=0; iPt<_ptrmesh->nPt(); iPt++)
  {
    _Matrix._pattern.n_zero=_Matrix._pattern.n_zero+labelPerNode[iPt].size();
  }
  _Matrix._pattern.J.resize(_Matrix._pattern.n_zero);
  _Matrix.K.resize(_Matrix._pattern.n_zero,0);
  // Here I  fill the pattern of the matrix
  _Matrix._pattern.I[0]=0;
  for(size_t iPt=0; iPt<_ptrmesh->nPt(); iPt++)
  {
    _Matrix._pattern.I[iPt+1]=_Matrix._pattern.I[iPt]+labelPerNode[iPt].size();
    
    unsigned counter=0;
    for(std::set<long int>::iterator it=labelPerNode[iPt].begin(); it!=labelPerNode[iPt].end(); ++it)
    {
      _Matrix._pattern.J[_Matrix._pattern.I[iPt]+counter]=*it;
      counter=counter+1;
    }
  }
}

std::vector<double> LaplaceSolver::localStiff(size_t iTet, double k)
{
  std::vector<double> lstif(16,0); //4X4 matrix
  double Volume=(_ptrmesh->VolTet(iTet));
  std::vector<double> invJt=(_ptrmesh->TetInvJacobianTransponse(iTet)); // 9X1
  std::vector<std::vector<double> > dphi;
  dphi.resize(4);
  double coef = k*Volume/6.0;
  for(short int jf=0; jf<4; jf++)
  {
    dphi[jf].resize(3,0);
    for(short int ic=0; ic<3; ic++)
    {
      for(short int jc=0; jc<3; jc++)
      {
        (dphi[jf])[ic]=(dphi[jf])[ic]+invJt[RMIndex(ic,jc,3)]  *(dphi0[jf])[jc];
      }
    }
  }
  // Now eval the matrix
  for(short int ir=0; ir<4; ir++)
  {
    for(short int jc=ir; jc<4; jc++)
    {
      double entry=0;
      for(short int icoor=0; icoor<3; icoor++)
      {
        entry=entry+((dphi[ir])[icoor])*((dphi[jc])[icoor]);
      }
      lstif[RMIndex(ir,jc,4)]=coef*entry;
      lstif[RMIndex(jc,ir,4)]=coef*entry;
    }
  }
  return(lstif);
}


LaplaceSolver::~LaplaceSolver()
{
  if(_ptrmesh)
  {
    _ptrmesh=NULL;
  }
  _Matrix.clear(true);
  dphi0.clear();
}


short int LaplaceSolver::RMIndex(short int irow, short int jcol, short int rank)
{
  short int index=irow*rank+jcol;
  return(index);
}

void LaplaceSolver::evaldphi0()
{
  dphi0.resize(4);
  for(short int jf=0; jf<4; jf++)
  {
      dphi0[jf].resize(3);
  }
  for(short int jc=0; jc<3; jc++)
  {
    dphi0[0][jc]=-1;
    dphi0[jc+1][jc]=1;
  }
  
  



}

